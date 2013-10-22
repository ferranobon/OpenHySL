#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <getopt.h>  /* For getopt_long() */
#include <time.h>    /* For time(), ctime() */

#include "Print_Messages.h"

#include "Algorithm_Aux.h"

#include "MatrixVector.h"
#include "Input_Load.h"
#include "New_State.h"
#include "EffK_Formulation.h"
#include "Error_Compensation.h"
#include "Rayleigh.h"
#include "GainMatrix.h"
#include "Substructure.h"
#include "Substructure_Experimental.h"
#include "Substructure_Auxiliary.h"  /* For Substructure_VectorXm(), Substructure_VectorXc(), ... */

#include "Definitions.h"

#if _HDF5_
#include "HDF5_Operations.h"
#endif

#if _MKL_
#include <mkl_pblas.h>
#include "Cblacs.h"
#include "Scalapack_Aux.h"
#else
#include "Netlib.h"
#endif

#include <mpi.h>

const char *Entry_Names[NUM_CHANNELS] = { "Sub-step",
					  "Control displacement actuator 1 [m]",
					  "Control displacement actuator 2 [m]",
					  "Measured displacement actuator 1 [m]",
					  "Measured displacement actuator 2 [m]",
					  "Control acceleration Actuator 1 [m/s^2]",
					  "Control acceleration Actuator 2 [m/s^2]",
					  "Measured acceleration actuator 1 [m/s^2]",
					  "Measured acceleration actuator 2 [m/s^2]",
					  "Acceleration TMD (y-direction) [m/s^2]",
					  "Acceleration TMD (x-direction) [m/s^2]",
					  "Acceleration Base-Frame (x-direction) [m/s^2]",
					  "Coupling Force 1 (y-direction) [N]",
					  "Coupling Force 2 (y-direction) [N]",
					  "Coupling Force 3 (x-direction) [N]",
					  "Displacement TMD (x-direction) [m]",
					  "Displacement TMD (y-direction) [m]",
					  "Control pressure [Pa]",
					  "Measured pressure [Pa]",
					  "Displacement at the coupling node [m]",
					  "Time spent doing the sub-step [ms]", 
					  "Sub-step time [ms]",
					  "Synchronisation time between PC and ADwin [ms]",
					  "Time between the first sub-step and arrival of the new update [ms]"
};

int main( int argc, char **argv ){


     /* MPI Variables */
     int rank, size;
     MPI_Status Status;

     /* BLACS Variables */
     int icontxt, bhandle;
     int nprow, npcol, myrow, mycol;

     /* Algorithm variables */
     unsigned int i, istep;
     AlgConst_t InitCnt;
     const char *FileConf;
     
#if _HDF5_
     hid_t hdf5_file;
#endif

     /* NETLIB Variables */
     int ione = 1;
     Scalars_t Constants;
     
     double *AccAll, *VelAll, *DispAll;

     PMatrixVector_t M, C, K;               /* Mass, Damping and Stiffness matrices */
     PMatrixVector_t Keinv;
     PMatrixVector_t Keinv_m;
     MatrixVector_t Keinv_c;

     PMatrixVector_t EffT;

     PMatrixVector_t DispT, DispTdT0, DispTdT;
     PMatrixVector_t DispTdT0_m;
     MatrixVector_t DispTdT0_c;
     PMatrixVector_t tempvec;

     PMatrixVector_t VelT, VelTdT;
     PMatrixVector_t AccT, AccTdT;

     PMatrixVector_t LoadVectorForm, LoadTdT, LoadTdT1;

     PMatrixVector_t fc, fcprevsub;
     PMatrixVector_t fu;

     PMatrixVector_t ErrCompT, ErrCompTdT;

     PMatrixVector_t Disp, Vel, Acc;

     CouplingNode_t CNodes;
     SaveTime_MPI_t Time;

     /* Options */
     int Selected_Option;
     struct option long_options[] = {
	  {"help", no_argument, 0, 'h'},
	  {"config-file", required_argument, 0, 'c'},
	  {0, 0, 0, 0}
     };

     MPI_Init( &argc, &argv );
     MPI_Comm_size( MPI_COMM_WORLD, &size );
     MPI_Comm_rank( MPI_COMM_WORLD, &rank );


     if( rank == 0 ){
	  /* Print Information */
	  printf( "\n\n" );
	  printf( "************************************************************\n" );
	  printf( "*                                                          *\n" );
	  printf( "*  This is Dorka's substructure algorithm as programed by  *\n" );
	  printf( "* Ferran Obón Santacana. Version alpha 0.5 'Heaven's Door' *\n" );
	  printf( "*                                                          *\n" );
	  printf( "************************************************************\n\n" );
	  /* Set de default value for the configuration file */
	  FileConf = "ConfFile_MPI.conf";

	  /* This is only used if there are no arguments */
	  if( argc == 1 ){
	       Print_Header( INFO );
	       printf( "Assuming the configuration file to be: %s\n", FileConf );

	  }

	  while ((Selected_Option = getopt_long( argc, argv, "c:h", long_options, NULL )) != -1 ){
	       switch( Selected_Option ){
	       case 'c':
		    FileConf = optarg;
		    break;
	       case 'h':
		    Algorithm_PrintHelp( argv[0] );
		    return EXIT_FAILURE;
		    break;
	       case '?':
		    /* Long options already prints an error message telling that there is an unrecognised
		     * option */
		    Algorithm_PrintHelp( argv[0] );
		    return EXIT_FAILURE;
	       case ':':
		    /* Long options already prints an error message telling that the option requieres an
		     * argument */
		    Algorithm_PrintHelp( argv[0] );
		    return EXIT_FAILURE;
	       }
	  }

	  /* Constants definitions. */
	  Algorithm_Init_MPI( FileConf, &InitCnt );

	  /* Read the coupling nodes from a file */
	  Substructure_ReadCouplingNodes( &InitCnt, &CNodes );
     }

     /* Send the Information read in process 0 to all the other processes using MPI_Bcast */
     Algorithm_BroadcastConfFile( &InitCnt );
     Substructure_BroadCastCouplingNodes( &CNodes );

     /* Allocate memory for saving the acceleration, displacement and velocity (input files) that will be used
      * during the test */
     if( InitCnt.Use_Absolute_Values ){
	  AccAll = NULL;
	  VelAll = (double *) calloc( (size_t) InitCnt.NStep, sizeof(double) );
	  DispAll = (double *) calloc( (size_t) InitCnt.NStep, sizeof(double) );
     } else {
	  AccAll = (double *) calloc( (size_t) InitCnt.NStep, sizeof(double) );
	  VelAll = NULL;
	  DispAll = NULL;
     }

     /* BLACS: Initialise BLACS */
     bhandle = Csys2blacs_handle( MPI_COMM_WORLD );
     icontxt = bhandle;

     /* BLACS: Initialise the process grid */
     Cblacs_gridinit( &icontxt, (char *) "Row", InitCnt.ProcessGrid.Rows, InitCnt.ProcessGrid.Cols );
     Cblacs_gridinfo( icontxt, &nprow, &npcol, &myrow, &mycol );

     /* Initialise the matrices and vectors that will be used in the Time Integration process */
     PMatrixVector_Create( icontxt, InitCnt.Order, InitCnt.Order, InitCnt.BlockSize.Rows,
			   InitCnt.BlockSize.Cols, &M );
     PMatrixVector_Create( icontxt, InitCnt.Order, InitCnt.Order, InitCnt.BlockSize.Rows,
			   InitCnt.BlockSize.Cols, &C );
     PMatrixVector_Create( icontxt, InitCnt.Order, InitCnt.Order, InitCnt.BlockSize.Rows,
			   InitCnt.BlockSize.Cols, &K );

     PMatrixVector_Create( icontxt, InitCnt.Order, InitCnt.Order, InitCnt.BlockSize.Rows,
			   InitCnt.BlockSize.Cols, &Keinv );

     if( CNodes.Order >= 1 ){
	  /* The coupling elements are not distributed across the processes and they are always located at
	   * process 0 */
	  if( rank == 0 ){
	       MatrixVector_Create( CNodes.Order, CNodes.Order, &Keinv_c );
	  }
	  PMatrixVector_Create( icontxt, InitCnt.Order - CNodes.Order, CNodes.Order, InitCnt.BlockSize.Rows,
				InitCnt.BlockSize.Cols, &Keinv_m );
     }

     PMatrixVector_Create( icontxt, InitCnt.Order, 1, InitCnt.BlockSize.Rows, InitCnt.BlockSize.Cols,
			   &tempvec );
     PMatrixVector_Create( icontxt, InitCnt.Order, 1, InitCnt.BlockSize.Rows, InitCnt.BlockSize.Cols,
			   &DispT );
     PMatrixVector_Create( icontxt, InitCnt.Order, 1, InitCnt.BlockSize.Rows, InitCnt.BlockSize.Cols,
			   &DispTdT0 );
     PMatrixVector_Create( icontxt, InitCnt.Order, 1, InitCnt.BlockSize.Rows, InitCnt.BlockSize.Cols,
			   &DispTdT );

     if( CNodes.Order >= 1 ){
	  /* The coupling nodes will be always located at rank 0 */
	  if( rank == 0 ){
	       MatrixVector_Create( CNodes.Order, 1, &DispTdT0_c );
	  }
	  PMatrixVector_Create( icontxt, InitCnt.Order - CNodes.Order, 1, InitCnt.BlockSize.Rows,
				InitCnt.BlockSize.Cols, &DispTdT0_m );
     }

     PMatrixVector_Create( icontxt, InitCnt.Order, 1, InitCnt.BlockSize.Rows, InitCnt.BlockSize.Cols, &VelT );
     PMatrixVector_Create( icontxt, InitCnt.Order, 1, InitCnt.BlockSize.Rows, InitCnt.BlockSize.Cols,
			   &VelTdT );

     PMatrixVector_Create( icontxt, InitCnt.Order, 1, InitCnt.BlockSize.Rows, InitCnt.BlockSize.Cols, &AccT );
     PMatrixVector_Create( icontxt, InitCnt.Order, 1, InitCnt.BlockSize.Rows, InitCnt.BlockSize.Cols,
			   &AccTdT );

     PMatrixVector_Create( icontxt, InitCnt.Order, 1, InitCnt.BlockSize.Rows, InitCnt.BlockSize.Cols, &EffT );

     PMatrixVector_Create( icontxt, InitCnt.Order, 1, InitCnt.BlockSize.Rows, InitCnt.BlockSize.Cols,
			   &LoadVectorForm );
     PMatrixVector_Create( icontxt, InitCnt.Order, 1, InitCnt.BlockSize.Rows, InitCnt.BlockSize.Cols,
			   &LoadTdT );
     PMatrixVector_Create( icontxt, InitCnt.Order, 1, InitCnt.BlockSize.Rows, InitCnt.BlockSize.Cols,
			   &LoadTdT1 );

     PMatrixVector_Create( icontxt, InitCnt.Order, 1, InitCnt.BlockSize.Rows, InitCnt.BlockSize.Cols, &fc );

     if( CNodes.Order >= 1 ){
	  PMatrixVector_Create( icontxt, CNodes.Order, 1, InitCnt.BlockSize.Rows, InitCnt.BlockSize.Cols,
				&fcprevsub );
     }

     PMatrixVector_Create( icontxt, InitCnt.Order, 1, InitCnt.BlockSize.Rows, InitCnt.BlockSize.Cols, &fu );


     PMatrixVector_Create( icontxt, InitCnt.Order, 1, InitCnt.BlockSize.Rows, InitCnt.BlockSize.Cols,
			   &ErrCompT );
     PMatrixVector_Create( icontxt, InitCnt.Order, 1, InitCnt.BlockSize.Rows, InitCnt.BlockSize.Cols,
			   &ErrCompTdT );

     PMatrixVector_Create( icontxt, InitCnt.Order, 1, InitCnt.BlockSize.Rows, InitCnt.BlockSize.Cols, &Disp );
     PMatrixVector_Create( icontxt, InitCnt.Order, 1, InitCnt.BlockSize.Rows, InitCnt.BlockSize.Cols, &Vel );
     PMatrixVector_Create( icontxt, InitCnt.Order, 1, InitCnt.BlockSize.Rows, InitCnt.BlockSize.Cols, &Acc );

     /* Read the matrices from a file */
     if( !InitCnt.Read_Sparse ){
	  PMatrixVector_FromFile( InitCnt.FileM, &M );
	  PMatrixVector_FromFile( InitCnt.FileK, &K );
     } else if ( InitCnt.Read_Sparse ){
	  PMatrixVector_FromFile_MM( InitCnt.FileM, &M );
	  PMatrixVector_FromFile_MM( InitCnt.FileK, &K );
     } else assert(0);

     if ( !InitCnt.Read_LVector ){
	  InputLoad_Generate_LoadVectorForm_MPI( InitCnt.ExcitedDOF, &LoadVectorForm );
     }  else {
	  PMatrixVector_FromFile_MM( InitCnt.FileLV, &LoadVectorForm );
     }

     /* Damping matrix */
     Rayleigh_Damping_MPI( &M, &K, &C, &InitCnt.Rayleigh );

     /* Calculate Matrix Keinv = 1.0*[K + a0*M + a1*C]^(-1) */
     Constants.Alpha = InitCnt.a0;    /* Mass matrix coefficient */
     Constants.Beta = InitCnt.a1;     /* Damping matrix coefficent */
     Constants.Gamma = 1.0;           /* Stiffness matrix coefficient */
     Constants.Lambda = 1.0;          /* Matrix inversion coefficient */
     
     /* Gain Matrix */
     IGainMatrix_MPI( &Keinv, &M, &C, &K, Constants );
     
     if( CNodes.Order >= 1 ){
	  Substructure_MatrixXc_MPI( MPI_COMM_WORLD, &CNodes, &Keinv, &Keinv_c );
	  Substructure_MatrixXcm_MPI( MPI_COMM_WORLD, &Keinv, &CNodes, &Keinv_m );  
	  if ( rank == 0 ){
	       /* Send the coupling part of the effective matrix if we are performing a distributed test */
	       for( i = 0; i < (unsigned int) CNodes.Order; i++ ){
		    if( CNodes.Sub[i].Type == EXP_ADWIN || CNodes.Sub[i].Type == REMOTE ){
			 Substructure_SendGainMatrix( Keinv_c.Array, (unsigned int) CNodes.Order, &CNodes.Sub[i] );
		    }
	       }
	  }
     }

     /* Read the earthquake data from a file */
     if( rank == 0 ){
	  if( InitCnt.Use_Absolute_Values ){
	       Algorithm_ReadDataEarthquake_AbsValues( InitCnt.NStep, InitCnt.FileData, InitCnt.Scale_Factor,
						  VelAll, DispAll );
	  } else {
	       Algorithm_ReadDataEarthquake_RelValues( InitCnt.NStep, InitCnt.FileData, InitCnt.Scale_Factor,
						       AccAll );
	  }
     }

     /* Broadcast the earthquake data */
     if( InitCnt.Use_Absolute_Values ){
	  MPI_Bcast( DispAll, (int) InitCnt.NStep, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	  MPI_Bcast( VelAll, (int) InitCnt.NStep, MPI_DOUBLE, 0, MPI_COMM_WORLD );
     } else {
	  MPI_Bcast( AccAll, (int) InitCnt.NStep, MPI_DOUBLE, 0, MPI_COMM_WORLD );
     }

     /* Open Output file. If the file cannot be opened, the program will exit, since the results cannot be
      * stored. */
#if _HDF5_
     HDF5_CreateFile_MPI( MPI_COMM_WORLD, InitCnt.FileOutput, &hdf5_file );

     HDF5_CreateGroup_Parameters( hdf5_file, &InitCnt, &CNodes, AccAll, VelAll, DispAll );
     HDF5_CreateGroup_TimeIntegration( hdf5_file, &InitCnt );
#endif

     Time.Start = MPI_Wtime();
     Time.Date_start = time( NULL );
     strcpy( Time.Date_time, ctime( &Time.Date_start ) );

     /* Calculate the input load */
     istep = 1;
     if( InitCnt.Use_Absolute_Values ){
	  InputLoad_Apply_LoadVectorForm_MPI( &LoadVectorForm, DispAll[istep - 1], &Disp );
	  InputLoad_Apply_LoadVectorForm_MPI( &LoadVectorForm, VelAll[istep - 1], &Vel );  
	  InputLoad_AbsValues_MPI( &K, &C, &Disp, &Vel, &LoadTdT );
     } else {
	  InputLoad_Apply_LoadVectorForm_MPI( &LoadVectorForm, AccAll[istep - 1], &Acc );
	  InputLoad_RelValues_MPI( &M, &Acc, &LoadTdT );
     }

     if( rank == 0 ){
	  Print_Header( INFO );
	  printf( "Starting stepping process.\n" );
     }

     while ( istep <= InitCnt.NStep ){

	  /* Calculate the effective force vector Fe = M*(a0*u + a2*v + a3*a) + C*(a1*u + a4*v + a5*a) */
	  EffK_EffectiveForce_MPI( &M, &C, &DispT, &VelT, &AccT, &tempvec, InitCnt.a0, InitCnt.a1, InitCnt.a2,
				   InitCnt.a3, InitCnt.a4, InitCnt.a5, &EffT );

	  /* Compute the new Displacement u0 */
	  Compute_NewState_MPI( &Keinv, &EffT, &LoadTdT, &fu, &tempvec, &DispTdT0 );

	  /* Split DispTdT into coupling and non-coupling part */
	  if( CNodes.Order >= 1 ){
	       Substructure_VectorXm_MPI( &DispTdT0, &CNodes, &DispTdT0_m );
	       Substructure_VectorXc_MPI( MPI_COMM_WORLD, &DispTdT0, &CNodes, &DispTdT0_c );
	  }

	  /* Perform substepping */
	  if( CNodes.Order >= 1 ){
/*	       Substructure_Substepping_MPI( &CNodes, Keinv_c.Array, DispTdT0_c.Array,
					     InitCnt.Delta_t*(double) istep, AccAll[istep - 1],
					     InitCnt.NSubstep, InitCnt.DeltaT_Sub, DispTdT.Array,
					     fcprevsub.Array, fc.Array);*/
	       Substructure_Substepping_MPI( &CNodes, Keinv_c.Array, DispTdT0_c.Array,
	                                     InitCnt.Delta_t*(double) istep, 0.0, InitCnt.NSubstep, 
					     InitCnt.DeltaT_Sub, DispTdT.Array, fcprevsub.Array, fc.Array );
	  }

	  if ( istep < InitCnt.NStep ){
	       /* Calculate the input load for the next step during the sub-stepping process. */
	       if( InitCnt.Use_Absolute_Values ){
		    InputLoad_Apply_LoadVectorForm_MPI( &LoadVectorForm, DispAll[istep], &Disp );
		    InputLoad_Apply_LoadVectorForm_MPI( &LoadVectorForm, VelAll[istep], &Vel );
		    InputLoad_AbsValues_MPI( &K, &C, &Disp, &Vel, &LoadTdT1 );
	       } else {
		    InputLoad_Apply_LoadVectorForm_MPI( &LoadVectorForm, AccAll[istep], &Acc );
		    InputLoad_RelValues_MPI( &M, &Acc, &LoadTdT1 );
	       }
	  }

	  /* Join the non-coupling part. DispTdT_m = Keinv_m*fc + DispTdT0_m. Although DispTdT0_m is what has
	   * been received from the other computer, it has the name of DispTdT_m to avoid further operations
	   * if using the NETLIB libraries. */
	  if( CNodes.Order >= 1 ){
	       Substructure_JoinNonCouplingPart_MPI( &DispTdT0_m, &Keinv_m, &fcprevsub, &CNodes,  &DispTdT );
	  } else {
	       hysl_pcopy( &DispTdT0.GlobalSize.Row, DispTdT0.Array, &ione, &ione, DispTdT0.Desc, &ione,
			DispTdT.Array, &ione, &ione, DispTdT.Desc, &ione );
	  }

	  /* Compute acceleration ai1 = a0*(ui1 -ui) - a2*vi -a3*ai */
	  EffK_ComputeAcceleration_MPI( &DispTdT, &DispT, &VelT, &AccT, InitCnt.a0, InitCnt.a2, InitCnt.a3,
					&AccTdT );

	  /* Compute Velocity. vi = vi + a6*ai +a7*ai */
	  EffK_ComputeVelocity_MPI( &VelT, &AccT, &AccTdT, InitCnt.a6, InitCnt.a7, &VelTdT );

	  /* Error Compensation. fu = LoatTdT + fc -(Mass*AccTdT + Damp*VelTdT + Stiff*DispTdT) */
	  if( InitCnt.PID.P != 0.0 || InitCnt.PID.I != 0.0 || InitCnt.PID.D != 0.0 ){
	       ErrorForce_PID_MPI( &M, &C, &K, &AccTdT, &VelTdT, &DispTdT, &fc, &LoadTdT, &InitCnt.PID, &fu );
	  }

#if _HDF5_
	  /* Save the result in a HDF5 file format */
	  HDF5_Store_TimeHistoryData_MPI( hdf5_file, &AccTdT, &VelTdT, &DispTdT, &fc, &fu, (int) istep, nprow, myrow, &InitCnt );
#endif

	  /* Backup vectors */
	  hysl_pcopy( &LoadTdT.GlobalSize.Row, LoadTdT1.Array, &ione, &ione, LoadTdT1.Desc, &ione, LoadTdT.Array,
		   &ione, &ione, LoadTdT.Desc, &ione ); /* ui = ui1 */
	  /* Backup vectors */
	  hysl_pcopy( &DispTdT.GlobalSize.Row, DispTdT.Array, &ione, &ione, DispTdT.Desc, &ione, DispT.Array,
		   &ione, &ione, DispT.Desc, &ione ); /* ui = ui1 */
	  hysl_pcopy( &VelTdT.GlobalSize.Row, VelTdT.Array, &ione, &ione, VelTdT.Desc, &ione, VelT.Array, &ione,
		   &ione, VelT.Desc, &ione ); /* vi = vi1 */
	  hysl_pcopy( &AccTdT.GlobalSize.Row, AccTdT.Array, &ione, &ione, AccTdT.Desc, &ione, AccT.Array, &ione,
		   &ione, AccT.Desc, &ione ); /* ai = ai1 */
	  istep = istep + 1;
     }

     Time.End = MPI_Wtime();
     /* Elapsed time in miliseconds */
     Time.Elapsed_time = (Time.End - Time.Start)*1000.0;

#if _HDF5_
     HDF5_Store_Time_MPI( hdf5_file, &Time );
     HDF5_CloseFile( &hdf5_file );
#endif

     if( rank == 0 ){
	  Print_Header( SUCCESS );
	  printf( "The time integration process has finished in %lf ms.\n",
		  Time.Elapsed_time );
     }

     /* Free initiation values */
     Algorithm_Destroy( &InitCnt );

     /* Free the memory */
     if( InitCnt.Use_Absolute_Values ){
	  free( VelAll );
	  free( DispAll );
     } else {
	  free( AccAll );
     }

     /* Destroy the data structures */
     PMatrixVector_Destroy( &M );
     PMatrixVector_Destroy( &K );
     PMatrixVector_Destroy( &C );

     PMatrixVector_Destroy( &Keinv );

     if( CNodes.Order >= 1 ){
	  if( rank == 0 ){
	       MatrixVector_Destroy( &Keinv_c );
	  }
	  PMatrixVector_Destroy( &Keinv_m );
     }

     PMatrixVector_Destroy( &tempvec );
 
     PMatrixVector_Destroy( &DispT );
     PMatrixVector_Destroy( &DispTdT0 );
     if( CNodes.Order >= 1 ){
	  if( rank == 0 ){
	       MatrixVector_Destroy( &DispTdT0_c );
	  }
	  PMatrixVector_Destroy( &DispTdT0_m );
     }
     PMatrixVector_Destroy( &DispTdT );

     PMatrixVector_Destroy( &VelT );
     PMatrixVector_Destroy( &VelTdT );

     PMatrixVector_Destroy( &AccT );
     PMatrixVector_Destroy( &AccTdT );

     PMatrixVector_Destroy( &LoadVectorForm );
     PMatrixVector_Destroy( &LoadTdT );
     PMatrixVector_Destroy( &LoadTdT1 );

     PMatrixVector_Destroy( &EffT );

     PMatrixVector_Destroy( &fc );
     if( CNodes.Order >= 1 ){
	  PMatrixVector_Destroy( &fcprevsub );
     }
     PMatrixVector_Destroy( &fu );

     PMatrixVector_Destroy( &ErrCompT );
     PMatrixVector_Destroy( &ErrCompTdT );

     PMatrixVector_Destroy( &Disp );
     PMatrixVector_Destroy( &Vel );
     PMatrixVector_Destroy( &Acc );

     /* Free the coupling nodes memory */
     if( CNodes.Order >= 1 ){
	  Substructure_DeleteCouplingNodes( &CNodes );
     }

     /* BLACS: Exit Grid and release bhandle */
     Cblacs_gridexit( icontxt );
     Cfree_blacs_system_handle( bhandle );

     MPI_Finalize();

     return 0;
}
