#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <getopt.h>  /* For getopt_long() */

#include "Print_Messages.h"

#include "Algorithm_Aux.h"

#include "MatrixVector.h"
#include "Common_Formulation.h"
#include "EffK_Formulation.h"
#include "Error_Compensation.h"
#include "Initiation.h"
#include "Substructure.h"
#include "Substructure_Auxiliary.h"  /* For Substructure_VectorXm(), Substructure_VectorXc(), ... */

#include "HDF5_Operations.h"

#if _SPARSE_
#include "MatrixVector_Sp.h"
#include "Common_Formulation_Sp.h"
#include "EffK_Formulation_Sp.h"
#include "Error_Compensation_Sp.h"
#include "Initiation_Sp.h"
#include <mkl_blas.h>
#else
#include "Netlib.h"
#endif

void Print_Help( const char *Program_Name );

int main( int argc, char **argv ){

     unsigned int istep;
     AlgConst_t InitCnt;
     const char *FileConf;
     
     int hdf5_file;
     
     /* NETLIB Variables */
     int incx, incy;
     Scalars_t Constants;
     
     double *AccAll, *VelAll, *DispAll;

     MatrixVector_t M, C, K;               /* Mass, Damping and Stiffness matrices */
     MatrixVector_t Keinv;
     MatrixVector_t Keinv_c, Keinv_m;

     MatrixVector_t EffT;

     MatrixVector_t DispT, DispTdT0, DispTdT;
     MatrixVector_t DispTdT0_c, DispTdT0_m;
     MatrixVector_t tempvec;

     MatrixVector_t VelT, VelTdT;
     MatrixVector_t AccT, AccTdT;

     MatrixVector_t LoadVectorForm, LoadTdT, LoadTdT1;

     MatrixVector_t fc, fcprevsub;
     MatrixVector_t fu;

     MatrixVector_t ErrCompT, ErrCompTdT;

     MatrixVector_t Disp, Vel, Acc;

#if _SPARSE_
     /* Sparse matrices */
     MatrixVector_Sp_t Sp_M, Sp_C, Sp_K;     /* Sparse representation of the M, C and K matrices */
#endif

     CouplingNode_t CNodes;
     HDF5time_t     Time;
     /* Options */
     int Selected_Option;
     struct option long_options[] = {
	  {"help", no_argument, 0, 'h'},
	  {"config-file", required_argument, 0, 'c'},
	  {0, 0, 0, 0}
     };

     /* Print Information */
     printf( "\n\n" );
     printf( "************************************************************\n" );
     printf( "*                                                          *\n" );
     printf( "*  This is Dorka's substructure algorithm as programed by  *\n" );
     printf( "* Ferran Obón Santacana. Version alpha 0.9 'Heaven's Door' *\n" );
     printf( "*                                                          *\n" );
     printf( "************************************************************\n\n" );
     /* Set de default value for the configuration file */
     FileConf = "ConfFile.conf";

     /* This is only used if there are no arguments */
     if( argc == 1 ){
	  Print_Header( INFO );
	  printf( "Assuming the configuration file to be: ConfFile.conf.\n" );

     }

     while ((Selected_Option = getopt_long( argc, argv, "c:h", long_options, NULL )) != -1 ){
	  switch( Selected_Option ){
	  case 'c':
	       FileConf = optarg;
	       break;
	  case 'h':
	       Print_Help( argv[0] );
	       return EXIT_FAILURE;
	       break;
	  case '?':
	       /* Long options already prints an error message telling that there is an unrecognised option */
	       Print_Help( argv[0] );
	       return EXIT_FAILURE;
	  case ':':
	       /* Long options already prints an error message telling that the option requieres an argument */
	       Print_Help( argv[0] );
	       return EXIT_FAILURE;
	  }
     }

     /* Constants definitions. */
     Algorithm_Init( FileConf, &InitCnt );

     /* Read the coupling nodes from a file */
     Substructure_ReadCouplingNodes( &CNodes, InitCnt.NStep, InitCnt.NSubstep, InitCnt.OrderSub, InitCnt.DeltaT_Sub, InitCnt.FileCNodes );

     /* Allocate memory for saving the acceleration, displacement and velocity (input files) that will
      * be used during the test */
     if( InitCnt.Use_Absolute_Values ){
	  AccAll = NULL;
	  VelAll = (double *) calloc( (size_t) InitCnt.NStep, sizeof(double) );
	  DispAll = (double *) calloc( (size_t) InitCnt.NStep, sizeof(double) );
     } else {
	  AccAll = (double *) calloc( (size_t) InitCnt.NStep, sizeof(double) );
	  VelAll = NULL;
	  DispAll = NULL;
     }

     /* Initialise the matrices and vectors that will be used in the Time Integration process */
     if( ((!InitCnt.Use_Sparse && !InitCnt.Read_Sparse) || (InitCnt.Use_Sparse && !InitCnt.Read_Sparse) ||
	  (!InitCnt.Use_Sparse && InitCnt.Read_Sparse)) && !InitCnt.Use_Packed ){
	  MatrixVector_Create( InitCnt.Order, InitCnt.Order, &M );
	  MatrixVector_Create( InitCnt.Order, InitCnt.Order, &K );
     } else if( ((!InitCnt.Use_Sparse && !InitCnt.Read_Sparse) || (InitCnt.Use_Sparse && !InitCnt.Read_Sparse) ||
	  (!InitCnt.Use_Sparse && InitCnt.Read_Sparse)) && InitCnt.Use_Packed ){
	  MatrixVector_Create_PS( InitCnt.Order, InitCnt.Order, &M );
	  MatrixVector_Create_PS( InitCnt.Order, InitCnt.Order, &K );
     } else if ( InitCnt.Use_Sparse && InitCnt.Read_Sparse ){
	  MatrixVector_SetRowsCols_Sp( InitCnt.Order, InitCnt.Order, &Sp_M );
	  MatrixVector_SetRowsCols_Sp( InitCnt.Order, InitCnt.Order, &Sp_K );
     } else assert(0);

     if( !InitCnt.Use_Packed ){
	  MatrixVector_Create( InitCnt.Order, InitCnt.Order, &Keinv );
     } else {
	  MatrixVector_Create_PS( InitCnt.Order, InitCnt.Order, &Keinv );
     }

     if( CNodes.Order >= 1 ){
	  MatrixVector_Create( CNodes.Order, CNodes.Order, &Keinv_c );
	  MatrixVector_Create( InitCnt.Order - CNodes.Order, CNodes.Order, &Keinv_m );
     }

     MatrixVector_Create( InitCnt.Order, 1, &tempvec );

     MatrixVector_Create( InitCnt.Order, 1, &DispT );
     MatrixVector_Create( InitCnt.Order, 1, &DispTdT0 );
     MatrixVector_Create( InitCnt.Order, 1, &DispTdT );

     if( CNodes.Order >= 1 ){
	  MatrixVector_Create( CNodes.Order, 1, &DispTdT0_c );
	  MatrixVector_Create( InitCnt.Order-CNodes.Order, 1, &DispTdT0_m );
     }

     MatrixVector_Create( InitCnt.Order, 1, &VelT );
     MatrixVector_Create( InitCnt.Order, 1, &VelTdT );

     MatrixVector_Create( InitCnt.Order, 1, &AccT );
     MatrixVector_Create( InitCnt.Order, 1, &AccTdT );

     MatrixVector_Create( InitCnt.Order, 1, &EffT );

     MatrixVector_Create( InitCnt.Order, 1, &LoadVectorForm );
     MatrixVector_Create( InitCnt.Order, 1, &LoadTdT );
     MatrixVector_Create( InitCnt.Order, 1, &LoadTdT1 );

     MatrixVector_Create( InitCnt.Order, 1, &fc );
     if( CNodes.Order >= 1 ){
	  MatrixVector_Create( CNodes.Order, 1, &fcprevsub );
     }
     MatrixVector_Create( InitCnt.Order, 1, &fu );

     MatrixVector_Create( InitCnt.Order, 1, &ErrCompT );
     MatrixVector_Create( InitCnt.Order, 1, &ErrCompTdT );

     MatrixVector_Create( InitCnt.Order, 1, &Disp );
     MatrixVector_Create( InitCnt.Order, 1, &Vel );
     MatrixVector_Create( InitCnt.Order, 1, &Acc );

     /* Read the matrices from a file */
     if( !InitCnt.Read_Sparse && !InitCnt.Use_Sparse && !InitCnt.Use_Packed ){
	  MatrixVector_FromFile( InitCnt.FileM, &M );
	  MatrixVector_FromFile( InitCnt.FileK, &K );
     } else if ( !InitCnt.Read_Sparse && !InitCnt.Use_Sparse && InitCnt.Use_Packed ){
	  MatrixVector_FromFile_GE2PS( InitCnt.FileM, &M );
	  MatrixVector_FromFile_GE2PS( InitCnt.FileK, &K );
     } else if ( !InitCnt.Read_Sparse && InitCnt.Use_Sparse ){
	  MatrixVector_FromFile( InitCnt.FileM, &M );
	  MatrixVector_Dense2CSR( &M, 0, &Sp_M );        /* Transform into CSR format */
	  MatrixVector_Destroy( &M );                    /* Destroy the dense matrix */

	  MatrixVector_FromFile( InitCnt.FileK, &K );
	  MatrixVector_Dense2CSR( &K, 0, &Sp_K );        /* Transform into CSR format */
	  MatrixVector_Destroy( &K );                    /* Destroy the dense matrix */
     } else if ( InitCnt.Read_Sparse && !InitCnt.Use_Sparse && !InitCnt.Use_Packed ){
	  MatrixVector_FromFile_MM( InitCnt.FileM, &M );
	  MatrixVector_FromFile_MM( InitCnt.FileK, &K );
     } else if ( InitCnt.Read_Sparse && !InitCnt.Use_Sparse && InitCnt.Use_Packed ){
	  MatrixVector_FromFile_MM_PS( InitCnt.FileM, &M );
	  MatrixVector_FromFile_MM_PS( InitCnt.FileK, &K );
     } else if ( InitCnt.Read_Sparse && InitCnt.Use_Sparse ){
#if _SPARSE_
	  MatrixVector_FromFile_MM_Sp( InitCnt.FileM, &Sp_M );
	  MatrixVector_FromFile_MM_Sp( InitCnt.FileK, &Sp_K );
#endif
     } else assert(0);

     if ( !InitCnt.Read_LVector ){
	  InputLoad_Generate_LoadVectorForm( InitCnt.ExcitedDOF, &LoadVectorForm );
     }  else {
	  MatrixVector_FromFile_MM( InitCnt.FileLV, &LoadVectorForm );
     }

     /* Calculate damping matrix using Rayleigh. C = alpha*M + beta*K */
     if ( InitCnt.Use_Sparse ) {
#if _SPARSE_
	  MatrixVector_Create_Sp( InitCnt.Order, InitCnt.Order, Sp_K.Num_Nonzero, &Sp_C );
	  Rayleigh_Damping_Sp( &Sp_M, &Sp_K, &Sp_C, &InitCnt.Rayleigh );
#endif
     } else if ( !InitCnt.Use_Packed && !InitCnt.Use_Sparse ){
	  MatrixVector_Create( InitCnt.Order, InitCnt.Order, &C );
	  Rayleigh_Damping( &M, &K, &C, &InitCnt.Rayleigh );
     } else if ( InitCnt.Use_Packed && !InitCnt.Use_Sparse ){
	  MatrixVector_Create_PS( InitCnt.Order, InitCnt.Order, &C );
	  Rayleigh_Damping_PS( &M, &K, &C, &InitCnt.Rayleigh );
     } else assert(0);

     /* Calculate Matrix Keinv = 1.0*[K + a0*M + a1*C]^(-1) */
     Constants.Alpha = InitCnt.a0;    /* Mass matrix coefficient */
     Constants.Beta = InitCnt.a1;     /* Damping matrix coefficent */
     Constants.Gamma = 1.0;           /* Stiffness matrix coefficient */
     Constants.Lambda = 1.0;          /* Matrix inversion coefficient */

     if( !InitCnt.Use_Sparse && !InitCnt.Use_Packed ){
	  IGainMatrix( &Keinv, &M, &C, &K, Constants );
     } else if( !InitCnt.Use_Sparse && InitCnt.Use_Packed ){
	  IGainMatrix_PS( &Keinv, &M, &C, &K, Constants );
     } else if ( InitCnt.Use_Sparse && !InitCnt.Use_Packed ) {
#if _SPARSE_
	  IGainMatrix_Sp( &Keinv, &Sp_M, &Sp_C, &Sp_K, Constants );
#endif
     } else if ( InitCnt.Use_Sparse && InitCnt.Use_Packed ) {
#if _SPARSE_
	  IGainMatrix_Sp_PS( &Keinv, &Sp_M, &Sp_C, &Sp_K, Constants );
#endif
     } else assert(0);

     if( CNodes.Order >= 1 ){
	  if( !InitCnt.Use_Packed ){
	       Substructure_MatrixXc( &Keinv, &CNodes, &Keinv_c );
	       Substructure_MatrixXcm( &Keinv, &CNodes, &Keinv_m );
	  } else {
	       Substructure_MatrixXc_PS( &Keinv, &CNodes, &Keinv_c );
	       Substructure_MatrixXcm_PS( &Keinv, &CNodes, &Keinv_m );
	  }
  
	  /* Send the coupling part of the effective matrix if we are performing a distributed test */
//	  Send_Effective_Matrix( Keinv_c.Array, (unsigned int) CNodes.Order, &Socket, InitCnt.Remote );
     }

     /* Read the earthquake data from a file */
     if( InitCnt.Use_Absolute_Values ){
	  Algorithm_ReadDataEarthquake_AbsValues( InitCnt.NStep, InitCnt.FileData, VelAll, DispAll );
     } else {
	  Algorithm_ReadDataEarthquake_RelValues( InitCnt.NStep, InitCnt.FileData, AccAll );
     }

     /* Open Output file. If the file cannot be opened, the program will exit, since the results cannot be stored. */

     hdf5_file = HDF5_CreateFile( InitCnt.FileOutput );
     HDF5_CreateGroup_Parameters( hdf5_file, &InitCnt, &CNodes );
     HDF5_CreateGroup_TimeIntegration( hdf5_file, &InitCnt );

     Time.Date_start = time( NULL );
     Time.Date_time = strdup( ctime( &Time.Date_start) );
     gettimeofday( &Time.start, NULL );        

     /* Calculate the input load */
     istep = 1;
     if( InitCnt.Use_Absolute_Values ){
	  /* Copy the diagonal elements of M */
	  InputLoad_Apply_LoadVectorForm( &LoadVectorForm, DispAll[istep - 1], &Disp );

	  InputLoad_Apply_LoadVectorForm( &LoadVectorForm, VelAll[istep - 1], &Vel );

	  if( !InitCnt.Use_Sparse && !InitCnt.Use_Packed ){
	       InputLoad_AbsValues( &K, &C, &Disp, &Vel, &LoadTdT );
	  } else if ( !InitCnt.Use_Sparse && InitCnt.Use_Packed ){
	       InputLoad_AbsValues_PS( &K, &C, &Disp, &Vel, &LoadTdT );
	  } else if ( InitCnt.Use_Sparse ){
#if _SPARSE_
	       InputLoad_AbsValues_Sp( &Sp_K, &Sp_C, &Disp, &Vel, &LoadTdT );
#endif
	  } else assert(0);
     } else {
	  InputLoad_Apply_LoadVectorForm( &LoadVectorForm, AccAll[istep - 1], &Acc );
	  if( !InitCnt.Use_Sparse && !InitCnt.Use_Packed ){
	       InputLoad_RelValues( &M, &Acc, &LoadTdT );
	  } else if ( !InitCnt.Use_Sparse && InitCnt.Use_Packed ){
	       InputLoad_RelValues_PS( &M, &Acc, &LoadTdT );
	  } else if ( InitCnt.Use_Sparse ){
#if _SPARSE_
	       InputLoad_RelValues_Sp( &Sp_M, &Acc, &LoadTdT );
#endif
	  } else assert(0);
     }

     incx = 1; incy = 1;
     Print_Header( INFO );
     printf( "Starting stepping process.\n" );
     while ( istep <= InitCnt.NStep ){

	  /* Calculate the effective force vector
	     Fe = M*(a0*u + a2*v + a3*a) + C*(a1*u + a4*v + a5*a) */
	  if( !InitCnt.Use_Sparse && !InitCnt.Use_Packed ){
	       EffK_EffectiveForce( &M, &C, &DispT, &VelT, &AccT, &tempvec, InitCnt.a0, InitCnt.a1, InitCnt.a2,
				    InitCnt.a3, InitCnt.a4, InitCnt.a5, &EffT );
	  } else if( !InitCnt.Use_Sparse && InitCnt.Use_Packed ){
	       EffK_EffectiveForce_PS( &M, &C, &DispT, &VelT, &AccT, &tempvec, InitCnt.a0, InitCnt.a1, InitCnt.a2,
				       InitCnt.a3, InitCnt.a4, InitCnt.a5, &EffT );
	  } else if( InitCnt.Use_Sparse ) {
#if _SPARSE_
	       EffK_EffectiveForce_Sp( &Sp_M, &Sp_C, &DispT, &VelT, &AccT, &tempvec, InitCnt.a0, InitCnt.a1,
				       InitCnt.a2, InitCnt.a3, InitCnt.a4, InitCnt.a5, &EffT );
#endif
	  }

	  /* Compute the new Displacement u0 */
	  if( !InitCnt.Use_Packed ){
	       Compute_NewState( &Keinv, &EffT, &LoadTdT, &fu, &tempvec, &DispTdT0 );
	  } else {
	       Compute_NewState_PS( &Keinv, &EffT, &LoadTdT, &fu, &tempvec, &DispTdT0 );
	  }

	  /* Split DispTdT into coupling and non-coupling part */
	  if( CNodes.Order >= 1 ){
	       Substructure_VectorXm( &DispTdT0, &CNodes, &DispTdT0_m );
	       Substructure_VectorXc( &DispTdT0, &CNodes, &DispTdT0_c );
	  }

	  /* Perform substepping */
	  if( CNodes.Order >= 1 ){
	       Substructure_Substepping( Keinv_c.Array, DispTdT0_c.Array, InitCnt.Delta_t*(double) istep, InitCnt.NSubstep,
					 InitCnt.DeltaT_Sub, &CNodes, DispTdT.Array, fcprevsub.Array, fc.Array );
	  }

	  if ( istep < InitCnt.NStep ){
	       /* Calculate the input load for the next step during the
		  sub-stepping process. */
	       if( InitCnt.Use_Absolute_Values ){
		    /* Copy the diagonal elements of M */
		    InputLoad_Apply_LoadVectorForm( &LoadVectorForm, DispAll[istep], &Disp );		    
		    InputLoad_Apply_LoadVectorForm( &LoadVectorForm, VelAll[istep], &Vel );

		    if( !InitCnt.Use_Sparse && !InitCnt.Use_Packed ){
			 InputLoad_AbsValues( &K, &C, &Disp, &Vel, &LoadTdT1 );
		    } else if ( !InitCnt.Use_Sparse && InitCnt.Use_Packed ){
			 InputLoad_AbsValues_PS( &K, &C, &Disp, &Vel, &LoadTdT1 );
		    } else if ( InitCnt.Use_Sparse ){
#if _SPARSE_
			 InputLoad_AbsValues_Sp( &Sp_K, &Sp_C, &Disp, &Vel, &LoadTdT1 );
#endif
		    } else assert(0);
	       } else {
		    InputLoad_Apply_LoadVectorForm( &LoadVectorForm, AccAll[istep], &Acc );
		    if( !InitCnt.Use_Sparse && !InitCnt.Use_Packed ){
			 InputLoad_RelValues( &M, &Acc, &LoadTdT1 );
		    } else if ( !InitCnt.Use_Sparse && InitCnt.Use_Packed ){
			 InputLoad_RelValues_PS( &M, &Acc, &LoadTdT1 );
		    } else if ( InitCnt.Use_Sparse ){
#if _SPARSE_
			 InputLoad_RelValues_Sp( &Sp_M, &Acc, &LoadTdT1 );
#endif
		    } else assert(0);
	       }
	  }

	  /* Join the non-coupling part. DispTdT_m = Keinv_m*fc + DispTdT0_m. Although DispTdT0_m is what has been received from the other computer,
	     it has the name of DispTdT_m to avoid further operations if using the NETLIB libraries. */
	  if( CNodes.Order >= 1 ){
	       Substructure_JoinNonCouplingPart( &DispTdT0_m, &Keinv_m, &fcprevsub, &CNodes,  &DispTdT );
	  } else {
	       dcopy( &DispTdT0.Rows, DispTdT0.Array, &incx, DispTdT.Array, &incy ); /* ui = ui1 */
	  }

	  /* Compute acceleration ai1 = a0*(ui1 -ui) - a2*vi -a3*ai */
	  EffK_ComputeAcceleration( &DispTdT, &DispT, &VelT, &AccT, InitCnt.a0, InitCnt.a2, InitCnt.a3, &AccTdT );

	  /* Compute Velocity. vi = vi + a6*ai +a7*ai */
	  EffK_ComputeVelocity( &VelT, &AccT, &AccTdT, InitCnt.a6, InitCnt.a7, &VelTdT );

	  /* Error Compensation. fu = LoatTdT + fc -(Mass*AccTdT + Damp*VelTdT + Stiff*DispTdT) */
	  if( InitCnt.PID.P != 0.0 || InitCnt.PID.I != 0.0 || InitCnt.PID.D != 0.0 ){
	       if( !InitCnt.Use_Sparse && !InitCnt.Use_Packed ){
		    ErrorForce_PID( &M, &C, &K, &AccTdT, &VelTdT, &DispTdT, &fc, &LoadTdT, &InitCnt.PID, &fu );
	       } else if( !InitCnt.Use_Sparse && InitCnt.Use_Packed ) {
		    ErrorForce_PID_PS( &M, &C, &K, &AccTdT, &VelTdT, &DispTdT, &fc, &LoadTdT, &InitCnt.PID, &fu );
	       } else if( InitCnt.Use_Sparse ){
#if _SPARSE_
		    ErrorForce_PID_Sp( &Sp_M, &Sp_C, &Sp_K, &AccTdT, &VelTdT, &DispTdT, &fc, &LoadTdT, &InitCnt.PID, &fu );
#endif
	       } else assert(0);
	  }

	  /* Save the result in a HDF5 file format */
	  HDF5_Store_TimeHistoryData( hdf5_file, &AccTdT, &VelTdT, &DispTdT, &LoadTdT, &fc, &fu, (int) istep, &InitCnt );

	  /* Backup vectors */
	  dcopy( &LoadTdT1.Rows, LoadTdT1.Array, &incx, LoadTdT.Array, &incy ); /* li = li1 */
	  dcopy( &DispTdT.Rows, DispTdT.Array, &incx, DispT.Array, &incy ); /* ui = ui1 */
	  dcopy( &VelTdT.Rows, VelTdT.Array, &incx, VelT.Array, &incy ); /* vi = vi1 */
	  dcopy( &AccTdT.Rows, AccTdT.Array, &incx, AccT.Array, &incy ); /* ai = ai1 */
	  istep = istep + 1;
     }

     gettimeofday( &Time.end, NULL );
     Time.Elapsed_time = (double) (Time.end.tv_sec - Time.start.tv_sec)*1000.0;
     Time.Elapsed_time += (double) (Time.end.tv_usec - Time.start.tv_usec)/1000.0;
     HDF5_StoreTime( hdf5_file, &Time );
     Print_Header( SUCCESS );
     printf( "The time integration process has finished in %lf ms.\n", Time.Elapsed_time );

     /* Close the Connection */
//     if( InitCnt.Remote.Type != NO_PROTOCOL ){
//	  Close_Connection( &Socket, hdf5_file, InitCnt.Remote.Type, (unsigned int) CNodes.Order, InitCnt.NStep, 4 );
//     }

     HDF5_CloseFile( hdf5_file );

     /* Free initiation values */
     Algorithm_Destroy( &InitCnt );

     /* Free the memory */
     if( InitCnt.Use_Absolute_Values ){
	  free( VelAll );
	  free( DispAll );
     } else {
	  free( AccAll );
     }

     /* Free Time string */
     free( Time.Date_time );

     /* Destroy the data structures */
     if( !InitCnt.Use_Sparse && !InitCnt.Use_Packed ){
	  MatrixVector_Destroy( &M );
	  MatrixVector_Destroy( &K );
	  MatrixVector_Destroy( &C );
     } else if ( !InitCnt.Use_Sparse && InitCnt.Use_Packed ){
	  MatrixVector_Destroy_PS( &M );
	  MatrixVector_Destroy_PS( &K );
	  MatrixVector_Destroy_PS( &C );
     } else if ( InitCnt.Use_Sparse ) {
#if _SPARSE_
	  MatrixVector_Destroy_Sp( &Sp_M );
	  MatrixVector_Destroy_Sp( &Sp_C );
	  MatrixVector_Destroy_Sp( &Sp_K );

#endif
     }

     if ( InitCnt.Use_Packed ){	  
	  MatrixVector_Destroy( &Keinv );
     } else if ( InitCnt.Use_Packed ){	 
	  MatrixVector_Destroy_PS( &Keinv );
     }
 
     if( CNodes.Order >= 1 ){
	  MatrixVector_Destroy( &Keinv_c );
	  MatrixVector_Destroy( &Keinv_m );
     }

     MatrixVector_Destroy( &tempvec );
 
     MatrixVector_Destroy( &DispT );
     MatrixVector_Destroy( &DispTdT0 );
     if( CNodes.Order >= 1 ){
	  MatrixVector_Destroy( &DispTdT0_c );
	  MatrixVector_Destroy( &DispTdT0_m );
     }
     MatrixVector_Destroy( &DispTdT );

     MatrixVector_Destroy( &VelT );
     MatrixVector_Destroy( &VelTdT );

     MatrixVector_Destroy( &AccT );
     MatrixVector_Destroy( &AccTdT );

     MatrixVector_Destroy( &LoadVectorForm );
     MatrixVector_Destroy( &LoadTdT );
     MatrixVector_Destroy( &LoadTdT1 );

     MatrixVector_Destroy( &EffT );

     MatrixVector_Destroy( &fc );
     if( CNodes.Order >= 1 ){
	  MatrixVector_Destroy( &fcprevsub );
     }
     MatrixVector_Destroy( &fu );

     MatrixVector_Destroy( &ErrCompT );
     MatrixVector_Destroy( &ErrCompTdT );

     MatrixVector_Destroy( &Disp );
     MatrixVector_Destroy( &Vel );
     MatrixVector_Destroy( &Acc );

     /* Free the coupling nodes memory */
     Substructure_DeleteCouplingNodes( &CNodes );

     return 0;
}

void Print_Help( const char *Program_Name )
{

     fprintf( stderr, "Usage: %s [-h] -c <ConfFile", Program_Name );
     fprintf( stderr,
	      "  -h  --help           This help text.\n"
	      "  -c  --config-file    The name of the configuration file. Default value: ConfFile.conf\n" );
}
