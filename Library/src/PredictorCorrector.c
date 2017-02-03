#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <getopt.h>
#include <math.h>
#include <time.h>

#include "Print_Messages.h"
#include "Algorithm_Aux.h"

#include "MatrixVector.h"
#include "MatrixVector_PS.h"
#include "Input_Load.h"
#include "GainMatrix.h"
#include "PCMethods.h"

#include "Substructure.h"

#include "Definitions.h"

#if _HDF5_
#include "HDF5_Operations.h"
#endif

#include "ADwin_Routines.h"

#if _SPARSE_
#include "MatrixVector_Sp.h"
#include <mkl_blas.h>
#else
#include "Netlib.h"
#endif

int main (int argc, char **argv){
     
     unsigned int istep;
     int incx, incy;
     
     AlgConst_t InitCnt;
     const char *FileConf;
     
#if _HDF5_
     hid_t hdf5_file;
#endif

     Scalars_t Constants;

     HYSL_FLOAT *AccAll1, *VelAll1, *DispAll1;
     HYSL_FLOAT *AccAll2, *VelAll2, *DispAll2;
     HYSL_FLOAT *AccAll3, *VelAll3, *DispAll3;
     
     MatrixVector_t M, C, K;               /* Mass, Damping and Stiffness matrices */
     MatrixVector_t Meinv;

     MatrixVector_t DispTdT, VelTdT, AccTdT;
     MatrixVector_t DispTdT_Pred, VelTdT_Pred;

     MatrixVector_t DispT, VelT, AccT;
     MatrixVector_t VelT_Pred;

     MatrixVector_t RForceTdT, RForceT;
     
     MatrixVector_t LoadVectorForm1, LoadVectorForm2, LoadVectorForm3, LoadT, LoadTdT;

     MatrixVector_t Disp, Vel, Acc;
     
     /* Options */
     int Selected_Option;
     struct option long_options[] = {
	  {"help", no_argument, 0, 'h'},
	  {"config-file", required_argument, 0, 'c'},
	  {0, 0, 0, 0}
     };

#if _SPARSE_
     /* Sparse matrices */
     MatrixVector_Sp_t Sp_M, Sp_C, Sp_K;     /* Sparse representation of the M, C and K matrices */
#endif
     
     CouplingNode_t CNodes;
     SaveTime_t     Time;

     /* Print Information */
     printf( "\n\n" );
     printf( "************************************************************\n" );
     printf( "*                                                          *\n" );
     printf( "*  This is Dorka's substructure algorithm as programed by  *\n" );
     printf( "* Ferran Obón Santacana. Version alpha 0.5 'Heaven's Door' *\n" );
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
	       Algorithm_PrintHelp( argv[0] );
	       return EXIT_FAILURE;
	       break;
	  case '?':
	       /* Long options already prints an error message telling that there is an unrecognised option */
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
     Algorithm_Init( FileConf, &InitCnt );
     
     /* Read the coupling nodes from a file */
     Substructure_ReadCouplingNodes( &InitCnt, &CNodes );
     
     /* Allocate memory for saving the acceleration, displacement and velocity (input files) that will be used
      * during the test */
     AccAll1 = (HYSL_FLOAT *) calloc( (size_t) InitCnt.NStep, sizeof(HYSL_FLOAT) );
     VelAll1 = (HYSL_FLOAT *) calloc( (size_t) InitCnt.NStep, sizeof(HYSL_FLOAT) );
     DispAll1 = (HYSL_FLOAT *) calloc( (size_t) InitCnt.NStep, sizeof(HYSL_FLOAT) );
     AccAll2 = (HYSL_FLOAT *) calloc( (size_t) InitCnt.NStep, sizeof(HYSL_FLOAT) );
     VelAll2 = (HYSL_FLOAT *) calloc( (size_t) InitCnt.NStep, sizeof(HYSL_FLOAT) );
     DispAll2 = (HYSL_FLOAT *) calloc( (size_t) InitCnt.NStep, sizeof(HYSL_FLOAT) );
     AccAll3 = (HYSL_FLOAT *) calloc( (size_t) InitCnt.NStep, sizeof(HYSL_FLOAT) );
     VelAll3 = (HYSL_FLOAT *) calloc( (size_t) InitCnt.NStep, sizeof(HYSL_FLOAT) );
     DispAll3 = (HYSL_FLOAT *) calloc( (size_t) InitCnt.NStep, sizeof(HYSL_FLOAT) );

     /* Initialise the matrices and vectors that will be used in the Time Integration process */
     if( ((!InitCnt.Use_Sparse && !InitCnt.Read_Sparse) || (InitCnt.Use_Sparse && !InitCnt.Read_Sparse) ||
	  (!InitCnt.Use_Sparse && InitCnt.Read_Sparse)) && !InitCnt.Use_Packed ){
	  MatrixVector_Create( InitCnt.Order, InitCnt.Order, &M );
	  MatrixVector_Create( InitCnt.Order, InitCnt.Order, &K );
	  if( InitCnt.Read_CMatrix ){
	       MatrixVector_Create( InitCnt.Order, InitCnt.Order, &C );
	  }
     } else if( ((!InitCnt.Use_Sparse && !InitCnt.Read_Sparse) || (InitCnt.Use_Sparse && !InitCnt.Read_Sparse)
		 || (!InitCnt.Use_Sparse && InitCnt.Read_Sparse)) && InitCnt.Use_Packed ){
	  MatrixVector_Create_PS( InitCnt.Order, InitCnt.Order, &M );
	  MatrixVector_Create_PS( InitCnt.Order, InitCnt.Order, &K );
	  if( InitCnt.Read_CMatrix ){
	       MatrixVector_Create_PS( InitCnt.Order, InitCnt.Order, &C );
	  }
     } else if ( InitCnt.Use_Sparse && InitCnt.Read_Sparse ){
#if _SPARSE_
	  MatrixVector_SetRowsCols_Sp( InitCnt.Order, InitCnt.Order, &Sp_M );
	  MatrixVector_SetRowsCols_Sp( InitCnt.Order, InitCnt.Order, &Sp_K );
	  if( InitCnt.Read_CMatrix ){
	       MatrixVector_SetRowsCols_Sp( InitCnt.Order, InitCnt.Order, &Sp_C );
	  }
#endif
     } else assert(0);

     if ( !InitCnt.Use_Packed ){
	  MatrixVector_Create( InitCnt.Order, InitCnt.Order, &Meinv );
     } else {
	  MatrixVector_Create_PS( InitCnt.Order, InitCnt.Order, &Meinv );
     }


     MatrixVector_Create( InitCnt.Order, 1, &DispT );
     MatrixVector_Create( InitCnt.Order, 1, &DispTdT );
     MatrixVector_Create( InitCnt.Order, 1, &DispTdT_Pred );

     MatrixVector_Create( InitCnt.Order, 1, &VelT );
     MatrixVector_Create( InitCnt.Order, 1, &VelT_Pred );
     MatrixVector_Create( InitCnt.Order, 1, &VelTdT );
     MatrixVector_Create( InitCnt.Order, 1, &VelTdT_Pred );

     MatrixVector_Create( InitCnt.Order, 1, &AccT );
     MatrixVector_Create( InitCnt.Order, 1, &AccTdT );

     MatrixVector_Create( InitCnt.Order, 1, &RForceT );
     MatrixVector_Create( InitCnt.Order, 1, &RForceTdT );
     
     MatrixVector_Create( InitCnt.Order, 1, &LoadVectorForm1 );
     MatrixVector_Create( InitCnt.Order, 1, &LoadVectorForm2 );
     MatrixVector_Create( InitCnt.Order, 1, &LoadVectorForm3 );
     MatrixVector_Create( InitCnt.Order, 1, &LoadT );
     MatrixVector_Create( InitCnt.Order, 1, &LoadTdT );

     MatrixVector_Create( InitCnt.Order, 1, &Disp );
     MatrixVector_Create( InitCnt.Order, 1, &Vel );
     MatrixVector_Create( InitCnt.Order, 1, &Acc );
     
          /* Read the matrices from a file */
     if( !InitCnt.Read_Sparse && !InitCnt.Use_Sparse && !InitCnt.Use_Packed ){
	  MatrixVector_FromFile( InitCnt.FileM, &M );
	  MatrixVector_FromFile( InitCnt.FileK, &K );
	  if( InitCnt.Read_CMatrix ){
	       MatrixVector_FromFile( InitCnt.FileC, &C );
	  }
     } else if ( !InitCnt.Read_Sparse && !InitCnt.Use_Sparse && InitCnt.Use_Packed ){
	  MatrixVector_FromFile_GE2PS( InitCnt.FileM, &M );
	  MatrixVector_FromFile_GE2PS( InitCnt.FileK, &K );
	  if( InitCnt.Read_CMatrix ){
	       MatrixVector_FromFile_GE2PS( InitCnt.FileC, &C );
	  }
     } else if ( !InitCnt.Read_Sparse && InitCnt.Use_Sparse ){
#if _SPARSE_
	  MatrixVector_FromFile( InitCnt.FileM, &M );
	  MatrixVector_Dense2CSR( &M, 0, &Sp_M );        /* Transform into CSR format */
	  MatrixVector_Destroy( &M );                    /* Destroy the dense matrix */

	  MatrixVector_FromFile( InitCnt.FileK, &K );
	  MatrixVector_Dense2CSR( &K, 0, &Sp_K );        /* Transform into CSR format */
	  MatrixVector_Destroy( &K );                    /* Destroy the dense matrix */

	  if( InitCnt.Read_CMatrix ){
	       MatrixVector_FromFile( InitCnt.FileK, &C );
	       MatrixVector_Dense2CSR( &C, 0, &Sp_C );        /* Transform into CSR format */
	       MatrixVector_Destroy( &C );                    /* Destroy the dense matrix */	  
	  }     
#endif
     } else if ( InitCnt.Read_Sparse && !InitCnt.Use_Sparse && !InitCnt.Use_Packed ){
	  MatrixVector_FromFile_MM( InitCnt.FileM, &M );
	  MatrixVector_FromFile_MM( InitCnt.FileK, &K );
	  if( InitCnt.Read_CMatrix ){
	       MatrixVector_FromFile_MM( InitCnt.FileC, &C );
	  }
     } else if ( InitCnt.Read_Sparse && !InitCnt.Use_Sparse && InitCnt.Use_Packed ){
	  MatrixVector_FromFile_MM_PS( InitCnt.FileM, &M );
	  MatrixVector_FromFile_MM_PS( InitCnt.FileK, &K );
	  if( InitCnt.Read_CMatrix ){
	       MatrixVector_FromFile_MM_PS( InitCnt.FileC, &C );
	  }
     } else if ( InitCnt.Read_Sparse && InitCnt.Use_Sparse ){
#if _SPARSE_
	  MatrixVector_FromFile_MM_Sp( InitCnt.FileM, &Sp_M );
	  MatrixVector_FromFile_MM_Sp( InitCnt.FileK, &Sp_K );
	  if( InitCnt.Read_CMatrix ){
	       MatrixVector_FromFile_MM_Sp( InitCnt.FileC, &Sp_C );
	  }
#endif
     } else assert(0);

     if ( !InitCnt.Read_LVector ){
	  InputLoad_Generate_LoadVectorForm( InitCnt.ExcitedDOF, &LoadVectorForm1 );
     }  else {
	  MatrixVector_FromFile_MM( InitCnt.FileLV1, &LoadVectorForm1 );
	  MatrixVector_FromFile_MM( InitCnt.FileLV2, &LoadVectorForm2 );
	  MatrixVector_FromFile_MM( InitCnt.FileLV3, &LoadVectorForm3 );
     }

     /* Calculate damping matrix using Rayleigh. C = alpha*M + beta*K */
     if( !InitCnt.Read_CMatrix ){
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
     }

     /* Calculate Matrix Meinv = 1.0*[M + (1 + alpha_H)*a7*C + (1 + alpha_H)*a8*K]^(-1) */
     Constants.Alpha = 1.0;                          /* Mass matrix coefficient */
     Constants.Beta = InitCnt.a7*(1.0 + InitCnt.TIntConst.HilberAlpha);  /* Damping matrix coefficent */
     Constants.Gamma = (1.0 + InitCnt.TIntConst.HilberAlpha)*InitCnt.a8; /* Stiffness matrix coefficient */
     Constants.Lambda = 1.0;                         /* Matrix inversion coefficient */

     if( !InitCnt.Use_Sparse && !InitCnt.Use_Packed ){
	  IGainMatrix( &Meinv, &M, &C, &K, Constants );
     } else if( !InitCnt.Use_Sparse && InitCnt.Use_Packed ){
	  IGainMatrix_PS( &Meinv, &M, &C, &K, Constants );
     } else if ( InitCnt.Use_Sparse && !InitCnt.Use_Packed ) {
#if _SPARSE_
	  IGainMatrix_Sp( &Meinv, &Sp_M, &Sp_C, &Sp_K, Constants );
#endif
     } else if ( InitCnt.Use_Sparse && InitCnt.Use_Packed ) {
#if _SPARSE_
	  IGainMatrix_Sp_PS( &Meinv, &Sp_M, &Sp_C, &Sp_K, Constants );
#endif
     } else assert(0);

     /* Read the earthquake data from a file */
     Algorithm_ReadDataEarthquake( InitCnt.NStep, InitCnt.FileData1, InitCnt.Scale_Factor, AccAll1,
				   VelAll1, DispAll1 );
     Algorithm_ReadDataEarthquake( InitCnt.NStep, InitCnt.FileData2, InitCnt.Scale_Factor, AccAll2,
				   VelAll2, DispAll2 );
     Algorithm_ReadDataEarthquake( InitCnt.NStep, InitCnt.FileData3, InitCnt.Scale_Factor, AccAll3,
				   VelAll3, DispAll3 );

     /* Open Output file. If the file cannot be opened, the program will exit, since the results cannot be
      * stored. */
#if _HDF5_
     hdf5_file = HDF5_CreateFile( Concatenate_Strings( 2, InitCnt.FileOutput, ".h5" ) );
     HDF5_CreateGroup_Parameters( hdf5_file, &InitCnt, &CNodes, AccAll1, VelAll1, DispAll1, AccAll2, VelAll2, DispAll2, AccAll3, VelAll3, DispAll3 );
     HDF5_CreateGroup_TimeIntegration( hdf5_file, &InitCnt );
#endif

     istep = 1;

     Time.Date_start = time( NULL );
     Time.Date_time = strdup( ctime( &Time.Date_start) );
     gettimeofday( &Time.Start, NULL );

     Print_Header( INFO );
     printf( "Starting stepping process.\n" );

     /* Initialise BLAS variables incx and incy */
     incx = 1; incy = 1;
     
     while ( istep <= InitCnt.NStep ){


	  /* Calculate input load */
	  if( InitCnt.Use_Absolute_Values ){
	       InputLoad_Apply_LoadVectorForm( &LoadVectorForm1, &LoadVectorForm2, &LoadVectorForm3, DispAll1[istep - 1], DispAll2[istep - 1], DispAll3[istep - 1], &Disp );
	       InputLoad_Apply_LoadVectorForm( &LoadVectorForm1, &LoadVectorForm2, &LoadVectorForm3, VelAll1[istep - 1], VelAll2[istep - 1], VelAll3[istep - 1], &Vel );
	       
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
	       InputLoad_Apply_LoadVectorForm( &LoadVectorForm1, &LoadVectorForm2, &LoadVectorForm3, AccAll1[istep - 1], AccAll2[istep - 1], AccAll3[istep - 1], &Acc );
	       
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

	  if( istep == 1 ){
	       hysl_copy( &AccTdT.Rows, Acc.Array, &incx, AccT.Array, &incy ); /* ai = ai1 */
	  }
	  
	  /* Calculate predictor step */
	  PC_PredictorStep_Displacement ( &DispT, &VelT, &AccT, InitCnt.a9, InitCnt.a10, &DispTdT_Pred );
	  PC_PredictorStep_Velocity ( &VelT, &AccT, InitCnt.a6, &VelTdT_Pred );

	  printf("DispTdT_Pred %lE\n", DispTdT_Pred.Array[0]);
	  printf("VelTdT_Pred %lE\n", VelTdT_Pred.Array[0]);

	  if( !InitCnt.Use_Sparse && !InitCnt.Use_Packed ){
	       /* Calculate linear reaction forces */
	       PC_ReactionForces_Numerical ( &DispTdT_Pred, &K, &RForceTdT );

	       printf("RForceTdT %lE\n", RForceTdT.Array[0]);
	       
	       /* Calculate accelerations at n + 1 */
	       PC_Calculate_Acceleration ( &LoadTdT, &LoadT, &RForceTdT, &RForceT,
					   &VelTdT_Pred, &VelT_Pred, &AccT, &K, &C, &Meinv,
					   InitCnt.TIntConst.HilberAlpha, InitCnt.a7, InitCnt.a8, &AccTdT );
	  } else if ( !InitCnt.Use_Sparse && InitCnt.Use_Packed ){
	       /* Calculate linear reaction forces */
	       PC_ReactionForces_Numerical_PS ( &DispTdT_Pred, &K, &RForceTdT );
	       
	       /* Calculate accelerations at n + 1 */
	       PC_Calculate_Acceleration_PS ( &LoadTdT, &LoadT, &RForceTdT, &RForceT,
					      &VelTdT_Pred, &VelT_Pred, &AccT, &K, &C, &Meinv,
					      InitCnt.TIntConst.HilberAlpha, InitCnt.a7, InitCnt.a8, &AccTdT );
	  } else if ( InitCnt.Use_Sparse && !InitCnt.Use_Packed){
	       Print_Header( ERROR );
	       fprintf( stderr, "Sparse Not supported" );
	       exit( EXIT_FAILURE );
	  } else if ( InitCnt.Use_Sparse && InitCnt.Use_Packed){
	       Print_Header( ERROR );
	       fprintf( stderr, "Sparse Not supported" );
	       exit( EXIT_FAILURE );
	  } else assert(0);
	  	  
	  /* Calculate corrector step */
	  PC_CorrectorStep_Displacement ( &DispTdT_Pred, &AccTdT, InitCnt.a8, &DispTdT );
	  PC_CorrectorStep_Velocity ( &VelTdT_Pred, &AccTdT, InitCnt.a7, &VelTdT );


	  /* Save the result in a HDF5 file format */
#if _HDF5_
	  HDF5_Store_TimeHistoryData( hdf5_file, &AccTdT, &VelTdT, &DispTdT, &LoadTdT, &LoadTdT, (int) istep, &InitCnt );
#else
	  
#endif
	  printf("%d\t %lE\n", istep, DispTdT.Array[0]);
	  /* Backup vectors */
	  hysl_copy( &DispTdT.Rows, DispTdT.Array, &incx, DispT.Array, &incy ); /* ui = ui1 */
	  hysl_copy( &VelTdT.Rows, VelTdT.Array, &incx, VelT.Array, &incy ); /* vi = vi1 */
	  hysl_copy( &AccTdT.Rows, AccTdT.Array, &incx, AccT.Array, &incy ); /* ai = ai1 */
	  hysl_copy( &LoadTdT.Rows, LoadTdT.Array, &incx, LoadT.Array, &incy ); /* ai = ai1 */
	  hysl_copy( &VelTdT_Pred.Rows, VelTdT_Pred.Array, &incx, VelT_Pred.Array, &incy ); /* vi' = vi1' */
	  hysl_copy( &RForceTdT.Rows, RForceTdT.Array, &incx, RForceT.Array, &incy ); /* ri = ri1 */
	  
	  istep = istep + 1;
     }


     gettimeofday( &Time.End, NULL );
     Time.Elapsed_time = (double) (Time.End.tv_sec - Time.Start.tv_sec)*1000.0;
     Time.Elapsed_time += (double) (Time.End.tv_usec - Time.Start.tv_usec)/1000.0;
#if _HDF5_
     HDF5_Store_Time( hdf5_file, &Time );
#endif

     Print_Header( SUCCESS );
     printf( "The time integration process has finished in %lf ms.\n", Time.Elapsed_time );

#if _HDF5_
     HDF5_CloseFile( &hdf5_file );
#endif

     /* Free initiation values */
     Algorithm_Destroy( &InitCnt );

     /* Free the memory */
     free( AccAll1 );
     free( VelAll1 );
     free( DispAll1 );
     free( AccAll2 );
     free( VelAll2 );
     free( DispAll2 );
     free( AccAll3 );
     free( VelAll3 );
     free( DispAll3 );
     
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

     if ( !InitCnt.Use_Packed ){	  
	  MatrixVector_Destroy( &Meinv );
     } else if ( InitCnt.Use_Packed ){	 
	  MatrixVector_Destroy_PS( &Meinv );
     }

     MatrixVector_Destroy( &DispT );
     MatrixVector_Destroy( &DispTdT );
     MatrixVector_Destroy( &DispTdT_Pred );

     MatrixVector_Destroy( &VelT );
     MatrixVector_Destroy( &VelT_Pred );
     MatrixVector_Destroy( &VelTdT );
     MatrixVector_Destroy( &VelTdT_Pred );

     MatrixVector_Destroy( &AccT );
     MatrixVector_Destroy( &AccTdT );

     MatrixVector_Destroy( &RForceT );
     MatrixVector_Destroy( &RForceTdT );

     MatrixVector_Destroy( &LoadVectorForm1 );
     MatrixVector_Destroy( &LoadVectorForm2 );
     MatrixVector_Destroy( &LoadVectorForm3 );
     MatrixVector_Destroy( &LoadTdT );
     MatrixVector_Destroy( &LoadT );

     MatrixVector_Destroy( &Disp );
     MatrixVector_Destroy( &Vel );
     MatrixVector_Destroy( &Acc );
     
     /* Free the coupling nodes memory and close sockets if appropiate */
     Substructure_DeleteCouplingNodes( &CNodes );
     
     return 0;
}
