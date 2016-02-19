#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <getopt.h>  /* For getopt_long() */
#include <math.h>
#include <time.h>

#include "Print_Messages.h"

#include "Algorithm_Aux.h"

#include "MatrixVector.h"
#include "Input_Load.h"
#include "New_State.h"
#include "EffK_Formulation.h"
#include "Error_Compensation.h"
#include "Rayleigh.h"
#include "Modal_Damping.h"
#include "GainMatrix.h"
#include "Substructure.h"
#include "Substructure_Experimental.h"
#include "Substructure_Exact.h"
#include "Substructure_Newmark.h"
#include "Substructure_Auxiliary.h"  /* For Substructure_VectorXm(), Substructure_VectorXc(), ... */

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


FILE *GlobalOutput;

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

     unsigned int istep, i;
     double temp1 = 0.0, temp2 = 0.0, temp3 = 0.0, done = 1.0, dzero = 0.0;
     AlgConst_t InitCnt;
     const char *FileConf;

#if _HDF5_
     hid_t hdf5_file;
#endif

     /* NETLIB Variables */
     int incx = 1, incy = 1;
     Scalars_t Constants;
     
     HYSL_FLOAT *AccAll, *VelAll, *DispAll;
     HYSL_FLOAT Alpha;
     MatrixVector_t M, C, K;               /* Mass, Damping and Stiffness matrices */
     MatrixVector_t Meinv;
     MatrixVector_t MatA, MatB;
     MatrixVector_t G_c, G_m;

     MatrixVector_t EffT;

     MatrixVector_t DispT0, DispT, DispTdT0, DispTdT;
     MatrixVector_t DispTdT0_c, DispTdT0_m;
     MatrixVector_t tempvec;

     MatrixVector_t VelT, VelTdT;
     MatrixVector_t AccT, AccTdT;

     MatrixVector_t LoadVectorForm, LoadTdT, LoadTdT1;

     MatrixVector_t ForceT0, ForceT;
     MatrixVector_t fc, fcprevsub;
     MatrixVector_t fu;

     MatrixVector_t Disp, Vel, Acc;

#if _SPARSE_
     /* Sparse matrices */
     MatrixVector_Sp_t Sp_M, Sp_C, Sp_K;     /* Sparse representation of the M, C and K matrices */
     MatrixVector_Sp_t Sp_MatA, Sp_MatB;     /* Constant matrices */
#endif

     CouplingNode_t CNodes;
     SaveTime_t     Time;

     NewmarkSim_t  *TMD;

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
     AccAll = (HYSL_FLOAT *) calloc( (size_t) InitCnt.NStep, sizeof(HYSL_FLOAT) );
     VelAll = (HYSL_FLOAT *) calloc( (size_t) InitCnt.NStep, sizeof(HYSL_FLOAT) );
     DispAll = (HYSL_FLOAT *) calloc( (size_t) InitCnt.NStep, sizeof(HYSL_FLOAT) );

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

     if( !InitCnt.Use_Packed ){
	  MatrixVector_Create( InitCnt.Order, InitCnt.Order, &Meinv );
     } else {
	  MatrixVector_Create_PS( InitCnt.Order, InitCnt.Order, &Meinv );
     }

     if( CNodes.Order >= 1 ){
	  MatrixVector_Create( CNodes.Order, CNodes.Order, &G_c );
	  MatrixVector_Create( InitCnt.Order - CNodes.Order, CNodes.Order, &G_m );
     }

     MatrixVector_Create( InitCnt.Order, 1, &tempvec );

     MatrixVector_Create( InitCnt.Order, 1, &DispT0 );
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
     MatrixVector_Create( InitCnt.Order, 1, &ForceT );
     MatrixVector_Create( InitCnt.Order, 1, &ForceT0 );

     MatrixVector_Create( InitCnt.Order, 1, &fc );
     if( CNodes.Order >= 1 ){
	  MatrixVector_Create( CNodes.Order, 1, &fcprevsub );
     }
     MatrixVector_Create( InitCnt.Order, 1, &fu );

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
	  InputLoad_Generate_LoadVectorForm( InitCnt.ExcitedDOF, &LoadVectorForm );
     }  else {
	  MatrixVector_FromFile_MM( InitCnt.FileLV, &LoadVectorForm );
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

     if ( InitCnt.Use_Sparse ){
	  MatrixVector_Create_Sp( InitCnt.Order, InitCnt.Order, Sp_K.Num_Nonzero, &Sp_MatA );
	  MatrixVector_Create_Sp( InitCnt.Order, InitCnt.Order, Sp_K.Num_Nonzero, &Sp_MatB );
     } else if ( InitCnt.Use_Packed ){
	  MatrixVector_Create_PS( InitCnt.Order, InitCnt.Order, &MatA );
	  MatrixVector_Create_PS( InitCnt.Order, InitCnt.Order, &MatB );
     } else if ( !InitCnt.Use_Packed && !InitCnt.Use_Sparse ){
	  MatrixVector_Create( InitCnt.Order, InitCnt.Order, &MatA );
	  MatrixVector_Create( InitCnt.Order, InitCnt.Order, &MatB );
     }

     Constants.Alpha = 2.0;           /* Mass matrix coefficient */
     Constants.Beta = -InitCnt.a16;     /* Damping matrix coefficent */
     Constants.Gamma = -InitCnt.a17;    /* Stiffness matrix coefficient */
     Constants.Lambda = 1.0;          /* Matrix inversion coefficient */
     if( !InitCnt.Use_Sparse && !InitCnt.Use_Packed ){
	  MatrixVector_Add3Mat( &M, &C, &K, Constants, &MatA );
     } else if ( InitCnt.Use_Sparse ){
	  MatrixVector_Add3Mat_Sp( &Sp_M, &Sp_C, &Sp_K, Constants, &Sp_MatA );
     } else if ( InitCnt.Use_Packed ){
	  MatrixVector_Add3Mat_PS( &M, &C, &K, Constants, &MatA );	  
     }


     Constants.Alpha = -1.0;           /* Mass matrix coefficient */
     Constants.Beta = InitCnt.a6;     /* Damping matrix coefficent */
     Constants.Gamma = -InitCnt.a18;    /* Stiffness matrix coefficient */
     Constants.Lambda = 1.0;          /* Matrix inversion coefficient */
     if( !InitCnt.Use_Sparse && !InitCnt.Use_Packed ){
	  MatrixVector_Add3Mat( &M, &C, &K, Constants, &MatB );
     } else if ( InitCnt.Use_Sparse ){
	  MatrixVector_Add3Mat_Sp( &Sp_M, &Sp_C, &Sp_K, Constants, &Sp_MatB );
     } else if ( InitCnt.Use_Packed ){
	  MatrixVector_Add3Mat_PS( &M, &C, &K, Constants, &MatB );	  
     }

     /* Calculate Matrix G = 1.0*[M + a7*C + a8*K]^(-1) */
     Constants.Alpha = 1.0;           /* Mass matrix coefficient */
     Constants.Beta = InitCnt.a7;     /* Damping matrix coefficent */
     Constants.Gamma = InitCnt.a8;    /* Stiffness matrix coefficient */
     Constants.Lambda = 1.0;          /* Matrix inversion coefficient */

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

     if( CNodes.Order >= 1 ){
	  if( !InitCnt.Use_Packed ){
	       Substructure_MatrixXc( &Meinv, &CNodes, &G_c );
	       Substructure_MatrixXcm( &Meinv, &CNodes, &G_m );
	  } else {
	       Substructure_MatrixXc_PS( &Meinv, &CNodes, &G_c );
	       Substructure_MatrixXcm_PS( &Meinv, &CNodes, &G_m );
	  }
	  hysl_scal( &G_c.Rows, &InitCnt.a8, G_c.Array, &incx );
	  hysl_scal( &G_m.Rows, &InitCnt.a8, G_m.Array, &incx );

	  /* Send the coupling part of the effective matrix if we are performing a distributed test */
	  for( i = 0; i < (unsigned int) CNodes.Order; i++ ){
	       if( CNodes.Sub[i].Type == EXP_ADWIN || CNodes.Sub[i].Type == REMOTE ){
		    Substructure_SendGainMatrix( G_c.Array, (unsigned int) CNodes.Order, &CNodes.Sub[i] );
	       }
	  }
     }
     /* Read the earthquake data from a file */
     Algorithm_ReadDataEarthquake( InitCnt.NStep, InitCnt.FileData, InitCnt.Scale_Factor, AccAll,
				   VelAll, DispAll );

     /* Open Output file. If the file cannot be opened, the program will exit, since the results cannot be
      * stored. */
#if _HDF5_
     hdf5_file = HDF5_CreateFile( Concatenate_Strings( 2, InitCnt.FileOutput, ".h5" ) );
     HDF5_CreateGroup_Parameters( hdf5_file, &InitCnt, &CNodes, AccAll, VelAll, DispAll );
     HDF5_CreateGroup_TimeIntegration( hdf5_file, &InitCnt );
#endif

     GlobalOutput = fopen( Concatenate_Strings( 2, InitCnt.FileOutput, ".txt" ), "w" );
     fprintf( GlobalOutput, "Step\tAcc1Right [m/s^2]\tAcc1Left [m/s^2]\tAcc2Right [m/s^2]\tAcc2Left [m/s^2]\tAcc2Right [m/s^2]\tAcc2Left [m/s^2]\t" );
     fprintf( GlobalOutput, "Disp1 [m]\tDisp2 [m]\tDisp3Mid [m]\tDisp3Left [m]\tDisp3Right [m]\t" );
     fprintf( GlobalOutput, "Coupling Force [N]\tError Force[N]\t" );
     fprintf( GlobalOutput, "Acc TMD [m/s^2]\tVel TMD [m/s]\tDisp TMD [m]\n" );
     TMD = CNodes.Sub[0].SimStruct;
     
     /* Calculate the input load */
     istep = 1;
     if( InitCnt.Use_Absolute_Values ){
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

     Time.Date_start = time( NULL );
     Time.Date_time = strdup( ctime( &Time.Date_start) );
     gettimeofday( &Time.Start, NULL );

     Print_Header( INFO );
     printf( "Starting stepping process.\n" );
//     double Ffs = 60.0, Ffst = 325.0, epsilon = 0.000001;
     double Ffs = 0.0, Ffst = 0.0, epsilon = 0.0;

     while ( istep <= InitCnt.NStep ){

	  /* Calculate the effective force vector Fe = M*(a0*u + a2*v + a3*a) + C*(a1*u + a4*v + a5*a) */
	  if( !InitCnt.Use_Sparse && !InitCnt.Use_Packed ){
	       Compute_NewState_Zienkiewicz( &Meinv, &MatA, &MatB, &DispT, &DispT0, &LoadTdT, &fu, &ForceT,
					     &ForceT0, InitCnt.a8, InitCnt.a17, InitCnt.a18, &tempvec, &DispTdT0 );
	  } else if( !InitCnt.Use_Sparse && InitCnt.Use_Packed ){
	       Compute_NewState_Zienkiewicz_PS( &Meinv, &MatA, &MatB, &DispT, &DispT0, &LoadTdT, &fu, &ForceT,
						&ForceT0, InitCnt.a8, InitCnt.a17, InitCnt.a18, &tempvec, &DispTdT0 );
	  } else if( InitCnt.Use_Sparse ) {
#if _SPARSE_
	       Compute_NewState_Zienkiewicz_Sp( &Meinv, &Sp_MatA, &Sp_MatB, &DispT, &DispT0, &LoadTdT, &fu, &ForceT,
						&ForceT0, InitCnt.a8, InitCnt.a17, InitCnt.a18, &tempvec, &DispTdT0 );
#endif
	  }

	  /* Split DispTdT into coupling and non-coupling part */
	  if( CNodes.Order >= 1 ){
	       Substructure_VectorXm( &DispTdT0, &CNodes, &DispTdT0_m );
	       Substructure_VectorXc( &DispTdT0, &CNodes, &DispTdT0_c );
	  }

	  /* Perform substepping */
	  if( CNodes.Order >= 1 ){
	       if( InitCnt.Use_Absolute_Values ){
		    Substructure_Substepping( G_c.Array, DispTdT0_c.Array, InitCnt.Delta_t*(HYSL_FLOAT) istep, 0.0,
					      InitCnt.NSubstep, InitCnt.DeltaT_Sub, &CNodes, DispTdT.Array, fcprevsub.Array,
					      fc.Array );
	       } else {
		    Substructure_Substepping( G_c.Array, DispTdT0_c.Array, InitCnt.Delta_t*(HYSL_FLOAT) istep,
					      AccAll[istep-1], InitCnt.NSubstep, InitCnt.DeltaT_Sub, &CNodes,
					      DispTdT.Array, fcprevsub.Array, fc.Array );
	       }
	  }

	  if ( istep < InitCnt.NStep ){
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

	  /* Join the non-coupling part. DispTdT_m = G_m*fc + DispTdT0_m. Although DispTdT0_m is what has
	   * been received from the other computer, it has the name of DispTdT_m to avoid further operations
	   * if using the NETLIB libraries. */
	  if( CNodes.Order >= 1 ){
	       Substructure_JoinNonCouplingPart( &DispTdT0_m, &G_m, &fcprevsub, &CNodes,  &DispTdT );
	  } else {
	       hysl_copy( &DispTdT0.Rows, DispTdT0.Array, &incx, DispTdT.Array, &incy ); /* ui = ui1 */
	  }

	  /* Compute acceleration ai1 = a0*(ui1 -ui) - a2*vi -a3*ai */
	  EffK_ComputeAcceleration( &DispTdT, &DispT, &VelT, &AccT, InitCnt.a0, InitCnt.a2, InitCnt.a3, &AccTdT );

	  /* Compute Velocity. vi = a1*(ui1 - ui) - a4*vi - a5*ai */
	  EffK_ComputeVelocity( &DispTdT, &DispT, &VelT, &AccT, InitCnt.a1, InitCnt.a4, InitCnt.a5, &VelTdT );

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

	  /* Print the rest of the information */
	  temp1 = TMD->aTdT + AccTdT.Array[47-1] + AccAll[istep-1];
	  temp2 = TMD->vTdT;// - VelTdT.Array[CNodes.Array[0]-1];
	  temp3 = TMD->dTdT;// - DispTdT.Array[CNodes.Array[0]-1];
	  fprintf( GlobalOutput, "%d\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t", istep, AccAll[istep-1] + AccTdT.Array[811-1], AccAll[istep-1] + AccTdT.Array[799-1], AccAll[istep-1]+ AccTdT.Array[929-1], AccAll[istep-1] + AccTdT.Array[942-1], AccAll[istep-1] + AccTdT.Array[304-1], AccAll[istep-1] + AccTdT.Array[258-1] );
	  fprintf( GlobalOutput, "%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t", DispAll[istep-1] + DispTdT.Array[720-1], DispAll[istep-1] + DispTdT.Array[227-1], DispAll[istep-1] + DispTdT.Array[179-1], DispAll[istep-1] + DispTdT.Array[258-1], DispAll[istep-1] + DispTdT.Array[304-1], DispAll[istep-1] + DispTdT.Array[CNodes.Array[0]-1] );
	  fprintf( GlobalOutput, "%lE\t%lE\t", fc.Array[CNodes.Array[0]-1], fu.Array[CNodes.Array[0]-1] );
	  fprintf( GlobalOutput, "%lE\t%lE\t%lE\n", temp1, temp2, temp3 );

	  /* Save the result in a HDF5 file format */
#if _HDF5_
	  HDF5_Store_TimeHistoryData( hdf5_file, &AccTdT, &VelTdT, &DispTdT, &fc, &fu, (int) istep, &InitCnt );
	  /* Store the relative displacements of the TMD */
//	  HDF5_Store_TMD( hdf5_file, &temp1, &temp2, &temp3, (int) istep );
#else
	  
#endif
	  /* Backup vectors */
	  hysl_copy( &ForceT0.Rows, ForceT.Array, &incx, ForceT0.Array, &incy ); /* ForceT0 = ForceT */
	  hysl_copy( &ForceT.Rows, LoadTdT.Array, &incx, ForceT.Array, &incy ); /* ForceT = li */
	  Alpha = 1.0;
	  hysl_axpy( &ForceT.Rows, &Alpha, fc.Array, &incx, ForceT.Array, &incy ); /* ForceT = li + fc */

	  hysl_copy( &LoadTdT1.Rows, LoadTdT1.Array, &incx, LoadTdT.Array, &incy ); /* li = li1 */
	  hysl_copy( &DispT.Rows, DispT.Array, &incx, DispT0.Array, &incy ); /* ui0 = ui */
	  hysl_copy( &DispTdT.Rows, DispTdT.Array, &incx, DispT.Array, &incy ); /* ui = ui1 */
	  hysl_copy( &VelTdT.Rows, VelTdT.Array, &incx, VelT.Array, &incy ); /* vi = vi1 */
	  hysl_copy( &AccTdT.Rows, AccTdT.Array, &incx, AccT.Array, &incy ); /* ai = ai1 */
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

#if _ADWIN_
     /* Save the data in ADwin into the HDF5 File */
     if( CNodes.Order >= 1 ){
	  for( i = 0; i < (unsigned int) CNodes.Order; i++ ){
	       if( CNodes.Sub[i].Type == EXP_ADWIN ){
#if _HDF5_
		    ADwin_SaveData_HDF5( hdf5_file, (int) InitCnt.NStep, (int) InitCnt.NSubstep,
					NUM_CHANNELS, Entry_Names, 97 );
#endif
	       }
	  }
     }
#endif

#if _HDF5_
     HDF5_CloseFile( &hdf5_file );
#endif

     fclose( GlobalOutput );

     /* Free initiation values */
     Algorithm_Destroy( &InitCnt );

     /* Free the memory */
     free( AccAll );
     free( VelAll );
     free( DispAll );

     /* Free Time string */
     free( Time.Date_time );

     /* Destroy the data structures */
     if( !InitCnt.Use_Sparse && !InitCnt.Use_Packed ){
	  MatrixVector_Destroy( &M );
	  MatrixVector_Destroy( &K );
	  MatrixVector_Destroy( &C );
	  MatrixVector_Destroy( &MatA );
	  MatrixVector_Destroy( &MatB );
     } else if ( !InitCnt.Use_Sparse && InitCnt.Use_Packed ){
	  MatrixVector_Destroy_PS( &M );
	  MatrixVector_Destroy_PS( &K );
	  MatrixVector_Destroy_PS( &C );
	  MatrixVector_Destroy_PS( &MatA );
	  MatrixVector_Destroy_PS( &MatB );
     } else if ( InitCnt.Use_Sparse ) {
#if _SPARSE_
	  MatrixVector_Destroy_Sp( &Sp_M );
	  MatrixVector_Destroy_Sp( &Sp_C );
	  MatrixVector_Destroy_Sp( &Sp_K );
	  MatrixVector_Destroy_Sp( &Sp_MatA );
	  MatrixVector_Destroy_Sp( &Sp_MatB );

#endif
     }

     if ( !InitCnt.Use_Packed ){	  
	  MatrixVector_Destroy( &Meinv );
     } else if ( InitCnt.Use_Packed ){	 
	  MatrixVector_Destroy_PS( &Meinv );
     }
 
     if( CNodes.Order >= 1 ){
	  MatrixVector_Destroy( &G_c );
	  MatrixVector_Destroy( &G_m );
     }

     MatrixVector_Destroy( &tempvec );
 
     MatrixVector_Destroy( &DispT0 );
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
     MatrixVector_Destroy( &ForceT );
     MatrixVector_Destroy( &ForceT0 );

     MatrixVector_Destroy( &EffT );

     MatrixVector_Destroy( &fc );
     if( CNodes.Order >= 1 ){
	  MatrixVector_Destroy( &fcprevsub );
     }
     MatrixVector_Destroy( &fu );

     MatrixVector_Destroy( &Disp );
     MatrixVector_Destroy( &Vel );
     MatrixVector_Destroy( &Acc );

     /* Free the coupling nodes memory and close sockets if appropiate */
     Substructure_DeleteCouplingNodes( &CNodes );

     return 0;
}
