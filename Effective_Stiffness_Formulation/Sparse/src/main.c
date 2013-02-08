#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>  /* For getopt_long() */
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <assert.h>

#include "MatrixVector.h"
#include "Initiation.h"
#include "Precalculations.h"
#include "ErrorHandling.h"
#include "Netlib.h"
#include "ComputeU0.h"
#include "EndingStep.h"
#include "Send_Receive_Data.h"

#include "HDF5_Operations.h"

#if REAL_TIME_
#include <sched.h>    /* For sched_setscheduler( ) */
#include <malloc.h>
#include <sys/mman.h> /* For mlockall( ) */
#include <unistd.h>   /* For sysconf( ) */
#define SOMESIZE (150*1024*1024) // 150MB
#endif

void Print_Help( const char *Program_Name );

int main( int argc, char **argv )
{
     
     /* Output file */
     int hdf5_file;

     unsigned int istep;		/* Counters */

     AlgConst InitCnt;
     const char *FileConf;

     /* NETLIB Variables */
     int incx, incy;
     Scalars Constants;
     
     double *AccAll, *VelAll, *DispAll;

     MatrixVector M, C, K;               /* Mass, Damping and Stiffness matrices */
     MatrixVector Keinv;
     MatrixVector Keinv_c, Keinv_m;

     MatrixVector EffT;

     MatrixVector DispT, DispTdT0, DispTdT;
     MatrixVector DispTdT0_c, DispTdT0_m;
     MatrixVector tempvec;

     MatrixVector VelT, VelTdT;
     MatrixVector AccT, AccTdT;

     MatrixVector LoadVectorForm, LoadTdT, LoadTdT1;

     MatrixVector fc, fcprevsub;
     MatrixVector fu;

     MatrixVector ErrCompT, ErrCompTdT;

     MatrixVector Disp, Vel, Acc;

     Coupling_Node CNodes;

     HDF5_time_t Time;

     /* TCP socket connection Variables */
     int Socket;

     /* Options */
     int Selected_Option;
     struct option long_options[] = {
	  {"help", no_argument, 0, 'h'},
	  {"config-file", required_argument, 0, 'c'},
	  {0, 0, 0, 0}
     };

#if _SPARSE_
     /* Sparse matrices */
     Sp_MatrixVector Sp_M, Sp_C, Sp_K;     /* Sparse representation of the M, C and K matrices */
#endif

#if REAL_TIME_
     struct sched_param sp;
     int policy;

     if((policy = sched_getscheduler(0) == -1)) {
	  fprintf(stderr, "Error getting scheduler" );
     }
     if( policy != SCHED_FIFO ) {
	  printf("Setting priority: SCHED_FIFO.\n");
	  sp.sched_priority = sched_get_priority_max(SCHED_FIFO);
	  if(sched_setscheduler(0, SCHED_FIFO, &sp) == -1 ){
	    PrintErrorAndExit("Filed real time\n");
	  }
     }

     // Allocate some memory
     int page_size;
     char* buffer;       

     // Now lock all current and future pages from preventing of being paged
     if (mlockall(MCL_CURRENT | MCL_FUTURE ))
     {
	  perror("mlockall failed:");
     }
       

     // Turn off malloc trimming.
     mallopt (M_TRIM_THRESHOLD, -1);
       

     // Turn off mmap usage.
     mallopt (M_MMAP_MAX, 0);
       

     page_size = sysconf(_SC_PAGESIZE);
     buffer = malloc(SOMESIZE);
       

     // Touch each page in this piece of memory to get it mapped into RAM
     for (i=0; i < SOMESIZE; i+=page_size){
	  // Each write to this buffer will generate a pagefault.
	  // Once the pagefault is handled a page will be locked in memory and never
	  // given back to the system.
	  buffer[i] = 0;
     }
     free(buffer);
     // buffer is now released. As glibc is configured such that it never gives back memory to
     // the kernel, the memory allocated above is locked for this process. All malloc() and new()
     // calls come from the memory pool reserved and locked above. Issuing free() and delete()
     // does NOT make this locking undone. So, with this locking mechanism we can build C++ applications
     // that will never run into a major/minor pagefault, even with swapping enabled.
       

     //<do your RT-thing>
#endif

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
	  PrintInfo( "Assuming the configuration file to be: ConfFile.conf\n" );

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
     InitConstants( &InitCnt, FileConf );

     /* Read the coupling nodes from a file */
     Read_Coupling_Nodes( &CNodes, InitCnt.OrderSub, InitCnt.DeltaT_Sub, InitCnt.FileCNodes );

     /* Allocate memory for saving the acceleration, displacement and velocity (input files) that will
      * be used during the test */
     if( InitCnt.Use_Absolute_Values ){
	  AccAll = NULL;
	  VelAll = (double *) calloc( (size_t) InitCnt.Nstep, sizeof(double) );
	  DispAll = (double *) calloc( (size_t) InitCnt.Nstep, sizeof(double) );
     } else {
	  AccAll = (double *) calloc( (size_t) InitCnt.Nstep, sizeof(double) );
	  VelAll = NULL;
	  DispAll = NULL;
     }

     /* Initialise the matrices and vectors that will be used in the Time Integration process */
     if( (!InitCnt.Use_Sparse && !InitCnt.Read_Sparse) || (InitCnt.Use_Sparse && !InitCnt.Read_Sparse) || (!InitCnt.Use_Sparse && InitCnt.Read_Sparse) ){
	  Init_MatrixVector( &M, InitCnt.Order, InitCnt.Order );
	  Init_MatrixVector( &K, InitCnt.Order, InitCnt.Order );
     } else if ( InitCnt.Use_Sparse && InitCnt.Read_Sparse ){
	  Set_RowsCols_Sp( &Sp_M, InitCnt.Order, InitCnt.Order );
	  Set_RowsCols_Sp( &Sp_K, InitCnt.Order, InitCnt.Order );
     }

     Init_MatrixVector( &Keinv, InitCnt.Order, InitCnt.Order );
 
     Init_MatrixVector( &Keinv_c, CNodes.Order, CNodes.Order );
     Init_MatrixVector( &Keinv_m, InitCnt.Order - CNodes.Order, CNodes.Order );

     Init_MatrixVector( &tempvec, InitCnt.Order, 1 );

     Init_MatrixVector( &DispT, InitCnt.Order, 1 );
     Init_MatrixVector( &DispTdT0, InitCnt.Order, 1 );
     Init_MatrixVector( &DispTdT, InitCnt.Order, 1 );

     Init_MatrixVector( &DispTdT0_c, CNodes.Order, 1 );
     Init_MatrixVector( &DispTdT0_m, InitCnt.Order-CNodes.Order, 1 );

     Init_MatrixVector( &VelT, InitCnt.Order, 1 );
     Init_MatrixVector( &VelTdT, InitCnt.Order, 1 );

     Init_MatrixVector( &AccT, InitCnt.Order, 1 );
     Init_MatrixVector( &AccTdT, InitCnt.Order, 1 );

     Init_MatrixVector( &EffT, InitCnt.Order, 1 );

     Init_MatrixVector( &LoadVectorForm, InitCnt.Order, 1 );
     Init_MatrixVector( &LoadTdT, InitCnt.Order, 1 );
     Init_MatrixVector( &LoadTdT1, InitCnt.Order, 1 );

     Init_MatrixVector( &fc, InitCnt.Order, 1 );
     Init_MatrixVector( &fcprevsub, CNodes.Order, 1 );
     Init_MatrixVector( &fu, InitCnt.Order, 1 );

     Init_MatrixVector( &ErrCompT, InitCnt.Order, 1 );
     Init_MatrixVector( &ErrCompTdT, InitCnt.Order, 1 );

     Init_MatrixVector( &Disp, InitCnt.Order, 1 );
     Init_MatrixVector( &Vel, InitCnt.Order, 1 );
     Init_MatrixVector( &Acc, InitCnt.Order, 1 );

     /* Read the matrices from a file */
     if( !InitCnt.Read_Sparse ){
	  MatrixVector_From_File( &M, InitCnt.FileM );
	  MatrixVector_From_File( &K, InitCnt.FileK );
     } else if ( InitCnt.Read_Sparse && !InitCnt.Use_Sparse ){
	  MatrixVector_From_File_Sp2Dense( &M, InitCnt.FileM );
	  Dense_to_CSR( &M, &Sp_M, 0 );            /* Transform into CSR format */
	  MatrixVector_From_File_Sp2Dense( &K, InitCnt.FileK );
	  Dense_to_CSR( &K, &Sp_K, 0 );            /* Transform into CSR format */
     } else if ( InitCnt.Read_Sparse && InitCnt.Use_Sparse ){
#if _SPARSE_
	  MatrixVector_From_File_Sp( &Sp_M, InitCnt.FileM );
	  MatrixVector_From_File_Sp( &Sp_K, InitCnt.FileK );
#endif
     } else {
	  assert( InitCnt.Read_Sparse && !InitCnt.Use_Sparse );
     }

     if ( !InitCnt.Read_LVector ){
	  Generate_LoadVectorForm( &LoadVectorForm, InitCnt.ExcitedDOF );
     }  else {
	  MatrixVector_From_File_Sp2Dense( &LoadVectorForm, InitCnt.FileLV );
     }

     /* Calculate damping matrix using Rayleigh. C = alpha*M + beta*K */
     if ( InitCnt.Use_Sparse && InitCnt.Read_Sparse ) {

#if _SPARSE_
	  Init_MatrixVector_Sp( &Sp_C, InitCnt.Order, InitCnt.Order, Sp_K.Num_Nonzero );
	  CalculateMatrixC_Sp( &Sp_M, &Sp_K, &Sp_C, &InitCnt.Rayleigh );
#endif
     } else {
	  Init_MatrixVector( &C, InitCnt.Order, InitCnt.Order );
	  CalculateMatrixC( &M, &K, &C, &InitCnt.Rayleigh );
     }

     /* Calculate Matrix Keinv = [K + a0*M + a1*C]^(-1) */
     Constants.Alpha = 1.0;
     Constants.Beta = InitCnt.a0;
     Constants.Gamma = InitCnt.a1;

     if( !InitCnt.Use_Pardiso ){
	  CalculateMatrixKeinv( &Keinv, &M, &C, &K, Constants );
     } else if ( InitCnt.Use_Pardiso && !InitCnt.Use_Sparse ){
#if _SPARSE_
	  CalculateMatrixKeinv_Pardiso( &Keinv, &M, &C, &K, Constants );
#endif
     } else if ( InitCnt.Use_Pardiso && InitCnt.Use_Sparse && InitCnt.Read_Sparse ){
#if _SPARSE_
	  CalculateMatrixKeinv_Pardiso_Sparse( &Keinv, &Sp_M, &Sp_C, &Sp_K,
					       Constants );
#endif
     }

     BuildMatrixXc( &Keinv, Keinv_c.Array, &CNodes );
     BuildMatrixXcm( &Keinv, &Keinv_m, &CNodes );
  
#if _SPARSE_
     /* Transform the matrices into CSR format */
     if( !InitCnt.Read_Sparse && InitCnt.Use_Sparse ){
	  Dense_to_CSR( &M, &Sp_M, 0 );            /* Transform into CSR format */
	  Destroy_MatrixVector( &M );              /* Destroy the dense matrix */

	  Dense_to_CSR( &K, &Sp_K, 0 );            /* Transform into CSR format */
	  Destroy_MatrixVector( &K );              /* Destroy the dense matrix */

	  Dense_to_CSR( &C, &Sp_C, 0 );            /* Transform into CSR format */
	  Destroy_MatrixVector( &C );              /* Destroy the dense matrix */
     }
#endif

     /* Send the coupling part of the effective matrix if we are performing a distributed test */
     Send_Effective_Matrix( Keinv_c.Array, (unsigned int) CNodes.Order, &Socket, InitCnt.Remote );

     /* Read the earthquake data from a file */
     if( InitCnt.Use_Absolute_Values ){
	  ReadDataEarthquake_AbsValues( VelAll, DispAll, InitCnt.Nstep, InitCnt.FileData );
     } else {
	  ReadDataEarthquake_RelValues( AccAll, InitCnt.Nstep, InitCnt.FileData );
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
	  Apply_LoadVectorForm( &Disp, &LoadVectorForm, DispAll[istep - 1] );

	  Apply_LoadVectorForm( &Vel, &LoadVectorForm, VelAll[istep - 1] );

	  if( !InitCnt.Use_Sparse ){
	       Calc_Input_Load_AbsValues( &LoadTdT, &K, &C, &Disp, &Vel );
	  } else {
#if _SPARSE_
	       Calc_Input_Load_AbsValues_Sparse( &LoadTdT, &Sp_K, &Sp_C, &Disp, &Vel );
#endif
	  }
     } else {
	  Apply_LoadVectorForm( &Acc, &LoadVectorForm, AccAll[istep - 1] );
	  if( !InitCnt.Use_Sparse ){
	       Calc_Input_Load_RelValues( &LoadTdT, &M, &Acc );
	  } else {
#if _SPARSE_
	       Calc_Input_Load_RelValues_Sparse( &LoadTdT, &Sp_M, &Acc );
#endif
	  }
     }

     incx = 1; incy = 1;
     PrintInfo( "Starting stepping process.\n" );
     while ( istep <= InitCnt.Nstep ){

	  /* Calculate the effective force vector
	     Fe = M*(a0*u + a2*v + a3*a) + C*(a1*u + a4*v + a5*a) */
	  if( !InitCnt.Use_Sparse ){
	       EffK_Calc_Effective_Force( &M, &C, &DispT, &VelT, &AccT, &tempvec,
					  InitCnt.a0, InitCnt.a1, InitCnt.a2,
					  InitCnt.a3, InitCnt.a4, InitCnt.a5,
					  &EffT );
	  } else {
#if _SPARSE_
	       EffK_Calc_Effective_Force_Sparse( &Sp_M, &Sp_C, &DispT, &VelT, &AccT,
						 &tempvec, InitCnt.a0, InitCnt.a1,
						 InitCnt.a2, InitCnt.a3, InitCnt.a4,
						 InitCnt.a5, &EffT );
#endif
	  }

	  /* Compute the new Displacement u0 */
	  EffK_ComputeU0( &EffT, &LoadTdT, &fu, InitCnt.PID.P, &Keinv, &tempvec, &DispTdT0 );

	  /* Split DispTdT into coupling and non-coupling part */
	  CreateVectorXm( &DispTdT0, &DispTdT0_m, &CNodes );
	  CreateVectorXc( &DispTdT0, DispTdT0_c.Array, &CNodes );

	  /* Perform substepping */
	  Do_Substepping( Keinv_c.Array, DispTdT0_c.Array, DispTdT.Array, fcprevsub.Array, fc.Array, InitCnt.Remote.Type,
			  InitCnt.Delta_t*(double) istep, Socket, &CNodes, InitCnt.NSubstep, InitCnt.DeltaT_Sub  );

	  if ( istep < InitCnt.Nstep ){
	       /* Calculate the input load for the next step during the
		  sub-stepping process. */
	       if( InitCnt.Use_Absolute_Values ){
		    Apply_LoadVectorForm( &Disp, &LoadVectorForm, DispAll[istep] );
		    Apply_LoadVectorForm( &Vel, &LoadVectorForm, VelAll[istep] );

		    if( !InitCnt.Use_Sparse ){
			 Calc_Input_Load_AbsValues( &LoadTdT1, &K, &C, &Disp, &Vel );
		    } else {
#if _SPARSE_
			 Calc_Input_Load_AbsValues_Sparse( &LoadTdT1, &Sp_K, &Sp_C, &Disp, &Vel );
#endif
		    }
	       } else {
		    Apply_LoadVectorForm( &Acc, &LoadVectorForm, AccAll[istep] );

		    if( !InitCnt.Use_Sparse ){
			 Calc_Input_Load_RelValues( &LoadTdT1, &M, &Acc );
		    } else {
#if _SPARSE_
			 Calc_Input_Load_RelValues_Sparse( &LoadTdT1, &Sp_M, &Acc );
#endif
		    }
	       }
	  }

	  /* Join the non-coupling part. DispTdT_m = Keinv_m*fc + DispTdT0_m. Although DispTdT0_m is what has been received from the other computer,
	     it has the name of DispTdT_m to avoid further operations if using the NETLIB libraries. */
	  JoinNonCouplingPart( &DispTdT0_m, &Keinv_m, &fcprevsub, &DispTdT, &CNodes );

	  /* Compute acceleration ai1 = a0*(ui1 -ui) - a2*vi -a3*ai */
	  Compute_Acceleration( &DispTdT, &DispT, &VelT, &AccT, InitCnt.a0, InitCnt.a2, InitCnt.a3, &AccTdT );

	  /* Compute Velocity. vi = vi + a6*ai +a7*ai */
	  Compute_Velocity( &VelT, &AccT, &AccTdT, InitCnt.a6, InitCnt.a7, &VelTdT );

	  /* Error Compensation. fu = LoatTdT + fc -(Mass*AccTdT + Damp*VelTdT + Stiff*DispTdT) */
	  if( InitCnt.PID.P != 0.0 || InitCnt.PID.I != 0.0 || InitCnt.PID.D != 0.0 ){
	       if( !InitCnt.Use_Sparse ){
		    Compute_Force_Error( &M, &C, &K, &AccTdT, &VelTdT, &DispTdT, &fc, &LoadTdT, &fu );
	       } else {
#if _SPARSE_
		    Compute_Force_Error_Sparse( &Sp_M, &Sp_C, &Sp_K, &AccTdT, &VelTdT, &DispTdT, &fc, &LoadTdT, &fu );
#endif
	       }
	  }

	  /* Save the result in a HDF5 file format */
	  HDF5_Store_TimeHistoryData( hdf5_file, &AccTdT, &VelTdT, &DispTdT, &LoadTdT, &fc, &fu, (int) istep, &InitCnt );

	  /* Backup vectors */
	  dcopy_( &LoadTdT1.Rows, LoadTdT1.Array, &incx, LoadTdT.Array, &incy ); /* li = li1 */
	  dcopy_( &DispTdT.Rows, DispTdT.Array, &incx, DispT.Array, &incy ); /* ui = ui1 */
	  dcopy_( &VelTdT.Rows, VelTdT.Array, &incx, VelT.Array, &incy ); /* vi = vi1 */
	  dcopy_( &AccTdT.Rows, AccTdT.Array, &incx, AccT.Array, &incy ); /* ai = ai1 */
	  istep = istep + 1;
     }

     gettimeofday( &Time.end, NULL );
     Time.Elapsed_time = (double) (Time.end.tv_sec - Time.start.tv_sec)*1000.0;
     Time.Elapsed_time += (double) (Time.end.tv_usec - Time.start.tv_usec)/1000.0;
     HDF5_StoreTime( hdf5_file, &Time );
     printf( "The stepping process has finished in %lf ms\n", Time.Elapsed_time );


     /* Close the Connection */
     if( InitCnt.Remote.Type != NO_PROTOCOL ){
	  Close_Connection( &Socket, hdf5_file, InitCnt.Remote.Type, (unsigned int) CNodes.Order, InitCnt.Nstep, 4 );
     }

     HDF5_CloseFile( hdf5_file );

     /* Free initiation values */
     Delete_InitConstants( &InitCnt );

     /* Free the memory */
     if( InitCnt.Use_Absolute_Values ){
	  free( VelAll );
	  free( DispAll );
     } else {
	  free( AccAll );
     }

     /* Free the coupling nodes memory */
     Delete_CouplingNodes( &CNodes );

     /* Free Time string */
     free( Time.Date_time );

     /* Destroy the data structures */
     if( (!InitCnt.Use_Sparse && !InitCnt.Read_Sparse) || (!InitCnt.Use_Sparse && InitCnt.Read_Sparse) ){
	  Destroy_MatrixVector( &M );
	  Destroy_MatrixVector( &K );
	  Destroy_MatrixVector( &C );
     } else if ( InitCnt.Use_Sparse ){
#if _SPARSE_
	  Destroy_MatrixVector_Sparse( &Sp_M );
	  Destroy_MatrixVector_Sparse( &Sp_C );
	  Destroy_MatrixVector_Sparse( &Sp_K );
#endif
     }
	  
     Destroy_MatrixVector( &Keinv );
     Destroy_MatrixVector( &Keinv_c );
     Destroy_MatrixVector( &Keinv_m );

     Destroy_MatrixVector( &tempvec );

     Destroy_MatrixVector( &DispT );
     Destroy_MatrixVector( &DispTdT0 );
     Destroy_MatrixVector( &DispTdT0_c );
     Destroy_MatrixVector( &DispTdT0_m );
     Destroy_MatrixVector( &DispTdT );

     Destroy_MatrixVector( &VelT );
     Destroy_MatrixVector( &VelTdT );

     Destroy_MatrixVector( &AccT );
     Destroy_MatrixVector( &AccTdT );

     Destroy_MatrixVector( &LoadVectorForm );
     Destroy_MatrixVector( &LoadTdT );
     Destroy_MatrixVector( &LoadTdT1 );

     Destroy_MatrixVector( &EffT );

     Destroy_MatrixVector( &fc );
     Destroy_MatrixVector( &fcprevsub );
     Destroy_MatrixVector( &fu );

     Destroy_MatrixVector( &ErrCompT );
     Destroy_MatrixVector( &ErrCompTdT );

     Destroy_MatrixVector( &Disp );
     Destroy_MatrixVector( &Vel );
     Destroy_MatrixVector( &Acc );

     return 0;
}

void Print_Help( const char *Program_Name )
{

     fprintf( stderr, "Usage: %s [-h] -c <ConfFile", Program_Name );
     fprintf( stderr,
	      "  -h  --help           This help text.\n"
	      "  -c  --config-file    The name of the configuration file. Default value: ConfFile.conf\n" );
}
