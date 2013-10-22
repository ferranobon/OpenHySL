#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>  /* For getopt_long() */
#include <time.h>

#include "Print_Messages.h"

#include "MatrixVector.h"
#include "Input_Load.h"

#include "Auxiliary_Math.h"
#include "Algorithm_Aux.h"
#include "Substructure_Exact.h"

#include "HDF5_Operations.h"

#include "Definitions.h"

#if _MKL_
#include "mkl_blas.h"
#else
#include "Netlib.h"
#endif

int main( int argc, char **argv )
{
     unsigned int istep;
     int incx, incy, i;

     AlgConst_t InitCnt;

     hid_t hdf5_file;

     const char *FileConf;

     MatrixVector_t Mass, Stiff, Damping_Ratios;
     MatrixVector_t EValues, EVectors;

     MatrixVector_t Init_Disp, Init_Vel;
     MatrixVector_t Load;
     MatrixVector_t End_Disp, End_Vel, End_Acc;

     MatrixVector_t LoadVectorForm, Acc;

     HYSL_FLOAT *AccAll;

     SaveTime_t     Time;
     /* Options */
     int Selected_Option;
     struct option long_options[] = {
	  {"help", no_argument, 0, 'h'},
	  {"config-file", required_argument, 0, 'c'},
	  {0, 0, 0, 0}
     };

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

     /* Allocate memory for saving the acceleration (input files) that will be used during the test */
     AccAll = (HYSL_FLOAT *) calloc( (size_t) InitCnt.NStep, sizeof(HYSL_FLOAT) );

     /* Read the matrices */
     MatrixVector_Create( InitCnt.Order, InitCnt.Order, &Mass );
     MatrixVector_Create( InitCnt.Order, InitCnt.Order, &Stiff );
     MatrixVector_Create( InitCnt.Order, 1, &Damping_Ratios );

     MatrixVector_Create( InitCnt.Order, InitCnt.Order, &EVectors );
     MatrixVector_Create( InitCnt.Order, 1, &EValues );

     MatrixVector_Create( InitCnt.Order, 1, &Init_Disp );
     MatrixVector_Create( InitCnt.Order, 1, &Init_Vel );

     MatrixVector_Create( InitCnt.Order, 1, &End_Disp );
     MatrixVector_Create( InitCnt.Order, 1, &End_Vel );
     MatrixVector_Create( InitCnt.Order, 1, &End_Acc );

     MatrixVector_Create( InitCnt.Order, 1, &Acc );
     MatrixVector_Create( InitCnt.Order, 1, &Load );
     MatrixVector_Create( InitCnt.Order, 1, &LoadVectorForm );

     /* Read Mass and Stiffness matrices in MatrixMarket format */
     MatrixVector_FromFile_MM( InitCnt.FileM, &Mass );
     MatrixVector_FromFile_MM( InitCnt.FileK, &Stiff );

     for( i = 0; i < InitCnt.Order; i++ ){
	  if( Mass.Array[i*InitCnt.Order + i] == 0.0 ){
#if _FLOAT_
	       Mass.Array[i*InitCnt.Order + i] = 1E-12f;
#else
	       Mass.Array[i*InitCnt.Order + i] = 1E-12;
#endif
	  }
     }

     /* Read the load vector form */
     if ( !InitCnt.Read_LVector ){
	  InputLoad_Generate_LoadVectorForm( InitCnt.ExcitedDOF, &LoadVectorForm );
     }  else {
	  MatrixVector_FromFile_MM( InitCnt.FileLV, &LoadVectorForm );
     }

     /* Read the earthquake data from a file */
     Algorithm_ReadDataEarthquake_RelValues( InitCnt.NStep, InitCnt.FileData, InitCnt.Scale_Factor,
					     AccAll );

     /* Compute Eigenvalues and eigenvectors */
     Compute_Eigenvalues_Eigenvectors( &Stiff, &Mass, &EValues, &EVectors );

     for( i = 0; i < 12; i++ ){
	  printf("%lE\t", EValues.Array[i] );
     }
     printf( "\n" );

     /* Compute the damping ratios */
     Compute_DampingRatios_Rayleigh( InitCnt.Rayleigh.Alpha, InitCnt.Rayleigh.Beta, EValues.Rows, EValues.Array,
				     Damping_Ratios.Array );
     /* Open Output file. If the file cannot be opened, the program will exit, since the results cannot be
      * stored. */
#if _HDF5_
     hdf5_file = HDF5_CreateFile( InitCnt.FileOutput );
     HDF5_CreateGroup_TimeIntegration( hdf5_file, &InitCnt );
#endif

     Time.Date_start = time( NULL );
     Time.Date_time = strdup( ctime( &Time.Date_start) );
     gettimeofday( &Time.Start, NULL );        

     istep = 1;
     incx = 1; incy = 1;
     while ( istep <= InitCnt.NStep ){

	  InputLoad_Apply_LoadVectorForm( &LoadVectorForm, AccAll[istep - 1], &Acc );
	  InputLoad_RelValues( &Mass, &Acc, &Load );

	  Duhamel_Integral( &Mass, &EValues, &EVectors, &Damping_Ratios, &Init_Disp,
			    &Init_Vel, &Load, &End_Disp, &End_Vel, &End_Acc, InitCnt.Delta_t );
#if _HDF5_
	  HDF5_Store_TimeHistoryData( hdf5_file, &End_Acc, &End_Vel, &End_Disp, NULL, NULL, (int) istep, &InitCnt );
#endif
	  /* Backup vectors */
	  hysl_copy( &End_Vel.Rows, End_Vel.Array, &incx, Init_Vel.Array, &incy );
	  hysl_copy( &End_Disp.Rows, End_Disp.Array, &incx, Init_Disp.Array, &incy );

	  istep = istep + 1;
     }

     gettimeofday( &Time.End, NULL );
     Time.Elapsed_time = (double) (Time.End.tv_sec - Time.Start.tv_sec)*1000.0;
     Time.Elapsed_time += (double) (Time.End.tv_usec - Time.Start.tv_usec)/1000.0;
#if _HDF5_
     HDF5_Store_Time( hdf5_file, &Time );
     HDF5_CloseFile( &hdf5_file );
#endif

     Print_Header( SUCCESS );
     printf( "The time integration process has finished in %lf ms.\n", Time.Elapsed_time );

     /* Free initiation values */
     Algorithm_Destroy( &InitCnt );

     free( AccAll );

     MatrixVector_Destroy( &Mass );
     MatrixVector_Destroy( &Stiff );
     MatrixVector_Destroy( &Damping_Ratios );

     MatrixVector_Destroy( &EVectors );
     MatrixVector_Destroy( &EValues );

     MatrixVector_Destroy( &Init_Disp );
     MatrixVector_Destroy( &Init_Vel );

     MatrixVector_Destroy( &End_Disp );
     MatrixVector_Destroy( &End_Vel );
     MatrixVector_Destroy( &End_Acc );

     MatrixVector_Destroy( &Load );
     MatrixVector_Destroy( &LoadVectorForm );
     MatrixVector_Destroy( &Acc );

     return 0;
}
