/**
 * \file Full_to_Sparse.c
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 4th of November 2012
 *
 * \brief Dense to sparse matrix conversion.
 *
 * Converts a ASCII file containing a full representation of a dense matrix into a 
 * sparse matrix following the format:
 *
 *     NumRows NumCols Num_Non_Zeros
 *
 *     Row;Col;Value_1
 *     Row;Col;Value_2;
 *     Row;Col;Value_3;
 *     ...
 *     Row;Col;Value_Num_Non_Zeros
 *
 * The entries of rows and columns follow the C naming convention and start at 0.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>  /* For strdup */
#include <getopt.h>  /* For getopt_long */

/**
 * \brief Prints the help for the main routine
 */
void print_help( char **argv, const struct option *long_opts );

int main( int argc, char** argv )
{

     unsigned int i, j;  /* Counters */
     int Rows, Cols;     /* Number of Rows and Columns */
     unsigned int nnz;   /* Number of non-zero elements */
     int rc, idx;
     int help_flag = 0;
     char *InputFileName;          /* Name of the input file */
     double Value;                  /* Read/write variable */
     FILE *InFile, *OutFile;        /* Input and Output files */

     struct option long_options[] = 
     {
	  {"input-file", required_argument, 0, 'i'},
	  {"output-file", required_argument, 0, 'o'},
	  {"num-rows", required_argument, 0, 'r'},
	  {"num-cols", required_argument, 0, 'c'},
	  {"help", no_argument, &help_flag, 'h'},
//	  {0, 0, 0, 0}
     };

     if( argc <= 1 && argc != 4 ){
	  fprintf( stderr, "Invalid number of arguments.\n" );
	  print_help( argv, long_options );
	  return( EXIT_FAILURE );
     }	  
     

     while( (rc = getopt_long( argc, argv, "i:o:r:c:h", long_options, &idx )) != -1 ){

	  switch( rc ){
	  case 'i':
	       InFile = fopen( optarg, "r" );
	       InputFileName = strdup( optarg );
	       if( InFile == NULL ){
		    fprintf( stderr, "Could not open file %s.\n", optarg );
		    fprintf( stderr, "Exiting.\n" );
		    exit( EXIT_FAILURE );
	       }
	       break;
	  case 'o':
	       OutFile = fopen( optarg, "w" );
	       if( OutFile == NULL ){
		    fprintf( stderr, "Could not open file %s.\n", optarg );
		    fprintf( stderr, "Exiting.\n" );
		    exit( EXIT_FAILURE );
	       }
	       break;
	  case 'r':
	       Rows = atoi( optarg );
	       if( Rows <= 0 ){
		    fprintf( stderr, "The number of rows cannot be less or equal than 0\n" );
		    exit( EXIT_FAILURE );
	       }
	       break;
	  case 'c':
	       Cols = atoi( optarg );
	       if( Cols <= 0 ){
		    fprintf( stderr, "The number of columns cannot be less or equal than 0\n" );
		    exit( EXIT_FAILURE );
	       }
	       break;
	  case 'h':
	       help_flag = 1;
	       break;
	  case '?':
	       /* Long options already prints an error message telling that there is an unrecognised option */
	       break;
	  case ':':
	       /* Long options already prints an error message telling that the option requieres an argument */
	       break;
	  }
     }
     /* Print the help output if requested by the user */
     if( help_flag ){
	  print_help( argv, long_options );
	  return ( EXIT_FAILURE );
     }

     /* The input file is read through to determine the number of non-zeros elements
      * in the dense matrix.
      */
     nnz = 0;
     for( i = 0; i < (unsigned int ) Rows; i++ ){
	  for (j = 0; j < (unsigned int) Cols; j++ ){
	       fscanf( InFile, "%lf", &Value );
	       if( j >= i ){
		    if ( Value != 0.0 ){
			 nnz = nnz + 1;
		    }
	       }
	  }
     }
     fclose( InFile );


     /* Open the file containing the dense matrix again in order to copy the non-zero
      * values to the file with the sparse representation
      */
     InFile = fopen( InputFileName, "r" );

     /* Print the header: Number of Rows Number of Columns and number of non-zero elements
      * in the sparse matrix */

     fprintf( OutFile, "%d %d %d\n", Rows, Cols, nnz );
     for( i = 0; i < (unsigned int ) Rows; i++ ){
	  for (j = 0; j < (unsigned int ) Cols; j++ ){
	       fscanf( InFile, "%lf", &Value );
	       if( j >= i ){		    
		    if ( Value != 0.0 ){
			 /* If the value is different than 0, save it in the sparse matrix
			  * representation
			  */
			 fprintf( OutFile, "%d;%d;%le\n", i, j, Value );
		    }
	       }
	  }
     }

     /* Close input and output files */
     fclose( InFile );
     fclose( OutFile );

     return 0;
}

void print_help( char **argv, const struct option *long_opts )
{
     int i;
     struct option *Current_option = NULL;
     const char *Explanation[5] = {"Input data file",
			     "Output data File",
			     "Number of rows in the matrix",
			     "Number of rows in the matrix",
			     "This help text"};

     Current_option= (struct option *) &long_opts[0];

     printf("Usage: %s [-h] -i <InputFile> -o <OutputFile> -r <NumRows> -c <NumCols>\n", argv[0] );
     printf("\n");
    

     printf("Summary of arguments\n");
     i = 0;
     while( (Current_option != NULL) && (Current_option->name != NULL ) ){
	  printf("  -%c  --%s    \t %s.\n", (char) Current_option->val, Current_option->name, Explanation[i] );
	  Current_option = Current_option + 1;
	  i = i + 1;
     }
}
