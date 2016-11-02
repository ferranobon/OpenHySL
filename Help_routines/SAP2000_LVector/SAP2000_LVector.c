#include <stdio.h>
#include <stdlib.h>
#include <string.h>  /* For strlen() */
#include <getopt.h>  /* For getopt_long */

void print_help( char **argv, const struct option *long_opts );

int main ( int argc, char **argv )
{

     int i, Num_Eqs, nnz;
     int Info[7];
     double DOF_Table[7];
     double *LVector;

     char Header[100], *FullString;
     FILE *InFile, *OutFile;
     int rc, idx;

     struct option long_options[] = 
     {
	  {"input-file", required_argument, 0, 'i'},
	  {"output-file", required_argument, 0, 'o'},
	  {"desired-DOF", required_argument, 0, 'd'},
	  {"num-equations", required_argument, 0, 'n'},
	  {"help", no_argument, 0, 'h'},
	  {0, 0, 0, 0}
     };

     if( argc <= 1 && argc != 5 ){
	  fprintf( stderr, "Invalid number of arguments.\n" );
	  print_help( argv, long_options );
	  return( EXIT_FAILURE );
     }	  

     while( (rc = getopt_long( argc, argv, "i:o:d:n:h", long_options, &idx )) != -1 ){

	  switch( rc ){
	  case 'i':
	       InFile = fopen( optarg, "r" );
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
	  case 'n':
	       Num_Eqs = atoi( optarg );
	       if( Num_Eqs <= 0 ){
		    fprintf( stderr, "The number of equations is equal or less than 0" );
		    fprintf( stderr, "Exiting.\n" );
		    exit( EXIT_FAILURE );
	       }
	       break;
	  case 'd':
	       FullString = strdup( optarg );

	       sscanf( FullString, "%lE %lE %lE %lE %lE %lE %lE", &DOF_Table[0], &DOF_Table[1], &DOF_Table[2], &DOF_Table[3], &DOF_Table[4], &DOF_Table[5], &DOF_Table[6]);
	       if( DOF_Table[0] <= 0.0 || DOF_Table[0] > 6.0 ){
		    fprintf( stderr, "The input number of degrees of freedom must be between 1 and 6 and not %lE\n", DOF_Table[0] );
		    fprintf( stderr, "Exiting.\n" );
	       }

	       /* Release the memory */
	       free( FullString );
	       break;
	  case 'h':
	       print_help( argv, long_options );
	       exit( EXIT_FAILURE );
	       break;
	  case '?':
	       /* Long options already prints an error message telling that there is an unrecognised option */
	       break;
	  case ':':
	       /* Long options already prints an error message telling that the option requieres an argument */
	       break;
	  }
     }

     fgets( Header, (size_t) 100, InFile );

     LVector = (double *) calloc( (size_t) Num_Eqs, sizeof(double) );
     nnz = 0;
     while( fscanf( InFile, "%i %i %i %i %i %i %i", &Info[0], &Info[1], &Info[2], &Info[3], &Info[4], &Info[5], &Info[6] ) != EOF ){
	  for( i = 1; i <= DOF_Table[0]; i++ ){
	       if( Info[i] > 0 && DOF_Table[i] > 0.0){
		    LVector[Info[i] - 1] = DOF_Table[i];
		    nnz = nnz + 1;
	       }
	  }
     }

     /* Print the vector in the MatrixMarket format */
     fprintf( OutFile, "%%%%MatrixMarket matrix coordinate real general\n" );
     fprintf( OutFile, "%i 1 %i\n", Num_Eqs, nnz );

     for( i = 0; i < Num_Eqs; i++ ){
	  if( LVector[i] > 0.0 ){
	       fprintf( OutFile, "%i 1 %lf\n", i + 1, LVector[i] );
	  }
     }

     free( LVector );

     return 0;
}

void print_help( char **argv, const struct option *long_opts )
{
     int i;
     struct option *Current_option = NULL;
     const char *Explanation[5] = {"Input matrix file",
				   "Output matrix file in MatrixMarket format",
				   "Desired Degrees of Freedom",
				   "Number of equations. This number is provided in the *.TXA file",
				   "This help text"};

     Current_option= (struct option *) &long_opts[0];

     printf("Usage: %s [-h] -i <InputFile> -o <OutputFile> -d <desiredDOF> -n <Num_Eqs>", argv[0] );
     printf("\n");
    

     printf("Summary of arguments\n");
     i = 0;
     while( (Current_option != NULL) && (Current_option->name != NULL ) ){
	  printf("  -%c  --%s    \t %s.\n", (char) Current_option->val, Current_option->name, Explanation[i] );
	  Current_option = Current_option + 1;
	  i = i + 1;
     }
}
