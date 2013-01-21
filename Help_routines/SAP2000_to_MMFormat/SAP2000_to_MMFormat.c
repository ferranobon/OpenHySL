#include <stdio.h>
#include <stdlib.h>
#include <string.h>  /* For strdup */
#include <getopt.h>  /* For getopt_long */

void print_help( char **argv, const struct option *long_opts );

int main ( int argc, char **argv )
{
     int i;
     int row, column;
     int NRows, NCols, nnz;
     double Value;

     FILE *InFile, *OutFile;
     char Header[100];
     char *InputFileName;
     int rc, idx;

     struct option long_options[] = 
     {
	  {"input-file", required_argument, 0, 'i'},
	  {"output-file", required_argument, 0, 'o'},
	  {"help", no_argument, 0, 'h'},
	  {0, 0, 0, 0}
     };

     if( argc <= 1 && argc != 3 ){
	  fprintf( stderr, "Invalid number of arguments.\n" );
	  print_help( argv, long_options );
	  return( EXIT_FAILURE );
     }	  

     while( (rc = getopt_long( argc, argv, "i:o:h", long_options, &idx )) != -1 ){

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
	  case 'h':
	       print_help( argv, long_options );
	       break;
	  case '?':
	       /* Long options already prints an error message telling that there is an unrecognised option */
	       break;
	  case ':':
	       /* Long options already prints an error message telling that the option requieres an argument */
	       break;
	  }
     }

     nnz = 0;
     row = 0; column = 0;
     NRows = 0; NCols = 0;
     fgets( Header, (size_t) 100, InFile );
     while( fscanf( InFile, "%i %i %lf", &row, &column, &Value ) != EOF ){
	  if( Value != 0.0 ){
	       nnz = nnz + 1;	       
	  }
	  if( row > NRows ){
	       NRows = row;
	       
	  }
	  if( column > NCols ){
	       NCols = column;
	  }
	  
	  if( NCols > NRows ){
	       NRows = NCols;
	  } else if ( NCols < NRows ){
	       NCols = NRows;
	  }
     }
     fclose( InFile );

     InFile = fopen( InputFileName, "r" );
     fgets( Header, (size_t) 100, InFile );

     /* Save in MatrixMarket format */
     fprintf( OutFile, "%%%%MatrixMarket matrix coordinate real symmetric\n" );
     fprintf( OutFile, "%i %i %i\n", NRows, NCols, nnz );

     i = 0;
     while ( i < nnz ){
	  fscanf( InFile, "%i %i %lf", &row, &column, &Value );
	  if( Value != 0.0 ){
	       fprintf( OutFile, "%i %i %.8lf\n", row, column, Value );
	       i = i + 1;
	  }
     }
     fclose( InFile );
     fclose( OutFile );

     return 0;
}

void print_help( char **argv, const struct option *long_opts )
{
     int i;
     struct option *Current_option = NULL;
     const char *Explanation[3] = {"Input matrix file",
			     "Output matrix file in MatrixMarket format",
			     "This help text"};

     Current_option= (struct option *) &long_opts[0];

     printf("Usage: %s [-h] -i <InputFile> -o <OutputFile>", argv[0] );
     printf("\n");
    

     printf("Summary of arguments\n");
     i = 0;
     while( (Current_option != NULL) && (Current_option->name != NULL ) ){
	  printf("  -%c  --%s    \t %s.\n", (char) Current_option->val, Current_option->name, Explanation[i] );
	  Current_option = Current_option + 1;
	  i = i + 1;
     }
}
