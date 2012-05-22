#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

void print_help( char **argv, const struct option *long_opts );

int main( int argc, char **argv )
{

     int i, Num_Elements;
     int idx, rc;
     int help_flag = 0;

     double Disp, Vel, Acc, factor;
     int step;

     FILE *InputFile, *OutputFile;

     struct option long_options[] =
     {
	  {"factor", required_argument, 0, 'f'},
	  {"input-file", required_argument, 0, 'i'},
	  {"help", no_argument, &help_flag, 'h'},
	  {"length", required_argument, 0, 'l'},
	  {"output-file", required_argument, 0, 'o'}
     };

     if( argc <= 1 ){
	  fprintf( stderr, "Invalid number of arguments.\n" );
	  print_help( argv, long_options );
	  return( EXIT_FAILURE );
     }	  

     while( (rc = getopt_long( argc, argv, "f:i:o:l:h", long_options, &idx )) != -1 ){
	  switch( rc ){
	  case 'f':
	       factor = atof( optarg );
	       break;
	  case 'i':
	       InputFile = fopen( optarg, "r" );
	       if( InputFile == NULL ){
		    fprintf( stderr, "Could not open file %s.\n", optarg );
		    fprintf( stderr, "Exiting.\n" );
		    exit( EXIT_FAILURE );
	       }
	       break;
	  case 'o':
	       OutputFile = fopen( optarg, "w" );
	       if( InputFile == NULL ){
		    fprintf( stderr, "Could not open file %s.\n", optarg );
		    fprintf( stderr, "Exiting.\n" );
		    exit( EXIT_FAILURE );
	       }
	       break;
	  case 'l':
	       Num_Elements= atoi( optarg );
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

     for ( i = 0; i < Num_Elements; i++ ){
	  fscanf( InputFile,"%d %lf %lf %lf", &step, &Acc, &Vel, &Disp );
	  fprintf( OutputFile, "%d\t%lf\t%lf\t%lf\n", step, Acc*factor, Vel*factor, Disp*factor );
     }

     fclose( InputFile );
     fclose( OutputFile );
     return 0;
}


void print_help( char **argv, const struct option *long_opts )
{
     int i;
     struct option *Current_option = NULL;
     const char *Explanation[5] = {"Factor to apply",
				   "File containing the input data",
				   "This help text",
				   "Number of elements in a file",
				   "File containing the output data"};

     Current_option= (struct option *) &long_opts[0];

     printf("Usage: %s [-h] -f <Factor> -i <InputFile> -l <TotalLength> -o <OutputFile>\n", argv[0] );
     printf("\n");

     printf("Summary of arguments\n");
     i = 0;
     while( (Current_option != NULL) && (Current_option->name != NULL ) ){
	  printf("  -%c  --%s      \t %s.\n", (char) Current_option->val, Current_option->name, Explanation[i] );
	  Current_option = Current_option + 1;
	  i = i + 1;
     }
}
