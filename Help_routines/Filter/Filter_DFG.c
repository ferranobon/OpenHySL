#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>

void print_help( char **argv, const struct option *long_opts );

int main( int argc, char **argv )
{

     int i;             /* A counter */
     int Max_Elements;  /* Maximum number of elements in the file */
     int Range_Max;     /* Limit of the filter */
     int idx, rc;
     int help_flag = 0;

     double InValue, OutValue; /* Input and Output value respectively */
     double Factor;     /* Factor applied to the input values */

     FILE *InFile, *OutFile;   /* Input and Output file respectively */

     struct option long_options[] = 
     {
	  {"input-file", required_argument, 0, 'i'},
	  {"output-file", required_argument, 0, 'o'},
	  {"length", required_argument, 0, 'l'},
	  {"range", required_argument, 0, 'r'},
	  {"help", no_argument, &help_flag, 'h'}
     };

     if( argc <= 1 ){  
	  fprintf( stderr, "Invalid number of arguments.\n" );
	  print_help( argv, long_options );
	  return( EXIT_FAILURE );
     }

     while( (rc = getopt_long( argc, argv, "i:o:l:r:h", long_options, &idx )) != -1 ){

	  switch( rc ){
	  case 'i':
	       InFile = fopen( optarg, "r" );
	       break;
	  case 'o':
	       OutFile = fopen( optarg, "w" );
	       break;
	  case 'l':
	       Max_Elements= atoi( optarg );
	       break;
	  case 'r':
	       Range_Max = atoi( optarg );
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

     if( help_flag != 0 ){
	  print_help( argv, long_options );
	  return ( EXIT_FAILURE );
     }

     for ( i = 0; i < Max_Elements; i++ ){
	  
	  if ( i < Range_Max ){
	       Factor = sin(((double)i/(double)Range_Max)*(3.14159625/2.0)) * sin(((double)i/(double)Range_Max)*(3.14159625/2.0));
	  } else {
	       Factor = 1.0;
	  }

	  fscanf( InFile, "%lf", &InValue );
	  OutValue = InValue*Factor;

	  fprintf( OutFile, "%lf\n", OutValue );
     }

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
			     "Number of elements in a file",
			     "Range of elements to be processed since the first position",
			     "This help text"};

     Current_option= (struct option *) &long_opts[0];

     printf("Usage: %s -i <InputFile> -o <OutputFile> -l <TotalLength> -r <Range>\n", argv[0] );
     printf("   or: %s -h\n", argv[0] );
     printf("\n");
    

     printf("Summary of arguments\n");
     i = 0;
     while( (Current_option != NULL) && (Current_option->name != NULL ) ){
	  printf("  -%c  --%s    \t %s.\n", (char) Current_option->val, Current_option->name, Explanation[i] );
	  Current_option = Current_option + 1;
	  i = i + 1;
     }
}
