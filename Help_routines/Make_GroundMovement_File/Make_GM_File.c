#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

void print_help( char **argv, const struct option *long_opts );

int main( int argc, char **argv )
{

     int i;  /* A counter */
     int Num_Elements;   /* Number of elements in the file */
     int idx, rc;
     int help_flag = 0;

     double Disp, Vel, Acc; /* Variables to store displacement, velocity
			     * and acceleration respectively
			     */

     FILE *Acc_File, *Vel_File, *Disp_File; /* Files containing the acceleration
					     * velocity and displacement data
					     * respectively 
					     */
     FILE *GM_File;      /* Output file of the format:
			  * #Row   Acceleration  Velocity  Displacement
			  */

     struct option long_options[] =
     {
	  {"acceleration-file", required_argument, 0, 'a'},
	  {"displacement-file", required_argument, 0, 'd'},
	  {"help", no_argument, &help_flag, 'h'},
	  {"length", required_argument, 0, 'l'},
	  {"velocity-file", required_argument, 0, 'v'}
     };

     if( argc <= 1 ){  
	  fprintf( stderr, "Invalid number of arguments.\n" );
	  print_help( argv, long_options );
	  return( EXIT_FAILURE );
     }

     while( (rc = getopt_long( argc, argv, "a:v:d:l:h", long_options, &idx )) != -1 ){

	  switch( rc ){
	  case 'a':
	       Acc_File = fopen( optarg, "r" );
	       if( Acc_File == NULL ){
		    fprintf( stderr, "Could not open file %s.\n", optarg );
		    fprintf( stderr, "Exiting.\n" );
		    exit( EXIT_FAILURE );
	       }
	       break;
	  case 'v':
	       Vel_File = fopen( optarg, "r" );
	       if( Vel_File == NULL ){
		    fprintf( stderr, "Could not open file: %s.\n", optarg );
		    fprintf( stderr, "Exiting.\n" );
		    exit( EXIT_FAILURE );
	       }
	       break;
	  case 'd':
	       Disp_File = fopen( optarg, "r" );
	       if( Disp_File == NULL ){
		    fprintf( stderr, "Could not open file: %s.\n", optarg );
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

     if( help_flag != 0 ){
	  print_help( argv, long_options );
	  return ( EXIT_FAILURE );
     }

     GM_File = fopen( "GroundMovement.txt", "w" );

     for ( i = 0; i < Num_Elements; i++ ){

	  fscanf( Acc_File, "%lf", &Acc );
	  fscanf( Vel_File, "%lf", &Vel );
	  fscanf( Disp_File, "%lf", &Disp );

	  fprintf( GM_File, "%i\t%lf\t%lf\t%lf\n", i+1, Acc, Vel, Disp );
     }

     fclose( Acc_File );
     fclose( Vel_File );
     fclose( Disp_File );
     fclose( GM_File );

     return 0;
}

void print_help( char **argv, const struct option *long_opts )
{
     int i;
     struct option *Current_option = NULL;
     const char *Explanation[5] = {"File containing acceleration data",
				   "File containing displacement data",
				   "This help text",
				   "Number of elements in a file",
				   "File containing velocity data"};

     Current_option= (struct option *) &long_opts[0];

     printf("Usage: %s [-h] -a <AccelerationData> -v <VelocityData> -d <DisplacementData> -l <TotalLength>\n", argv[0] );
     printf("\n");

     printf("Summary of arguments\n");
     i = 0;
     while( (Current_option != NULL) && (Current_option->name != NULL ) ){
	  printf("  -%c  --%s      \t %s.\n", (char) Current_option->val, Current_option->name, Explanation[i] );
	  Current_option = Current_option + 1;
	  i = i + 1;
     }
}
