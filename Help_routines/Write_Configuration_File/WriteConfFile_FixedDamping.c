#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>  /* For va_arg( ), ... */  
#include <string.h>
#include <getopt.h>

void print_help( char **argv, const struct option *long_opts );
void Calculate_Rayleigh( const double *const Freq, const double *const Damp, double *const Rayleigh );
char* concat(int count, ...);

int main( int argc, char **argv )
{

     FILE *TheFile;

     double Freq[2];       /* Frequencies */
     double Damp[2];       /* Damping ratios */
     double Rayleigh[2];   /* Rayleigh: Mass and stiffness proportional coefficients */

     char temp[2];
     char *MMatrix; char *KMatrix;
     char *OutputFile; char *ConfFile;
     int rc, idx, i;

     double factor[3] = {0.02, 0.05, 0.1};

     struct option long_options[] = 
	  {
	       {"conf-file", required_argument, 0, 'c'},
	       {"mass-matrix", required_argument, 0, 'm'},
	       {"stiffness-matrix", required_argument, 0, 'k'},
	       {"first-frequency", required_argument, 0, 'F'},
	       {"second-frequency", required_argument, 0, 'f'},
	       {"first-damping", required_argument, 0, 'D'},
	       {"second-damping", required_argument, 0, 'd'},
	       {"OutFile", required_argument, 0, 'o'},
	       {"help", no_argument, 0, 'h'},
	       {0, 0, 0, 0}
	  };

     if( argc <= 1 || argc != 17 ){
	  fprintf( stderr, "Invalid number of arguments.\n" );
	  print_help( argv, long_options );
	  return( EXIT_FAILURE );
     }	  

     while( (rc = getopt_long( argc, argv, "c:m:k:F:f:D:d:o:h", long_options, &idx )) != -1 ){

	  switch( rc ){
	  case 'c':
	       ConfFile = strdup( optarg );
	       break;
	  case 'm':
	       MMatrix = strdup( optarg );
	       break;
	  case 'k':
	       KMatrix = strdup( optarg );
	       break;
	  case 'F':
	       Freq[0] = atof( optarg );
	       break;
	  case 'f':
	       Freq[1] = atof( optarg );
	       break;
	  case 'D':
	       Damp[0] = atof( optarg );
	       break;
	  case 'd':
	       Damp[1] = atof( optarg );
	       break;
	  case 'o':
	       OutputFile = strdup( optarg );
	       break;
	  case 'h':
	       print_help( argv, long_options );
	       return( EXIT_FAILURE );
	       break;
	  case '?':
	       /* Long options already prints an error message telling that there is an unrecognised option */
	       break;
	  case ':':
	       /* Long options already prints an error message telling that the option requieres an argument */
	       break;
	  }
     }

     i = 0;
     Rayleigh[0] = 0.0; Rayleigh[1] = 0.0;


     while( i < 3 ){
	  Damp[0] = factor[i];
	  Damp[1] = factor[i];
	  Calculate_Rayleigh( Freq, Damp, Rayleigh );
	  sprintf( temp, "%.2lf", factor[i] );

	  TheFile = fopen( concat( 4, ConfFile, "_Damping" , temp , ".conf" ), "w" );
	  if( TheFile == NULL ){
	       fprintf( stderr, "Could not open file %s.\n", ConfFile );
	       fprintf( stderr, "Exiting.\n" );
	       exit( EXIT_FAILURE );
	  }
     
	  /* Print the contents of the configuration file */
	  fprintf( TheFile, "[General]\n" );
	  fprintf( TheFile, "Use_Absolute_Values = 0\n" );
	  fprintf( TheFile, "Read_Sparse = 1\n" );
	  fprintf( TheFile, "Use_Sparse = 1\n" );
	  fprintf( TheFile, "Use_Packed = 0\n" );
	  fprintf( TheFile, "Order = 960				# Number of Degrees of freedom. Order of the matrices.\n" );
	  fprintf( TheFile, "Read_LVector = 1\n" );
	  fprintf( TheFile, "Excited_DOF = \"6 1 1 1 0 0 0\"		# Excited degrees of Freedom x y z Rx Ry Rz. Only used if Read_LVector = 0\n" );
	  fprintf( TheFile, "Read_CMatrix = 0\n" );
	  fprintf( TheFile, "Num_Steps = 4096			# Number of steps.\n" );
	  fprintf( TheFile, "Delta = 0.01				# Time increment between two consecutive steps.\n" );
	  fprintf( TheFile, "Scale_Factor = 0.001\n" );
	  fprintf( TheFile, "\n" );

	  fprintf( TheFile, "[Grid]\n" );
	  fprintf( TheFile, "Rows = 2\n" );
	  fprintf( TheFile, "Cols = 2\n" );
	  fprintf( TheFile, "\n" );

	  fprintf( TheFile, "[Block_Size]\n" );
	  fprintf( TheFile, "Rows = 2\n" );
	  fprintf( TheFile, "Cols = 2\n" );
	  fprintf( TheFile, "\n" );

	  fprintf( TheFile, "# This section defines the values of the Rayleigh damping.\n" );
	  fprintf( TheFile, "[Rayleigh]\n" );
	  fprintf( TheFile, "Alpha = %.6lE\n", Rayleigh[0] );
	  fprintf( TheFile, "Beta = %.6lE\n", Rayleigh[1] );
	  fprintf( TheFile, "\n" );

	  fprintf( TheFile, "# Newmark Alpha and Beta values\n" );
	  fprintf( TheFile, "[Newmark]\n" );
	  fprintf( TheFile, "Gamma = 0.5\n" );
	  fprintf( TheFile, "Beta = 0.25\n" );
	  fprintf( TheFile, "\n" );

	  fprintf( TheFile, "# Values for the PID parameters. Currently only P is implemented in the algorithm\n" );
	  fprintf( TheFile, "[PID]\n" );
	  fprintf( TheFile, "P = 0.95\n" );
	  fprintf( TheFile, "I = 0.0\n" );
	  fprintf( TheFile, "D = 0.0\n" );
	  fprintf( TheFile, "\n" );

	  fprintf( TheFile, "# File names for the different matrices required by the algorith.\n" );
	  fprintf( TheFile, "[FileNames]\n" );
	  fprintf( TheFile, "Mass_Matrix = \"%s\"\n", MMatrix );
	  fprintf( TheFile, "Stiffness_Matrix = \"%s\"\n", KMatrix );
	  fprintf( TheFile, "Damping_Matrix = \"None\"\n" );
	  fprintf( TheFile, "Load_Vector = \"LV960_Test_MM.txt\"\n" );
	  fprintf( TheFile, "Coupling_Nodes = \"Couple_Nodes.txt\"	# File containing the coupling nodes.\n" );
	  fprintf( TheFile, "Ground_Motion = \"GroundMovement_Sinesweep.txt\"	# Ground movement.\n" );
	  fprintf( TheFile, "OutputFile = \"%s_Damping%.2lf\"			# Output file\n", OutputFile, factor[i] );
	  fprintf( TheFile, "\n" );

	  fprintf( TheFile, "[Substructure]\n" );
	  fprintf( TheFile, "Num_Substeps = 4			# Number of sub-steps.\n" );
	  fprintf( TheFile, "Order = 1\n" );

	  i = i + 1;
     }

     fclose( TheFile );
     free( MMatrix );
     free( KMatrix );
     free( ConfFile );
     free( OutputFile );	  

     return 0;
}

void Calculate_Rayleigh( const double *const Freq, const double *const Damp, double *const Rayleigh )
{

     double omega[2];
     int i;

     for ( i = 0; i < 2; i++ ){
	  omega[i] = 2.0*3.141592*Freq[i];
     }

     Rayleigh[1] = (2.0*omega[0]*Damp[0] - 2.0*omega[1]*Damp[1])/(omega[0]*omega[0] - omega[1]*omega[1]);
     Rayleigh[0] = 2.0*omega[0]*Damp[0] - Rayleigh[1]*omega[0]*omega[0];

     if( Rayleigh[0] < 0.0 || Rayleigh[1] < 0.0 ){
	  fprintf( stderr, "Negative rayleigh damping coefficients for frequencies (%lf,%lf) and damping factors (%.6lE, %.6lE).\n", Freq[0], Freq[1], Damp[0], Damp[1] );
	  exit( EXIT_FAILURE );
     }
}

char* concat(int count, ...)
{
     va_list ap;
     int i, null_pos;
     char *merged, *s;

     /* Find required length to store merged string */
     int len = 1; /* Space for NULL */
     va_start( ap, count );

     for( i = 0; i < count; i++ )
	  len += strlen( va_arg( ap, char* ) );
     va_end(ap);

     /* Allocate memory to concat strings */
     merged = calloc( sizeof(char) ,len );
     null_pos = 0;

     /* Actually concatenate strings */
     va_start( ap, count );

     for( i = 0; i < count; i++ )
     {
	  s = va_arg( ap, char* );
	  strcpy( merged+null_pos, s );
	  null_pos += strlen( s );
     }
     va_end( ap );

     return merged;
}

void print_help( char **argv, const struct option *long_opts )
{
     int i;
     struct option *Current_option = NULL;
     const char *Explanation[9] = {"Name of the configuration file (without extension)",
				   "Mass matrix file",
				   "Stiffness matrix file",
				   "First frequency for Rayleigh coefficients",
				   "Second frequency for Rayleigh coefficients",
				   "Damping ratio for the first frequency",
				   "Damping ratio for the second frequency",
				   "Name of the output file (without extension)",
				   "This help text"};

     Current_option= (struct option *) &long_opts[0];

     printf("Usage: %s [-h] -c <ConfFile> -m <MassMatrix> -k <StiffnessMatrix> -F <FirstFreq.> -f <SecondFreq.> -D <DampFirstFreq.> -d <DampSecondFreq.> -o <OutputFile>", argv[0] );
     printf("\n");
    

     printf("Summary of arguments\n");
     i = 0;
     while( (Current_option != NULL) && (Current_option->name != NULL ) ){
	  printf("  -%c  --%s    \t %s.\n", (char) Current_option->val, Current_option->name, Explanation[i] );
	  Current_option = Current_option + 1;
	  i = i + 1;
     }
}
