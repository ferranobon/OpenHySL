#include <stdio.h>
#include <stdlib.h>      /* For atoi( ) and exit() */
#include <string.h>      /* For memset( ) */
#include <getopt.h>      /* For getopt_long() */

#include "ErrorHandling.h"
#include "Substructure.h"
#include "Send_Receive_Data.h"
#include "Custom_Server.h"

#if ADWIN_
#include "RoutinesADwin.h"
#endif

int main( int argc, char **argv )
{
     int i;

     /* Variables concerning data communication */
     int Server_Socket;                /* Socket for the server */
     int Client_Socket;                /* Socket for the client */
     int Port;

     /* Variables for the sub-structure testing/simulation */
     int Is_Not_Finished;
     int Length;

     ConstSub Cnst;

     float *Gc;
     float *u0c, *uc;
     float *fcprev, *fc;
     float *Send;

     TMD_Sim Num_TMD;
     UHYDE_Sim Num_UHYDE;

#if ADWIN_
     /* Array where the data from ADwin will be stored */
     float *ADWIN_DATA;
#endif

     /* Variables to deal with arguments */
     int Mode, Selected_Option;

     struct option long_options[] = {
	  {"help", no_argument, 0, 'h'},
	  {"mode", required_argument, 0, 'm'},
	  {"server-port", required_argument, 0, 'p'},
	  {0, 0, 0, 0}
     };

     /* Set the default value for Port and Mode before the user input. */
     Port = 3333;  /* Default port */
     Mode = 1;     /* Simulate the sub-structure using an exact solution. */

    /* This is only used if there are no arguments */
     if ( argc == 1 ){	  
	  printf("Defaulting on mode 1.\n");
     }

     /* Assign each argument to the correct variable */
     while( (Selected_Option = getopt_long( argc, argv, "m:p:h", long_options, NULL )) != -1 ){
	  switch( Selected_Option ){
	  case 'm':
	       Mode = atoi( optarg );

	       if ( Mode < 0 || Mode > 3 ){
		    fprintf( stderr, "Mode %d is not a valid mode value.\n", Mode );
		    Print_Help( argv[0] );
		    return EXIT_FAILURE;
	       }
	       break;
	  case 'p':
	       Port = atoi( optarg );
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

     /* Create a TCP/IP socket for the server */
     Server_Socket = Init_TCP_Server_Socket( Port );
     printf( "TCP Server created.\nListening on port %d for incoming client connections...\n", Port );

     /* Accept the connection from the client */
     Client_Socket = Accept_TCP_Client_Connection( Server_Socket );

     /* Initialise the constants of the substructure */
     Init_Constants_Substructure( &Cnst );

     /* Dynamically allocate memory */
     Gc = calloc( Cnst.Order_Couple*Cnst.Order_Couple, sizeof( float ) );
 
     u0c = calloc( Cnst.Order_Couple, sizeof( float ) );
     uc = calloc( Cnst.Order_Couple, sizeof( float ) );

     fcprev = calloc( Cnst.Order_Couple, sizeof( float ) );
     fc = calloc( Cnst.Order_Couple, sizeof( float ) );

     Send = calloc( 3*Cnst.Order_Couple, sizeof( float ) );

     /* Receive matrix Gc */
     Length = Cnst.Order_Couple*Cnst.Order_Couple;
     Receive_Data_TCP( Gc, Length, Client_Socket );

     if ( Mode == USE_ADWIN ){
#if ADWIN_
	  /* Run with ADwin */
	  ADWIN_SetGc( Gc, Cnst.Order_Couple*Cnst.Order_Couple );
	  printf( "Using ADwin to perform the sub-stepping process.\n" );
	  ADWIN_DATA = calloc( Cnst.Num_Sub*Cnst.Num_Steps*NUM_CHANNELS, sizeof( float ) );
#else
	  fprintf(stderr, "The program was not compiled with ADwin support.\n");
	  exit( EXIT_FAILURE );
#endif
     } else if ( Mode == USE_EXACT ){
	  /* Simulate the substructure numerically */
	  printf( "Simulating the sub-structure using an exact integration method.\n");
	  ExactSolution_Init( 285, 352.18177, 68000, Cnst.DeltaT_Sub, &Num_TMD );
     } else if ( Mode == USE_UHYDE ){
	  printf( "Simulating the friction device UHYDE-fbr.\n" );
	  Simulate_UHYDE_1D_Init( 0.0002, 0.9, 500.0, &Num_UHYDE );
     } else {
	  printf( "Simulating the sub-structure using measured values as an input.\n");
	  /* Do nothing for the moment */
     }

     Is_Not_Finished = 1;
     while( Is_Not_Finished ){

	  /* Receive the displacement */
	  Length = Cnst.Order_Couple;
	  Receive_Data_TCP( u0c, Length, Client_Socket );

	  if ( u0c[0] == -9999.0 ){
	       Is_Not_Finished = 0;
	  } else {
	       /* Perform the substepping process */
	       if ( Mode == USE_ADWIN ){
#if ADWIN_
		    /* Run using ADwin */
		    ADWIN_Substep( u0c, uc, fcprev, fc, Cnst.Order_Couple, Cnst.Num_Sub, Cnst.DeltaT_Sub );
#endif
	       } else if ( Mode == USE_EXACT ){
		    /* Run without ADwin and simulating the substructure using an exact
		     * solution.
		     */
		    Simulate_Substructure( &Num_TMD, Mode, Gc, u0c, uc, fcprev, fc, Cnst.Order_Couple, Cnst.Num_Sub, Cnst.DeltaT_Sub );
	       } else if ( Mode == USE_UHYDE ){
		    /*
		     * Simulate the UHYDE-fbr without ADwin
		     */
		    Simulate_Substructure( &Num_UHYDE, Mode, Gc, u0c, uc, fcprev, fc, Cnst.Order_Couple, Cnst.Num_Sub, Cnst.DeltaT_Sub );
	       } else {
		    /* Run without ADwin and simulating the substructure using measured
		     * values of the coupling force.
		     */
		    Simulate_Substructure_Measured_Values( "fc.txt", Gc, u0c, uc, fcprev, fc, Cnst.Order_Couple, Cnst.Num_Sub, Cnst.DeltaT_Sub );
	       }

	       /* Compose the data to send */
	       for (i = 0; i < Cnst.Order_Couple; i++) {
		    Send[i] = uc[i];
		    Send[i+Cnst.Order_Couple] = fcprev[i];
		    Send[i+2*Cnst.Order_Couple] = fc[i];
	       }

	       /* Send the response */
	       Length = 3*Cnst.Order_Couple;
	       Send_Data_TCP( Send, Length, Client_Socket );
	  }
     }

     /* Close the connection with the Client */
     Close_Socket( &Client_Socket );

     if ( Mode == USE_ADWIN ){
#if ADWIN_
	  /* Get the Data from ADwin */
	  printf("Getting the data from ADwin...");
	  GetDataADwin( Cnst.Num_Steps, Cnst.Num_Sub, ADWIN_DATA );
	  printf(" DONE!\n");
     
	  free( ADWIN_DATA );
#endif
     } else {
	  printf("The simulatiovn has finished\n");
     }
     
     /* Free the dinamically allocated memory */
     free( Gc );

     free( u0c );
     free( uc );

     free( fcprev );
     free( fc );

     free( Send );

     return 0;
}



void Print_Help( const char *Program_Name )
{

     fprintf( stderr, "Usage: %s [-h] -m <Mode> -p <Port>", Program_Name );
     fprintf( stderr,
	      "  -h  --help    This help text.\n"
	      "  -m  --mode    The mode used by the program. Default value 1.\n"
	      "                  0 - ADwin will be used to perform the sub-stepping process.\n"
	      "                  1 - The Substructure will be simulated using an exact solution.\n"
	      "                  2 - The Substructure will be simulated using measured values.\n"
	      "  -p  --port    Port used for communication with the client program. Default value 3333.\n" );
}
