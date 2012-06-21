#include <stdio.h>
#include <stdlib.h>      /* For atoi( ) and exit() */

#include <string.h>      /* For memset( ) */
#include <sys/socket.h>  /* For socket( ), bind( ) and connect( ) */
#include <arpa/inet.h>   /* For sockaddr_in and inet_ntoa( ) */
#include <unistd.h>      /* For close( ) */

#include <getopt.h>      /* For getopt_long() */

#include "ErrorHandling.h"
#include "Substructure.h"
#include "Send_Receive_Data.h"
#include "RoutinesADwin.h"
#include "Custom_Server.h"

#define MAXPENDING   5    /* Maximum outstanding connection requests */

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

     /* Array where the data from ADwin will be stored */
     float *ADWIN_DATA;

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
	  printf("Defaulting on mode 1: ");
     }

     /* Assign each argument to the correct variable */
     while( (Selected_Option = getopt_long( argc, argv, "m:p:h", long_options, NULL )) != -1 ){
	  switch( Selected_Option ){
	  case 'm':
	       Mode = atoi( optarg );

	       if ( Mode < 0 || Mode > 2 ){
		    fprintf( stderr, "Mode %d is not a valid mode value.\n", Mode );
		    Print_Help( argv[0] );
		    return EXIT_FAILURE;
	       } else {
		    printf( "Using mode %d: ", Mode );
		    break;
	       }
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

     if ( Mode == USE_ADWIN ){
	  /* Run with ADwin */
	  printf( "using ADwin to perform the sub-stepping process.\n" );
	  ADWIN_DATA = calloc( Cnst.Num_Sub*Cnst.Num_Steps*NUM_CHANNELS, sizeof( float ) );
     } else if ( Mode == USE_EXACT ){
	  /* Simulate the substructure numerically */
	  printf( "simulating the sub-structure using an exact integration method.\n");
	  ExactSolution_Init( 285, 352.18177, 68000, Cnst.DeltaT_Sub, &Num_TMD );
     } else {
	  printf( "simulating the sub-structure using measured values as an input.\n");
	  /* Do nothing for the moment */
     }


     /* Create a TCP/IP socket for the server */
     Server_Socket = Init_TCP_Server_Socket( Port );
     printf( "TCP Server created. Listening on port %d for incoming client connections...\n", Port );

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
     Receive_Data( Gc, Length, Client_Socket );

     if ( Mode == 0 ){
	  ADWIN_SetGc( Gc, Cnst.Order_Couple*Cnst.Order_Couple );
     } else {
	  printf("Simulating the substructure\n");
     }

     Is_Not_Finished = 1;
     while( Is_Not_Finished ){

	  /* Receive the displacement */
	  Length = Cnst.Order_Couple;
	  Receive_Data( u0c, Length, Client_Socket );

	  if ( u0c[0] == -9999.0 ){
	       Is_Not_Finished = 0;
	  } else {
	       /* Perform the substepping process */

	       if ( Mode == USE_ADWIN ){
		    /* Run using ADwin */
		    ADWIN_Substep( u0c, uc, fcprev, fc, Cnst.Order_Couple, Cnst.Num_Sub, Cnst.DeltaT_Sub );
	       } else if ( Mode == USE_EXACT ){
		    /* Run without ADwin and simulating the substructure using an exact
		     * solution.
		     */
		    Simulate_Substructure( &Num_TMD, Gc, u0c, uc, fcprev, fc, Cnst.Order_Couple, Cnst.Num_Sub, Cnst.DeltaT_Sub );
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
	       Send_Data( Send, Length, Client_Socket );
	  }
     }

     /* Close the connection with the Client */
     close( Client_Socket );

     if ( Mode == USE_ADWIN ){
	  /* Get the Data from ADwin */
	  printf("Getting the data from ADwin...");
	  GetDataADwin( Cnst.Num_Steps, Cnst.Num_Sub, ADWIN_DATA );
	  printf(" DONE!\n");
     
	  free( ADWIN_DATA );
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

int Init_TCP_Server_Socket( const unsigned short int Server_Port )
{
     
     int Socket;
     struct sockaddr_in Server_Addr;

     /* Create socket for incoming connections */
     if ( (Socket = socket( PF_INET, SOCK_STREAM, IPPROTO_TCP ) ) < 0){
	  PrintErrorAndExit( "socket() failed." );
	  exit( EXIT_FAILURE );
     }
      
     /* Construct local address structure */
     memset( &Server_Addr, 0, sizeof(Server_Addr) );    /* Zero out structure */
     Server_Addr.sin_family = AF_INET;                  /* Internet address family */
     Server_Addr.sin_addr.s_addr = htonl( INADDR_ANY ); /* Any incoming interface */
     Server_Addr.sin_port = htons( Server_Port);        /* Local port */

     /* Bind to the local address */
     if ( bind( Socket, (struct sockaddr *) &Server_Addr, sizeof(Server_Addr)) < 0){
	  PrintErrorAndExit( "bind() failed." );
     }

     /* Mark the socket so it will listen for incoming connections */
     if ( listen( Socket, MAXPENDING) < 0){
	  PrintErrorAndExit( "listen() failed." );
     }

     return Socket;

}

int Accept_TCP_Client_Connection( int Server_Socket )
{
     
     int Client_Socket;                 /* Socket descriptor for client */
     struct sockaddr_in Client_Addr;    /* Client address */
     unsigned int Client_Length;        /* Length of client address data structure */

     /* Set the size of the in-out parameter */
     Client_Length = sizeof( Client_Addr );
    
     /* Wait for a client to connect */
     if ( (Client_Socket = accept( Server_Socket, (struct sockaddr *) &Client_Addr, &Client_Length)) < 0 ){
	  PrintErrorAndExit( "accept() failed.\n" );
     }
    
     /* The Client is connected. Display a proper message */
     printf("Handling client %s\n", inet_ntoa( Client_Addr.sin_addr ) );

     return Client_Socket;
}

void Print_Help( const char *Program_Name )
{

     fprintf( stderr, "Usage: %s [-h] -m <Mode> -p <Port>, argv[0]" );
     fprintf( stderr,
	      "  -h  --help    This help text.\n"
	      "  -m  --mode    The mode used by the program.\n"
	      "                  0 - ADwin will be used to perform the sub-stepping process.\n"
	      "                  1 - The Substructure will be simulated using an exact solution.\n"
	      "                  2 - The Substructure will be simulated using measured values.\n"
	      "  -p  --port    Port used for communication with the client program.\n" );
}
