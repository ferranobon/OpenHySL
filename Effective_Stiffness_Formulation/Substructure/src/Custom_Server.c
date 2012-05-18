#include <stdio.h>
#include <stdlib.h>      /* For atoi( ) and exit() */

#include <string.h>      /* For memset( ) */
#include <sys/socket.h>  /* For socket( ), bind( ) and connect( ) */
#include <arpa/inet.h>   /* For sockaddr_in and inet_ntoa( ) */
#include <unistd.h>      /* For close( ) */


#include "ErrorHandling.h"
#include "Substructure.h"
#include "Send_Receive_Data.h"
#include "RoutinesADwin.h"
#include "Custom_Server.h"

#define MAXPENDING 5    /* Maximum outstanding connection requests */

int main( int argc, char **argv )
{
     int i;
     /* Variables concerning data communication */
     int Server_Socket;                /* Socket for the server */
     int Client_Socket;                /* Socket for the client */

     int Is_Not_Finished;
     int Length;

     ConstSub Cnst;

     float *Gc;
     float *u0c, *uc;
     float *fcprev, *fc;
     float *Send;

#if SIMULATE_SUB_
     TMD_Sim Num_TMD;
#endif

     /* Array where the data from ADwin will be stored */
     float *ADWIN_DATA;

     /* Test the correct number of arguments */
     if (argc != 2){
	  fprintf( stderr, "Usage: %s <Server Port>\n", argv[0] );
	  exit( EXIT_FAILURE );
     }

     /* Create a TCP/IP socket for the server */
     Server_Socket = Init_TCP_Server_Socket( atoi( argv[1] ) );

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

#if SIMULATE_SUB_
     /* Do nothing */
     ExactSolution_Init( 285, 352.18177, 68000, Cnst.DeltaT_Sub, &Num_TMD );
#else
     ADWIN_DATA = calloc( Cnst.Num_Sub*Cnst.Num_Steps*NUM_CHANNELS, sizeof( float ) );
#endif

     Send = calloc( 3*Cnst.Order_Couple, sizeof( float ) );

     /* Receive matrix Gc */
     Length = Cnst.Order_Couple*Cnst.Order_Couple;
     Receive_Data( Gc, Length, Client_Socket );

#if SIMULATE_SUB_  /* Run this without ADwin */
     /* Do nothing */
     printf("Simulating the substructure\n");
#else
     ADWIN_SetGc( Gc, Cnst.Order_Couple*Cnst.Order_Couple );
#endif

     Is_Not_Finished = 1;
     while( Is_Not_Finished ){

	  /* Receive the displacement */
	  Length = Cnst.Order_Couple;
	  Receive_Data( u0c, Length, Client_Socket );

	  if ( u0c[0] == -9999.0 ){
	       Is_Not_Finished = 0;
	  } else {
	       /* Perform the substepping process */

#if SIMULATE_SUB_  /* Run this without ADwin */
	       Simulate_Substructure( &Num_TMD, Gc, u0c, uc, fcprev, fc, Cnst.Order_Couple, Cnst.Num_Sub, Cnst.DeltaT_Sub );
#else              /* Run using ADwin */
	       ADWIN_Substep( u0c, uc, fcprev, fc, Cnst.Order_Couple, Cnst.Num_Sub, Cnst.DeltaT_Sub );
#endif	  
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

#if SIMULATE_SUB_  /* Run this without ADwin */
     printf("The simulatiovn has finished\n");
#else
     /* Get the Data from ADwin */
     printf("Getting the data from ADwin...");
     GetDataADwin( Cnst.Num_Steps, Cnst.Num_Sub, ADWIN_DATA );
     printf(" DONE!\n");
     
     free( ADWIN_DATA );
#endif
     
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
