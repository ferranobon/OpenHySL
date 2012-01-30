#include <stdio.h>
#include <stdlib.h>       /* For exit( ) and atoi( ) */

#include <string.h>       /* For memset( ) */
#include <sys/socket.h>   /* For socket( ), bind( ) and connect( ) */
#include <arpa/inet.h>    /* For sockaddr_in and inet_ntoa( ) */
#include <unistd.h>       /* For close( ) */

#include "RoutinesADwin.h"
#include "Substructure.h" /* For NUM_CHANNELS, ConstSub and Init_Constants_Substructure( ) */
#include "OPSocket.h"

int main ( int argc, char **argv )
{
     int i;
     /* Variables concerning data communication */
     int Server_Socket;      /* Socket for the server */
     int Client_Socket;      /* Socket for the client */
     unsigned int Port;      /* Port */

     int DataTypeSize;       /* Size of the data type to be transfered */
     char *Data;             /* To Send and receive data */
     int Length;             /* Length of the data to be transfered */
     int ierr;               /* Error. OpenFresco routines */
     int Is_Not_Finished;    /* To check if the process is finished or not */

     /* Variables required by OpenFresco */
     int iData[11];

     ConstSub Cnst;

     float *Gc;
     float *u0c, *uc;
     float *fcprev, *fc;
     float *Send, *Recv;

     /* Array where the data from ADwin will be stored */
     float *ADWIN_DATA;

     /* Test the correct number of arguments */
     if (argc != 2){
	  fprintf( stderr, "Usage: %s <Server Port>\n", argv[0] );
	  exit( EXIT_FAILURE );
     }

     /* Store the port to be used */
     Port = (unsigned int) atoi( argv[1] );

     /* Create a TCP/IP socket for the server using OpenFresco routine */
     setupconnectionserver( &Port, &Server_Socket );
     if ( Server_Socket >= 0 ){
	  printf("Server connection successfully configured\n" );
     } else {
	  exit( EXIT_FAILURE );
     }

     /* Receive an ID from OpenFresco with sizeCtrl, sizeDaq and dataSize */
     Data = (char *) iData;
     DataTypeSize = sizeof( int );
     Length = 11;
     recvdata( &Server_Socket, &DataTypeSize, Data, &Length, &ierr );

     /* TODO: Implement receive data as a file */

     /* Initialise the constants of the substructure */
     Init_Constants_Substructure( &Cnst );

     Gc = calloc( Cnst.Order_Couple*Cnst.Order_Couple, sizeof( float ) );

     u0c = calloc( Cnst.Order_Couple, sizeof( float ) );
     uc = calloc( Cnst.Order_Couple, sizeof( float ) );

     fcprev = calloc( Cnst.Order_Couple, sizeof( float ) );
     fc = calloc( Cnst.Order_Couple, sizeof( float ) );

     /* Array where the data from ADwin will be stored */
     ADWIN_DATA = calloc( Cnst.Num_Sub*Cnst.Num_Steps*NUM_CHANNELS, sizeof( float ) );

     /* The size of the data to be exchanged is given by the last element of iData */
     Length = iData[10];
     Send = calloc( Length, sizeof(float) );
     Recv = calloc( Length, sizeof(float) );

     /* The process is not finished. */
     Is_Not_Finished = 1;  
     while( Is_Not_Finished ){

	  /* Receive the control values. */
	  Data = (char *) Recv;
	  recvdata( &Server_Socket, &DataTypeSize, Data, &Length, &ierr );
	  
	  /* Check if what has been received are the control values */
	  if ( Recv[0] == 3.0 ){
	       for ( i = 0; i < Cnst.Order_Couple; i++ ){
		    u0c[i] = Recv[1+i];
	       }

	       /* Perform the substepping process in ADwin */
	       ADWIN_Substep( u0c, uc, fcprev, fc, Cnst.Order_Couple, Cnst.Num_Sub, Cnst.DeltaT_Sub );
	  }

	  /* Check if the Daq values are being asked */
	  if ( Recv[0] == 6.0 ){
	       Send[0] = 10.0;  /* What is the first value ????? */

	       /* Compose the data */
	       for ( i = 0; i < Cnst.Order_Couple; i++ ){
		    Send[i + 1] = uc[i];
		    Send[i + Cnst.Order_Couple + 1] = fcprev[i];
		    Send[i + 2*Cnst.Order_Couple + 1] = fc[i];
	       }

	       /* Send the requested data */
	       Data = (char *) Send;
	       senddata( &Server_Socket, &DataTypeSize, Data, &Length, &ierr );
	  }

	  /* Check if the process is finished */
	  if ( Recv[0] == 99.0 ){
	       Is_Not_Finished = 0;
	  }

     }

     /* End the connection */
     closeconnection( &Server_Socket, &ierr );

     /* Get the Data from ADwin */
     GetDataADwin( Cnst.Num_Steps, Cnst.Num_Sub, ADWIN_DATA );

     /* Free the dinamically allocated memory */
     free( ADWIN_DATA );

     free( Gc );

     free( u0c );
     free( uc );

     free( fcprev );
     free( fc );

     return 0;
}
