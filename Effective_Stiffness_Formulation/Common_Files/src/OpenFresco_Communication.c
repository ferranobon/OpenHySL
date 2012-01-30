#include <stdio.h>
#include <stdlib.h>    /* For calloc( ) and free( ) */
#include <string.h>    /* For strlen( ) */

#include "Send_Receive_Data.h"
#include "OPSocket.h"

/* Global Variables */
int SocketID;
int DataSize = 256;

int Communicate_With_OpenFresco( const float *const Data_To_Send, float *const Data_To_Receive, int Size, int WhatToDo )
{

     /* Local Variables */
     static int i, iData[11];

     static float *sData, *rData;

     int ierr, nleft, DataTypeSize;
     char *gMsg;

     if( WhatToDo == 1 ){
	  int Size_Machine_Inet;
	  Remote_Machine_Info RMachine;
	  unsigned int Port;

	  /* Allocate memory to send and receive vectors */
	  DataSize = ( 2*Size > DataSize ) ? 2*Size : DataSize;
	  DataSize = ( Size*Size > DataSize ) ? Size*Size : DataSize;
     
	  sData = (float *) calloc( DataSize, sizeof (float) );
	  rData = (float *) calloc( DataSize, sizeof (float) );

	  /* Setup the connection */
	  GetServerInformation( &RMachine );
	  Size_Machine_Inet = strlen( RMachine.IP );
	  Port = (unsigned int) RMachine.Port;
	  setupconnectionclient( &Port, RMachine.IP, &Size_Machine_Inet, &SocketID );


	  if ( SocketID < 0 ){
	       fprintf( stderr, "Cannot connect to the OpenFresco Server.\n" );
	       return -1;
	  }

	  /* Set the Data size for the experimental element */
	  /* SizeCtrl */
	  iData[0] = Size;   /* Displacement */
	  iData[1] = 0;      /* Velocity */
	  iData[2] = 0;      /* Acceleration */
	  iData[3] = 0;      /* Force */
	  iData[4] = 0;      /* Time */
	  /* sizeDaq */
	  iData[5] = Size;   /* Displacement */
	  iData[6] = 0;      /* Velocity */
	  iData[7] = 0;      /* Acceleration */
	  iData[8] = 2*Size; /* Force */
	  iData[9] = 0;      /* Time */
	  /* DataSize */
	  iData[10] = DataSize;

	  gMsg = ( char *) iData;
	  DataTypeSize = sizeof( int );
	  nleft = 11;
	  senddata( &SocketID, &DataTypeSize, gMsg, &nleft, &ierr );

     } else if ( WhatToDo == 2 ) {
	  /* Send Trial response to the experimental site */
	  sData[0] = 3;
	  for ( i = 0; i < Size; i++ ){
	       sData[i + 1] = Data_To_Send[i];
	  }
	  gMsg = (char *) sData;
	  DataTypeSize = sizeof(float);
	  nleft = DataSize;
	  senddata( &SocketID, &DataTypeSize, gMsg, &nleft, &ierr );

	  /* Ask for displacement. */
	  sData[0] = 7.0;

	  gMsg = (char *) sData;
	  nleft = DataSize;
	  senddata( &SocketID, &DataTypeSize, gMsg, &nleft, &ierr );

	  /* Receive the displacement from the experimental site */
	  gMsg = (char *) rData;
	  nleft = DataSize;
	  recvdata( &SocketID, &DataTypeSize, gMsg, &nleft, &ierr );

	  for ( i = 0; i < Size; i++ ){
	       Data_To_Receive[i]= rData[i];
	  }

	  /* Ask for force. */
	  sData[0] = 10.0;

	  gMsg = (char *) sData;
	  nleft = DataSize;
	  senddata( &SocketID, &DataTypeSize, gMsg, &nleft, &ierr );

	  /* Get the force from the experimental site. First the force from the previous
	   * substep and afterwards the force from the last substep
	   */
	  gMsg = (char *) rData;
	  nleft = DataSize;
	  recvdata( &SocketID, &DataTypeSize, gMsg, &nleft, &ierr );

	  for ( i = 0; i < 2*Size; i++ ){
	       Data_To_Receive[Size + i]= rData[i];
	  }
	  
     } else if ( WhatToDo == 10 ){ /* Disconnect from the experimental site */

	  sData[0] = 99.0;
	  gMsg = (char *) sData;
	  DataTypeSize = sizeof(float);
	  nleft = DataSize;
	  senddata( &SocketID, &DataTypeSize, gMsg, &nleft, &ierr );

	  /* Close the connection and release dinamically allocated memory */
	  closeconnection( &SocketID, &ierr );
	  free( sData );
	  free( rData );
     }

     return 0;
}
