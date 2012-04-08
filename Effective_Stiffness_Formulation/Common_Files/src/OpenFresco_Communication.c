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

     static double *sData, *rData;

     int ierr, nleft, DataTypeSize;
     char *gMsg;

     if( WhatToDo == 1 ){
	  int Size_Machine_Inet;
	  Remote_Machine_Info RMachine;
	  unsigned int Port;

	  /* Allocate memory to send and receive vectors */
	  DataSize = ( 2*Size > DataSize ) ? 2*Size : DataSize;
	  DataSize = ( Size*Size > DataSize ) ? Size*Size : DataSize;
     
	  sData = (double *) calloc( DataSize, sizeof (double) );
	  rData = (double *) calloc( DataSize, sizeof (double) );

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

     } else if ( WhatToDo == 3 ) {
	  /* Send Trial response to the experimental site */
	  sData[0] = 3;
	  for ( i = 0; i < Size; i++ ){
	       sData[i + 1] = (double) Data_To_Send[i];
	  }
	  gMsg = (char *) sData;
	  DataTypeSize = sizeof(double);
	  nleft = DataSize;
	  senddata( &SocketID, &DataTypeSize, gMsg, &nleft, &ierr );

	  /* Ask for DAQ values. */
	  sData[0] = 6.0;

	  gMsg = (char *) sData;
	  nleft = DataSize;
	  senddata( &SocketID, &DataTypeSize, gMsg, &nleft, &ierr );

	  /* Get the force from the experimental site. First the force from the previous
	   * substep and afterwards the force from the last substep
	   */
	  gMsg = (char *) rData;
	  nleft = DataSize;
	  recvdata( &SocketID, &DataTypeSize, gMsg, &nleft, &ierr );

	  for ( i = 0; i < 3*Size; i++ ){
	       Data_To_Receive[i]= (float) rData[i];
	       Data_To_Receive[i + Size]= (float) rData[i + Size];
	       Data_To_Receive[i + 2*Size]= (float) rData[i + 2*Size];
	  }
	  
     } else if ( WhatToDo == 10 ){ /* Disconnect from the experimental site */

	  sData[0] = 99.0;
	  gMsg = (char *) sData;
	  DataTypeSize = sizeof(double);
	  nleft = DataSize;
	  senddata( &SocketID, &DataTypeSize, gMsg, &nleft, &ierr );

	  /* Close the connection and release dinamically allocated memory */
	  closeconnection( &SocketID, &ierr );
	  free( sData );
	  free( rData );
     }

     return 0;
}
