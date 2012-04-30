/**
 * \file OpenFresco_Communication.c
 * \author Ferraon Obón Santacana
 * \version 1.0
 * \date 9th of April 2012
 *
 * \brief Source code of the interface between the sub-structure algorithm and
 * OpenFresco. This is in fact a \c GenericClient adaptation and must be run in
 * conjunction with \c ECGeneric experimental control.
 *
 * This file contains the source code of the communication routines used to
 * interface OpenFresco and the Dorka's sub-structure algorithm. This interface
 * is in fact an adaptation of the \c GenericClient of OpenFresco and, in order
 * to run successfully it has to be used together with the \c ECGeneric
 * experimental control of OpenFresco. It uses the NEES routines to handle the
 * TCP/IP communication. 
 */
#include <stdio.h>
#include <stdlib.h>            /* For calloc( ) and free( ) */
#include <string.h>            /* For strlen( ) */

#include "Send_Receive_Data.h" /* Common routines to handle TCP/IP or UDP
				* communication */
#include "OPSocket.h"          /* OpenFresco routines setupconnectionclient()
				* senddata(), recvdata() and closeconnection() */

/* Global variables */
int DataSize;
int SocketID;

int Communicate_With_OpenFresco( const float *const Data_To_Send, float *const Data_To_Receive, int Size, int WhatToDo )
{

     /* Local Variables */
     static int i, iData[11];

     /* OpenFresco uses type double to send/receive messages */
     static double *sData, *rData;

     /* Data communication variables */
     int ierr, nleft, DataTypeSize;
     char *gMsg;

     /* Setup the connection with the OpenFresco's SimAppServer */
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
	  setupconnectionclient( &Port, RMachine.IP, &Size_Machine_Inet,
				 &SocketID );

	  /* Check if the connection could be established */
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
	       /* ADwin can only handle float correctly. Therefore we must stick
		* with it and make the conversion from float to double so that
		* the messages can be sent correctly using OpenFresco. */
	       sData[i + 1] = (double) Data_To_Send[i];
	  }
	  gMsg = (char *) sData;
	  DataTypeSize = sizeof(double);
	  nleft = DataSize;
	  senddata( &SocketID, &DataTypeSize, gMsg, &nleft, &ierr );

	  /* Ask for DAQ values.*/
	  sData[0] = 6.0;
	  gMsg = (char *) sData;
	  nleft = DataSize;
	  senddata( &SocketID, &DataTypeSize, gMsg, &nleft, &ierr );

	  /* The DAQ values are received. The data received is bounded to a
	   * particular order in order to work properly: the new displacement
	   * in the first position, followed by the coupling force from the
	   * previous sub-step and the coupling force at the last sub-step. */
	  gMsg = (char *) rData;
	  nleft = DataSize;
	  recvdata( &SocketID, &DataTypeSize, gMsg, &nleft, &ierr );

	  for ( i = 0; i < 3*Size; i++ ){
	       Data_To_Receive[i]= (float) rData[i];
	       Data_To_Receive[i + Size]= (float) rData[i + Size];
	       Data_To_Receive[i + 2*Size]= (float) rData[i + 2*Size];
	  }
	  
     } else if ( WhatToDo == 10 ){
	  /* Disconnect from the experimental site */
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
