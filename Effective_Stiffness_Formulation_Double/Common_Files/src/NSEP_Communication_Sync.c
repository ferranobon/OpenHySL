/**
 * \file NSEP_Communication_Sync.c
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 4th of November 2011
 *
 * \brief Routines for the algorithm of Prof. Dorka to work together with NSEP protocol
 * 
 * This file contains the communication routines that are required in order to be compliant with
 * the NSEP protocol.
 *
 * \todo Implement the asynchronous version of the CGM to send/receive data from the server without a
 * request from the CGM.
 */

#include <stdio.h>       /* For printf() and fprintf() */
#include <stdlib.h>      /* For atoi() and exit( ) */
#include <string.h>      /* For memset() */

#if WIN32
#include <winsock2.h>
#else
#include <sys/types.h>
#include <sys/socket.h>
#endif

#include "NSEP_Communication_Sync.h"    /* Prototypes of several functions involved in the CGM */
#include "NSEP_Definitions.h"           /* Definition of various NSEP constants */
#include "Send_Receive_Data.h"          /* Send and receive data routines */
#include "ErrorHandling.h"              /* Error handling routines */
#include "Conf_Parser.h"

void Communicate_With_PNSE( const int WhatToDo, double Time,
			    const double *const Data_To_Send, double *Data_To_Receive, const unsigned int Order )
{
     static unsigned int Length;
     static int NSEP_Type;
     static int Socket;
     static int Error;  /* To help identifying the error source */
     static int What;
     static int ExpState;

     char Send_Buffer[MAXBUFLEN], Receive_Buffer[MAXBUFLEN], *pos;


     /* Connect and login to the server and set the state to running */
     if ( WhatToDo == 0 ){
	  Remote_Machine_Info Server;
	  ConfFile *Config;

	  /* Get the information regarding the PNSE server, Account Name and Password */
     	  Config = ConfFile_Create( 5 );
	  ConfFile_ReadFile( Config, "ConfFile.conf" );
	  GetNetworkInformation( &Server, Config );
	  /* Free the configuration File */
	  ConfFile_Free( Config );

	  /* Connect to the PNSE server */
	  Socket = Setup_Client_Socket( Server, PROTOCOL_TCP );
     
	  /* Login to the server */
	  Login_Server( Socket, Server, Send_Buffer, Receive_Buffer, &Error );
	  if ( Error == 1 ){
	       Close_Socket( &Socket );
	       PrintErrorAndExit( "CGM: send() sent a different number of bytes than expected. Exiting." );
	  } else if ( Error == 2 ){
	       Close_Socket( &Socket );
	       PrintErrorAndExit( "CGM: recv() failed or connection closed prematurely. Exiting.\n" );
	  } else if ( Error == 3 ){
	       Close_Socket( &Socket );
	       PrintErrorAndExit( "CGM: Exiting." );
	  }

          /* Tell the PNSE Server to change the state to Running.
	   * Run State = 2
	   */
	  Send_Client_State( Socket, Send_Buffer, &Error, NSEP_CS_RUNNING );
	  if ( Error == 1 ){
	       Close_Socket( &Socket );
	       PrintErrorAndExit( "CGM: send() sent a different number of bytes than expected. Exiting." );
	  }
	  ConfFile_Free( Config );
	  Delete_NetworkInformation( &Server );
	  
     } else if ( WhatToDo == 1 ) { /* Send the displacement to the server */

	  /* Compose the NSEP_CMD packet and send it to the server */
	  Length = 12 + 4*Order;
	  NSEP_Type = NSEP_CMD;

	  pos = Send_Buffer;
	  memcpy( pos, &Length, sizeof (int) );	

	  pos = pos + sizeof (int);
	  memcpy( pos , &NSEP_Type, sizeof (int) );

	  pos = pos + sizeof (int);
	  memcpy( pos, &Time, sizeof (double) );

	  pos = pos + sizeof (double);
	  memcpy( pos, Data_To_Send, sizeof(double)*Order );
	  if ( send( Socket, Send_Buffer, Length, 0 ) != (int) sizeof (char)*Length ){
	       Close_Socket( &Socket );
	       PrintErrorAndExit( "CGM: send() sent a different number of bytes than expected. Exiting." );
	  }
     } else if ( WhatToDo == 2 ) { /* Request a message from the server */
	  /* Send the NSEP_QUERY packet to query for the critical signal */
	  Length = 12;
	  NSEP_Type = NSEP_QUERY;
	  What = 9;

	  pos = Send_Buffer;
	  memcpy( pos, &Length, sizeof (int) );

	  pos = pos + sizeof (int);
	  memcpy( pos,& NSEP_Type, sizeof (int) );

	  pos = pos + sizeof (int);
	  memcpy( pos ,&What, sizeof (int) );
	  if ( send( Socket, Send_Buffer, Length, 0 ) != (int) sizeof (char)*Length ){
	       Close_Socket( &Socket );
	       PrintErrorAndExit( "CGM: send() sent a different number of bytes than expected. Exiting." );
	  }

	  /* Wait for the NSEP_CSIG to arrive */
	  printf("CGM: Waiting for NSEP_CSIG...\n" );
	  
	  if ( recv( Socket, Receive_Buffer, MAXBUFLEN, 0 ) <= 0 ){
	       Close_Socket( &Socket );
	       PrintErrorAndExit( "CGM: recv() failed or connection closed prematurely. Exiting.\n" );
	  } else {
	       /* Get the restoring force */
	       memcpy( Data_To_Receive, Receive_Buffer + 2*sizeof(int), sizeof (double)*Order );
	  }

     } else if ( WhatToDo == 3 ) { /* Request a message from the server (FCM) */
	  /* Send the NSEP_QUERY packet to query for the critical signal */
	  Length = 12;
	  NSEP_Type = NSEP_QUERY;
	  What = NSEP_CMD; /* NSEP_CMD = 8 */

	  pos = Send_Buffer;
	  memcpy( pos, &Length, sizeof (int) );

	  pos = pos + sizeof (int);
	  memcpy( pos,& NSEP_Type, sizeof (int) );

	  pos = pos + sizeof (int);
	  memcpy( pos ,&What, sizeof (int) );
	  if ( send( Socket, Send_Buffer, Length, 0 ) != (int) sizeof (char)*Length ){
	       Close_Socket( &Socket );
	       PrintErrorAndExit( "FCM: send() sent a different number of bytes than expected. Exiting." );
	  }

	  /* Wait for the NSEP_CMD to arrive */
	  printf("FCM: Waiting for NSEP_CMD...\n" );
	  
	  if ( recv( Socket, Receive_Buffer, MAXBUFLEN, 0 ) <= 0 ){
	       Close_Socket( &Socket );
	       PrintErrorAndExit( "FCM: recv() failed or connection closed prematurely. Exiting.\n" );
	  } else {
	       /* Get the displacement */
	       pos = Receive_Buffer + 2*sizeof( int );
	       memcpy( &Time, pos, sizeof (double) );
	       pos = pos + sizeof( double );
	       memcpy( Data_To_Receive, pos, sizeof (double)*Order );
	       printf("Time %e\t Data Received:%e\n", Time, Data_To_Receive[0]);
	  }
     } else if ( WhatToDo == 4 ) { /* Send the force to the server */

	  /* Compose the NSEP_CSIG packet and send it to the server */
	  Length = 8 + 4*Order;
	  NSEP_Type = NSEP_CSIG; /*NSEP_CSIG;*/

	  pos = Send_Buffer;
	  memcpy( pos, &Length, sizeof (int) );	

	  pos = pos + sizeof (int);
	  memcpy( pos , &NSEP_Type, sizeof (int) );

	  pos = pos + sizeof (double);
	  memcpy( pos, Data_To_Send, sizeof(double)*Order );
	  if ( send( Socket, Send_Buffer, Length, 0 ) != (int) sizeof (char)*Length ){
	       Close_Socket( &Socket );
	       PrintErrorAndExit( "FCM: send() sent a different number of bytes than expected. Exiting." );
	  }
	  printf("%i\t%e\n", Send_Buffer[0], Data_To_Send[0] );

     } else if ( WhatToDo == 5 ) {
	  /* Ask the PNSE server for the experiment state */
	  Length = 12;
	  NSEP_Type = NSEP_QUERY;
	  What = NSEP_EXPSTATE; /* NSEP_EXPSTATE */

	  pos = Send_Buffer;
	  memcpy( pos, &Length, sizeof (int) );

	  pos = pos + sizeof (int);
	  memcpy( pos,& NSEP_Type, sizeof (int) );

	  pos = pos + sizeof (int);
	  memcpy( pos ,&What, sizeof (int) );
	  if ( send( Socket, Send_Buffer, Length, 0 ) != (int) sizeof (char)*Length ){
	       Close_Socket( &Socket );
	       PrintErrorAndExit( "FCM: send() sent a different number of bytes than expected. Exiting." );
	  }

	  if ( recv( Socket, Receive_Buffer, MAXBUFLEN, 0 ) <= 0 ){
	       Close_Socket( &Socket );
	       PrintErrorAndExit( "FCM: recv() failed or connection closed prematurely. Exiting.\n" );
	  } else {
	       /* Get the Experimental state */
	       memcpy( &ExpState, Receive_Buffer + 2*sizeof(int), sizeof (int) );
	  }

	  if ( ExpState == NSEP_ES_FINISHED ){
	       Data_To_Receive[0] = -9999.0;
	  }
	  
     } else if ( WhatToDo == 6 ) { /* Terminate Connection with the server */
	  
          /* Tell the PNSE Server to change the state to finished
	   * Run State = 3
	   */
	  Send_Client_State( Socket, Send_Buffer, &Error, NSEP_CS_FINISHED );
	  if ( Error == 1 ){
	       Close_Socket( &Socket );
	       PrintErrorAndExit( "CGM: send() sent a different number of bytes than expected. Exiting." );
	  }
	  Close_Socket( &Socket );
     } else {
	  Close_Socket( &Socket );
	  fprintf( stderr, "CGM: Invalid WhatToDo value %i. WhatToDo should be between 0 and 3.\n", WhatToDo );
     }
} 

void Login_Server( const int Socket, const Remote_Machine_Info Login, char *const Send_Buffer, char *const Receive_Buffer, int *Error )
{
     int Error_TCP;  /* To control and identify TCP communication problems */
     int Error_Code; /* To identify the answer from the server */
     int NSEP_Type;  /* Type of message received from the PNSE server */

     /* Send the account name and password */
     Send_Login_Information( Socket, Login, Send_Buffer, &Error_TCP );
     if ( Error_TCP ){
	  (*Error) = 1;  /* Error during send( ). */
     }

     /* Receive the login result */
     Receive_Login_Information( Socket, Receive_Buffer, &Error_TCP, &NSEP_Type, &Error_Code );
     if ( Error_TCP ){
	  (*Error) = 2; /* Error during recv( ). */ 
     } else {
	  if( NSEP_Type == __AP_LOGIN_REPLY ) {
	       switch ( Error_Code ) {
	       case __AP_LOGIN_REPLY_OK:		    
		    printf( "PNSE: login successful.\n" );
		    break;
	       case __AP_LOGIN_REPLY_EXCEEDTRIALTIME:
		    fprintf( stderr, "PNSE: Exceeded trial time to login.\n" );
		    (*Error) = 3; /* Login information incorrect*/
		    break;
	       case __AP_LOGIN_REPLY_WRONGPASSWORD:
		    fprintf( stderr, "PNSE: Unrecognised Password.\n" );
		    (*Error) = 3; /* Login information incorrect*/
		    break;
	       case __AP_LOGIN_REPLY_UNRECOGNIZEDACCOUNT:
		    fprintf( stderr, "PNSE: Unrecognised Account.\n" );
		    (*Error) = 3; /* Login information incorrect*/
		    break;
	       case __AP_LOGIN_REPLY_NOTLOGINPACKET:
		    fprintf( stderr, "PNSE: The sent packet to the server was not a login packet.\n" );
		    (*Error) = 3; /* Login information incorrect*/
		    break;
	       }
	  }
     }
}

void Send_Login_Information( const int Socket, const Remote_Machine_Info Login, char *const Send_Buffer, int *Error_TCP )
{
     size_t Length;    /* Length of the message */
     int NSEP_Type; /* NSEP type of message */
     char *pos;     /* Pointer to identify the position within the message */

     Length = 8 + strlen( Login.Account_Name ) + 1 + strlen( Login.Account_Password ) + 1;
     NSEP_Type = __AP_LOGIN_REQUEST;

     pos = Send_Buffer;
     memcpy( pos, &Length, sizeof (int) );
     
     pos = pos + sizeof (int);
     memcpy( pos, &NSEP_Type, sizeof (int) );

     pos = pos + sizeof (int);
     memcpy( pos, Login.Account_Name, (strlen( Login.Account_Name ) + 1)*sizeof (char) );
     
     pos = pos + (strlen( Login.Account_Name ) + 1)*sizeof (char);
     memcpy( pos, Login.Account_Password, (strlen( Login.Account_Password ) + 1)*sizeof (char) );

     /* Send the message through TCP/IP and check if there were any problems with
      * the communication
      */
     if ( send( Socket, Send_Buffer, Length, 0) != (ssize_t) Length ){
	  (*Error_TCP) = 1;
     } else {
	  (*Error_TCP) = 0;
     }
}

void Receive_Login_Information( const int Socket, char *const Receive_Buffer,
				int *Error_TCP, int *const NSEP_Type, int *const Error_Code ){

     unsigned int Length; /* Length of the message */

     /* Receive the message from the server and check if there were any problems
      * with the communication
      */
     if ( recv( Socket, Receive_Buffer, MAXBUFLEN, 0 ) <= 0 ){
	  (*Error_TCP) = 1;
     } else {
	  (*Error_TCP) = 0;
     }
     
     /* Put the contents of the received message into the correct variables */
     memcpy( &Length, Receive_Buffer, sizeof (int) );
     memcpy( &(*NSEP_Type), Receive_Buffer + sizeof (int), sizeof (int) );
     memcpy( &(*Error_Code), Receive_Buffer + 2*sizeof (int), sizeof (int) );
}

void Send_Client_State( const int Socket, char *const Send_Buffer,
			int *Error_TCP, const int Run_State )
{
     size_t Length;    /* Length of the message */
     int NSEP_Type;    /* NSEP type of message */
     char *pos;        /* Pointer to identify the position within the message */

     Length = 12;
     NSEP_Type = NSEP_CLNSTATE;
     	  
     pos = Send_Buffer;
     memcpy( pos, &Length, sizeof (int) );
     
     pos = pos + sizeof (int);
     memcpy( pos, &NSEP_Type, sizeof (int) );
     
     pos = pos + sizeof (int);    
     memcpy( pos, &Run_State, sizeof (int) );
     if( send( Socket, Send_Buffer, Length, 0 ) != (ssize_t) Length ){
	  (*Error_TCP) = 1;
     }

}
