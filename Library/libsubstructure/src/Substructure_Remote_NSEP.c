#include <stdio.h>
#include <stdlib.h>  /* For exit() */
#include <string.h>  /* For memcpy() */

#include "Substructure_Remote.h"
#include "Print_Messages.h"
#include "Substructure_Remote_NSEP.h"

#include "Definitions.h"

void Substructure_Remote_NSEP( const Remote_t *const Remote, const int WhatToDo, const hysl_float_t Time,
			       const unsigned int Size, const hysl_float_t *const Data_To_Send, hysl_float_t *const Data_To_Receive )
{
     static unsigned int Length;
     static int NSEP_Type;
     static int What;
     static int ExpState;

     char Send_Buffer[MAXBUFLEN_NSEP], Receive_Buffer[MAXBUFLEN_NSEP], *pos;


     /* Connect and login to the server and set the state to running */
     if ( WhatToDo == 0 ){

	  /* Login to the server */
	  Substructure_Remote_NSEP_LoginServer( Remote, Send_Buffer, Receive_Buffer );

	  /*
	   * Tell the PNSE Server to change the state to Running.
	   * Run State = 2
	   */
	  Substructure_Remote_NSEP_SetClientState( Remote->Socket, NSEP_CS_RUNNING, Send_Buffer );

     } else if ( WhatToDo == NSEP_LOG ) { /* Send the displacement to the server */

	  /* Compose the NSEP_CMD packet and send it to the server */
	  Length = 12 + 4*Size;
	  NSEP_Type = NSEP_CMD;

	  pos = Send_Buffer;
	  memcpy( pos, &Length, sizeof (int) );	

	  pos = pos + sizeof (int);
	  memcpy( pos , &NSEP_Type, sizeof (int) );

	  pos = pos + sizeof (int);
	  memcpy( pos, &Time, sizeof (hysl_float_t) );

	  pos = pos + sizeof (hysl_float_t);
	  memcpy( pos, Data_To_Send, sizeof(hysl_float_t)*Size );
	  Substructure_Remote_Send( Remote->Socket, Length, sizeof(char), Send_Buffer );
     } else if ( WhatToDo == NSEP_REQUEST_CSIG ) { /* Request a message from the server */
	  /* Send the NSEP_QUERY packet to query for the critical signal */
	  Length = 12;
	  NSEP_Type = NSEP_QUERY;
	  What = NSEP_CSIG;

	  pos = Send_Buffer;
	  memcpy( pos, &Length, sizeof (int) );

	  pos = pos + sizeof (int);
	  memcpy( pos,& NSEP_Type, sizeof (int) );

	  pos = pos + sizeof (int);
	  memcpy( pos ,&What, sizeof (int) );
	  Substructure_Remote_Send( Remote->Socket, Length, sizeof(char), Send_Buffer );
	  
	  Substructure_Remote_Receive( Remote->Socket, MAXBUFLEN_NSEP, sizeof(char), Receive_Buffer );
	  memcpy( Data_To_Receive, Receive_Buffer + 2*sizeof(int), sizeof (hysl_float_t)*Size );

     } else if ( WhatToDo == NSEP_REQUEST_CMD ) { /* Request a message from the server (FCM) */
	  /* Send the NSEP_QUERY packet to query for the command */
	  Length = 12;
	  NSEP_Type = NSEP_QUERY;
	  What = NSEP_CMD; /* NSEP_CMD = 8 */

	  pos = Send_Buffer;
	  memcpy( pos, &Length, sizeof (int) );

	  pos = pos + sizeof (int);
	  memcpy( pos,& NSEP_Type, sizeof (int) );

	  pos = pos + sizeof (int);
	  memcpy( pos ,&What, sizeof (int) );
	  Substructure_Remote_Send( Remote->Socket, Length, sizeof(char), Send_Buffer );

	  Substructure_Remote_Receive( Remote->Socket, MAXBUFLEN_NSEP, sizeof(char), Receive_Buffer );
	  /* Get the displacement */
	  pos = Receive_Buffer + 2*sizeof( int );
	  memcpy( &Time, pos, sizeof (hysl_float_t) );
	  pos = pos + sizeof( hysl_float_t );
	  memcpy( Data_To_Receive, pos, sizeof (hysl_float_t)*Size );

     } else if ( WhatToDo == NSEP_SEND_CSIG ) { /* Send the force to the server */

	  /* Compose the NSEP_CSIG packet and send it to the server */
	  Length = 8 + 4*Size;
	  NSEP_Type = NSEP_CSIG; /*NSEP_CSIG;*/

	  pos = Send_Buffer;
	  memcpy( pos, &Length, sizeof (int) );	

	  pos = pos + sizeof (int);
	  memcpy( pos , &NSEP_Type, sizeof (int) );

	  pos = pos + sizeof (hysl_float_t);
	  memcpy( pos, Data_To_Send, sizeof(hysl_float_t)*Size );
	  Substructure_Remote_Send( Remote->Socket, Length, sizeof(char), Send_Buffer );

     } else if ( WhatToDo == NSEP_REQUEST_STATE ) {
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
	  Substructure_Remote_Send( Remote->Socket, Length, sizeof(char), Send_Buffer );

	  Substructure_Remote_Receive( Remote->Socket, MAXBUFLEN_NSEP, sizeof(char), Receive_Buffer );
	  memcpy( &ExpState, Receive_Buffer + 2*sizeof(int), sizeof (int) );

	  if ( ExpState == NSEP_ES_FINISHED ){
	       Data_To_Receive[0] = -9999.0;
	  }
  
     } else if ( WhatToDo == NSEP_SET_TO_FINISHED ) { /* Terminate Connection with the server */
	  
          /* Tell the PNSE Server to change the state to finished
	   * Run State = 3
	   */
	  Substructure_Remote_NSEP_SetClientState( Remote->Socket, NSEP_CS_FINISHED, Send_Buffer );
     } else {
	  Print_Header( ERROR );
	  fprintf( stderr, "CGM: Invalid WhatToDo value %i. WhatToDo should be between 0 and 3.\n",
		   WhatToDo );
	  exit( EXIT_FAILURE );
     }
} 

void Substructure_Remote_NSEP_LoginServer( const Remote_t *const Remote, const char *const Send_Buffer,
					    char *const Receive_Buffer )
{
     int Error_Code; /* To identify the answer from the server */
     int NSEP_Type;  /* Type of message received from the PNSE server */

     /* Send the account name and password */
     Substructure_Remote_NSEP_SendLoginInformation( Remote->Socket, Remote->Account_Name, Remote->Account_Password, Send_Buffer);

     /* Receive the login result */
     Substructure_Remote_NSEP_ReceiveLoginInformation( Remote->Socket, &NSEP_Type, &Error_Code, Receive_Buffer );

     if( NSEP_Type == __AP_LOGIN_REPLY ) {
	  switch ( Error_Code ) {
	  case __AP_LOGIN_REPLY_OK:
	       Print_Header( SUCCESS );		    
	       printf( "PNSE: login successful.\n" );
	       break;
	  case __AP_LOGIN_REPLY_EXCEEDTRIALTIME:
	       Print_Header( ERROR );
	       fprintf( stderr, "PNSE: Exceeded trial time to login.\n" );
	       exit( EXIT_FAILURE );
	       break;
	  case __AP_LOGIN_REPLY_WRONGPASSWORD:
	       Print_Header( ERROR );
	       fprintf( stderr, "PNSE: Unrecognised Password.\n" );
	       exit( EXIT_FAILURE );
	       break;
	  case __AP_LOGIN_REPLY_UNRECOGNIZEDACCOUNT:
	       Print_Header( ERROR );
	       fprintf( stderr, "PNSE: Unrecognised Account.\n" );
	       exit( EXIT_FAILURE );
	       break;
	  case __AP_LOGIN_REPLY_NOTLOGINPACKET:
	       Print_Header( ERROR );
	       fprintf( stderr, "PNSE: The sent packet to the server was not a login packet.\n" );
	       exit( EXIT_FAILURE );
	       break;
	  }
     }
}

void Substructure_Remote_NSEP_SendLoginInformation( const int Socket, const char *const Name,
						     const char *const Password, const char *const Send_Buffer )
{
     size_t Length;    /* Length of the message */
     int NSEP_Type; /* NSEP type of message */
     char *pos;     /* Pointer to identify the position within the message */

     Length = 8 + strlen( Name ) + 1 + strlen( Password ) + 1;
     NSEP_Type = __AP_LOGIN_REQUEST;

     pos = Send_Buffer;
     memcpy( pos, &Length, sizeof (int) );
     
     pos = pos + sizeof (int);
     memcpy( pos, &NSEP_Type, sizeof (int) );

     pos = pos + sizeof (int);
     memcpy( pos, Name, (strlen( Name ) + 1)*sizeof (char) );
     
     pos = pos + (strlen( Name ) + 1)*sizeof (char);
     memcpy( pos, Password, (strlen( Password ) + 1)*sizeof (char) );

     /* Send the message through TCP/IP and check if there were any problems with the communication */
     Substructure_Remote_Send( Socket, (unsigned int) Length, sizeof(char), Send_Buffer );
}

void Substructure_Remote_NSEP_ReceiveLoginInformation( const int Socket, int *const NSEP_Type,
							int *const Error_Code, char *const Receive_Buffer )
{

     /* Receive the message from the server and check if there were any problems
      * with the communication
      */
     Substructure_Remote_Receive( Socket, MAXBUFLEN_NSEP, sizeof(char), Receive_Buffer );
     
     /* Put the contents of the received message into the correct variables */
     /* The first value of length sizeof(int) contains the Length of the message */
     memcpy( &(*NSEP_Type), Receive_Buffer + sizeof (int), sizeof (int) );
     memcpy( &(*Error_Code), Receive_Buffer + 2*sizeof (int), sizeof (int) );
}

void Substructure_Remote_NSEP_SetClientState( const int Socket, const int Run_State, char *const Send_Buffer )
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
     Substructure_Remote_Send( Socket, (unsigned int) Length, sizeof(char), Send_Buffer );
}
