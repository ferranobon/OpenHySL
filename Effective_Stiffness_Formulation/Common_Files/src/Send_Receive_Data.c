/**
 * \file Send_Receive_Data.c
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 10th of August 2011
 *
 *  \brief Source code of the communication subroutines.
 *
 *  This file contains the source code of the communication subroutines used in the substructure algorithm. This includes opening
 *  and closing sockets and sending and receiving data. For the moment, only the TCP/IP protocol has been considered.
 *
 *  \todo Implement UDP protocol.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if _WIN32_
#include <winsock2.h>
#else
#include <sys/socket.h>
#include <sys/types.h>
#include <arpa/inet.h>
#include <unistd.h>
#endif

#include "Send_Receive_Data.h"         /* Definition of struct Remote_Machine_Info. */
#include "ErrorHandling.h"             /* Headers for Error Handling functions. */
#include "OpenFresco_Communication.h"  /* OpenFresco header files */
#include "NSEP_Definitions.h"          /* NSEP Definitions and constants */
#include "NSEP_Communication_Sync.h"   /* NSEP header files */

#if ADWIN_
#include "RoutinesADwin.h"             /* Communicate with ADwin */
#endif


void GetServerInformation( Remote_Machine_Info *const Server )
{
	FILE *InFile;

	/* Read the IP address and the port from file */
	InFile = fopen( "Connection.txt", "r" );

	if ( InFile != NULL ){

	     /* See the type of client */	     
	     fscanf( InFile, "%s", (*Server).Type );
	     /* Server IP address */
	     fscanf( InFile, "%s", (*Server).IP );

	     /* Server port */
	     fscanf( InFile, "%hu", &(*Server).Port );

	     /* If the Server is of type PNSE, then read the account name and the
	      *	password associated to it.
	      */
	     if ( !strcmp( (*Server).Type, "PNSE" ) ){
		  /* Account name */
		  fscanf( InFile, "%s", (*Server).Account_Name );
		  /* Password */
		  fscanf( InFile, "%s", (*Server).Account_Password );
	     }
	     fclose( InFile );		  
	} else {
		ErrorFileAndExit( "Error while setting up the connection because it was not possible to open ", "Connection.txt" );
	}

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

void OpenSocket( const Remote_Machine_Info Server, int *Socket )
{
  
     struct sockaddr_in Server_Addr;  /* Type for handling the connections (IP, port) */

     if ( (*Socket = socket( AF_INET, SOCK_STREAM, 0 )) < 0 ){
	  PrintErrorAndExit( "socket() failed" );
     }

     memset( &Server_Addr, 0, sizeof( Server_Addr ) );
     /* Specify the address family as internet */
     Server_Addr.sin_family = AF_INET;

     /* Set the server address */
#if _WIN32_
     Server_Addr.sin_addr.S_un.S_addr = inet_addr( Server.IP );
#else     
     Server_Addr.sin_addr.s_addr = inet_addr( Server.IP );
#endif

     /* Set the server port */     
     Server_Addr.sin_port = htons( Server.Port );

     if ( connect( *Socket, (struct sockaddr *) &Server_Addr, sizeof( Server_Addr ) ) < 0 ){
	  Close_Socket( Socket );
	  PrintErrorAndExit( "connect() failed." );
     } else {
	  printf("Successfully connected to the %s Server: %s on port %hi.\n", Server.Type, Server.IP, Server.Port );	  
     }
}

void Send_Data( const float *Data, const int DATA_LENGTH, const int sock )
{

  if ( send(sock, Data, sizeof Data  * DATA_LENGTH, 0) != (int)sizeof Data * DATA_LENGTH ){
	  PrintErrorAndExit( "send() sent a different number of bytes than expected" );
  }

}

void Receive_Data( float *Data, const int DATA_LENGTH, const int sock )
{
  int bytesRcvd, totalBytesRcvd;

  totalBytesRcvd = 0;

  while (totalBytesRcvd < (int)sizeof Data * DATA_LENGTH)
    {
	  if ((bytesRcvd = recv(sock, Data, sizeof Data * DATA_LENGTH,0)) <= 0){
		  PrintErrorAndExit( "recv() failed or connection closed prematurely" );
	  }
      totalBytesRcvd += bytesRcvd;
    }
}

void Send_Effective_Matrix( const float *const Eff_Mat, const int Protocol_Type, const int OrderC, int *const Socket )
{
     int i, j; /* Counters */
     Remote_Machine_Info Server;
     float *Send, *Recv;

     Send = calloc( OrderC, sizeof(float) );
     Recv = calloc( 3*OrderC, sizeof(float) );

     switch( Protocol_Type ){

#if ADWIN_
     case PROTOCOL_ADWIN:
	  printf( "Running without TCP communication.\n" );

	  /* Send matrix Gc to ADwin */
	  ADWIN_SetGc( Eff_Mat, OrderC*OrderC );
	  break;
#endif
     case PROTOCOL_CUSTOM:  
	  /* Using custom communication protocol */
	  GetServerInformation( &Server );
	  printf( "Establishing connection with the Server.\n" );
	  
	  OpenSocket( Server, &(*Socket) );
	  
	  Send_Data( Eff_Mat, OrderC*OrderC, (*Socket) );
	  break;
     case PROTOCOL_NSEP:
	  /* Using NSEP Protocol */
	  /* Open the Socket */
	  printf( "Establishing connection with PNSE server.\n" );
	  Communicate_With_PNSE( 0, 0.0, Send, Recv, 0 );
	  
	  /* Send the matrix Gc to the PNSE server in order to reach the FCM */
	  /*
	   * Note that the CGM can only send a NSEP_CMD message with the size of order, therefore
	   * it is needed to send the matrix Gc per rows to fullfill this requisite. This also means, that
	   * the FCM must send as many NSEP_CSIG packets as the number of rows to keep everything synchronised
	   */

	  for ( i = 0; i < OrderC; i++ ){
	       for ( j = 0; j < OrderC; j++ ){
		    Send[j] = Eff_Mat[i*OrderC + j];
	       }
	       /* Send the matrix Keinv_c to PNSE Server */
	       Communicate_With_PNSE( 1, 0.0, Send, Recv, OrderC );
	       /* This is done so that PNSE do not overtake the first step */
	       Communicate_With_PNSE( 2, 0.0,  Send, Recv, OrderC );
	  }
	  break;
     case PROTOCOL_OF:
	  /* Using OpenFresco */
	  printf( "Establishing connection with OpenFresco.\n" );
	  Communicate_With_OpenFresco( Eff_Mat, Recv, OrderC*OrderC, 1 );
	  /* TODO Implement Send the Matrix G in OpenFresco. Wait for the answer from Andreas */
	  break;
     }

     free( Send );
     free( Recv );
}

void Do_Substepping( const float *const DispTdT0_c, float *const DispTdT, float *const fcprevsub, float *const fc, const int Protocol_Type, const float Time, const float DeltaT, const int Num_Sub, const int Socket, const int OrderC, const int Pos_Couple )
{

     static int i;
     float *Recv;

     Recv = calloc( 3*OrderC, sizeof(float) );


     switch ( Protocol_Type ){
#if ADWIN_
     case PROTOCOL_ADWIN:
	  /* Tell ADwin to perform the substepping process */
	  ADWIN_Substep( DispTdT0_c, &Recv[0], &Recv[1], &Recv[2], OrderC, Num_Sub, DeltaT/(float)Num_Sub );
	  //Recv[0] = uc[0]; Recv[1] = 0.0; Recv[2] = 0.0;
	  break;
#endif
     case PROTOCOL_CUSTOM:
	  /* Using custom communication protocol */
	  Send_Data( DispTdT0_c, OrderC, Socket );

	  Receive_Data( Recv, 3*OrderC, Socket );
	  break;
     case PROTOCOL_NSEP:
	  /* Using NSEP Protocol */
	  Communicate_With_PNSE( 1, Time, DispTdT0_c, Recv, OrderC );
	  /* Receive the force from the PNSE server. WhatToDo = 2 */
	  Communicate_With_PNSE( 2, Time, DispTdT0_c, Recv, 3*OrderC );
	  break;
     case PROTOCOL_OF:
	  /* Using OpenFresco */
	  Communicate_With_OpenFresco( DispTdT0_c, Recv, OrderC, 3 ); 
	  break;
     }

     for ( i = 0; i < OrderC; i++ ){
	  DispTdT[(Pos_Couple - 1) + i] = Recv[i];
	  fcprevsub[(Pos_Couple - 1) + i] = Recv[OrderC + i];
	  fc[(Pos_Couple - 1 ) + i] = Recv[2*OrderC + i];
     }

     free( Recv );

}
/*
void Send_Data( const int Socket, const int Data_Type_Size, char* const To_Send, const int Data_Length )
{
     int Bytes_Left, Bytes_Sent;

     char *Data = To_Send;

     Bytes_Left = Data_Length * Data_Type_Size;

     while ( Bytes_Left > 0 ){
	  Bytes_Sent = send( Socket, Data, Bytes_Left, 0);
	  Bytes_Left -= Bytes_Sent;
	  Data += Bytes_Sent;
     }

}

void Receive_Data( const int Socket, const int Data_Type_Size, char *const To_Receive, const int Data_Length )
{
     char *Data = To_Receive;
     int Bytes_Left, Bytes_Read;

     Bytes_Left = Data_Length * Data_Type_Size;

     while ( Bytes_Left > 0 ){
	  Bytes_Read = recv( Socket, Data, Bytes_Left,0 );
	  Bytes_Left -= Bytes_Read;
	  Data += Bytes_Read;
    }
}
*/

void Close_Connection( int *Socket, const int Protocol_Type, const int OrderC, const int Num_Steps, const int Num_Sub )
{

     float *Send;

#if ADWIN_
     float *ADWIN_DATA;
#endif

     Send = calloc( OrderC, sizeof(float) );

     switch ( Protocol_Type ){
#if ADWIN_
     case PROTOCOL_ADWIN:
	  /* Connect directly to ADwin */
	  ADWIN_DATA = calloc( Num_Sub*Num_Steps*22, sizeof( float ) );
	  GetDataADwin( Num_Steps, Num_Sub, ADWIN_DATA );
	  free( ADWIN_DATA );
	  break;
#endif
     case PROTOCOL_CUSTOM:
	  /* Using custom communication protocol */
	  Send[0] = -9999.0;
	  Send_Data( Send, OrderC, *Socket );
	  /* Close the socket */
	  Close_Socket( Socket );
	  break;
     case PROTOCOL_NSEP:
	  /* Using NSEP Protocol */
	  /* Say to PNSE server that the process has finished. WhatToDo = 6 */
	  Communicate_With_PNSE( 6, 0.0,  Send, Send, 0 );
	  break;
     case PROTOCOL_OF:
	  /* Using OpenFresco */
	  Communicate_With_OpenFresco( Send, Send, 1, 10 );
	  break;
     }

     free( Send );
}

void Close_Socket( int *Socket )
{

#if _WIN32_
     closesocket( *Socket );
#else
     close( *Socket );
#endif
}
