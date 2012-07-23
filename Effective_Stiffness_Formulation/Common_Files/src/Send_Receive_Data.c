/**
 * \file Send_Receive_Data.c
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 10th of August 2011
 *
 *  \brief Source code of the communication subroutines.
 *
 *  This file contains the source code of the communication subroutines used in the substructure algorithm. This includes opening
 *  and closing sockets and sending and receiving data. It supports both TCP and UDP protocols.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if WIN32
#include <winsock2.h> /* For send(), recv(), sockadrr... */
#include <WS2tcpip.h> /* For socklen_t */
#else
#include <sys/socket.h>
#include <sys/types.h>
#include <arpa/inet.h>
#include <netdb.h>
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
	  fscanf( InFile, "%s", (*Server).Port );

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

int Setup_Server_Socket( const char* Port, const int Socket_Type )
{

     
     struct addrinfo addrCriteria; /* Create a generic address storage to handle both IPv6 and IPv4 addresses */
     struct addrinfo *Server_Addr; /* List of server addresses */
     struct addrinfo *addr;
     struct sockaddr_storage Local_Addr;
     socklen_t addrSize;
     int rtnVal;
     int Socket = -1;


     memset( &addrCriteria, 0, sizeof(addrCriteria) ); /* Initialisee the structure */
     addrCriteria.ai_family = AF_UNSPEC;               /* Set the value to accept any type of addresses */
     addrCriteria.ai_flags = AI_PASSIVE;               /* Acccept on any address/port */
     if ( Socket_Type == PROTOCOL_TCP ){
	  addrCriteria.ai_socktype = SOCK_STREAM;           /* Only streaming sockets */
	  addrCriteria.ai_protocol = IPPROTO_TCP;           /* Use TCP protocol */
     } else {
	  addrCriteria.ai_socktype = SOCK_DGRAM;            /* Only datagrams sockets will be used */
	  addrCriteria.ai_protocol = IPPROTO_UDP;           /* Use UDP protocol */
     }
     /* The rest of the values of the structure have already been set to "0", meaning that the fields are set
      * to "don't care" */

     /* Get the addresses */
     rtnVal = getaddrinfo( NULL, Port, &addrCriteria, &Server_Addr );
     if ( rtnVal != 0 ){
	  PrintErrorDetailAndExit( "getaddrinfo() failed", gai_strerror( rtnVal ) );
     }
     
     for( addr = Server_Addr; addr != NULL; addr = addr->ai_next ){
	  /* Create a TCP or UDP socket */
	  Socket = socket( Server_Addr->ai_family, Server_Addr->ai_socktype, Server_Addr->ai_protocol );
	  if( Socket < 0 ){
	       continue;
	  }

	  /* Bind to the local address */
	  if( Socket_Type == PROTOCOL_TCP ){
	       if ( (bind( Socket, Server_Addr->ai_addr, Server_Addr->ai_addrlen ) == 0 ) && (listen( Socket, MAXPENDING) == 0) ){
		    /* Print local address of socket */
		    addrSize = sizeof( Local_Addr );
		    if (getsockname( Socket, (struct sockaddr *) &Local_Addr, &addrSize ) < 0 ){
			 PrintErrorAndExit( "getsockname() failed" );
		    }
		    printf( "Binding to " );
		    PrintSocketAddress( (struct sockaddr *) &Local_Addr );
		    printf( "\n" );
		    break;
	       }
	  } else {
	       if ( bind( Socket, Server_Addr->ai_addr, Server_Addr->ai_addrlen ) == 0 ){
		       /* Print local address of socket */
		    addrSize = sizeof( Local_Addr );
		    if (getsockname( Socket, (struct sockaddr *) &Local_Addr, &addrSize ) < 0 ){
			 PrintErrorAndExit( "getsockname() failed" );
		    }
		    printf( "Binding to " );
		    PrintSocketAddress( (struct sockaddr *) &Local_Addr );
		    printf( "\n" );
		    break;
	       }
	  }
	  
	  close( Socket );
	  Socket = -1;
     }


     /* Free the address list allocated by getaddrinfo */
     freeaddrinfo( Server_Addr );


     return Socket; 
}

void PrintSocketAddress( struct sockaddr *const address )
{
     if (address == NULL ){
	  return;
     }

     void *numericAddress;  /* Pointer to binary address */
     char addrBuffer[INET6_ADDRSTRLEN]; /* buffer to contain result */
     in_port_t port;    /* Port to print */

     /* Set pointer to address based on address family */
     switch ( address->sa_family ){
     case AF_INET:
	  numericAddress = &((struct sockaddr_in *) address)->sin_addr;
	  port = ntohs(((struct sockaddr_in *) address)->sin_port);
	  break;
     case AF_INET6:
	  numericAddress = &((struct sockaddr_in6 *) address)->sin6_addr;
	  port = ntohs((( struct sockaddr_in6 *) address)->sin6_port);
	  break;
     default:
	  fprintf( stderr, "[unknown type]\n" );
	  return;
     }

     /* Convert binary to printable address */
     if( inet_ntop(address->sa_family, numericAddress, addrBuffer, sizeof(addrBuffer)) == NULL){
	  fprintf( stderr, "[Invalid address]\n" ); /* Unable to convert */
     } else {
	  printf( "%s", addrBuffer );
	  if ( port != 0 ){
	       printf( "-%u", port );
	  }
     }

}

int Accept_TCP_Client_Connection( int Server_Socket )
{
     
     int Client_Socket;                      /* Socket descriptor for client */
     struct sockaddr_storage Client_Addr;    /* Client address */
     socklen_t Client_Addr_Length;           /* Length of client address data structure */

     /* Set the size of the client address structure (in-out parameter) */
     Client_Addr_Length = sizeof( Client_Addr );
    
     /* Wait for a client to connect */
     Client_Socket = accept( Server_Socket, (struct sockaddr *) &Client_Addr, &Client_Addr_Length);
     if ( Client_Socket < 0 ){
	  PrintErrorAndExit( "accept() failed.\n" );
     }
    
     /* The Client is connected. Display a proper message */
     printf("Handling client: ");
     PrintSocketAddress( (struct sockaddr *) &Client_Addr );
     printf( "\n" );

     return Client_Socket;
}

int Setup_Client_Socket( const Remote_Machine_Info Server, const int Socket_Type )
{

     
     struct addrinfo addrCriteria; /* Create a generic address storage to handle both IPv6 and IPv4 addresses */
     struct addrinfo *Server_Addr; /* List of server addresses */
     struct addrinfo *addr;
     int rtnVal;
     int Socket;

     memset( &addrCriteria, 0, sizeof(addrCriteria) ); /* Initialisee the structure */
     addrCriteria.ai_family = AF_UNSPEC;               /* Set the value to accept any type of addresses */
     if ( Socket_Type == PROTOCOL_TCP ){
	  addrCriteria.ai_socktype = SOCK_STREAM;           /* Only streaming sockets */
	  addrCriteria.ai_protocol = IPPROTO_TCP;           /* Use TCP protocol */
     } else {
	  addrCriteria.ai_socktype = SOCK_DGRAM;            /* Only datagrams sockets will be used */
	  addrCriteria.ai_protocol = IPPROTO_UDP;           /* Use UDP protocol */
     }
     /* The rest of the values of the structure have already been set to "0", meaning that the fields are set
      * to "don't care" */

     /* Get the addresses */
     rtnVal = getaddrinfo( Server.IP, Server.Port, &addrCriteria, &Server_Addr );
     if ( rtnVal != 0 ){
	  PrintErrorDetailAndExit( "getaddrinfo() failed", gai_strerror( rtnVal ) );
     }
     
     Socket = -1;
     for ( addr = Server_Addr; addr != NULL; addr = addr->ai_next ) {
	  /* Create the TCP/UDP socket */
	  Socket = socket( addr->ai_family, addr->ai_socktype, addr->ai_protocol );
	  if ( Socket < 0 ){
	       continue;  /* Socket creation failed; try next address */
	  }

	  /* Establish the connection with the Server to start the simulation */
	  if ( connect( Socket, addr->ai_addr, addr->ai_addrlen ) == 0 ){
	       printf( "Successfully connected to the %s Server: %s on port %s.\n", Server.Type, Server.IP, Server.Port );
	       break;    /* The socket has been successfully created, break and return Socket */
	  }
	  
	  close( Socket ); /* The creation of the socket has failed. Try the next address */
	  Socket = -1;
     }
     
     /* Free the address list allocated by getaddrinfo() */
     freeaddrinfo(Server_Addr);

     return Socket;
}

void Send_Data( float *const Data, const unsigned int Data_Length, const int sock )
{
     char *Msg;
     size_t Length;

     Msg = (char *) Data;
     Length = Data_Length*sizeof(float);

     if ( send(sock, Msg, Length, 0) != (ssize_t) Length ){
	  PrintErrorAndExit( "send() sent a different number of bytes than expected" );
  }

}

void Receive_Data( float *const Data, const unsigned int Data_Length, const int sock )
{
     char *Msg;
     ssize_t bytesRcvd, totalBytesRcvd;
     size_t Length;

     totalBytesRcvd = 0;
     Length = sizeof(float) * Data_Length;
     Msg = (char *) Data;

     while (totalBytesRcvd < (ssize_t) Length )
     {
	  if ((bytesRcvd = recv(sock, Msg, Length,0)) <= 0){
	       PrintErrorAndExit( "recv() failed or connection closed prematurely" );
	  }
	  totalBytesRcvd += bytesRcvd;
     }
}

void Send_Effective_Matrix( float *const Eff_Mat, const int Protocol_Type, const unsigned int OrderC, int *const Socket )
{
     unsigned int i, j; /* Counters */
     Remote_Machine_Info Server;
     float *Send = NULL, *Recv = NULL;

     Send = (float *) calloc( (size_t) OrderC, sizeof(float) );
     Recv = (float *) calloc( (size_t) 3*OrderC, sizeof(float) );

     switch( Protocol_Type ){

#if ADWIN_
     case PROTOCOL_ADWIN:
	  printf( "Running without TCP communication.\n" );

	  /* Send matrix Gc to ADwin */
	  ADWIN_SetGc( Eff_Mat, OrderC*OrderC );
	  break;
#endif
     case PROTOCOL_TCP:  
	  /* Using TCP communication protocol */
	  GetServerInformation( &Server );
	  printf( "Establishing connection with the TCP Server.\n" );
	  
	  *Socket = Setup_Client_Socket( Server, PROTOCOL_TCP );
	  if ( *Socket < 0 ){
	       PrintErrorAndExit( "Setup_Client_Socket() failed. Unable to connect" );
	  }
	  
	  Send_Data( Eff_Mat, OrderC*OrderC, (*Socket) );
	  break;
     case PROTOCOL_UDP:
	  /* Using the UDP protocol */
	  GetServerInformation( &Server );
	  printf( "Establishing connection with the UDP Server.\n" );


	  *Socket = Setup_Client_Socket( Server, PROTOCOL_UDP );
	  if ( *Socket < -1 ){
	       PrintErrorAndExit( "Setup_Client_Socket() failed. Unable to connect" );
	  }

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
	  /* Using OpenFresco. Exit with failure if the connection could not be established */
	  if ( Communicate_With_OpenFresco( Eff_Mat, Recv, OrderC*OrderC, 1 ) < 0 ){
	       exit( EXIT_FAILURE );
	  }
	  if ( Communicate_With_OpenFresco( Eff_Mat, Recv, OrderC*OrderC, 3 ) < 0 ){
	       exit( EXIT_FAILURE );
	  }
	  /* TODO Implement Send the Matrix G in OpenFresco. Wait for the answer from Andreas */
	  break;
     }

     free( Send );
     free( Recv );
}

void Do_Substepping( float *const DispTdT0_c, float *const DispTdT, float *const fcprevsub, float *const fc, const int Protocol_Type, const float Time, const int Socket, const unsigned int OrderC, const unsigned int *Pos_Couple )
{

     unsigned int i;
     float *Recv = NULL;

     Recv = (float *) calloc( (size_t) 3*OrderC, sizeof(float) );


     switch ( Protocol_Type ){
#if ADWIN_
     case PROTOCOL_ADWIN:
	  /* Tell ADwin to perform the substepping process */
	  ADWIN_Substep( DispTdT0_c, &Recv[0], &Recv[1], &Recv[2], OrderC );
	  break;
#endif
     case PROTOCOL_TCP:
	  /* Using TCP communication protocol */
	  Send_Data( DispTdT0_c, OrderC, Socket );

	  Receive_Data( Recv, 3*OrderC, Socket );
	  break;
     case PROTOCOL_UDP:
	  /* Using UDP communication protocol */

	  Send_Data( DispTdT0_c, OrderC, Socket );
	  if ( recv( Socket, Recv, sizeof(float)*3*OrderC,0) != (int) sizeof(float)*3*OrderC ){    /* sizeof returns an unsigned integer ? */
	       PrintErrorAndExit( "recv() failed in connected UDP mode" );
	  }
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
	  DispTdT[(Pos_Couple[i] - 1) + i] = Recv[i];
	  fcprevsub[(Pos_Couple[i] - 1) + i] = Recv[OrderC + i];
	  fc[(Pos_Couple[i] - 1 ) + i] = Recv[2*OrderC + i];
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

void Receive_Data_TCP( const int Socket, const int Data_Type_Size, char *const To_Receive, const int Data_Length )
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

void Close_Connection( int *Socket, const int Protocol_Type, const unsigned int OrderC, const unsigned int Num_Steps, const unsigned int Num_Sub )
{

     float *Send = NULL;

     float *ADWIN_DATA = NULL;

     Send = (float *) calloc( (size_t) OrderC, sizeof(float) );

     switch ( Protocol_Type ){
#if ADWIN_
     case PROTOCOL_ADWIN:
	  /* Connect directly to ADwin */
	  ADWIN_DATA = (float *) calloc( (size_t) Num_Sub*Num_Steps*NUM_CHANNELS, sizeof( float ) );
	  GetDataADwin( Num_Steps, Num_Sub, ADWIN_DATA );
	  free( ADWIN_DATA );
	  break;
#endif
     case PROTOCOL_TCP:
	  /* Using custom communication protocol */
	  Send[0] = -9999.0;
	  Send_Data( Send, OrderC, *Socket );
	  /* Close the socket */
	  Close_Socket( Socket );
	  break;
     case PROTOCOL_UDP:
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

#if WIN32
     closesocket( *Socket );
#else
     close( *Socket );
#endif
}
