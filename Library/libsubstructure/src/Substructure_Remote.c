#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>

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

#include "Definitions.h"
#include "Print_Messages.h"
#include "Substructure_Remote.h"
#include "Substructure_Remote_NSEP.h"
#include "Substructure_Remote_OpenFresco.h"

void Substructure_Remote_Init( const char *RemoteType, const char *IPAddress, const char *Port, const char *Acc_Name, const char *Acc_Password, const unsigned int NSub,
			       const int *const DOF, const char *Description, Remote_t *const Remote )
{

     unsigned int i;   /* A counter */

     Remote->Type = Substructure_Remote_Identify( RemoteType );
     
     Remote->IP = strdup( IPAddress );
     Remote->Port = strdup( Port );

     Remote->Account_Name = strdup( Acc_Name );
     Remote->Account_Password = strdup( Acc_Password );

     Remote->NSub = NSub;

     Remote->DOFs = (int *) calloc( (size_t) NSub, sizeof(int) );
     if( Remote->DOFs == NULL ){
	  Print_Header( ERROR );
	  fprintf( stderr, "Substructure_Remote_Init(): Out of memory.\n" );
	  exit( EXIT_FAILURE );
     }

     for( i = 0; i < NSub; i++ ){
	  Remote->DOFs[i] = DOF[i];
     }

     if( Remote->Type >= 0 || Remote->Type < NUM_REMOTE_TYPE ){
	  Substructure_Remote_SetupClientSocket( Remote );
     }

     if ( Remote->Type == REMOTE_NSEP ){
	  /* Log into the PNSE server and set the client state to running */
	  Substructure_Remote_NSEP( Remote, NSEP_LOG, 0.0f, 0, NULL, NULL );
     } else if ( Remote->Type == REMOTE_OF ){
	  /* Configure the connection with OpenFresco */
	  Substructure_Remote_OpenFresco( Remote->Socket, OF_REMOTE_SETUP, Remote->NSub, NULL, NULL );
     }

     Remote->Description = strdup( Description );

}

int Substructure_Remote_Identify( const char *RemoteType )
{

     int ID;   /* A counter */
     bool Found = false;

     ID = 0;
     while ( ID < NUM_REMOTE_TYPE && !Found ){
	  if ( strcmp( Substructure_Remote_Type[ID], RemoteType ) == 0 ){
	       Found = true;
	  } else {
	       ID = ID + 1;
	  }
     }

     /* The substructure in Type is not supported.*/
     if ( !Found ){
	  Print_Header( ERROR );
	  fprintf( stderr, "Substructure_Remote_Identify: The remote type '%s' is not supported. Valid remote options are:\n", RemoteType );
	  for( ID = 0; ID < NUM_REMOTE_TYPE; ID++ ){
	       fprintf( stderr, "[......] %d) %s.\n", ID+1, Substructure_Remote_Type[ID] );
	  }
	  exit( EXIT_FAILURE );
     } else {
	  /* Assign the identity */
	  return ID;
     }
}

int Substructure_Remote_SetupServer( const char* Port, const int Socket_Type )
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
     if ( Socket_Type == REMOTE_TCP ){
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
	  Print_Header( ERROR );
	  fprintf( stderr, "getaddrinfo() failed: %s.\n", gai_strerror( rtnVal ) );
	  exit( EXIT_FAILURE );
     }
     
     for( addr = Server_Addr; addr != NULL; addr = addr->ai_next ){
	  /* Create a TCP or UDP socket */
	  Socket = socket( Server_Addr->ai_family, Server_Addr->ai_socktype, Server_Addr->ai_protocol );
	  if( Socket < 0 ){
	       continue;
	  }

	  /* Bind to the local address */
	  if( Socket_Type == REMOTE_TCP ){
	       if ( (bind( Socket, Server_Addr->ai_addr, Server_Addr->ai_addrlen ) == 0 ) && (listen( Socket, MAXPENDING) == 0) ){
		    /* Print local address of socket */
		    addrSize = sizeof( Local_Addr );
		    if (getsockname( Socket, (struct sockaddr *) &Local_Addr, &addrSize ) < 0 ){
			 Print_Header( ERROR );
			 fprintf( stderr, "getsockname() failed.\n" );
			 exit( EXIT_FAILURE );
		    }
		    Print_Header( INFO );
		    printf( "Binding to " );
		    Substructure_Remote_PrintSocketAddress( (struct sockaddr *) &Local_Addr );
		    printf( "\n" );
		    break;
	       }
	  } else {
	       if ( bind( Socket, Server_Addr->ai_addr, Server_Addr->ai_addrlen ) == 0 ){
		       /* Print local address of socket */
		    addrSize = sizeof( Local_Addr );
		    if (getsockname( Socket, (struct sockaddr *) &Local_Addr, &addrSize ) < 0 ){
			 Print_Header( ERROR );
			 fprintf( stderr, "getsockname() failed.\n" );
		    }
		    Print_Header( INFO );
		    printf( "Binding to " );
		    Substructure_Remote_PrintSocketAddress( (struct sockaddr *) &Local_Addr );
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

void Substructure_Remote_PrintSocketAddress( struct sockaddr *const address )
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
	  Print_Header( WARNING );
	  fprintf( stderr, "[unknown type].\n" );
	  return;
     }

     /* Convert binary to printable address */
     if( inet_ntop(address->sa_family, numericAddress, addrBuffer, sizeof(addrBuffer)) == NULL){
	  Print_Header( WARNING );
	  fprintf( stderr, "[Invalid address].\n" ); /* Unable to convert */
     } else {
	  printf( "%s", addrBuffer );
	  if ( port != 0 ){
	       printf( "-%u", port );
	  }
     }

}

int Substructure_Remote_AcceptTCPClientConnection( const int Server_Socket )
{
     
     int Client_Socket;                      /* Socket descriptor for client */
     struct sockaddr_storage Client_Addr;    /* Client address */
     socklen_t Client_Addr_Length;           /* Length of client address data structure */

     /* Set the size of the client address structure (in-out parameter) */
     Client_Addr_Length = sizeof( Client_Addr );
    
     /* Wait for a client to connect */
     Client_Socket = accept( Server_Socket, (struct sockaddr *) &Client_Addr, &Client_Addr_Length);
     if ( Client_Socket < 0 ){
	  fprintf( stderr, "accept() failed.\n" );
	  exit( EXIT_FAILURE );
     }
    
     /* The Client is connected. Display a proper message */
     Print_Header( INFO );
     printf("Handling client: ");
     Substructure_Remote_PrintSocketAddress( (struct sockaddr *) &Client_Addr );
     printf( "\n" );

     return Client_Socket;
}

void Substructure_Remote_SetupClientSocket( Remote_t *const RemoteNode )
{

     struct addrinfo addrCriteria;    /* Create a generic address storage to handle both IPv6 and IPv4
				       * addresses */
     struct addrinfo *Server_Addr;    /* List of server addresses */
     struct addrinfo *addr;
     int rtnVal;

     memset( &addrCriteria, 0, sizeof(addrCriteria) );   /* Initialise the structure */
     addrCriteria.ai_family = AF_UNSPEC;                 /* Set the value to accept any type of addresses */
     if ( RemoteNode->Type == REMOTE_TCP || RemoteNode->Type == REMOTE_NSEP ||
	  RemoteNode->Type == REMOTE_OF || RemoteNode->Type == REMOTE_CELESTINA){

	  addrCriteria.ai_socktype = SOCK_STREAM;        /* Only streaming sockets */
	  addrCriteria.ai_protocol = IPPROTO_TCP;        /* Use TCP protocol */

     } else if ( RemoteNode->Type == REMOTE_UDP ){

	  addrCriteria.ai_socktype = SOCK_DGRAM;         /* Only datagrams sockets will be used */
	  addrCriteria.ai_protocol = IPPROTO_UDP;        /* Use UDP protocol */

     } else { 
	  assert( RemoteNode->Type >= 0 || RemoteNode->Type < NUM_REMOTE_TYPE );
     }

     /* The rest of the values of the structure have already been set to "0", meaning that the fields are set
      * to "don't care" */

     /* Get the addresses */
     rtnVal = getaddrinfo( RemoteNode->IP, RemoteNode->Port, &addrCriteria, &Server_Addr );
     if ( rtnVal != 0 ){
	  Print_Header( ERROR );
	  fprintf( stderr, "getaddrinfo() failed: %s", gai_strerror( rtnVal ) );
	  exit( EXIT_FAILURE );
     }

     RemoteNode->Socket = -1;
     for ( addr = Server_Addr; addr != NULL; addr = addr->ai_next ) {
	  /* Create the TCP/UDP socket */
	  RemoteNode->Socket = socket( addr->ai_family, addr->ai_socktype, addr->ai_protocol );
	  if ( RemoteNode->Socket < 0 ){
	       continue;               /* Socket creation failed; try next address */
	  }

	  /* Establish the connection with the Server to start the simulation */
	  if ( connect( RemoteNode->Socket, addr->ai_addr, addr->ai_addrlen ) == 0 ){
	       Print_Header( SUCCESS );
	       if( RemoteNode->Type == REMOTE_TCP ){
		    printf( "Successfully connected to the TCP Server: %s on port %s.\n", RemoteNode->IP,
			    RemoteNode->Port );
	       } else if ( RemoteNode->Type == REMOTE_UDP ) {
		    printf( "Successfully connected to the UDP Server: %s on port %s.\n", RemoteNode->IP,
			    RemoteNode->Port );
	       } else if( RemoteNode->Type == REMOTE_NSEP ){
		    printf( "Successfully connected to the NSEP Server: %s on port %s.\n", RemoteNode->IP,
			    RemoteNode->Port );
	       } else if( RemoteNode->Type == REMOTE_OF ){
		    printf( "Successfully connected to the OpenFresco Server: %s on port %s.\n", RemoteNode->IP,
			    RemoteNode->Port );
	       } else if( RemoteNode->Type == REMOTE_OF ){
		    printf( "Successfully connected to the Celestina Server: %s on port %s.\n", RemoteNode->IP,
			    RemoteNode->Port );
	       }
	       break;                  /* The socket has been successfully created, break and return Socket */
	  }
	  
	  /* The creation of the socket has failed. Try the next address */
#if WIN32
	  closesocket( RemoteNode->Socket );
#else
	  close( RemoteNode->Socket );
#endif
	  RemoteNode->Socket = -1;
     }
     
     if ( RemoteNode->Socket < 0 ){
	  Print_Header( ERROR );
	  fprintf( stderr, "The connection to the TCP Server %s on port %s could not be established.\n", RemoteNode->IP, RemoteNode->Port );
	  exit( EXIT_FAILURE );
     }

     /* Free the address list allocated by getaddrinfo() */
     freeaddrinfo(Server_Addr);
}	  

void Substructure_Remote_Send( const int Socket, const unsigned int Data_Length, const size_t Datatype_Size, const char *const Data )
{
     size_t Length;

     Length = Data_Length*Datatype_Size;

     if ( send( Socket, Data, Length, 0) != (ssize_t) Length ){
	  Print_Header( ERROR );
	  fprintf( stderr, "send() sent a different number of bytes than expected.\n" );
	  exit( EXIT_FAILURE );
     }

}

void Substructure_Remote_Receive( const int Socket, const unsigned int Data_Length, const size_t Datatype_Size, char *const Data )
{
     ssize_t bytesRcvd, totalBytesRcvd;
     size_t Length;

     totalBytesRcvd = 0;
     Length = Datatype_Size * Data_Length;

     while (totalBytesRcvd < (ssize_t) Length )
     {
	  if ((bytesRcvd = recv( Socket, Data, Length, 0 )) <= 0){
	       Print_Header( ERROR );
	       fprintf( stderr, "recv() failed or connection closed prematurely.\n" );
	       exit( EXIT_FAILURE );
	  }
	  totalBytesRcvd += bytesRcvd;
     }
}

void Substructure_Remote_Destroy( Remote_t *const Remote )
{
     HYSL_FLOAT *Send = NULL;
     
     Send = (HYSL_FLOAT *) calloc( (size_t) Remote->NSub + 1, sizeof(HYSL_FLOAT) );


     free( Remote->DOFs );

     if( Remote->Type == REMOTE_TCP || Remote->Type == REMOTE_UDP || Remote->Type == REMOTE_CELESTINA ){
	  Send[0] = -9999.0;
	  Substructure_Remote_Send( Remote->Socket, (unsigned int) Remote->NSub + 1, sizeof(HYSL_FLOAT), (char *const) Send );
     } else if( Remote->Type == REMOTE_NSEP ){
	  /* Using NSEP Protocol */
	  /* Say to PNSE server that the process has finished. WhatToDo = 6 */
#if _FLOAT_
	  Substructure_Remote_NSEP( Remote, NSEP_SET_TO_FINISHED, 0.0f, 0, Send, NULL );
#else
	  Substructure_Remote_NSEP( Remote, NSEP_SET_TO_FINISHED, 0.0, 0, Send, NULL );
#endif

	  free( Remote->Account_Name );
	  free( Remote->Account_Password );
     } else if( Remote->Type == REMOTE_OF ){
	  /* Using OpenFresco */
	  Substructure_Remote_OpenFresco( Remote->Socket, OF_REMOTE_DIE, Remote->NSub, NULL, NULL );
     } else assert( Remote->Type >= 0 || Remote->Type < NUM_REMOTE_TYPE );

     Remote->NSub = 0;
     free( Send );

     Substructure_Remote_CloseSocket( &Remote->Socket );
     Print_Header( SUCCESS );
     printf( "Socket with %s on port %s successfully closed.\n", Remote->IP, Remote->Port );

     free( Remote->IP );
     free( Remote->Port );
     free( Remote->Description );
}

void Substructure_Remote_CloseSocket( int *const Socket )
{
#if WIN32
     closesocket( *Socket );
#else
     close( *Socket );
#endif

}