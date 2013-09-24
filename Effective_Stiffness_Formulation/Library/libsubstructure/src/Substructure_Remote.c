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

#include "Substructure.h"
#include "Substructure_Remote.h"
#include "Print_Messages.h"


void Substructure_Remote_Init( const char *RemoteType, const char *IPAddress, const char *Port, const int NSub,
			       const int *const DOF, const char *Description, Remote_t *const Remote )
{

     int i;   /* A counter */

     Substructure_Remote_Identify( RemoteType, &Remote->Type );
     
     Remote->IP = strdup( IPAddress );
     Remote->Port = strdup( Port );
     Remote->NSub = NSub;

     Remote->DOFs = (int *) calloc( (size_t) NSub, sizeof(int) );
     for( i = 0; i < NSub; i++ ){
	  Remote->DOFs[i] = DOF[i];
     }

     if( Remote->Type == REMOTE_TCP || Remote->Type == REMOTE_UDP ){
	  Substructure_Remote_SetupClientConnection( Remote );
     } else if ( Remote->Type == REMOTE_NCREE ){
     } else if ( Remote->Type == REMOTE_OF ){
     } else if ( Remote->Type == REMOTE_CELESTINA ){
     } else assert(0);

     Remote->Description = strdup( Description );

}

void Substructure_Remote_Identify( const char *RemoteType, int *const Type )
{

     int ID;   /* A counter */
     bool Found = false;

     ID = 0;
     while ( ID < NUM_REMOTE_TYPE && !Found ){
	  if ( strcmp( Substructure_RemoteType[ID], RemoteType ) == 0 ){
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
	       fprintf( stderr, "[......] %d) %s.\n", ID+1, Substructure_RemoteType[ID] );
	  }
	  exit( EXIT_FAILURE );
     } else {
	  /* Assign the identity */
	  (*Type) = ID;
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

void Substructure_Remote_SetupClientConnection( Remote_t *const RemoteNode )
{

     struct addrinfo addrCriteria;    /* Create a generic address storage to handle both IPv6 and IPv4
				       * addresses */
     struct addrinfo *Server_Addr;    /* List of server addresses */
     struct addrinfo *addr;
     int rtnVal;

     memset( &addrCriteria, 0, sizeof(addrCriteria) );   /* Initialise the structure */
     addrCriteria.ai_family = AF_UNSPEC;                 /* Set the value to accept any type of addresses */
     if ( RemoteNode->Type == REMOTE_TCP ){
	  addrCriteria.ai_socktype = SOCK_STREAM;        /* Only streaming sockets */
	  addrCriteria.ai_protocol = IPPROTO_TCP;        /* Use TCP protocol */
     } else if ( RemoteNode->Type == REMOTE_UDP ){
	  addrCriteria.ai_socktype = SOCK_DGRAM;         /* Only datagrams sockets will be used */
	  addrCriteria.ai_protocol = IPPROTO_UDP;        /* Use UDP protocol */
     } else { 
	  assert ( RemoteNode->Type == REMOTE_TCP || RemoteNode->Type == REMOTE_UDP );
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
	       } else {
		    printf( "Successfully connected to the UDP Server: %s on port %s.\n", RemoteNode->IP,
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

void Substructure_Remote_Send( const double *const Data, const unsigned int Data_Length, const int Socket )
{
     char *Msg;
     size_t Length;

     Msg = (char *) Data;
     Length = Data_Length*sizeof(double);

     if ( send( Socket, Msg, Length, 0) != (ssize_t) Length ){
	  Print_Header( ERROR );
	  fprintf( stderr, "send() sent a different number of bytes than expected.\n" );
     }

}

void Substructure_Remote_Receive( double *const Data, const unsigned int Data_Length, const int Socket )
{
     char *Msg;
     ssize_t bytesRcvd, totalBytesRcvd;
     size_t Length;

     totalBytesRcvd = 0;
     Length = sizeof(double) * Data_Length;
     Msg = (char *) Data;

     while (totalBytesRcvd < (ssize_t) Length )
     {
	  if ((bytesRcvd = recv( Socket, Msg, Length,0)) <= 0){
	       Print_Header( ERROR );
	       fprintf( stderr, "recv() failed or connection closed prematurely.\n" );
	       exit( EXIT_FAILURE );
	  }
	  totalBytesRcvd += bytesRcvd;
     }
}

void Substructure_Remote_Destroy( Remote_t *const Remote, const int Order )
{
     double *Send = NULL;
     
     Send = (double *) calloc( (size_t) Order + 1, sizeof(double) );

     Remote->NSub = 0;
     free( Remote->DOFs );

     if( Remote->Type == REMOTE_TCP || Remote->Type == REMOTE_UDP || Remote->Type == REMOTE_CELESTINA ){
	  Send[0] = -9999.0;
	  Substructure_Remote_Send( Send, (unsigned int) Order + 1, Remote->Socket );
     } else if( Remote->Type == REMOTE_NCREE ){
	  /* Using NSEP Protocol */
	  /* Say to PNSE server that the process has finished. WhatToDo = 6 */
	  //Communicate_With_PNSE( 6, 0.0,  Send, Send, 0 );
     } else if( Remote->Type == REMOTE_OF ){
	  /* Using OpenFresco */
	  //Communicate_With_OpenFresco( Send, Send, 1, 10 );
     } else assert( Remote->Type >= 0 || Remote->Type < NUM_REMOTE_TYPE );

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
