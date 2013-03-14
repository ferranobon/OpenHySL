#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Substructure_Remote.h"

void Substructure_Remote_Init( const char *IPAddress, const char *Port, const char *Description,
			       Remote_t *const Sub )
{
     
     Sub->IP = strdup( IPAddress );
     Sub->Port = strdup( Port );

     Sub->Description = strdup( Description );

}

void Substructure_Remote_Connect( Remote_t *const RemoteNode, const int Type )
{

     struct addrinfo addrCriteria;    /* Create a generic address storage to handle both IPv6 and IPv4
				       * addresses */
     struct addrinfo *Server_Addr;    /* List of server addresses */
     struct addrinfo *addr;
     int rtnVal;

     memset( &addrCriteria, 0, sizeof(addrCriteria) );   /* Initialise the structure */
     addrCriteria.ai_family = AF_UNSPEC;                 /* Set the value to accept any type of addresses */
     if ( Type == REMOTE_TCP ){
	  addrCriteria.ai_socktype = SOCK_STREAM;        /* Only streaming sockets */
	  addrCriteria.ai_protocol = IPPROTO_TCP;        /* Use TCP protocol */
     } else if ( Type == REMOTE_UDP ){
	  addrCriteria.ai_socktype = SOCK_DGRAM;         /* Only datagrams sockets will be used */
	  addrCriteria.ai_protocol = IPPROTO_UDP;        /* Use UDP protocol */
     } else { 
	  assert ( Type == REMOTE_TCP || Type == REMOTE_UDP );
     }

     /* The rest of the values of the structure have already been set to "0", meaning that the fields are set
      * to "don't care" */

     /* Get the addresses */
     rtnVal = getaddrinfo( RemoteNode->IP, RemoteNode->Port, &addrCriteria, &Server_Addr );
     if ( rtnVal != 0 ){
	  Print_Header( ERROR );
	  fprintf( stderr, "getaddrinfo() failed", gai_strerror( rtnVal ) );
     }
     
     RemoteNode->Socket = -1;
     for ( addr = Server_Addr; addr != NULL; addr = addr->ai_next ) {
	  /* Create the TCP/UDP socket */
	  RemoteNode->Socket = socket( addr->ai_family, addr->ai_socktype, addr->ai_protocol );
	  if ( Socket < 0 ){
	       continue;           /* Socket creation failed; try next address */
	  }

	  /* Establish the connection with the Server to start the simulation */
	  if ( connect( RemoteNode->Socket, addr->ai_addr, addr->ai_addrlen ) == 0 ){
	       Print_Hdeader( SUCCESS );
	       printf( "Successfully connected to Server: %s on port %s.\n", Type, RemoteNode->IP,
		       RemoteNode->Port );
	       break;              /* The socket has been successfully created, break and return Socket */
	  }
	  
	  close( Remote->Socket ); /* The creation of the socket has failed. Try the next address */
	  Remote->Socket = -1;
     }
     
     /* Free the address list allocated by getaddrinfo() */
     freeaddrinfo(Server_Addr);
}

void Substructure_Remote_Destroy( Remote_t *const Sub )
{
     free( Sub->IP );
     free( Sub->Port );
     free( Sub->Description );
}
