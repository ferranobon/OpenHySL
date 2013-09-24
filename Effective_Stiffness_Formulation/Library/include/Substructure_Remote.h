/**
 * \file Substructure_Remote.h
 * \author Ferran Obón Santacana
 * \version 1.0 
 * \date 24th of August 2013
 *
 * \brief Routines to deal with remote substructures.
 *
 * These routines deal with the creation and destriction of remote substructures.
 */
#ifndef SUBSTRUCTURE_REMOTE_H_
#define SUBSTRUCTURE_REMOTE_H_

#include <sys/socket.h>  /* For struct sockaddr */

#define NUM_REMOTE_TYPE 5   /*!< \brief Number of recognised remote sub-structure types. */
#define MAXPENDING      5   /*!< \brief Maximum outstanding connection requests */

/**
 * \brief Supported sub-structure types
*/
enum Substructure_Remote_ID { REMOTE_TCP,       /*!< \brief Standard TCP/IP connection */
			      REMOTE_UDP,       /*!< \brief Standard UDP connection */
			      REMOTE_NCREE,     /*!< \brief NSEP protocol */
			      REMOTE_OF,        /*!< \brief OpenFresco protocol */
			      REMOTE_CELESTINA  /*!< \brief Celestina protocol */
};

static const char *Substructure_RemoteType[] = {"TCP",
						"UDP",
						"NCREE",
						"OpenFresco",
						"Celestina"};

/**
 * \brief Structure to handle the Server IP and the communication port.
 *
 * This structure is used to store the IP address of the Server (IPv4) where the sub-stepping part of the
 * algorithm is performed. The port that will be used throughout the TCP/IP communication is also stored.
 */
typedef struct Remote{
     char *IP;               /*!< IP address (IPv4) of the server. */
     char *Port;             /*!< Port that will be used for TCP/IP communication. */
     char *Account_Name;     /*!< Account information */
     char *Account_Password; /*!< Account password */
     int   Socket;           /*!< Socket */
     int   NSub;             /*!< Number of substructures to be simulated in the remote site */
     int  *DOFs;             /*!< Degrees of freedom that are affected by the substructures. Position of the
			      *   substructures in the global matrix */
     char *Description;      /*!< \brief Description of the experimental sub-structure.*/
     int Type;               /* Type of remote substructures */
} Remote_t;

int Substructure_Remote_AcceptTCPClientConnection( int Server_Socket );
void Substructure_Remote_CloseSocket( int *const Socket );
void Substructure_Remote_Destroy( Remote_t *const Remote, const int Order );
void Substructure_Remote_Init( const char *RemoteType, const char *IPAddress, const char *Port, const int NSub,
			       const int *const DOF, const char *Description, Remote_t *const Remote );
void Substructure_Remote_Identify( const char *RemoteType, int *const Type );
void Substructure_Remote_Send( double *const Data, const unsigned int Data_Length, const int Socket );
void Substructure_Remote_SetupClientConnection( Remote_t *const RemoteNode );
int Substructure_Remote_SetupServer( const char* Port, const int Socket_Type );
void Substructure_Remote_Receive( double *const Data, const unsigned int Data_Length, const int Socket );
void Substructure_Remote_PrintSocketAddress( struct sockaddr *const address );


/**
 * \brief Frees the memory allocated during GetServerInformation.
 *
 * The memory allocated during the GetServerInformation() routine by the strdup() function is deallocated.
 * 
 * \param[out] Server Struct storing the IP address, port and login information.
 *
 */
void Delete_NetworkInformation( Remote_t *const Server );


#endif /* SUBSTRUCTURE_REMOTE_H_ */
