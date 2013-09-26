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
#define OF_INFO_DATA_LENGTH 11
#define OF_DATA_SIZE 256

/**
 * \brief Supported sub-structure types
*/
enum Substructure_Remote_ID { REMOTE_TCP,       /*!< \brief Standard TCP/IP connection */
			      REMOTE_UDP,       /*!< \brief Standard UDP connection */
			      REMOTE_NSEP,     /*!< \brief NSEP protocol */
			      REMOTE_OF,        /*!< \brief OpenFresco protocol */
			      REMOTE_CELESTINA  /*!< \brief Celestina protocol */
};

static const char *Substructure_RemoteType[] = {"TCP",
						"UDP",
						"NSEP",
						"OpenFresco",
						"Celestina"};

/**
 * \brief Structure to handle the Server IP and the communication port.
 *
 * This structure is used to store the IP address of the Server (IPv4) where the sub-stepping part of the
 * algorithm is performed. The port that will be used throughout the TCP/IP communication is also stored.
 */
typedef struct Remote {
     char *IP;               /*!< \brief IP address (IPv4) of the server. */
     char *Port;             /*!< \brief Port that will be used for TCP/IP communication. */
     char *Account_Name;     /*!< \brief Account information (NSEP protocol only). */
     char *Account_Password; /*!< \brief Account password (NSEP protocol only). */
     int   Socket;           /*!< \brief Socket. */
     int   NSub;             /*!< \brief Number of substructures to be simulated at the remote site. */
     int  *DOFs;             /*!< \brief Degrees of freedom that are affected by the substructures. Position
			      * of the substructures in the global matrix. */
     char *Description;      /*!< \brief Description of the experimental substructure.*/
     int Type;               /*!< \brief Type of remote substructures. */
} Remote_t;

/**
 * \brief Accepts an incoming connection from a client.
 *
 * Accepts an incoming connectrion from a client as long as it uses the TCP/IP protocol. If for whatever
 * reasons the connection was not accepted, the routine exits with \c EXIT_FAILURE.
 *
 * \pre \c Server_Socket must be a valid TCP socket created with Substructure_Remote_SetupServer().
 *
 * \param[in] Server_Socket TCP Socket to be used (server side).
 *
 * \return A TCP socket from the client whose connection has been accepted.
 */
int Substructure_Remote_AcceptTCPClientConnection( const int Server_Socket );

/**
 * \brief Closes a TCP or a UDP socket. This routine works on both Unix/Linux and Windows operating systems.
 *
 * \pre \c Socket must contain a valid TCP/IP or UDP socket.
 *
 * \param[out] Socket The TCP/UDP socket to be closed.
 *
 * \post The desired socket is closed and is no longer valid.
 */
void Substructure_Remote_CloseSocket( int *const Socket );

/**
 *
 * \param[in]  Socket
 * \param[in]  WhatToDo
 * \param[in]  Size
 * \param[in]  Data_To_Send
 * \param[out] Data_To_Receive
 */
void Substructure_Remote_OpenFresco( const int Socket, const int WhatToDo, const int Size, const double *const Data_To_Send, double *const Data_To_Receive );

/**
 * \brief Finalises the connection with the remote facility and frees dynamically allocated memory.
 *
 * The routine informs the remote facility that the connection is going to be stopped. The way this is
 * performed depends on the different protocols:
 *
 * - \c REMOTE_TCP, \c REMOTE_UDP or \c REMOTE_CELESTINA: A message is sent with a value of -9999.0 
 * - \c REMOTE_NSEP: A message of type \c XXXX is sent to the NSEP server.
 * - \c REMOTE_OF: A message is sent with a 99 identifier.
 *
 * 
 * \param[in,out] Remote
 * \param[in] Order
 */
void Substructure_Remote_Destroy( Remote_t *const Remote, const int Order );

/**
 * \param[in]  RemoteType
 * \param[in]  IPAddress
 * \param[in]  Port
 * \param[in]  NSub
 * \param[in]  DOF
 * \param[in]  Description
 * \param[out] Remote
 */
void Substructure_Remote_Init( const char *RemoteType, const char *IPAddress, const char *Port, const int NSub,
			       const int *const DOF, const char *Description, Remote_t *const Remote );

/**
 * \param[in]  RemoteType
 * \param[out] Type
 */
void Substructure_Remote_Identify( const char *RemoteType, int *const Type );

/**
 * \param[in] Socket
 * \param[in] Data_Length
 * \param[in] Datatype_Size
 * \param[in] Data
 */
void Substructure_Remote_Send( const int Socket, const unsigned int Data_Length, const size_t Datatype_Size, const char *const Data );

/**
 * \param[in,out] RemoteNode
 */
void Substructure_Remote_SetupClientSocket( Remote_t *const RemoteNode );

/**
 * \param[in] Port
 * \param[in] Socket_Type
 *
 * \return
 */
int Substructure_Remote_SetupServer( const char* Port, const int Socket_Type );

/**
 * \param[in]  Socket
 * \param[in]  Data_Length
 * \param[in]  Datatype_Size
 * \param[out] Data
 */
void Substructure_Remote_Receive( const int Socket, const unsigned int Data_Length, const size_t Datatype_Size, char *const Data );

/**
 * \brief Prints the IP address and port to \c stdout.
 *
 * \param[in] address Address (IP and port) to be printed on \c stdout.
 *
 * \post The routine would provide warnings if the given address cannot be converted/is invalid.
 */
void Substructure_Remote_PrintSocketAddress( struct sockaddr *const address );


/**
 * \brief Frees the memory allocated during GetServerInformation().
 *
 * The memory allocated during the GetServerInformation() routine by the strdup() function is deallocated.
 * 
 * \param[in,out] Server Struct storing the IP address, port and login information.
 *
 */
void Delete_NetworkInformation( Remote_t *const Server );

#endif /* SUBSTRUCTURE_REMOTE_H_ */
