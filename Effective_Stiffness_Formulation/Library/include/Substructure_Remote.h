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
			      REMOTE_NSEP,      /*!< \brief NSEP protocol */
			      REMOTE_OF,        /*!< \brief OpenFresco protocol */
			      REMOTE_CELESTINA  /*!< \brief Celestina protocol */
};


static const char *Substructure_Remote_Type[NUM_REMOTE_TYPE] = {"TCP",    
								"UDP",  
								"NSEP",
								"OpenFresco",
								"Celestina"
};

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
     unsigned int   NSub;    /*!< \brief Number of substructures to be simulated at the remote site. */
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
 * \brief Finalises the connection with the remote facility and frees dynamically allocated memory.
 *
 * The routine informs the remote facility that the connection is going to be stopped. The way this is
 * performed depends on the different protocols:
 *
 * - \c REMOTE_TCP, \c REMOTE_UDP or \c REMOTE_CELESTINA: A message is sent with a value of -9999.0 
 * - \c REMOTE_NSEP: A message of type \c NSEP_CS_FINISHED is sent to the NSEP server in order to set the
 *   client state to finished and end the test
 * - \c REMOTE_OF: A message is sent of type OF_REMOTE_DIE in order to finalise the connection with the
 *   OpenFresco server.
 *
 * The memory that has been dynamically allocated during Substructure_Remote_Init() is freed and returned to
 * the system.
 *
 * \pre \c Remote must be properly initialised through Substructure_Remote_Init().
 * 
 * \param[in,out] Remote Substrucute of type \c Remote_t. On input only the Remote.NSub is referenced. On
 * output, memory is deallocated and all the variables are reseted.
 *
 * \sa Remote_t.
 */
void Substructure_Remote_Destroy( Remote_t *const Remote );

/**
 * \brief Initialises a substructure that runs on a remote facility.

 * \param[in]  RemoteType   Type of remote facility. It can be one of the following:
 *                          - TCP: The exchange of data is performed using TCP/IP.
 *			    - UDP: The exchange of data is performed using UDP.
 *			    - NSEP: The exchange of data is performed according to the NSEP specifications
 *                            \cite Wang_2013.
 *			    - OpenFresco: The exchange of data is performed according to the OpenFresco
 *                            specifications.
 *			    - Celestina: The exchange of data is performed through Celestina.
 * \param[in]  IPAddress    IP address of the remote facility.
 * \param[in]  Port         Port used in exchanging data.
 * \param[in]  Acc_Name     Account name to login into a PNSE server. It is ignored otherwise.
 * \param[in]  Acc_Password Account password to login into a PNSE server. It is ignored otherwise.
 * \param[in]  NSub         Number of degrees of freedom that the remote facility has to handle.
 * \param[in]  DOF          Degrees of freedom that the remote facility has to handle. Its length must be \f$l
 *                          \geq N_{Sub}\f$.
 * \param[in]  Description  Brief description of the substructure.
 * \param[out] Remote       Substructure of type \c REMOTE. On output it stores all the information presented
 *                         above.
 *
 * \post
 * - \c Remote.Description is a duplicate of the string in \c Description.
 * - If there is not enough memory or \c RemoteType contains an invalid entry, the program exits with \c
 *   EXIT_FAILURE.
 * - The memory should be deallocated through Substructure_DeleteCouplingNodes() routine.
 *
 * \sa Remote_t, MAX_DESCRIPTION, Substructure_Id, Substructure_Remote_Id.
 */
void Substructure_Remote_Init( const char *RemoteType, const char *IPAddress, const char *Port, const char *Acc_Name,
			       const char *Acc_Password, const unsigned int NSub, const int *const DOF,
			       const char *Description, Remote_t *const Remote );

/**
 * \brief Given a remote type, the routine identifies if it is valid and returns the proper identifier.
 *
 * \param[in] RemoteType String with the remote type to be identified
 *
 * \return The type of remote substructure. It can be one of the following
 *         - REMOTE_TCP: The exchange of data is performed using TCP/IP.
 *	   - REMOTE_UDP: The exchange of data is performed using UDP.
 *	   - REMOTE_NSEP: The exchange of data is performed according to the NSEP specifications
 *           \cite Wang_2013.
 *	   - REMOTE_OF: The exchange of data is performed according to the OpenFresco specifications.
 *	   - REMOTE_CELESTINA: The exchange of data is performed through Celestina.
 *
 * \post The returned integer is a valid identifier for a remote substructure. If none is found, the routine
 * exits with \c EXIT_FAILURE.
 *
 * \sa Substructure_Remote_ID.
 */
int Substructure_Remote_Identify( const char *RemoteType );

/**
 * \brief Transmits a message to another socket.
 *
 * This routine transmits a message to another socket. The amount of bytes that are going to be sent is \f$
 * Bytes = Data\_Length\cdot Datatype\_Size\f$. The routine will exit with \c EXIT_FAILURE if an error during
 * the trasnmission is found.
 *
 * \pre 
 * - \c Socket must be a valid TCP socket created with Substructure_Remote_SetupServer() or
 *   Substructure_Remote_SetupClientSocket().
 * - \c Data must be of length \f$l \geq Data\_Length\f$.
 * - \c Data must be a tokenised array.
 *
 * \param[in] Socket         TCP/UDP Socket to be used.
 * \param[in] Data_Length    Amount of data to be transmitted through the socket.
 * \param[in] Datatype_Size  Size of the non-tokenised data. Helps in specifying the amount to bytes to be
 *                           transmited.
 * \param[in] Data           Tokenised array that contains the data to be sent. The size of the original data
 *                           type is specified in \c Datatype_Size.
 */
void Substructure_Remote_Send( const int Socket, const unsigned int Data_Length, const size_t Datatype_Size, const char *const Data );

/**
 * \brief Establishes a client/server connection using TCP or UDP protocols. Client side.
 *
 * This routine connects through, TCP or UDP, to a server. Upon success, a notification is issued but the
 * routine will exit with \c EXIT_FAILURE if the connection could not be established.
 *
 * \pre 
 * - A server of the specified type has to be waiting to accept incoming connections.
 * - \c RemoteNode has to be properly initialised through the Substructure_Remote_Init() routine.
 *
 * \param[in,out] RemoteNode On input, the IP address, port and the type of protocol are used. On output, \c
 *                           RemoteNode.Socket contains a valid socket.
 *
 * \sa Remote_t, Substructure_Remote_ID.
 */
void Substructure_Remote_SetupClientSocket( Remote_t *const RemoteNode );

/**
 *
 * \brief Configures a TCP or a UDP server.
 *
 * Configures a TCP or a UDP server so that it can accept incoming client connections. It binds to \c
 * AF_UNSPEC (any type of addresses) and to the given port. Upon success, a notification is issued but the
 * routine will exit with \c EXIT_FAILURE if the set-up is not successful.

 * \param[in] Port        Port where the server will listen to.
 * \param[in] Socket_Type Identifies the type of server (TCP/IP or UDP).
 *
 * \return The socket to be used for incoming client connections.
 */
int Substructure_Remote_SetupServer( const char* Port, const int Socket_Type );

/**
 * \brief Receives a message from a socket
 *
 * This routine receives a message from a socket. The amount of bytes that are expected to be recieved is \f$
 * Bytes = Data\_Length\cdot Datatype\_Size\f$. The routine will wait for a message to be received and it will
 * exit with \c EXIT_FAILURE if an error during the reception is found.
 *
 * \pre 
 * - \c Socket must be a valid TCP socket created with Substructure_Remote_SetupServer() or
 *   Substructure_Remote_SetupClientSocket().
 * - \c Data must be of length \f$l \geq Data\_Length\f$.
 * - \c Data must be a tokenised array.
 *
 * \param[in]  Socket         TCP/UDP Socket to be used.
 * \param[in]  Data_Length    Amount of data to be received through the socket.
 * \param[in]  Datatype_Size  Size of the non-tokenised data. Helps in specifying the amount to bytes to be
 *                            received.
 * \param[out] Data           Tokenised array that on output contains the received data. The size of the original data
 *                            type is specified in \c Datatype_Size.
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


#endif /* SUBSTRUCTURE_REMOTE_H_ */
