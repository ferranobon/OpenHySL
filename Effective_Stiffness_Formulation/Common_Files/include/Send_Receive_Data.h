/**
 * \file Send_Receive_Data.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 10th of August 2011
 *
 *  \brief Prototypes of the communication subroutines.
 *
 *  This file contains the prototypes of the communication subroutines used in the substructure algorithm. This includes opening
 *  and closing sockets and sending and receiving data. It supports both TCP and UDP protocols.
 *
 */

#ifndef SEND_RECEIVE_DATA_H_
#define SEND_RECEIVE_DATA_H_

#if WIN32
/* Do nothing */
#else 
#include <netdb.h>  /* For struct sockaddr */
#endif

#include "Conf_Parser.h"  /* For ConfFile struct */

#define MAXPENDING   5    /* Maximum outstanding connection requests */

/**
 * \brief Structure to handle the Server IP and the communication port
 *
 * This structure is used to store the IP address of the Server (IPv4) where the sub-stepping part
 * of the algorithm is performed. The port that will be used throughout the TCP/IP communication is also stored.
 */
typedef struct {
     char *IP;               /*!< IP address (IPv4) of the server. */
     char *Port;             /*!< Port that will be used for TCP/IP communication. */
     char *Account_Name;     /*!< Account information */
     char *Account_Password; /*!< Account password */
     int Type;             /*!< Type of Server to connect: OpenFresco, PNSE or Custom */
} Remote_Machine_Info;

#define NO_PROTOCOL     0
#define PROTOCOL_ADWIN  1
#define PROTOCOL_TCP    2
#define PROTOCOL_UDP    3
#define PROTOCOL_NSEP   4
#define PROTOCOL_OF     5

/**
 * \brief Gets the information to establish the TCP/IP connection.
 *
 * This routine gets the IP address of the server and the port from a file named "Connection.txt". In case the case
 * of error opening/reading the file, and error is thrown and the program exits with \c EXIT_FAILURE.
 *
 * \pre
 * - The input text file should be in ASCII format and it must contain no header.
 *      - The first line must contain the type of the Server: OpenFresco, PNSE, Custom
 * 	- The second line must contain the IP address (IPv4).
 * 	- The third line must specify the port.
 *      - The fourth line must specify the account name (PNSE only).
 *      - The fifth line must specify the Password (PNSE only).
 * - The file must be named "Connection.txt".
 *
 * \param[out] Server Struct to store the IP address of the server and the port that will be used.
 * \param[in] CFile Stores the variables specified in the configuration file.
 *
 * \post
 * - \c Server has the server IP address and the port that will be used for the TCP/IP communication.
 *
 * \sa Remote_Machine_Info.
 */
void GetServerInformation( Remote_Machine_Info *const Server, const ConfFile *const CFile );

/**
 * \brief Frees the memory allocated during GetServerInformation.
 *
 * The memory allocated during the GetServerInformation() routine by the strdup()
 * function is deallocated.
 * 
 * \param[out] Server Struct storing the IP address, port and login information.
 *
 */
void Delete_ServerInformation( Remote_Machine_Info *const Server );

/**
 * \brief Identify the communication protocol to be used
 *
 * The communication protocol to be used is identified, and a proper return value is given. It makes use of the function
 * Get_Server_Information().
 *
 * \pre The first line of the file \c Connection.txt must contain the desired protocol type
 *
 * \param[in] Type It contains a string with the names of the supported protocols.
 *
 * \return 
 * - 0 if the desired protocol is of type \c Custom.
 * - 1 if the desired protocol is of type \c PNSE.
 * - 2 if the desired protocol is of type \c OpenFresco.
 * - -1 if the desired protocol is not recognised.
 */
int Get_Type_Protocol( const char *Type );

int Setup_Server_Socket( const char* Port, const int Socket_Type );
void PrintSocketAddress( struct sockaddr *const address );
int Accept_TCP_Client_Connection( int Server_Socket );

int Setup_Client_Socket( const Remote_Machine_Info Server, const int Socket_Type );

/**
 * \brief Sends data to the server (blocking)
 *
 * This routine sends an array of data to the server using blocking communication. In case of errors during the
 * communication an error will be thrown and the program will exit with \c EXIT_FAILURE.
 *
 * \pre
 * - \c sock must be an open TCP/IP socket.
 * - The data array must be properly initialised.
 *
 * \param[in] Data The array to be sent through TCP/IP.
 * \param[in] DATA_LENGTH The length of the array to be sent.
 * \param[in] sock An open TCP/IP socket.
 *
 * \post
 * - The data inside \c Data is sent successfully through TCP/IP.
 */
void Send_Data( float *const Data, const unsigned int Data_Length, const int sock );

/**
 * \brief Receives data from the server (blocking)
 *
 * This routine receives an array of data from the server using blocking communication. In case of errors during the
 * communication an error will be thrown and the program will exit with \c EXIT_FAILURE.
 *
 * \pre
 * - \c sock must be an open TCP/IP socket.
 * - The data array must be properly initialised.
 *
 * \param[out] Data The array to be received through TCP/IP.
 * \param[in] DATA_LENGTH The length of the array to be sent.
 * \param[in] sock An open TCP/IP socket.
 *
 * \post
 * - The received data from the Server through TCP/IP is stored in the array \c Data.
 */
void Receive_Data( float *const Data, const unsigned int Data_Length, const int sock );

void Send_Effective_Matrix( float *const Eff_Mat, const unsigned int OrderC, int *const Socket, const Remote_Machine_Info Server );

void Do_Substepping( float *const DispTdT0_c, float *const DispTdT, float *const fcprevsub, float *const fc, const int Protocol_Type, const float Time, const int Socket, const unsigned int OrderC, const unsigned int *Pos_Couple );
void Close_Connection( int *Socket, const int Protocol_Type, const unsigned int OrderC, const unsigned int Num_Steps, const unsigned int Num_Sub );

void Close_Socket( int *Socket );

#endif /* SEND_RECEIVE_DATA_H_ */
