/**
 * \file Send_Receive_Data.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 10th of August 2011
 *
 *  \brief Prototypes of the communication subroutines.
 *
 *  This file contains the prototypes of the communication subroutines used in the substructure algorithm. This includes opening
 *  and closing sockets and sending and receiving data. For the moment, only the TCP/IP protocol has been considered.
 *
 *  \todo Implement UDP protocol.
 */

#ifndef SEND_RECEIVE_DATA_H_
#define SEND_RECEIVE_DATA_H_

/**
 * \brief Structure to handle the Server IP and the communication port
 *
 * This structure is used to store the IP address of the Server (IPv4) where the sub-stepping part
 * of the algorithm is performed. The port that will be used throughout the TCP/IP communication is also stored.
 */
typedef struct {
     char IP[20];  /*!< IP address (IPv4) of the server. */
     unsigned short int Port; /*!< Port that will be used for TCP/IP communication. */
     char Type[11];  /* Type of Server to connect: OpenFresco, PNSE or Custom */
     char Account_Name[10]; /*!< Account information */
     char Account_Password[10]; /*!< Account password */
} Remote_Machine_Info;

#define PROTOCOL_ADWIN  0
#define PROTOCOL_CUSTOM 1
#define PROTOCOL_NSEP   2
#define PROTOCOL_OF     3

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
 *
 * \post
 * - \c Server has the server IP address and the port that will be used for the TCP/IP communication.
 *
 * \sa Remote_Machine_Info.
 */
void GetServerInformation( Remote_Machine_Info *const Server );

/**
 * \brief Opens a TCP/IP socket.
 *
 * This routine opens a TCP/IP socket that will be open throughout the algorithm.
 *
 * \pre
 * - \c Server must contain the IP address of the server and the port.
 *
 * \param[in] Server Struct that contains the IP address of the server and the port that will be used.
 * \param[out] sock The socket.
 *
 * \post
 * - The variable \c sock contains an open TCP/IP socket between the Server and the Client.
 *
 * \sa Remote_Machine_Info.
 */
void OpenSocket( const Remote_Machine_Info Server, int *sock );

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
void Send_Data( const float *Data, const int DATA_LENGTH, const int sock );

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
void Receive_Data( float *Data, const int DATA_LENGTH, const int sock );

void Send_Effective_Matrix( const float *const Eff_Mat, const int Protocol_Type, const int OrderC, int *const Socket );

void Do_Substepping( const float *const DispTdT0_c, float *const DispTdT, float *const fcprevsub, float *const fc, const int Protocol_Type, const float Time, const float DeltT, const int Num_Sub, const int Socket, const int OrderC, const int Pos_Couple );
void Close_Connection( int *Socket, const int Protocol_Type, const int OrderC, const int Num_Steps, const int Num_Sub );

void Close_Socket( int *Socket );

#endif /* SEND_RECEIVE_DATA_H_ */
