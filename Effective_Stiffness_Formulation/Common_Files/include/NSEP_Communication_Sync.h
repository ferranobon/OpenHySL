/**
 * \file NSEP_Communication_Sync.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 4th of November 2011
 * 
 * \brief Prototypes for the communication subroutines compliant with NSEP protocol
 *
 * This file contains the prototypes of the communication subroutines that are used
 * when connecting to a PNSE server type.
 *
 * \todo 
 * - Implement the asynchronous version to send/receive data from the server.
 * - Document the several NSEP packages.
 */

#ifndef NSEP_COMMUNICATION_SYNC_H_
#define NSEP_COMMUNICATION_SYNC_H_

#include "Send_Receive_Data.h"

#define MAXBUFLEN  1000  /**< Maximum size of the buffer used in sending/receiving messages between the PNSE server and the CGM */

/**
 * \brief Main routine to communicate with PNSE using the NSEP protocol as a CGM.
 *
 * This routine handles several aspects of the communication of a CGM with the PNSE server. It covers the aspects of
 * opening a socket, login, send and receive messages and to define the running status of the CGM in a compliant way with
 * the NSEP protocol. This is achieved through the \c WhatToDo variable.
 *
 * \pre
 * - A PNSE type server must be waiting for connections.
 *
 * \param[in] WhatToDo Action that will be performed.
 * - 0     Connect to the PNSE, login and set the state to running.
 * - 1     Send something to the server ( can be whatever ) NSEP_CMD Packet.
 * - 2     Request something from the server ( can be whatever ) NSEP_QUERY Packet.
 * - 3     Terminate Connection (NSEP_CLNSTATE packet).
 * \param[in] Time Only used when \f$WhatToDo = 1\f$. It contains the pseudo-time (experiment time) when the data to be sent
 * \c Data_To_Send was aquired or calculated.
 * \param[in] Data_To_Send Only used when \f$WhatToDo = 1\f$. It contains the data to be sent to the PNSE server.
 * \param[out] Data_To_Receive Only used when \f$WhatToDo = 2\f$. On exit it contains the data received from the server.
 * \param[in] Order Only referenced when \f$WhatToDo = 1\f$ or \f$WhatToDo = 2\f$. It is the amount of data that has to be
 * sent or received from the PNSE server.
 *
 * \post
 * - If any error occurs during the execution of this routine, like communication problems, the program shows the appropiate
 * error message befor exiting with EXIT_FAILURE.
 */
void Communicate_With_PNSE( const int WhatToDo, float Time, const float *const Data_To_Send, float *const Data_To_Receive, const int Order );

/**
 * \brief Login to the server in order to proceed with the test.
 *
 * Routine that deals with the login process required by the NSEP protocol. It makes use of the routines Send_Login_Information()
 * and Receive_Login_Information() to deal with this process.
 *
 * \pre
 * - The TCP/IP connection with the PNSE server must be stablished.
 * - \c Login must be properly initialised through GetServerInformation and contain an account name and password.
 *
 * \param[in] Socket An open TCP/IP socket with the PNSE server.
 * \param[in] Send_Buffer The message to be sent to the server. In this case a __AP_LOGIN_REQUEST. Not referenced on entry,
 * although it must contain an initialised memory space.
 * \param[out] Receive_Buffer On exit this will contain the answer comming from the server regarding the login request.
 * \param[out] Error On exit it will contain the kind of error, if any, that ocurred within this routine.
 * - 0 No error ocurred.
 * - 1 Error during sending the message to the PNSE server.
 * - 2 Error while receiving the message from the PNSE server.
 * - 3 Error during the login process.
 *
 * \post A __AP_LOGIN_REQUEST is sent to the server and afterwards a __AP_LOGIN_REPLY is received. The routines inform through
 * \c Error the result of such operation.
 */
void Login_Server( const int Socket, const Remote_Machine_Info Login, char *const Send_Buffer, char *const Receive_Buffer, int *Error );

/**
 *
 * \brief Sends a __AP_LOGIN_REQUEST packet to the PNSE server.
 *
 * A __AP_LOGIN_REQUEST packet is sent to the server. It must be the next message sent by the client to the PNSE server. 
 * All the other messages will be rejected.
 *
 * \pre 
 * - The TCP/IP connection with the PNSE server must be stablished.
 * - \c Login must be properly initialised through GetServerInformation and contain an account name and password.
 * 
 * \param[in] Socket An open TCP/IP socket with the PNSE server.
 * \param[in] Login Contains the account name and password to be sent to the PNSE server.
 * \param[in] Send_Buffer The message to be sent to the server. In this case a __AP_LOGIN_REQUEST. Not referenced on entry,
 * although it must contain an initialised memory space.
 * \param[out] Error_TCP It can have two possible values:
 * - 0 if there was no problem sending the message.
 * - 1 if there was an error while sending a message to the PNSE server. In that case the an appropiate error message is shown
 * and the program exits with EXIT_FAILURE.
 *
 * \post A __AP_LOGIN_REQUEST packet contaning the account name and password to login has been sent to the PNSE server
 * using TCP/IP if no communication errors are found.
 */
void Send_Login_Information( const int Socket, const Remote_Machine_Info Login, char *const Send_Buffer, int *Error_TCP );

/**
 * \brief Receives a __AP_LOGIN_REPLY packet from the PNSE server after sending a __AP_LOGIN_REQUEST packet to the server.
 *
 * After sending a __AP_LOGIN_REQUEST to the server it send a packet back (__AP_LOGIN_REPLY) through TCP/IP to inform the client
 * about the request, wether it was successful or if there was a problem like invalid account or password.
 *
 * \pre
 * - The TCP/IP connection with the PNSE server must be stablished.
 * 
 * \param[in] Socket An open TCP/IP socket with the PNSE server.
 * \param[out] Receive_Buffer On exit this will contain the answer comming from the server regarding the login request.
 * \param[out] Error_TCP It can have two possible values:
 * - 0 if there was no problem receiving the message.
 * - 1 if there was an error while receiving a message from the PNSE server. In that case the an appropiate error message
 * is shown and the program exits with EXIT_FAILURE.
 * \param[out] NSEP_Type Part of the received message that contains what is received. In this case it should be __AP_LOGIN_REPLY.
 * \param[out] Error_Code Part of the received message that contains the answer from the PNSE server. It can be
 * __AP_LOGIN_REPLY_NOTLOGINPACKET, __AP_LOGIN_REPLY_UNRECOGNIZEDACCOUNT, __AP_LOGIN_REPLY_WRONGPASSWORD
 * __AP_LOGIN_REPLY_EXCEEDTRIALTIME and __AP_LOGIN_REPLY_OK
 *
 * \post A __AP_LOGIN_REPLY packet has been received from the server with an answer regarding the previous __AP_LOGIN_REQUEST if
 * no errors during communication are found.
 */
void Receive_Login_Information( const int Socket, char *const Receive_Buffer, int *Error_TCP, int *const NSEP_Type, int *const Error_Code );

/**
 * \brief Send a NSEP_CLNSTATE package to determine the current running state
 *
 * This routine sends a message to the PNSE server that contains the state of the client. That is, it sends a NSEP_CLNSTATE
 * package to the server.
 *
 * \pre 
 * - The TCP/IP connection with the PNSE server must be stablished.
 * - The Client must be successfully logged in.
 * 
 * \param[in] Socket An open TCP/IP socket with the PNSE server.
 * \param[in] Send_Buffer The message to be sent to the server. In this case a NSEP_CLNSTATE packet. Not referenced on entry,
 * although it must contain an initialised memory space.
 * \param[out] Error_TCP It can have two possible values:
 * - 0 if there was no problem sending the message.
 * - 1 if there was an error while sending a message to the PNSE server. In that case the an appropiate error message is shown
 * and the program exits with EXIT_FAILURE.
 * \param[in] Run_State It describes the running state of the client. It can have the values of NSEP_CS_DISCONNECTED,
 * NSEP_CS_NOTREADY, NSEP_CS_RUNNING or NSEP_CS_FINISHED.
 *
 * \post A NSEP_CLNSTATE package contaning the running state of the client has been sent to the PNSE server using TCP/IP
 *  if no communication errors are found.
 * 
 */
void Send_Client_State( const int Socket, char *const Send_Buffer, int *Error_TCP, const int Run_State );

#endif /* NSEP_COMMUNICATION_SYNC_H */
