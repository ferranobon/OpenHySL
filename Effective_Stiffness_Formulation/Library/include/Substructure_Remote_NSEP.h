/**
 * \file Substructure_Remote_NSEP.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 26th of September 2013
 * 
 * \brief Prototypes for the communication subroutines compliant with NSEP protocol and Definition of various
 * NSEP and PNSE constants
 *
 * This file contains the prototypes of the communication subroutines that are used when connecting to a PNSE
 * server type and the definition of the constants of the NSEP protocol. This definitions are used to identify
 * the type of message that is sent between NSEP compliant facilities.
 *
 * The information here presented regarding the different NSEP packages has been extracted from 
 * \cite Wang_2013. The reader will find more in depth information about the NSEP protocol specifications in
 * there.
 *
 * \todo 
 * - Implement the asynchronous version to send/receive data from the server.
 */

#ifndef SUBSTRUCTURE_REMOTE_NSEP_H
#define SUBSTRUCTURE_REMOTE_NSEP_H

#include "Substructure_Remote.h"
#include "Definitions.h"

#define MAXBUFLEN_NSEP 1000 /*!< \brief Maximum size of the buffer used in sending/receiving messages between the
			     * PNSE server and the CGM */

/* Login information related constants */
#define __AP_ERROR          0  /*!< \brief Identifies a packet that is an error. */
#define __AP_LOGIN_REQUEST  1  /*!< \brief Identifies a packet that contains information (client name and
				* password) for login procedure. */
#define __AP_LOGIN_REPLY    2  /*!< \brief Identifies a packet sent by the PNSE server to the client program
				* as a response of a previously sent \c __AP_LOGIN_REQUEST packet. This packet
				* contains information that indicates whether or not the login attempt is
				* accepted. */
#define NSEP_SETPUSH	    3  /*!< \brief Identifies a packet that can be sent by a PNSE client to the server
				* to turn on or off the server's push mode of data transmission
				* thenceforward. */
#define NSEP_QUERY          4  /*!< \brief Identifies a packet that is used to request relevant information
				* from the server.*/
#define NSEP_EXPINFO	    5  /*!< \brief Identifies a packet that contains descriptive information of the
				* collaborative experiment.*/
#define NSEP_EXPSTATE       6  /*!< \brief Identifies a packet that contains the current running state of the
				* whole experiment as well as the running states of all clients.*/
#define NSEP_CLNSTATE       7  /*!< \brief Identifies a packet that contains the client's current running
				* state.*/
#define NSEP_CMD            8  /*!< \brief Identifies a packet that contains command values to be generated
				* for the FCM's.*/
#define NSEP_CSIG           9  /*!< \brief Identifies a packet that contains critical signal values as the
				* result of command execution.*/
#define NSEP_IM             10 /*!< \brief Identifies a packet that is used used to implement instant
				* messaging. */

#define __AP_LOGIN_REPLY_NOTLOGINPACKET		0  /*!< \brief Indicates that the connected application failed
						    * to send a packet as its first packet to the server after
						    * the connection was established. Depending on the
						    * implementation of the PNSE server, after sending this \c
						    * __AP_LOGIN_REPLY packet, the server might immediately
						    * cut the connection; or it might allow the connected
						    * application further chances of login attempt. In the
						    * latter case, the connected application can initiate
						    * another login request by sending a new \c
						    * __AP_LOGIN_REQUEST packet.*/
#define __AP_LOGIN_REPLY_UNRECOGNIZEDACCOUNT	1  /*!< \brief Indicates that the \c CLNTNAME contained in the
						    * previously received \c __AP_LOGIN_REQUEST packet cannot
						    * can not be recognized by the PNSE server. Depending on
						    * the implementation of the PNSE server, after sending
						    * this \c __AP_LOGIN_REPLY packet, the server might
						    * immediately cut the connection; or it might allow the
						    * connected application further chances of login
						    * attempt. In the latter case, the connected application
						    * can initiate another login request by sending a new \c
						    * __AP_LOGIN_REQUEST packet.*/
#define __AP_LOGIN_REPLY_WRONGPASSWORD		2  /*!< \brief Indicates that in the previous \c
						    * __AP_LOGIN_REQUEST packet the connected application
						    * submitted a recognizable client name but along with a
						    * wrong password. Depending on the implementation of the
						    * PNSE server, after sending this \c __AP_LOGIN_REPLY
						    * packet, the server might immediately cut the connection;
						    * or it might allow the connected application further
						    * chances of login attempt. In the latter case, the
						    * connected application can initiate another login request
						    * by sending a new \c __AP_LOGIN_REQUEST packet.*/
#define __AP_LOGIN_REPLY_EXCEEDTRIALTIME	3  /*!< \brief Indicates that the connected application failed
						    * to successfully login to the server by providing correct
						    * pair of client name and password for certain amount of
						    * times. The maximal login trial times is determined by
						    * the implementation of the PNSE server program. After
						    * receiving the \c __AP_LOGIN_REPLY packet with this error
						    * code, the application will be disconnected immediately
						    * by the server and a new connection is required to be
						    * established for further login attempts.*/
#define __AP_LOGIN_REPLY_ALREADYLOGGEDIN	4  /*!< \brief Indicates that the connected application failed
						    * to successfully login to the server since the client has
						    * already loged in before. */
#define __AP_LOGIN_REPLY_OK			5  /*!< \brief Indicates that the connected application
						    * successfully logins to the PNSE server. This means that
						    * it can start to get services from the server after
						    * receiving this error code.*/

/* NSEP client state constants */
#define NSEP_CS_DISCONNECTED  0 /*!< \brief Indicates that the client has not finished the login procedure to
				 * the PNSE server.*/
#define NSEP_CS_NOTREADY      1 /*!< \brief Indicates that the client has already finished the login procedure
				 * to the PNSE server but is still not ready to take its responsibility to
				 * participate in the experiment or that the experiment has been startd but
				 * the client has withdrawn its participation in progressing the test for
				 * certain reasons and the withdrawal is temporary.*/
#define NSEP_CS_RUNNING       2 /*!< \brief Indicates that the client is in a state that it can participate in
				 * the collaborative experiment. For the CGM, this indicates that the CGM is
				 * either ready to or in a short later time will send the \c NSEP_CMD packet,
				 * or it has already sent the \c NSEP_CMD packet and is currently waiting for
				 * the \c NSEP_CSIG packet. For a FCM, this indicates that the FCM is either
				 * ready to receive a \c NSEP_CMD packet, or it is executing the command
				 * previously received and will send a \c NSEP_CSIG packet in a later time.*/
#define NSEP_CS_FINISHED      3 /*!< \brief indicates that the client has already stopped its participation in
				 * generating or executing the commands in the collaborative experiment.*/

/* Experimental State */
#define NSEP_ES_NOTREADY     0 /*!< \brief Initial value of the \c EXPSTATE parameter. It indicates that at
				* least one clients have the running state \c NSEP_CS_DISCONNECTED or \c
				* NSEP_CS_NOTREADY. When the experiment state is \c NSEP_ES_NOTREADY, it means
				* that the experiment is currently not running. That is, the experiment is
				* either "not started," or it is temporarily paused if it has already been
				* started.*/
#define NSEP_ES_RUNNING      1 /*!< \brief To indicate that the experiment is currently running. This implies
				* that all the clients have the running state \c NSEP_CS_RUNNING.*/
#define NSEP_ES_INTERRUPTED  2 /*!< \brief indicates that the experiemnt progress had been interrupted and the
				* test can not be resumed.*/
#define NSEP_ES_FINISHED     3 /*!< \brief Indicates that the experiment had been successfully completed.*/

#define NSEP_LOG              0 /*!< \brief Identifier to perform a log in action from a client. */
#define NSEP_SEND_CMD         1 /*!< \brief Identifier to send a \c NSEP_CMD packet. */
#define NSEP_REQUEST_CSIG     2 /*!< \brief Identifier to request a \c NSEP_CSIG packet. */
#define NSEP_REQUEST_CMD      3 /*!< \brief Identifier to request a \c NSEP_CMD packet. */
#define NSEP_SEND_CSIG        4 /*!< \brief Identifier to send a \c NSEP_CSIG packet.*/
#define NSEP_REQUEST_STATE    5 /*!< \brief Identifier to request a \c NSEP_EXPSTATE packet.*/
#define NSEP_SET_TO_FINISHED  6 /*!< \brief Identifier to set the status of a running client to \c
				   * NSEP_CS_FINISHED.*/

/**
 * \brief Main routine to communicate with PNSE using the NSEP protocol as a CGM or a FCM.
 *
 * This routine handles several aspects of the communication of a CGM or a FCM with the PNSE server. It covers
 * the aspects of login, send and receive messages and to define the running status of the CGM/FCM in a
 * compliant way with the NSEP protocol. This is achieved through the \c WhatToDo variable.
 *
 * \pre
 * - A PNSE type server must be waiting for connections.
 * - \c Remote must be properly initialised through Substructure_Remote_Init() and
 * Substructure_Remote_SetupClientSocket().
 *
 * \param[in]  Remote          Substructure of type \c REMOTE_NSEP.
 * \param[in]  WhatToDo        Action that will be performed.
 *                             - \c NSEP_LOG (0) Connect to the PNSE, login and set the state to running.
 *                             - \c NSEP_SEND_CMD (1) Send a \c NSEP_CMD package to the server.
 *                             - \c NSEP_REQUEST_CSIG (2) Request a \c NSEP_CSIG package from the server.
 *                             - \c NSEP_REQUEST_CMD (3) Request a \c NSEP_CMD package from the server.
 *                             - \c NSEP_SEND_CSIG (4) Send a \c NSEP_CSIG package to the server.
 *                             - \c NSEP_REQUEST_STATE (5) Ask the server fo the experiment state (\c
 *                               NSEP_EXPSTATE package).
 *                             - \c NSEP_SET_TO_FINISHED (6) Terminate the connection (\c NSEP_CLNSTATE packet
 *                               with \c NSEP_CS_FINISHED).
 * \param[in]  Time            Only used when \f$WhatToDo = 1\f$ or \f$WhatToDo = 3\f$. It contains the
 *                             pseudo-time (experiment time) when the data is sent by the CGM or received by
 *                             the FCM through a \c NSEP_CMD package.
 * \param[in]  Size            Referenced when \f$WhatToDo \in [1,4]\f$. It is the amount of data that has to
 *                             be sent or received from the PNSE server.
 * \param[in]  Data_To_Send    Only used when \f$WhatToDo = 1\f$ (\c NSEP_CMD) or \f$WhatToDo = 4\f$ (\c
 *                             NSEP_CSIG). It contains the data to be sent to the PNSE server.
 * \param[out] Data_To_Receive Only used when \f$WhatToDo = 2\f$ (\c NSEP_CSIG), \f$WhatToDo = 3\f$ (NSEP_CMD)
 *                             or \f$WhatToDo = 5\f$. On exit it contains the data received from the
 *                             server. This is not the case when \f$WhatToDo = 5\f$ since it will be used to
 *                             indicate the FCM that the process has finished <tt>Data_To_Receive[0] =
 *                             -9999.0</tt>.
 *
 * \post If any error occurs during the execution of this routine, like communication problems, the program
 * shows the appropiate error message befor exiting with \c EXIT_FAILURE.
 *
 * \sa Remote_t.
 */
void Substructure_Remote_NSEP( const Remote_t *const Remote, const int WhatToDo, const HYSL_FLOAT Time,
			       const unsigned int Size, const HYSL_FLOAT *const Data_To_Send, HYSL_FLOAT *const Data_To_Receive );

/**
 * \brief Login to the server in order to proceed with the test.
 *
 * Routine that deals with the login process required by the NSEP protocol. It makes use of the routines
 * Substructure_Remote_SendLoginInformation() and Substructure_Remote_ReceiveLoginInformation() to deal with
 * this process.
 *
 * \pre
 * - The TCP/IP connection with the PNSE server must be stablished.
 * - \c Remote must be properly initialised through Substructure_Remote_Init() and
 *   Substructure_Remote_SetupClientSocket().
 *
 * \param[in]  Remote         Substructure of type \c REMOTE_NSEP.
 * \param[in]  Send_Buffer    The message to be sent to the server. In this case a \c __AP_LOGIN_REQUEST. Not
 *                            referenced on entry, although it must contain an initialised memory space.
 * \param[out] Receive_Buffer On exit this will contain the answer comming from the server regarding the login
 *                            request.
 *
 * \post \c A __AP_LOGIN_REQUEST is sent to the server and afterwards a \c __AP_LOGIN_REPLY is received.
 *
 * \sa Remote_t.
 */
void Substructure_Remote_NSEP_LoginServer( const Remote_t *const Remote, const char *const Send_Buffer,
					   char *const Receive_Buffer );

/**
 * \brief Sends a \c __AP_LOGIN_REQUEST packet to the PNSE server.
 *
 * A \c __AP_LOGIN_REQUEST packet is sent to the server. It must be the next message sent by the client to the
 * PNSE server.  All the other messages will be rejected.
 *
 * \pre The TCP/IP connection with the PNSE server must be stablished.
 * 
 * \param[in] Socket      An open TCP/IP socket with the PNSE server.
 * \param[in] Name        Account name.
 * \param[in] Password    Account password.
 * \param[in] Send_Buffer The message to be sent to the server. In this case a \c __AP_LOGIN_REQUEST. Not
 *                        referenced on entry, although it must contain an initialised memory space.
 *
 * \post A \c __AP_LOGIN_REQUEST packet contaning the account name and password to login has been sent to the
 * PNSE server using TCP/IP.
 */
void Substructure_Remote_NSEP_SendLoginInformation( const int Socket, const char *const Name,
						    const char *const Password, const char *const Send_Buffer );

/**
 * \brief Receives a \c __AP_LOGIN_REPLY packet from the PNSE server.
 *
 * After sending a \c __AP_LOGIN_REQUEST to the server, the server ansers with a \c __AP_LOGIN_REPLY through
 * TCP/IP to inform the client about the request; whether it was successful or if there was a problem like
 * invalid account or password.
 *
 * \pre
 * - The TCP/IP connection with the PNSE server must be stablished.
 * - A \c __AP_LOGIN_REQUEST has to be already sent to the server.
 * 
 * \param[in]  Socket         An open TCP/IP socket with the PNSE server.
 * \param[out] NSEP_Type      Part of the received message that contains what is received. In this case it
 *                            should be \c __AP_LOGIN_REPLY.
 * \param[out] Error_Code     Part of the received message that contains the answer from the PNSE server. It
 *                            can be \c __AP_LOGIN_REPLY_NOTLOGINPACKET, \c
 *                            __AP_LOGIN_REPLY_UNRECOGNIZEDACCOUNT, \c __AP_LOGIN_REPLY_WRONGPASSWORD, \c
 *                            __AP_LOGIN_REPLY_EXCEEDTRIALTIME and \c __AP_LOGIN_REPLY_OK.
 * \param[out] Receive_Buffer On exit this will contain the answer comming from the server regarding the login request.

 *
 * \post \c A __AP_LOGIN_REPLY packet has been received from the server with an answer regarding the
 * previously sent \c __AP_LOGIN_REQUEST package.
 */
void Substructure_Remote_NSEP_ReceiveLoginInformation( const int Socket, int *const NSEP_Type,
						       int *const Error_Code, char *const Receive_Buffer );

/**
 * \brief Send a \c NSEP_CLNSTATE package to determine the current running state
 *
 * This routine sends a message to the PNSE server that contains the state of the client. That is, it sends a
 * NSEP_CLNSTATE package to the server.
 *
 * \pre 
 * - The TCP/IP connection with the PNSE server must be stablished.
 * - The Client must be successfully logged in.
 * 
 * \param[in] Socket      An open TCP/IP socket with the PNSE server.
 * \param[in] Run_State   It describes the running state of the client. It can have the values of \c
 *                        NSEP_CS_DISCONNECTED, \c NSEP_CS_NOTREADY, \c NSEP_CS_RUNNING or \c
 *                        NSEP_CS_FINISHED.
 * \param[in] Send_Buffer The message to be sent to the server. In this case a \c NSEP_CLNSTATE packet. Not
 *                        referenced on entry, although it must contain an initialised memory space.
 *
 * \post A \c NSEP_CLNSTATE package contaning the running state of the client has been sent to the PNSE server
 * using TCP/IP.
 */
void Substructure_Remote_NSEP_SetClientState( const int Socket, const int Run_State, char *const Send_Buffer );

#endif /* SUBSTRUCTURE_REMOTE_NSEP_H */
