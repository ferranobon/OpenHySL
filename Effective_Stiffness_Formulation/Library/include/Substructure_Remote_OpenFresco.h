/**
 * \file Substructure_Remote_OpenFresco.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 26th of September 2013
 * 
 * \brief Prototypes for the communication subroutines compliant with OpenFresco.
 *
 * This file contains the prototypes of the communication subroutines that are used when connecting to an OpenFresco
 * server type.
 */

#ifndef SUBSTRUCTURE_REMOTE_OPENFRESCO_H
#define SUBSTRUCTURE_REMOTE_OPENFRESCO_H

#include "Definitions.h"

#define OF_INFO_DATA_LENGTH 11 /*!< Data length of the information message sent at the beginning of the
				* connection in order to establish the amount of data to be exchanged. */
#define OF_DATA_SIZE 256  /*!< Default amount of data exchanged between the OpenFresco client/server each
			   * time */

#define OF_REMOTE_OPEN                1 /*!< \brief Open the OpenFresco socket.*/
#define OF_REMOTE_SETUP               2 /*!< \brief Set-up the exchange of data*/
#define OF_REMOTE_SET_TRIAL_RESPONSE  3 /*!< \brief Send the trial response to the remote facility.*/
#define OF_REMOTE_EXECUTE             4 /*!< \brief Executes a command.*/
#define OF_REMOTE_COMMIT_STATE        5 /*!< \brief Commits a state.*/
#define OF_REMOTE_GET_DAQ_RESPONSE    6 /*!< \brief Get the DAQ response from the remote facility.*/
#define OF_REMOTE_GET_DISP            7 /*!< \brief Get the displacement from the remote facility.*/
#define OF_REMOTE_GET_VEL             8 /*!< \brief Get the velocity from the remote facility.*/
#define OF_REMOTE_GET_ACCEL           9 /*!< \brief Get the acceleration from the remote facility.*/
#define OF_REMOTE_GET_FORCE          10 /*!< \brief Get the restoring force from the remote facility.*/
#define OF_REMOTE_GET_TIME           11 /*!< \brief Get the elapsed time from the remote facility.*/
#define OF_REMOTE_GET_INITIAL_STIFF  12 /*!< \brief Get the initial stiffness matrix.*/
#define OF_REMOTE_GET_TANGENT_STIFF  13 /*!< \brief Get the tangent stiffness matrix.*/
#define OF_REMOTE_GET_DAMP           14 /*!< \brief Get the damping matrix.*/
#define OF_REMOTE_GET_MASS           15 /*!< \brief Get the mass matrix.*/
#define OF_REMOTE_FINISH_SETUP       98 /*!< \brief Finish the set-up of the controller (ECGEneric only).*/
#define OF_REMOTE_DIE                99 /*!< \brief Finish the test.*/


/**
 * \brief Main routine to communicate with OpenFresco remote facilities.
 *
 * This routine handles several aspects of the communication of a OpenFresco. It covers the aspects of set-up
 * and sending/receiving messages. This is achieved through the \c WhatToDo variable.
 *
 * \pre
 * - An OpenFresco server must be waiting for connections.
 * - \c Remote must be properly initialised through Substructure_Remote_Init() and
 * Substructure_Remote_SetupClientSocket().
 *
 * \param[in]  Socket          An open TCP/IP socket with the PNSE server.
 * \param[in]  WhatToDo        Action that will be performed.
 *                             - OF_REMOTE_SETUP (2) Setup the connection with the OpenFresco server.
 *                             - OF_REMOTE_SET_TRIAL_RESPONSE (3) send the new state to the OpenFresco server.
 *                             - OF_REMOTE_GET_DAQ_RESPONSE (6) to request the DAQ response from the remote
 *                               site.
 *                             - OF_REMOTE_DIE (99) to tell OpenFresco that the test has finished. Frees
 *                               internal memory.
 * \param[in]  Size            Referenced when \f$WhatToDo = \{2,3,6\}\f$. It is the number of coupling
 *                             degrees of freedom that the remote site must handle.
 * \param[in]  Data_To_Send    Only used when \f$WhatToDo = 3\f$. It contains the data to be sent to the
 *                             remote site.
 * \param[out] Data_To_Receive Only used when \f$WhatToDo = 6\f$. On exit it contains the data received from
 *                             the remote facility.
 *
 * \post If any error occurs during the execution of this routine, like communication problems, the program
 * shows the appropiate error message befor exiting with \c EXIT_FAILURE.
 */
void Substructure_Remote_OpenFresco( const int Socket, const int WhatToDo, const unsigned int Size, const HYSL_FLOAT *const Data_To_Send, HYSL_FLOAT *const Data_To_Receive );

#endif /* SUBSTRUCTURE_REMOTE_OPENFRESCO_H */
