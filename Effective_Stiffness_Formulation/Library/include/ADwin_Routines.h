/**
 * \file ADwin_Routines.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 19th of February 2013
 * 
 * \brief Prototypes of the functions regarding ADwin communication routines.
 *
 * \warning The routines described in this page require the presence of the ADwin libraries and a properly
 * installed device to run smoothly.
 *
 * This file contains the prototypes of the routines that involve synchronisation, communication or process
 * management of ADwin.
 */

#ifndef ADWIN_ROUTINES_H_
#define ADWIN_ROUTINES_H_

/**
 * \brief Routine to boot ADwin.
 * 
 * \warning The graphical interface written in QT4 is recomended to manage ADwin. This function is deprecated.
 * 
 * \param[in] Device_Number Device number of ADwin.
 * \param[in] Boot_Path Full path of the file to be used in order to boot ADwin.
 *
 * \post ADwin is loaded and ready to be used. In any other case, an appropiate error is displayed.
 */
void ADwin_Boot( const int32_t Device_Number, const char *Boot_Path );

/**
 * \brief Checks the process status in ADwin.
 *
 * \warning The graphical interface written in QT4 is recomended to manage ADwin. This function is deprecated.
 *
 * \pre ADwin should be booted through the ADwin_Boot() routine.
 *
 * \param[in] ProcessNumber identifies the process to be checked.
 *
 * \post The status of the process is displayed and the program will terminate if the process is not running.
 *
 * \sa ADwin_Boot().
 */
void ADwin_CheckProcessStatus( const int ProcessNumber );

/**
 * \brief Loads, starts and stops an ADwin process.
 *
 * \warning The graphical interface written in QT4 is recomended to manage ADwin. This function is deprecated.
 *
 * \pre 
 * - ADwin should be booted through the ADwin_Boot() routine.
 * - The variable \c dowhat must fulfill the following expression \f$ 1 \leq dowhat \leq 3\f$.
 * - When loading a process, PName should point to the ADwin binary file specifically compiled for the ADwin
 *   processor that will be run on.
 *
 * Depending on the value of the \c dowhat variable the process will be loaded (\f$ dowhat = 1\f$), started
 * (\f$ dowhat = 2\f$) or stopped (\f$ dowhat = 3\f$). In the case of loading a process, only the \c PName
 * variable is referenced, while in the other cases only the \c PNum variable is used.
 *
 * \param[in] PName  Binary file to be loaded into ADwin. Only referenced if \f$dowhat = 1\f$.
 * \param[in] PNum   Process number to be started or stopped. Only referenced if \f$dowhat = 2\f$ or \f$dowhat
 *                   = 3\f$.
 * \param[in] dowhat Identifies the operation to be performed. Load (\f$ dowhat = 1\f$), started (\f$ dowhat =
 *                   2\f$) or stop (\f$ dowhat = 3\f$) a process.
 *
 * \post The specified operation is performed and checked if it has been successful.
 *
 * \sa ADwin_Boot().
 */
void ADwin_ManageProcess( const char* PName, const int PNum, const int dowhat );

/**
 * \brief Sends an array to ADwin.
 *
 * An array of the specified length is sent into the ADwin system and will overwrite any data that is present
 * into the pointed data array in ADwin.
 *
 * \pre
 * - ADwin must be properly booted.
 * - The array must be properly initialised and have as many elements as specified in \c Length.
 * - \c Index cannot exceed the number of data arrays in the ADwin system.
 *
 * \param[in] Index  Array number in ADwin.
 * \param[in] Array  It contains the data to be transfered (Array length \f$\geq Length\f$).
 * \param[in] Length Number of elements to be transfered.
 *
 * \post The number of variables contained in Array and transfered to ADwin is equal to \c Length.
 */
void ADwin_SendArray( const unsigned int Index, const double *const Array, const unsigned int Length );

/**
 * \brief Sub-stepping process in ADwin.
 *
 * This routine communicates with ADwin to perform the sub-stepping process in ADwin. It does so by sending
 * the explicit displacement to the ADwin system. It waits a certain amount of time before retrieving the data
 * from the ADwin system. This is done in order to avoid overloading the system with network requests when the
 * system is not yet ready. A good time to get the data from ADwin is:
 *
 * \f[Time_To_Wait = \frac{Round Trip Time}{0.9}\f]
 *
 * Note that this value has to be coherent also with when ADwin measures the last coupling force from the
 * speciment.
 *
 * This routine relies on a synchronization variable when communicating with ADwin.
 * - When sending the data to ADwin, the first value of the array should be \c 1.0. This indicates the running
 *   process in ADwin that a new displacement value is available.
 * - If \f$ADwinReady = 0\f$ the sub-stepping process is not finished and therefore the last coupling force is
 *   not yet measured. The routine will keep asking ADwin for the last coupling force until the first received
 *   value from ADwin is equal to -1.0. When this happens the substepping process is finished.

 *
 * \pre
 * - ADwin must be successfully booted and the substepping process must be running.
 * - The sub-stepping process in ADwin should expect the new displacement vectors to be set in array 2. The
 *   first value should be a 1.0 in order to indicate that a new value is available.
 * - The sub-stepping process in ADwin should save the new displacement and force values on array 3. When the
 *   values are reaady (the sub-stepping process is finished) the first value should be a -1.0.
 * - A valid gain matrix of size \$OrderC\$x\$OrderC\$ should be on ADwin.
 * - Only the coupling elements should be on the given vectors since the routine does not consider the global
 *   position of the coupling degrees of freedom.
 *
 * \param[in]  VecTdT_0c    Explicit displacement at the begining of the step to be sent to ADwin. It should
 *                          be at least \f$Length \geq OrderC\f$.
 * \param[in]  OrderC       Amount of data to be sent to the ADwin system.
 * \param[in]  Time_To_Wait Amount of time that the computer should wait before asking ADwin for the first
 *                          time for the new values.
 * \param[out] VecTdT_c     Displacement at the end of the sub-stepping process. It should be at least
 *                          \f$Length \geq OrderC\f$.
 * \param[out] fcprev_c     Coupling forces measured at the sub-step \f$n_{sub} - 1\f$, where \f$n_{sub}\f$ is
  *                         the total number of sub-steps. It should be at least \f$Length \geq OrderC\f$.
 * \param[out] fc_c         Coupling forces measured at the sub-step \f$n_{sub}\f$. It should be at least
 *                          \f$Length \geq OrderC\f$.
 * 
 * \post The vectors \c VecTdT_c, \c fcprev_c, \c fc_c contain the last control values and the coupling force
 * at the previous and last sub-steps.
 */
void ADwin_Substep( const double *const VecTdT_0c, const unsigned int OrderC, const double Time_To_Wait, double *const VecTdT_c,
		    double *const fcprev_c, double *const fc_c );

/**
 * \brief Checks if the loaded driver matches the ADwin system.
 * 
 * \warning The graphical interface written in QT4 is recomended to manage ADwin. This function is deprecated.
 *
 * \pre ADwin must be booted.
 */
void ADwin_TestVersion( void );

/**
 * \brief Reads the data from ADwin and stores it in a file.
 *
 * The data stored in the specified array ADwin is read and stored into the desired ASCII file with as many
 * columns as the specified number of chanels and as many rows as \f$N_{steps}\cdot N_{substeps}\f$.
 *
 * \pre
 * - ADwin must be properly booted and the substructure test must be finished.
 * - Data must be a properly initialised array of length \f$L \geq N_{steps}\cdot N_{substeps}\cdot
 *   N_{channels}\f$.
 *
 * \param[in] FileName Name of the output file.
 * \param[in] Num_Steps Number of steps.
 * \param[in] Num_Sub Number of sub-steps.
 * \param[in] Num_Channels Number of data acquisition channels.
 * \param[in] Chan_Names List of channel names.
 * \param[in] DataIndex points out which data array is to be accessed in ADwin.
 *
 * \post FileName will contain the data stored in the array identified by \c Data_Index in ADwin with:
 * - A row containing the name of the channels
 * - Number of columns equal to the number of channels.
 * - Number of rows equal to \f$N_{step}\cdot N_{substep}\f$.
 */
void ADwin_SaveData_ASCII( const char *FileName, const unsigned int Num_Steps, const unsigned int Num_Sub,
			   const unsigned short int Num_Channels, const char **Chan_Names,
			   const int DataIndex );


#endif /* ADWIN_ROUTINES_H_ */

