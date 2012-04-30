/**
 * \file RoutinesADwin.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 3rd of October 2011
 * 
 * \brief Prototypes of the functions regarding ADwin communication routines
 *
 * This file contains the prototypes of the routines that involve synchronisation
 * communication or process management of ADwin. 
 */

#ifndef ROUTINESADWIN_HPP_
#define ROUTINESADWIN_HPP_

/**
 * \brief Routine to boot ADwin.
 * 
 * This routine calls the ADwin routine Boot, using the Boot() routine specified
 * in the ADwin documentation.
 *
 * \pre The ADwin device and libraries must be properly installed in the system.
 * 
 * \param[in] Device_Number Device number of ADwin.
 * \param[in] Boot_Path Full path of the file to be used in order to boot ADwin.
 * \returns 1 if ADwin is loaded correctly and 0 otherwise.
 *
 * \post ADwin is loaded and ready to be used. In any other case, an appropiate
 * error is displayed.
 */
int BootADWIN( const int Device_Number, const char *Boot_Path );
int ADWIN_TestVersion( void );
void ADWIN_ManageProcess( const char* PName, const int PNum, const int dowhat );
void ADWIN_CheckProcessStatus( int ProcessNumber );

/**
 * Sets the Matrix Gc into ADwin
 *
 * \param[in] Gc
 * \param[in] length
 */
void ADWIN_SetGc( const float *const Gc, const int length );

/**
 * \brief Routine to perform the substepping process of Dorka's substructure
 * algorithm in ADwin.
 *
 * This routine performs the substepping process of Dorka's substructure
 * algorithm in ADwin. In this routine, displacement control is considered.
 * In order to avoid changing each time the code in ADwin, the input and output
 * variables are always stored in the same space in the ADwin memory.
 *
 * \pre ADwin must have, as a process running, the substepping part of Dorka's
 * substructure algorithm. The process in ADwin must look for the following
 * variables in:
 *
 * \param[in] u0c
 * \param[out] uc
 * \param[out] fcprev
 * \param[out] fc
 * \param[in] OrderC
 * \param[in] NSub Number of sub-steps in the sub-stepping process.
 * \param[in] Deltat_Sub 
 *
 * \post ADwin will perform the substepping process in displacement control.
 */
void ADWIN_Substep( const float *const u0c, float *const uc, float *const fcprev, float *const fc, const int OrderC, const int NSub, const float Deltat_Sub );

void GetDataADwin( const int Num_Steps, const int Num_Sub, float *const Data );

#endif /* ROUTINESADWIN_HPP */
