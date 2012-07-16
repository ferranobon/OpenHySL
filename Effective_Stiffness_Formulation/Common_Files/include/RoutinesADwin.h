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

#define NUM_CHANNELS 24  /* Number of channels to save from ADwin */

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
 * algorithm in ADwin considering displacement control.
 * In order to avoid changing each time the code in ADwin, the input and output
 * variables are always stored in the same space in the ADwin memory.
 *
 * \pre - ADwin must have, as a process running, the substepping part of Dorka's
 * substructure algorithm. The process in ADwin must look for the following
 * variables in:
 * - The vectors u0c, uc, fcprev and fc must be properly initialised and must be of length \f$L \geq OrderC\f$.
 *
 * \param[in] u0c      The implicit displacement resulting from the newmark integration process. This array
 *                       contains only the displacement of the coupling nodes.
 * \param[out] uc      Displacement containing the implicit and the explicit part at the final sub-step.
 * \param[out] fcprev  Measured force at the sub-step \f$ j \leq N_{sub\-step} - 1\f$. Used to compute the new
                         acceleration (Compute_Acceleration()) and velocity (Compute_Velocity()) in the main program. Only the coupling nodes are present.
 * \param[out] fc      Measured force at the final sub-step. Used in order to calculate the error force at the end
                         of the step (Compute_Force_Error()).
 * \param[in] OrderC   Order of the vectors. This is equal to the number of coupling degrees of freedom.
 *
 * \post ADwin will perform the substepping process in displacement control.
 */
void ADWIN_Substep( const float *const u0c, float *const uc, float *const fcprev, float *const fc, const int OrderC );

/**
 * \brief Reads the data from ADwin and stores it in a file.
 *
 * The data stored in the array 97 of ADwin is read and stored in a file named "Data.txt" with as many columns as
 * the number of channels (\c NUM_CHANNELS) and as many rows as \f$N_{steps}\cdot N_{substeps}\f$
 *
 * \pre Data must be a properly initialised array of length \f$L \geq N_{steps}\cdot N_{substeps}\cdot NUM\_CHANNELS\f$.
 *
 * \param[in] Num_Steps Number of steps. Variable involved in the amount of data to be read.
 * \param[in] Num_Sub Number of substeps. Variable involved in the amount of data to be read.
 * \param[out] Data Array to store the desired data
 *
 * \post Data.txt will contain the data stored in the array 97 of ADwin with:
 * - A header file with the names of the channels.
 * - Number of columns equal to \c NUM_CHANNELS.
 * - Number of rows equal to \f$N_{step}\cdot N_{substep}\f$.
 */
void GetDataADwin( const int Num_Steps, const int Num_Sub, float *const Data );

#endif /* ROUTINESADWIN_HPP */
