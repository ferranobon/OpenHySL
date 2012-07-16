/**
 * \file ErrorHandling.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 10th of August 2011
 *
 * \brief Error handling prototypes
 *
 * ErrorHandling.h contains the prototypes for the functions that throw an error to the standard output \c stderr.
 * They serve as error handling and they are designed to be generic.
 */

#ifndef ERRORHANDLING_H_
#define ERRORHANDLING_H_

/**
 * \brief Prints the error message and exits the program.
 *
 * The desired error message is printed using the standard error output stderr. Afterwards the return signal
 * \c EXIT_FAILURE is thrown and the program exits abnormally using function exit in stdlib.h.
 *
 * \pre A failure condition in a routine is found.
 *
 * \param[in] ErrorMessage It contains the desired error message to display.
 *
 * \post The program exits with \c EXIT_FAILURE.
 */
void PrintErrorAndExit ( const char *ErrorMessage );

/**
 * \brief Prints an error message with a more detailed information.
 *
 * The desired error message is printed using the standard error output stderr with some additional information
 * for those routines supporting it. Afterwards the program exits with \c EXIT_FAILURE.
 *
 * \pre A failure condition in a routine that supports extended output is found.
 *
 * \param[in] ErrorMessage It contains the desired error message to display.
 * \param[in] Detail It contains additional information regarding the error message.
 *
 * \post The program exits with \c EXIT_FAILURE.
 */
void PrintErrorDetailAndExit( const char *ErrorMessage, const char *Detail );

/**
 * \brief Prints the error message due to an I/O exception and exits the program.
 *
 * The desired error message is printed using the standard error output stderr. Afterwards the return signal
 * \c EXIT_FAILURE is thrown and the program exits abnormally using function exit in stdlib.h.
 *
 * \pre A failure condition due to an error in file operations is found.
 *
 * \param[in] ErrorMessage It contains the desired error message to display.
 * \param[in] Filename The name of the file that has caused the error.
 *
 * \post The program exits with \c EXIT_FAILURE.
 */
void ErrorFileAndExit( const char *ErrorMessage, const char *Filename );

/**
 * \brief Prints the error message of LAPACK routines.
 *
 * The error message comming from a LAPACK routine is displayed using the standard error output stderr. Afterwards
 * the return signal \c EXIT_FAILURE is thrown and the program exits abnormally using function exit in stdlib.h.
 *
 * \pre A failure condition due to an invalid output of a LAPACK routine is found.
 *
 * \param[in] ErrorMessage1 It contains the first part of the error message,
 * which informs the failing routine.
 * \param[in] info The output of the LAPACK function.
 * \param[in] ErrorMessage2 It contains the second part of the error message,
 * and gives more information concerning the type of error.
 *
 * \post The program exits with \c EXIT_FAILURE.
 */
void LAPACKPErrorAndExit( const char *ErrorMessage1, int info, const char* ErrorMessage2 );

#endif /* ERRORHANDLING_H_ */
