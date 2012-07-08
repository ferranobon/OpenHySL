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
 * The desired error message is printed using the standard error output std:cerr. Afterwards the return signal
 * \c EXIT_FAILURE is thrown and the program exits abnormally using function exit in stdlib.h.
 *
 * \pre An failure condition in a routine is found.
 *
 * \param[in] ErrorMessage an argument of type const char*. It contains the desired error message to display.
 *
 * \post The program must exit with \c EXIT_FAILURE.
 */
void PrintErrorAndExit ( const char *ErrorMessage );

void PrintErrorDetailAndExit( const char *ErrorMessage, const char *Detail );

/**
 * \brief Prints the error message due to an I/O exception and exits the program.
 *
 * The desired error message is printed using the standard error output std:cerr. Afterwards the return signal
 * \c EXIT_FAILURE is thrown and the program exits abnormally using function exit in stdlib.h.
 *
 * \pre An failure condition due to an error in file operations is found.
 *
 * \param[in] ErrorMessage an argument of type const char*. It contains the desired error message to display.
 * \param[in] Filename an argument of type const char* The name of the file with errors.
 *
 * \post The program must exit with \c EXIT_FAILURE.
 */
void ErrorFileAndExit( const char *ErrorMessage, const char *Filename );

/**
 * \brief Prints the error message of LAPACK routines.
 *
 * The error message comming from a LAPACK routine is displayed using the standard error output std:cerr. Afterwards
 * the return signal \c EXIT_FAILURE is thrown and the program exits abnormally using function exit in stdlib.h.
 *
 * \pre An failure condition due to an invalid output of a LAPACK routine is found.
 *
 * \param[in] ErrorMessage1 an argument of type const char*. It contains the first part of the error message,
 * which informs the failing routine.
 * \param[in] info an argument of type integer and the output of the LAPACK function.
 * \param[in] ErrorMessage2 an argument of type const char*. It contains the second part of the error message,
 * and gives more information concerning the type of error.
 *
 * \post The program must exit with \c EXIT_FAILURE.
 */
void LAPACKPErrorAndExit( const char *ErrorMessage1, int info, const char* ErrorMessage2 );

#endif /* ERRORHANDLING_H_ */
