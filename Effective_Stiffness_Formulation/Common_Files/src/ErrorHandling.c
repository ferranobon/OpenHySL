/**
 * \file ErrorHandling.c
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 10th of August 2011
 *
 * \brief Source code of Error handling functions
 *
 * ErrorHandling.c contains the source code of the functions that throw an error to the standard output \c stderr.
 * They serve as error handling and they are designed to be generic.
 */


#include <stdio.h>  /* For fprintf and stderr */
#include <stdlib.h>  /* For exit( ) */

#include "ErrorHandling.h" /* Prototypes of the Error Handling functions.*/

void PrintErrorAndExit ( const char *ErrorMessage )
{
	fprintf( stderr, "%s.\n", ErrorMessage );
	fprintf( stderr, "Exiting program." );
	exit( EXIT_FAILURE );
}

void ErrorFileAndExit( const char *ErrorMessage, const char *Filename )
{

	fprintf( stderr, "%s%s.\n", ErrorMessage, Filename );
	fprintf( stderr, "Exiting program." );
	exit( EXIT_FAILURE );
}

void LAPACKPErrorAndExit( const char *ErrorMessage1, int info, const char* ErrorMessage2 )
{
	fprintf( stderr, "%s%d%s.\n", ErrorMessage1, info, ErrorMessage2 );
	exit( EXIT_FAILURE );
}
