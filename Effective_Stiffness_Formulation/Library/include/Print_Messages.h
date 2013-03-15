/**
 * \file Print_Messages.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 19th of February 2013
 *
 * \brief Routines to print formated messages to the screen. Color definition.
 *
 * Functions to print formated messages to the standard ouptut, including stderr if the message in question is
 * an error. The ansi colors are also defined.
 */

#ifndef PRINT_MESSAGES_H_
#define PRINT_MESSAGES_H_

#define ERROR   0                      /*!< Error message type.*/
#define SUCCESS 1                      /*!< Success message type.*/
#define INFO    2                      /*!< Informative message type.*/
#define WARNING 3                      /*!< Warning message type.*/

#define RESET   "\033[0m"              /*!< Reset to default lettering.*/
#define BLACK   "\033[30m"             /*!< Black color.*/
#define RED     "\033[31m"             /*!< Red color.*/
#define GREEN   "\033[32m"             /*!< Green color.*/
#define YELLOW  "\033[33m"             /*!< Yellow color.*/
#define BLUE    "\033[34m"             /*!< Blue color.*/
#define MAGENTA "\033[35m"             /*!< Magenta color.*/
#define CYAN    "\033[36m"             /*!< Cyan color.*/
#define WHITE   "\033[37m"             /*!< White color.*/

/* Bold colors */
#define BOLDBLACK   "\033[1m\033[30m"  /*!< Bold Black color.*/
#define BOLDRED     "\033[1m\033[31m"  /*!< Bold Red color.*/
#define BOLDGREEN   "\033[1m\033[32m"  /*!< Bold Green color.*/
#define BOLDYELLOW  "\033[1m\033[33m"  /*!< Bold Yellow color.*/
#define BOLDBLUE    "\033[1m\033[34m"  /*!< Bold Blue color.*/
#define BOLDMAGENTA "\033[1m\033[35m"  /*!< Bold Magenta color.*/
#define BOLDCYAN    "\033[1m\033[36m"  /*!< Bold Cyan color.*/
#define BOLDWHITE   "\033[1m\033[37m"  /*!< Bold White color.*/

/**
 * \brief Prints a formated header to the screen.
 *
 * The desired formated header is printed either to \c stdout or \c stderr depending on the type of message \c
 * ERROR (stderr), \c WARNING (stderr), \c SUCCESS (stdout) or \c INFO (stdout). The resulting output is:
 *
 * - "[FAILED] " if the specified type is \c ERROR in \c stderr.
 * - "[ WARN ] " if the specified type is \c WARNING in \c stderr.
 * - "[ INFO ] " if the specified type is \c INFO in \c stdout.
 * - "[  OK  ] " if the specified type is \c SUCCESS in \c stdout.
 *
 * \pre \c Mess_Type should be : \c ERROR, \c WARNING, \c SUCCESS or \c INFO.
 *
 * \param[in] Mess_Type Type of the message to print.
 *
 * \post A formated header is printed to \c stdout or \c stderr depending on the message type.
 *
 */
void Print_Header( const int Mess_Type );

#endif /* PRINT_MESSAGES_H_ */
