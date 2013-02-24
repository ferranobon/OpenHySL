/**
 * \file Print_Messages.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 19th of February 2013
 *
 * \brief Routines to print formated messages to the screen. Color definition.
 *
 * Functions to print formated messages to the standard ouptut, including stderr if the message in question
 * is an error. The ansi colors are also defined.
 */

#ifndef PRINT_MESSAGES_H_
#define PRINT_MESSAGES_H_

/**
 * \brief Supported types for using in conjunction with the PrintMessage() function.
 */
enum MyTypes {INT,    /*!< Integer type. */
	      UINT,   /*!< unsigned integer type. */
	      FLOAT,  /*!< float type. */
	      DOUBLE, /*!< double type. */
	      STRING  /*!< string type. */
};

#define ERROR   0  /*!< Error message type.*/
#define SUCCESS 1  /*!< Success message type.*/
#define INFO    2  /*!< Informative message type.*/

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
 * \brief Prints a formated message to the screen.
 *
 * The desired message is printed to stdout in the form. The function uses a variable number of
 * arguments in order to handle different situations. Also, the message is redirected to \c stdout
<<<<<<< HEAD
 * or \c stderr depending on the type of message: \c ERROR (stderr), \c WARNING (stderr), \c SUCCESS
 * (stdout) or \c INFO (stdout). Also, depending on the type of message, a prefix will be added to
 * the line, being:
 *
 * - "[FAILED] " if the specified type is \c ERROR in \c stderr.
 * - "[ WARN ] " if the specified type is \c WARNING in \c stderr.
=======
 * or \c stderr depending on the type of message: \c ERROR (stderr), \c SUCCESS (stdout) or \c INFO
 * (stdout). Also, depending on the type of message, a prefix will be added to the line, being:
 *
 * - "[FAILED] " if the specified type is \c ERROR in \c stderr.
>>>>>>> parent of 73b846c... The Print_Message() routine now can handle warnings. The routines ADwin_SaveData* make use of this additional message type.
 * - "[ INFO ] " if the specified type is \c INFO in \c stdout.
 * - "[  OK  ] " if the specified type is \c SUCCESS in \c stdout.
 *
 * In all cases a new line escape sequence is added at the end
 * of the message. For example. calling the function as:
 *
 * PrintMessage( INFO, 4, INT, 2, DOUBLE, 2.5, DOUBLE, 5.0, STRING, "This is a message." );
 *
 * will yield in:
 *
 * <tt>(stdout):  [ INFO ] 2 2.5 5.0 This is a message.<tt>
 *
 * \pre
 * - The first argument should be indicative of the type of the message: \c ERROR, \c WARNING, \c SUCCESS
 * or \c INFO.
 * - The second argument indicates the number of pairs of arguments that follow.
 * - The next pairs of arguments must be always in the format: TYPE, Value. Where TYPE must \b match one of
 * the entries of enum MyType and \b match the type of the next argument.
 *
 * \param[in] Mess_Type Type of the message to print.
 * \param[in] Num_Pairs Number of pairs of arguments.
 *
 * \post A formated message is print to stdout or stderr depending on the message type.
 *
 * \sa MyTypes.
 */
void Print_Message( const int Mess_Type, const int Num_Pairs, ... );

#endif /* PRINT_MESSAGES_H_ */
