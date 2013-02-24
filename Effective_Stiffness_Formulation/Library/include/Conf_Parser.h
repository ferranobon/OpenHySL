/**
 * \file Conf_Parser.h
 * \author Ferran Obón Santacana
 * \version 1.0 
 * \date 21st of February 2013
 * 
 * \brief Routines for parsing INI type configuration files.
 *
 * Routines for reading INI type configuration files. It supports configuration files with section names
 * and therefore it can distinguish two variables of identical name but that are placed in different sections.
 * The supported types are: integer, float, double and string. A typical configuration file is as follows:
 *
 * \warning The routines are case sensitive.
 *
 * \verbatim
 * # Lines that begin with # are treated as comments.
 * [Rayleigh]                       # A section is identified by a word between [ and ]
 * Alpha = 0.197                    # This is a floating point variable, either double or float. It can also be treated like a string.
 * Beta = 0.0000929
 *
 * # Newmark Alpha and Beta values
 * [Newmark]
 * Gamma = 0.5
 * Beta = 0.25                      # This variable is not the same as the one in Rayleigh section
 *
 * [FileNames]
 * Mass_Matrix = "504DOF Mass.txt"   # The double quotes are ignored and only the text in between is stored. This is the safest way to
 *                                   # input strings.
 * \endverbatim
 */

#ifndef _CONF_PARSER_H
#define _CONF_PARSER_H

#include <string.h> /* For size_t */

/**
 * \brief Structure to save the contents of a configuration file.
 */
typedef struct ConfigFile {
     size_t NumEntries;    /*!< \brief Number of entries stored. This variable acts as a counter. */
     size_t Length;        /*!< \brief Maximum number of entries allowed*/
     char **Keys;          /*!< \brief Stores the string on the left of the equal sign with the section
			    * name appended: "Section_Name:Key_Name" */
     char **Values;        /*!< \brief Stores the string on the right of the equal sign without any comments. */
} ConfFile_t;

#define CONFFILE_LEMPTY   0   /*!< \brief Line is empty */  
#define CONFFILE_LCOMMENT 1   /*!< \brief Line is a comment */
#define CONFFILE_LSECTION 2   /*!< \brief The line contains a section definition */
#define CONFFILE_LVALUE   3   /*!< \brief Line is a entry value */
#define CONFFILE_ERROR   -1   /*!< \brief This line is not valid */

#define MAX_LENGTH 256  /*!< \brief Maximum allowed length of a line. Increase this number for supporting longer lines. */

/**
 * \brief Allocates space for storing the entries in the configuration file.
 *
 * \param[in] Size Number of Keys and Values entries allowed.
 * \returns NULL if the allocation of memory failed. A pointer to the allocated
 * memory otherwise.
 *
 * \sa ConfFile.
 */
ConfFile* ConfFile_Create( const size_t Size );

/**
 * \brief Reads the contents of a configuration file.
 *
 * \pre \c CFile must be properly initialised through ConfFile_Create().
 *
 * A configuration file is read and stored in memory for further processing. It assumes a maximum lenght of a line to be \c MAX_LENGTH.
 * It evaluates, process and finally stores every line in the configuration file and it only stops until either the maximum number of
 * entries is reached or the end of file identifier is found.
 *
 * \param[in] FileName Name of the configuration file.
 * \param[in,out] CFile Stores the entries of the configuration file. On entry only the maximum number of entries, \c Length, is referenced.
 *
 * \sa ConfFile_t.
 */
int ConfFile_ReadFile( const char *FileName, ConfFile_t *const CFile );

/**
 * \brief Process the information of a line.
 * 
 * \pre The line must have been evaluated through ConfFile_EvalLine().
 *
 * The contents of a line are processed and stored accordingly. The string on the left of
 * the equal sign '\b =' is stored on the \c Key variable while the right string on the 
 * \c Value variable. If a comment sign, '\b #' or '\b ;', is found, the rest of the line
 * is ignored. If a entry contains no value, a 0 is assumend.
 *
 * \param[in] Line String to be processed.
 * \param[out] Key On exit contains the string on the right side of the equal sign '\b ='.
 * \param[out] Value On exit contains the string, until a comment sign is found if any, on the left
 * side of the equal sign '\b ='.
 * \returns 1 if the operation was successful, 0 otherwise.
 */
int ConfFile_ProcessLine( const char *Line, char *Key, char *Value );

/**
 * \brief Adds a new entry to the structure.
 *
 * \pre
 * - A new \c Key and \c Value obtained by evaluating and processing the information of a line
 * by ConfFile_EvalLine() and ConfFile_ProcessLine() respectively.
 * - \c Entry must be formated as 'Section_Name:Key_Name'.
 *
 * \param[in] Entry A string that will serve as an identifier for further referencing.
 * \param[in] Value The string on the right side of the equal sign, '\b =', without comments.
 * \param[in,out] CFile Configuration file. On entry only the \c NumEntries variable and the maximum
 * number of entries (\c Length) is referenced.
 *
 * \post \c Entry and \c Value strings are saved in the position \c i of the \c Keys and \c Values arrays
 * respectively.
 *
 * \sa ConfFile_t
 */
void ConfFile_SetEntry( const char *Entry, const char *Value, ConfFile_t *const CFile );

/**
 * \brief Identifies what kind of line has been read.
 * 
 * \param[in] Line String to be evaluated.
 * \returns
 * - \c CONFFILE_LEMPTY if the line is empty.
 * - \c CONFFILE_LCOMMENT it the line starts by '\b #' or '\b ;'. The line will be treated as a comment.
 * - \c CONFFILE_LSECTION if the line starts by '\b [' and ends by '\b ]'. This line contains a section
 * identifier.
 * - \c CONFFILE_LVALUE if none of the others. This is a key value line.
 */
short int ConfFile_EvalLine( const char *Line );

/**
 * \brief Returns the string identified by \c Key.
 *
 * \pre The configuration file must be properly initialised through ConfFile_Create().
 *
 * \param[in] CFile Structure with the entries of the configuration file.
 * \param[in] Key Entry to be found in \c CFile. It should be formated as "SectionName:KeyName".
 * \returns The string value assigned to \c Key in the configuration file.
 *
 * \sa ConfFile_t.
 */
char* ConfFile_GetString( const ConfFile_t *const CFile, const char *Key );

/**
 * \brief Returns the integer identified by \c Key.
 *
 * \pre The configuration file must be properly initialised through ConfFile_Create().
 *
 * \param[in] CFile Structure with the entries of the configuration file.
 * \param[in] Key Entry to be found in \c CFile. It should be formated as "SectionName:KeyName".
 * \returns The integer value assigned to the same position of the \c Key array in the configuration file.
 *
 * \sa ConfFile_t.
 */
int ConfFile_GetInt( const ConfFile_t *const CFile, const char *Key );

/**
 * \brief Returns the float identified by \c Key.
 *
 * \pre The configuration file must be properly initialised through ConfFile_Create().
 *
 * \param[in] CFile Structure with the entries of the configuration file.
 * \param[in] Key Entry to be found in \c CFile. It should be formated as "SectionName:KeyName".
 * \returns The float value assigned to the same position of the \c Key array in the configuration file.
 *
 * \sa ConfFile_t.
 */
float ConfFile_GetFloat( const ConfFile_t *const CFile, const char *Key );

/**
 * \brief Returns the double identified by \c Key.
 *
 * \pre The configuration file must be properly initialised through ConfFile_Create().
 *
 * \param[in] CFile Structure with the entries of the configuration file.
 * \param[in] Key Entry to be found in \c CFile. It should be formated as "SectionName:KeyName".
 * \returns The double value assigned to the same position of the \c Key array in the configuration file.
 *
 * \sa ConfFile_t.
 */
double ConfFile_GetDouble( const ConfFile_t *const CFile, const char *Key );

/**
 * \brief Finds the position within ConFile with a specified entry.
 *
 * \pre The configuration file must be properly initialised through ConfFile_Create().
 *
 * \param[in] CFile Structure with the entries of the configuration file.
 * \param[in] Key Entry to be found in \c CFile.
 * \returns The string value assigned to the same position of the \c Key array in the configuration file.
 *
 * \sa ConfFile_t.
 */
size_t ConfFile_FindPosition( const ConfFile_t *const CFile, const char *Key );

/**
 * \brief Deallocates the dynamically allocated memory in the ConfFile_t struct.
 *
 * \pre \c CFile must be properly initialised through ConfFile_Create().
 *
 * \param[out] CFile Configuration file to be feed.
 *
 * \sa ConfFile_t.
 */
void ConfFile_Destroy( ConfFile_t *const CFile );

#endif
