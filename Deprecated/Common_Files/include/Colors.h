/**
 * \file Colors.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 16th of January 2013
 *
 * \brief Colors for output
 * 
 * Definition of ANSI colors to be used by printf and other routines.
 */

#ifndef COLORS_H
#define COLORS_H

#define RESET   "\033[0m"              /*!< Reset to default lettering */
#define BLACK   "\033[30m"             /*!< Black color */
#define RED     "\033[31m"             /*!< Red color */
#define GREEN   "\033[32m"             /*!< Green color */
#define YELLOW  "\033[33m"             /*!< Yellow color */
#define BLUE    "\033[34m"             /*!< Blue color */
#define MAGENTA "\033[35m"             /*!< Magenta color */
#define CYAN    "\033[36m"             /*!< Cyan color */
#define WHITE   "\033[37m"             /*!< White color */

/* Bold colors */
#define BOLDBLACK   "\033[1m\033[30m"  /* Bold Black color */
#define BOLDRED     "\033[1m\033[31m"  /* Bold Red color */
#define BOLDGREEN   "\033[1m\033[32m"  /* Bold Green color */
#define BOLDYELLOW  "\033[1m\033[33m"  /* Bold Yellow color */
#define BOLDBLUE    "\033[1m\033[34m"  /* Bold Blue color */
#define BOLDMAGENTA "\033[1m\033[35m"  /* Bold Magenta color */
#define BOLDCYAN    "\033[1m\033[36m"  /* Bold Cyan color */
#define BOLDWHITE   "\033[1m\033[37m"  /* Bold White color */

#endif
