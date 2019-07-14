/**
 * \file Custom_Server.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 19th of June 2013
 *
 * \brief Routines that depend on the implementation of the custom server.
 */
#ifndef _CUSTOM_SEVER_H_
#define _CUSTOM_SEVER_H_

#include "Algorithm_Aux.h" /* For TIntegration_t */

/**
 * \brief Definition of constant values and filenames that will be used during the Custom-Server program execution.
 *
 * This routine reads the values specified in a configuration file, such as the order of the sub-matrices, the
 * number of sub-steps and file names, and are stored for further use.
 *
 *
 * \pre The values of the Newmark integration must be coherent/feasible. The routine will not perform checks
 * on them.
 *
 * \param[in]  FileName  Name of the configuration file.
 * \param[out] InitConst A structure that comprises of several constants.
 *
 * \post
 * - The size of the matrices will determine the memory that will be allocated when defining a \c
 *   MatrixVector_t type and also also how many elements will be read/written from/to the files.
 *
 * \sa AlgConst_t, PID_t, Rayleigh_t and TIntegration_t.
 */
void CustomServer_Init( const char *FileName, AlgConst_t *const InitConst );

/**
 * \brief Frees the memory allocated during the CustomServer_Init() routine.
 *
 * \pre InitConst must be properly initialised through CustomServer_Init().
 * 
 * \param[out] InitConst Structure containing the data to be deallocated.
 *
 * \sa AlgConst_t.
 */
void CustomServer_Destroy( AlgConst_t *const InitConst );

/**
 * \brief Prints a help text with the different options available when launching the program.
 *
 * \param[in] Program_Name Name of the executable.
 */
void CustomServer_PrintHelp( const char *Program_Name );

#endif /* _CUSTOM_SEVER_H_ */
