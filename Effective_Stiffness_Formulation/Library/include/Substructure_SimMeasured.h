/**
 * \file Substructure_SimMeasured.h
 * \author Ferran Obón Santacana
 * \version 1.0 
 * \date 20th of June 2013
 *
 * \brief Routines to deal with substructures with recorded values.
 *
 * These routines deal with the creation, use and destruction of those sub-structures that have the response
 * stored in a file. This file must contain, therefore, a time history of the coupling forces at the sub-step
 * level.
 */

#ifndef SUBSTRUCTURE_SIMMEASURED_H_
#define SUBSTRUCTURE_SIMMEASURED_H_

#include "Substructure.h"

/**
 * \brief Stores the information necessary to deal with substructures that come in form of a time history of
 * the coupling force.
 */
typedef struct MeasuredSim {
     char *Description;    /*!< \brief Description of the sub-structure.*/
     double *Values;       /*!< \brief Time history of coupling force values. */
     unsigned int Length;  /*!< \brief Length of the \c Values array.*/
} MeasuredSim_t;

/**
 * \brief Initialises the parameters of a substructure that comes as a time history of coupling forces
 * recorded at the sub-step level.
 *
 * \pre
 * - FileName must point to a file with tab or space separated values of the coupling force in a single
 *   column.
 * - The number of entries in the file should be \f$N_{entries} \geq N_{steps}\cdot N_{sub\-steps}\f$.
 * 
 * \param[in]  FileName    Name of the file where the time history is located.
 * \param[in]  NSteps      Number of steps of the test.
 * \param[in]  NSubsteps   Number of sub-steps in the sub-stepping process.
 * \param[in]  Description Brief description of the sub-structure.
 * \param[out] Sub         Substructure of type \c SIM_MEASURED.
 *
 * \post
 * - \c Sub.Length is equal to \f$N_{steps}\cdot N_{sub\-steps}\f$.
 * - \c Sub.Values contains the first \f$N_{steps}\cdot N_{sub\-steps}\f$ entries in the file.
 * - \c Sub.Description is a duplicate of the string in \c Description.
 * - If there is not enough memory, the program exits with \c EXIT_FAILURE.
 * - The memory should be deallocated through Substructure_DeleteCouplingNodes() routine.
 *
 * \sa MeasuredSim_t, MAX_DESCRIPTION, Substructure_Id.
 */
void Substructure_SimMeasured_Init( const char *FileName, const unsigned int NSteps, const unsigned int NSubsteps, const char *Description, MeasuredSim_t *const Sub );

/**
 * Returns the coupling force in the \a i-th iteration.
 *
 * This routine returns the value stored at the \a i-th position of the \c Sub.Value array. The position is
 * increased through a \c static variable in increments of 1.
 *
 * \pre \c Sub must be properly initialised through Substructure_SimMeasured_Init().
 *
 * \param[in]  Sub Sub-structure of type \c SIM_MEASURED.
 * \param[out] fc  Coupling force at the \a i-th postion of the array \c Sub.Value.
 *
 * \sa MeasuredSim_t, Substructure_SimMeasured_Init().
 */
void Substructure_SimMeasured( const MeasuredSim_t *const Sub, double *const fc );

/**
 * \brief Frees the memory dynamically allocated in Substructure_SimMeasured_Init().
 *
 * \pre \c Sub must be properly initialised through Substructure_SimMeasured_Init().
 *
 * \param[out] Sub Substructure of type \c SIM_MEASURED to be destroyed.
 *
 * \post
 * - \c Sub.Length is set to 0.
 * - \c The memory allocated in \c Sub.Description and \c Sub.Values is freed.
 *
 * \sa MeasuredSim_t, Substructure_SimMeasured_Init().
 */
void Substructure_SimMeasured_Destroy( MeasuredSim_t *const Sub );

#endif /* SUBSTRUCTURE_SIMMEASURED_H_ */
