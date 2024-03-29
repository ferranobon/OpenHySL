/**
 * \file Substructure_Experimental.h
 * \author Ferran Obón Santacana
 * \version 1.0 
 * \date 21st of June 2013
 *
 * \brief Routines to deal with experimental substructures.
 *
 * These routines deal with the creation and destriction of experimental substructures. They do not depend on
 * specific hardware (eg. ADwin).
 */
#ifndef SUBSTRUCTURE_EXPERIMENTAL_H_
#define SUBSTRUCTURE_EXPERIMENTAL_H_

#include "Substructure.h"

#define NUM_CHANNELS 19  /*!< \brief Number of channels of recorded data in ADwin */

/**
 * \brief Definition of a experimental substructure.
 */
typedef struct ExpSub {
     char *Description;  /*!< \brief Description of the experimental sub-structure.*/
} ExpSub_t;

/**
 * \brief Initialises the experimental substructure.
 * 
 * \param[in]  Description Brief description of the sub-structure.
 * \param[out] Sub         Experimental substructure of type \c EXP_ADWIN.
 *
 * \post
 * - \c Sub.Description is a duplicate of the string in \c Description.
 * - The memory should be deallocated through Substructure_DeleteCouplingNodes() routine.
 *
 * \sa ExpSub_t, MAX_DESCRIPTION, Substructure_Id.
 */ 
void Substructure_Experimental_Init( const char *Description, ExpSub_t *const Sub );

/**
 * \brief Frees the memory dynamically allocated in Substructure_Experimental_Init().
 *
 * \pre \c Sub must be properly initialised through Substructure_Experimental_Init().
 *
 * \param[out] Sub Experimental substructure.
 *
 * \post
 * - The memory allocated in \c Sub.Description is freed.
 *
 * \sa ExpSub_t and Substructure_Experimental_Init().
 */
void Substructure_Experimental_Destroy( ExpSub_t *const Sub );

#endif /* SUBSTRUCTURE_EXACT_H_ */
