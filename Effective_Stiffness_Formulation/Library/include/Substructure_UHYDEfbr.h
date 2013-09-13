/**
 * \file Substructure_UHYDEfbr.h
 * \author Ferran Obón Santacana
 * \version 1.0 
 * \date 13th of September 2013
 *
 * \brief Routines to deal with substructures that represent the UHYDE-\em fbr device.
 *
 * These routines deal with the creation, destruction and simulation of those substructures that represent the
 * UHYDE-\em fbr device.
 */

#ifndef SUBSTRUCTURE_UHYDEFBR_H_
#define SUBSTRUCTURE_UHYDEFBR_H_

#include "Substructure.h" /* For MAX_DESCRIPTION */

#define UHYDE_NUMPARAM_INIT 3 /*!< Number of required parameters in order to initialise a substructure of type
			       * UHYDE-\em fbr */
/**
 * \brief To be used with UHYDE-\em fbr substructures and routines.
 */
typedef struct UHYDEfbrSim {
     double u0c_old;    /*!< Old displacement (the one from the previous sub-step. */
     double q;          /*!< Displacement in the device. */
     double qyield;     /*!< Yield displacement. */
     double qplastic;   /*!< Plastic displacement. */
     double k;          /*!< Initial stiffness. */
     char *Description; /*!< Optional description of the substructure. */
} UHYDEfbrSim_t;

void Substructure_SimUHYDE_1D( const double u0c, const double DeltaT, UHYDEfbrSim_t *const Sub, double *const Friction_Force );
void Substructure_SimUHYDE_1D_Init( const double qyield, const double yield_factor, const double Friction, const char *Description, UHYDEfbrSim_t *const Sub );

/**
 * \brief Frees the memory dynamically allocated in Substructure_SimUHYDE_1D_Init().
 *
 * \pre \c Sub must be properly initialised through Substructure_SimUHYDE_1D_Init().
 *
 * \param[in,out] Sub Substructure representing a UHYDE-\em fbr device.
 *
 * \post
 * - The memory allocated in \c Sub.Description is freed.
 *
 * \sa UHYDEfbrSim_t and Substructure_SimUHYDE_1D_Init().
 */
void Substructure_SimUHYDE_Destroy( UHYDEfbrSim_t *const Sub );

#endif /* SUBSTRUCTURE_UHYDEFBR_H_ */
