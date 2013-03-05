#ifndef SUBSTRUCTURE_UHYDEFBR_H_
#define SUBSTRUCTURE_UHYDEFBR_H_

#include "Substructure.h" /* For MAX_DESCRIPTION */

#define UHYDE_NUMPARAM_INIT 3 /*!< Number of required parameters in order to initialise a substructure of type
			       * UHYDE-\em fbr */
/**
 * \brief To b e used with UHYDE-\em fbr substructures.
 */
typedef struct UHYDEfbrSim {
     double u0c_old;                     /*!< Old displacement (the one from the previous sub-step. */
     double q;                           /*!< Displacement in the device. */
     double qyield;                      /*!< Yield displacement. */
     double qplastic;                    /*!< Plastic displacement. */
     double k;                           /*!< Initial stiffness. */
     char *Description;  /*!< Optional description of the substructure. */
} UHYDEfbrSim_t;

void Substructure_SimUHYDE_1D( const double u0c, const double DeltaT, UHYDEfbrSim_t *const Sub, double *const Friction_Force );
void Substructure_SimUHYDE_1D_Init( const double qyield, const double yield_factor, const double Friction, const char *Description, UHYDEfbrSim_t *const Sub );
void Substructure_SimUHYDE_Destroy( UHYDEfbrSim_t *const Sub );

#endif /* SUBSTRUCTURE_UHYDEFBR_H_ */
