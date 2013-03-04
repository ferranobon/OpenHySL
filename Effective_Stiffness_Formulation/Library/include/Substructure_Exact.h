#ifndef SUBSTRUCTURE_EXACT_H_
#define SUBSTRUCTURE_EXACT_H_

#include "Substructure.h"

#define EXACT_NUMPARAM_INIT 3  /*!< Number of required parameters in order to initialise a substructure of
			        * type Exact fbr */
typedef struct ExactSim {
     double Mass, Damp, Stiff;
     double Disp0, Vel0;
     double Disp, Vel;
     double Force_0, Force_1;
     double u0c_old;      /* Backup of the displacement at the coupling node */
     double v0c;          /* Velocity at the coupling node */

     double A,B,C,D,E,F,G,H;  /* Several constants */
     char *Description;  /*!< Optional description of the substructure. */
} ExactSim_t;


void Substructure_ExactSolution_Init( const double Mass, const double Damp, const double Stiff, const double DeltaT, const char *Description, ExactSim_t *const Num );
void Substructure_ExactSolution_SDOF( const double u0c, const double DeltaT, ExactSim_t *const Num, double *const fc );

#endif /* SUBSTRUCTURE_EXACT_H_ */
