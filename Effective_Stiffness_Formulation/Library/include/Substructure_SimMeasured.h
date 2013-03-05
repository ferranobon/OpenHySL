#ifndef SUBSTRUCTURE_SIMMEASURED_H_
#define SUBSTRUCTURE_SIMMEASURED_H_

#include "Substructure.h"

typedef struct MeasuredSim {
     char *Description;
     double *Values;
     unsigned int Length;
} MeasuredSim_t;


void Substructure_SimMeasured_Init( const char *FileName, const unsigned int NSteps, const unsigned int NSubsteps, const char *Description, MeasuredSim_t *const Sub );
void Substructure_SimMeasured( const MeasuredSim_t *const Sub, double *const fc );
void Substructure_SimMeasured_Destroy( MeasuredSim_t *const Sub );
#endif /* SUBSTRUCTURE_SIMMEASURED_H_ */
