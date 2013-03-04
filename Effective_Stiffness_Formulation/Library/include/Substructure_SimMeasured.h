#ifndef SUBSTRUCTURE_SIMMEASURED_H_
#define SUBSTRUCTURE_SIMMEASURED_H_

#include "Substructure.h"

typedef struct MeasuredSim {
     char Description[MAX_DESCRIPTION];
     double *Values;
     unsigned int Length;
} MeasuredSim_t;


void Substructure_SimMeasured_Init( const char *FileName, const unsigned int NSteps, const unsigned int NSubsteps, MeasuredSim_t *const Sub );
void Substructure_SimMeasured( const MeasuredSim_t *const Sub, double *const fc );

#endif /* SUBSTRUCTURE_SIMMEASURED_H_ */
