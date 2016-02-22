#ifndef _SUBSTRUCTURE_STONEDRUMS_H_
#define _SUBSTRUCTURE_STONEDRUMS_H_

#include "Substructure.h"
#include "Definitions.h"

#define STONEDRUMS_NUMPARAM_INIT 1

typedef struct StoneDrums {

     int PrevDOF;
     
     HYSL_FLOAT Result_Fc;
     HYSL_FLOAT Disp0, DispT;
     HYSL_FLOAT Vel0, VelT, VelTdT;
     HYSL_FLOAT Acc0, AccT, AccTdT;
     HYSL_FLOAT a0, a2, a3, a6, a7;
     
     char *Description;  /*!< \brief Optional description of the substructure. */
} StoneDrums_t;

void Substructure_StoneDrums_Init (const int PrevDOF, const char *Description, StoneDrums_t *const Sub );
void Substructure_StoneDrums ( const HYSL_FLOAT DispTdT, const HYSL_FLOAT ramp, StoneDrums_t * const Sub, HYSL_FLOAT *const fc );
void Substructure_StoneDrums_Destroy( StoneDrums_t *const Sub );

#endif /* _SUBSTRUCTURE_STONEDRUMS_H_ */
