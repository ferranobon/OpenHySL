#ifndef _SUBSTRUCTURE_STONEDRUMS_H_
#define _SUBSTRUCTURE_STONEDRUMS_H_

#include "Substructure_BoucWen.h"
#include "Definitions.h"

#define STONEDRUMS_NUMPARAM_INIT 40

typedef struct StoneDrums {
    
     char *Description;  /*!< \brief Optional description of the substructure. */
     int JointID, PrevJointID, NextJointID, Dir;

     BoucWenSurface_t BoucWen_Sliding;
     BoucWen_t BoucWen_Torsion;
     BoucWenSurface_t BoucWen_Bending;
} StoneDrums_t;

void Substructure_StoneDrums_Init (const int JointID, const int PrevJointID, const int NextJointID,
				   const int Dir, const HYSL_FLOAT *const BW,
				   const int *const BoucWen_Type, const char *Description,
				   StoneDrums_t *const Sub);

void Substructure_StoneDrums ( const HYSL_FLOAT *const DispTdT, const int Order, StoneDrums_t *const Sub, HYSL_FLOAT *force );

void Substructure_StoneDrums_Destroy( StoneDrums_t *const Sub );

#endif /* _SUBSTRUCTURE_STONEDRUMS_H_ */
