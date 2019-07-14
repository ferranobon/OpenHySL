#ifndef PRECALCULATIONS_H_
#define PRECALCULATIONS_H_

#include <petscmat.h>

void ReadDataEarthquake_AbsValues( PetscScalar *Velocity, PetscScalar *Displacement, PetscInt NumSteps, const char *Filename );
void ReadDataEarthquake_RelValues( PetscScalar *Acceleration, PetscInt NumSteps, const char *Filename );
void Calc_Input_Load_AbsValues( Vec InLoad, Mat Stif, const Mat Damp, Vec D, Vec V );
void Calc_Input_Load_RelValues( Vec InLoad, Mat Mass, Vec A );
void Apply_LoadVectorForm ( Vec Vector, Vec LoadForm, PetscScalar Value );

#endif
