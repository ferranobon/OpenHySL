/*
 * Precalculations.h
 *
 *  Created on: 22/07/2011
 *      Author: ferran
 */

#ifndef PRECALCULATIONS_H_
#define PRECALCULATIONS_H_

#include <mpi.h>

#include "PMatrixVector.h"

void ReadDataEarthquake_AbsValues( float *Velocity, float *Displacement, const unsigned int NumSteps, const char *Filename );
void ReadDataEarthquake_RelValues( float *Acceleration, const unsigned int NumSteps, const char *Filename );
void CopyDiagonalValues( MPI_Comm Comm, PMatrixVector *const Mat, PMatrixVector *const Vec );
void Calc_Input_Load_AbsValues( PMatrixVector *const InLoad, const PMatrixVector *const Stif, const PMatrixVector *const Damp, const PMatrixVector *const D, const PMatrixVector *const V );
void Calc_Input_Load_RelValues( PMatrixVector *const InLoad, const PMatrixVector *const Mass, const PMatrixVector *const A );
void Apply_LoadVectorForm( PMatrixVector *const Vector, const PMatrixVector *const LoadForm, const float Value );

#endif /* PRECALCULATIONS_H_ */
