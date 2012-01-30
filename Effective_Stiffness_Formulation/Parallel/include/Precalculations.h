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

void ReadDataEarthquake( float *displacement, float *velocity, float *acceleration, const int NumSteps, const char *Filename );
void CopyDiagonalValues( MPI_Comm Comm, PMatrixVector *const Mat, PMatrixVector *const Vec );
void Calc_Input_Load( PMatrixVector *const InLoad, PMatrixVector *const Stif, PMatrixVector *const Damp, PMatrixVector *const Mass, PMatrixVector *const DiagM, PMatrixVector *const D, PMatrixVector *const V, PMatrixVector *const A );

#endif /* PRECALCULATIONS_H_ */
