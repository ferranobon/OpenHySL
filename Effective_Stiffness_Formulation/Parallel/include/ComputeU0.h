/*
 * ComputeU0.h
 *
 *  Created on: 22/07/2011
 *      Author: ferran
 */

#ifndef COMPUTEU0_H_
#define COMPUTEU0_H_

#include <mpi.h>

#include "PMatrixVector.h"

void EffK_Calc_Effective_Force( PMatrixVector *const Mass, PMatrixVector *const Damp,
				PMatrixVector *const Disp, PMatrixVector *const Vel,
				PMatrixVector *const Acc, PMatrixVector *const Tempvec,
				const float a0, const float a1, const float a2,
				const float a3, const float a4, const float a5,
			       PMatrixVector *const Eff_Force );
void EffK_ComputeU0( PMatrixVector *const Eff_Force, PMatrixVector *const In_Load,
		     PMatrixVector *const Err_Force, const float PID_P, PMatrixVector *const Keinv,
		     PMatrixVector *const Tempvec, PMatrixVector *const Disp0 );
void CreateVectorXm( PMatrixVector *const VectorX, PMatrixVector *const VectorXm, const int PosCoupl, const int OrderC );
void CreateVectorXc( MPI_Comm Comm, PMatrixVector *const VecX, float *VecXc, int PosCoupl, int OrderC );

#endif /* COMPUTEU0_H_ */
