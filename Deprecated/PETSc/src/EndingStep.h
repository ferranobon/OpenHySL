#ifndef ENDINGSTEP_H_
#define ENDINGSTEP_H_

#include <petscmat.h>
#include "Initiation.h"

void JoinNonCouplingPart( MPI_Comm Comm, Vec VecXm, Mat MatrixXm, Vec fcprevsub, Vec VecX, Coupling_Node *const CNodes );
void Compute_Acceleration( Vec DispTdT, Vec DispT, Vec VelT, Vec AccT,
			   PetscScalar a0, PetscScalar a2, PetscScalar a3,
			   Vec AccTdT );
void Compute_Velocity( Vec VelT, Vec AccT, Vec AccTdT,
		       PetscScalar a6, PetscScalar a7, 
		       Vec VelTdT );
void Compute_Force_Error( Mat Mass, Mat Damp, Mat Stiff,
			  Vec AccTdT, Vec VelTdT, Vec DispTdT,
			  Vec fc, Vec LoadTdT, Vec fu );

#endif
