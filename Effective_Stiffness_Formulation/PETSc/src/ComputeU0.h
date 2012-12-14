#ifndef COMPUTEU0_H_
#define COMPUTEU0_H_

#include <petscmat.h>
#include "Initiation.h"

void EffK_Calc_Effective_Force( Mat Mass, Mat Damp, Vec Disp, Vec Vel, Vec Acc, Vec Tempvec, PetscScalar a0, PetscScalar a1, PetscScalar a2, PetscScalar a3, PetscScalar a4, PetscScalar a5, Vec Eff_Force );
void EffK_ComputeU0( Vec Eff_Force, Vec In_Load, Vec Err_Force, PetscScalar PID_P, Mat Keinv, Vec Tempvec, Vec Disp0 );

void CreateVectorXm( MPI_Comm Comm, Vec VecX, Vec VecXm, Coupling_Node *const CNodes );
void CreateVectorXc( MPI_Comm Comm, Vec VecX, PetscScalar *VecXc, Coupling_Node *const CNodes );
PetscInt GetOwner_Position_Vector( MPI_Comm Comm, Vec Vector, PetscInt Row );
PetscScalar* GetVec_Values( MPI_Comm Comm, Vec Vector, PetscInt NumRows, PetscInt *Rows );

#endif
