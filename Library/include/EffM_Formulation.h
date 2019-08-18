#ifndef EFFK_FORMULATION_H_
#define EFFK_FORMULATION_H_

#include "MatrixVector.h"
#include "MatrixVector_Sp.h"
#include "MatrixVector_MPI.h"

void EffM_EffectiveForce( const MatrixVector_t *const Damp, const MatrixVector_t *const Stiff,
			  const MatrixVector_t *const DispT, const MatrixVector_t *const VelT,
			  const MatrixVector_t *const AccT, MatrixVector_t *const Tempvec,
			  const hysl_float_t a6, const hysl_float_t a9, const hysl_float_t a10,
			  MatrixVector_t *const Eff_ForceT );
void EffM_EffectiveForce_PS( const MatrixVector_t *const Stiff, const MatrixVector_t *const Damp,
			     const MatrixVector_t *const DispT, const MatrixVector_t *const VelT,
			     const MatrixVector_t *const AccT, MatrixVector_t *const Tempvec,
			     const hysl_float_t a6, const hysl_float_t a9, const hysl_float_t a10,
			     MatrixVector_t *const Eff_ForceT );

void EffM_ComputeDisplacement( const MatrixVector_t *const DispT, const MatrixVector_t *const VelT,
			       const MatrixVector_t *const AccT, const MatrixVector_t *const AccTdT,
			       const hysl_float_t a8, const hysl_float_t a9, const hysl_float_t a10, MatrixVector_t *const DispTdT );


void EffM_ComputeVelocity( const MatrixVector_t *const VelT, const MatrixVector_t *const AccT,
			   const MatrixVector_t *const AccTdT, const hysl_float_t a6, const hysl_float_t a7,
			   MatrixVector_t *const VelTdT );


void EffM_EffectiveForce_Sp( const MatrixVector_Sp_t *const Stiff, const MatrixVector_Sp_t *const Damp, const MatrixVector_t *const DispT,
			     const MatrixVector_t *const VelT, const MatrixVector_t *const AccT, MatrixVector_t *const Tempvec,
			     const hysl_float_t a6, const hysl_float_t a9, const hysl_float_t a10,
			     MatrixVector_t *const Eff_ForceT );

void EffM_EffectiveForce_MPI( PMatrixVector_t *const Stiff, PMatrixVector_t *const Damp,
			      PMatrixVector_t *const DispT, PMatrixVector_t *const VelT,
			      PMatrixVector_t *const AccT, PMatrixVector_t *const Tempvec,
			      const hysl_float_t a6, const hysl_float_t a9, const hysl_float_t a10,
			      PMatrixVector_t *const Eff_ForceT );

void EffM_ComputeDisplacement_MPI( PMatrixVector_t *const DispT, PMatrixVector_t *const VelT,
				   PMatrixVector_t *const AccT, PMatrixVector_t *const AccTdT,
				   const hysl_float_t a8, const hysl_float_t a9, const hysl_float_t a10,
				   PMatrixVector_t *const DispTdT );

void EffM_ComputeVelocity_MPI( PMatrixVector_t *const VelT, PMatrixVector_t *const AccT,
			       PMatrixVector_t *const AccTdT, const hysl_float_t a6, const hysl_float_t a7,
			       PMatrixVector_t *const VelTdT );

#endif /* EFFM_FORMULATION_H_ */
