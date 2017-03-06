/**
 * \file New_State.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 9th of February 2013
 *
 * \todo Add single precision routines.
 * 
 * \brief Routines for computing the new state.
 *
 * Routines for calculating the new state. The routines make use of the BLAS library to perform the linear
 * algebra operations and they support both single and double precision. Both general and packed storage
 * versions are available.
 */

#ifndef NEW_STATE_H_
#define NEW_STATE_H_

#include "MatrixVector.h"
#include "MatrixVector_Sp.h"
#include "MatrixVector_MPI.h"

/**
 * \brief Calculates the new state. General storage version.
 *
 * The explicit part of the new state is computed. This routine does not depend on a specific formulation and
 * therefore it can be used in displacement, velocity and acceleration control formulations.
 * 
 * \f[\vec n_0^{t + \Delta t} = \mathcal{G}^{-1}(\vec f_{eff}^t + \vec l_i^t - \vec f_{err}^t)\f]
 *
 * where:
 * - \f$\vec n_0^{t + \Delta t}\f$ is the explicit part of the new state vector at time
 *   \f$t + \Delta t\f$,
 * - \f$\mathcal{G}\f$ is the gain matrix,
 * - \f$\vec f_{eff}^t\f$ is the effective force vector at time \f$t\f$,
 * - \f$\vec l_i^t\f$ is the input load vector at time \f$t\f$,
 * - and \f$\vec f_{err}^t\f$ is the error compensation force at time \f$t\f$.
 *
 * It makes use of BLAS routines to perform the lineal algebra operations. For the packed storage version, the
 * routine Compute_NewState_PS() should be used instead.
 * 
 * \pre
 * - All elements of type \c MatrixVector_t must be properly initialised through the  MatrixVector_Create()
 *   routine.
 * - \c IGain must be a symmetrical matrix in general storage. Only the upper part will be referenced (lower
 *   part in FORTRAN).
 * - The size of the vectors and matrices must be coherent since it will not be checked in the routine.
 * 
 * \param[in]     IGain      The inverted gain matrix \f$\mathcal{G}^{1}\f$.
 * \param[in]     Eff_ForceT The effective force vector \f$\vec f_{eff}^t\f$.
 * \param[in]     In_LoadT   The input load vector \f$\vec l_i^t\f$.
 * \param[in]     Err_ForceT The error force vector \f$\vec f_{err}^t\f$.
 * \param[in,out] Tempvec    Temporal vector that helps in the calculations. This is included in order to
 *                           avoid allocating and deallocating memory each step. On input only the number of
 *                           rows is used.
 * \param[out]    VecTdT_0   New explicit state (displacement, velocity or acceleration depending on the
 *                           selected formulation) \f$\vec n_0^{t + \Delta t}\f$.
 *
 * \post
 * - \c VecTdT_0 is the result of the operation:
 *
 * \f[\vec n_0^{t + \Delta t} = \mathcal{G}^{-1}(\vec f_{eff}^t + \vec l_i^t + \vec f_{err}^t)\f]
 *
 * \sa MatrixVector_t and EffK_EffectiveForce().
 */
void Compute_NewState( const MatrixVector_t *const IGain, const MatrixVector_t *const Eff_ForceT,
		       const MatrixVector_t *const In_LoadT, const MatrixVector_t *const Err_ForceT,
		       MatrixVector_t *const Tempvec, MatrixVector_t *const VecTdT_0 );

/**
 * \brief Calculates the new state. Packed storage version.
 *
 * The explicit part of the new state is computed. This routine does not depend on a specific formulation and
 * therefore it can be used in displacement, velocity and acceleration control formulations.
 * 
 * \f[\vec n_0^{t + \Delta t} = \mathcal{G}^{-1}(\vec f_{eff}^t + \vec l_i^t - \vec f_{err}^t)\f]
 *
 * where:
 * - \f$\vec n_0^{t + \Delta t}\f$ is the explicit part of the new state vector at time
 *   \f$t + \Delta t\f$,
 * - \f$\mathcal{G}\f$ is the gain matrix,
 * - \f$\vec f_{eff}^t\f$ is the effective force vector at time \f$t\f$,
 * - \f$\vec l_i^t\f$ is the input load vector at time \f$t\f$,
 * - and \f$\vec f_{err}^t\f$ is the error compensation force at time \f$t\f$.
 *
 * It makes use of BLAS routines to perform the lineal algebra operations. For the general storage version,
 * the routine Compute_NewState() should be used instead.
 * 
 * \pre
 * - \c IGain should be a symmetrical matrix in packed storage format with the upper triangular part
 *   referenced (lower part in FORTRAN). It should also be properly initialised through the
 *   MatrixVector_Create_PS() routine.
 * - The rest of the elements of type \c MatrixVector_t must be properly initialised through the
 *   MatrixVector_Create() routine.
 * - The size of the vectors and the matrix must be coherent since it will not be checked in the routine.
 * 
 * \param[in]     IGain      The inverted gain matrix \f$\mathcal{G}^{1}\f$.
 * \param[in]     Eff_ForceT The effective force vector \f$\vec f_{eff}^t\f$.
 * \param[in]     In_LoadT   The input load vector \f$\vec l_i^t\f$.
 * \param[in]     Err_ForceT the error force vector \f$\vec f_{err}^t\f$.
 * \param[in,out] Tempvec    Temporal vector that helps in the calculations. This is included in order to
 *                           avoid allocating and deallocating memory each step. On input only the number of
 *                           rows is used.
 * \param[out]    VecTdT_0   is new explicit state (displacement, velocity or acceleration depending on the
 *                           selected formulation) \f$\vec n_0^{t + \Delta t}\f$.
 *
 * \post
 * - \c VecTdT_0 is the result of the operation:
 *
 * \f[\vec n_0^{t + \Delta t} = \mathcal{G}^{-1}(\vec f_{eff}^t + \vec l_i^t + \vec f_{err}^t)\f]
 *
 * \sa MatrixVector_t and EffK_EffectiveForce_PS().
 */
void Compute_NewState_PS( const MatrixVector_t *const IGain, const MatrixVector_t *const Eff_ForceT,
			  const MatrixVector_t *const In_LoadT, const MatrixVector_t *const Err_ForceT,
			  MatrixVector_t *const Tempvec, MatrixVector_t *const VecTdT_0 );

void Compute_NewState_HHT( const MatrixVector_t *const IGain, const MatrixVector_t *const Eff_ForceT,
			   const MatrixVector_t *const In_LoadT, const MatrixVector_t *const In_LoadTdT,
			   const MatrixVector_t *const Err_ForceT, MatrixVector_t *const Tempvec,
			   const HYSL_FLOAT Alpha_H, MatrixVector_t *const VecTdT_0 );

void Compute_NewState_Zienkiewicz( const MatrixVector_t *const Meff, const MatrixVector_t *const MatA,
				   const MatrixVector_t *const MatB, const MatrixVector_t *const DispT,
				   const MatrixVector_t *const DispT0, const MatrixVector_t *const In_LoadT,
				   const MatrixVector_t *const Err_ForceT, const MatrixVector_t *const ForceT,
				   const MatrixVector_t *const ForceT0, const HYSL_FLOAT a8, const HYSL_FLOAT a17,
				   const HYSL_FLOAT a18, MatrixVector_t *const Tempvec, MatrixVector_t *const VecTdT_0 );
void Compute_NewState_Zienkiewicz_PS( const MatrixVector_t *const Meff, const MatrixVector_t *const MatA,
				      const MatrixVector_t *const MatB, const MatrixVector_t *const DispT,S
				      const MatrixVector_t *const DispT0, const MatrixVector_t *const In_LoadT,
				      const MatrixVector_t *const Err_ForceT, const MatrixVector_t *const ForceT,
				      const MatrixVector_t *const ForceT0, const HYSL_FLOAT a8, const HYSL_FLOAT a17,
				      const HYSL_FLOAT a18, MatrixVector_t *const Tempvec, MatrixVector_t *const VecTdT_0 );

void Compute_NewState_Zienkiewicz_Sp( const MatrixVector_t *const Meff, const MatrixVector_Sp_t *const MatA,
				      const MatrixVector_Sp_t *const MatB, const MatrixVector_t *const DispT,
				      const MatrixVector_t *const DispT0, const MatrixVector_t *const In_LoadT,
				      const MatrixVector_t *const Err_ForceT, const MatrixVector_t *const ForceT,
				      const MatrixVector_t *const ForceT0, const HYSL_FLOAT a8, const HYSL_FLOAT a17,
				      const HYSL_FLOAT a18, MatrixVector_t *const Tempvec, MatrixVector_t *const VecTdT_0 );

/**
 * \brief Calculates the new state. MPI version.
 *
 * The explicit part of the new state is computed. This routine does not depend on a specific formulation and
 * therefore it can be used in displacement, velocity and acceleration control formulations.
 * 
 * \f[\vec n_0^{t + \Delta t} = \mathcal{G}^{-1}(\vec f_{eff}^t + \vec l_i^t - \vec f_{err}^t)\f]
 *
 * where:
 * - \f$\vec n_0^{t + \Delta t}\f$ is the explicit part of the new state vector at time
 *   \f$t + \Delta t\f$,
 * - \f$\mathcal{G}\f$ is the gain matrix,
 * - \f$\vec f_{eff}^t\f$ is the effective force vector at time \f$t\f$,
 * - \f$\vec l_i^t\f$ is the input load vector at time \f$t\f$,
 * - and \f$\vec f_{err}^t\f$ is the error compensation force at time \f$t\f$.
 *
 * It makes use of PBLAS routines to perform the lineal algebra operations. For the general and packed storage
 * version, the routines Compute_NewState() or Compute_NewState_PS() should be used instead.
 * 
 * \pre
 * - All elements of type \c PMatrixVector_t must be properly initialised through the PMatrixVector_Create()
 *   routine.
 * - \c IGain must be a symmetrical matrix in general storage. Only the upper part will be referenced (lower
 *   part in FORTRAN).
 * - The size of the vectors and matrices must be coherent since it will not be checked in the routine.
 * 
 * \param[in]     IGain      The inverted gain matrix \f$\mathcal{G}^{1}\f$.
 * \param[in]     Eff_ForceT The effective force vector \f$\vec f_{eff}^t\f$.
 * \param[in]     In_LoadT   The input load vector \f$\vec l_i^t\f$.
 * \param[in]     Err_ForceT The error force vector \f$\vec f_{err}^t\f$.
 * \param[in,out] Tempvec    Temporal vector that helps in the calculations. This is included in order to
 *                           avoid allocating and deallocating memory each step. On input only the number of
 *                           rows is used.
 * \param[out]    VecTdT_0   New explicit state (displacement, velocity or acceleration depending on the
 *                           selected formulation) \f$\vec n_0^{t + \Delta t}\f$.
 *
 * \post
 * - \c VecTdT_0 is the result of the operation:
 *
 * \f[\vec n_0^{t + \Delta t} = \mathcal{G}^{-1}(\vec f_{eff}^t + \vec l_i^t + \vec f_{err}^t)\f]
 *
 * \sa PMatrixVector_t and EffK_EffectiveForce_MPI().
 */
void Compute_NewState_MPI( PMatrixVector_t *const IGain, PMatrixVector_t *const Eff_ForceT,
			   PMatrixVector_t *const In_LoadT, PMatrixVector_t *const Err_ForceT,
			   PMatrixVector_t *const Tempvec, PMatrixVector_t *const VecTdT_0 );

#endif /* NEW_STATE_H_ */
