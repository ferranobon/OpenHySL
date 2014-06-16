/**
 * \file EffK_Formulation.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 9th of February 2013
 *
 * \brief Routines specific to the effective stiffness matrix formulation of the algorithm.
 *
 * Routines specific to the effective stiffness matrix formulation of Dorka's substructure algorithm. The
 * routines are available for general, packed and sparse storage and they make use of the BLAS library to
 * perform the linear algebra operations and they support both single and double precision. Sparse BLAS
 * operations are supported through the Intel MKL library.
 */

#ifndef EFFK_FORMULATION_H_
#define EFFK_FORMULATION_H_

#include "MatrixVector.h"
#include "MatrixVector_Sp.h"
#include "MatrixVector_MPI.h"

/**
 * \brief Calculates the effective force vector according to the formulation using the effective
 * stiffness matrix. General storage version.
 * 
 * The effective force vector \f$\vec f_{eff}^t\f$ is calculated according to the formulation using the
 * effective stiffness matrix described in page 53 (\cite Dorka_1998):
 *
 * \f[\vec f_{eff}^t = \mathcal{M}(a_0\vec u^t + a_2\dot{\vec u}^t + a_3\ddot{\vec u}^t) + \mathcal{C}(a_1\vec
 * u^t + a_4\dot{\vec u}^t + a_5\ddot{\vec u}^t)\f]
 *
 * where:
 * - \f$\vec f_{eff}^t\f$ is the effective force vector at time \f$t\f$,
 * - \f$\mathcal{M}\f$ is the mass matrix,
 * - \f$\mathcal{C}\f$ is the proportional viscous damping matrix,
 * - \f$\vec u^t\f$ is the displacement vector at time \f$t\f$,
 * - \f$\dot{\vec u}^t\f$ is the velocity vector at time \f$t\f$,
 * - \f$\ddot{\vec u}^t\f$ is the acceleration vector at time \f$t\f$,
 * - and \f$a_0\f$, \f$a_1\f$, \f$a_2\f$, \f$a_3\f$, \f$a_4\f$ and \f$a_5\f$ are integration constants (see
 *   \cite Dorka_1998).
 *
 * It makes use of BLAS routines to perform the linear algebra operations. For sparse matrices or matrices in
 * packed storage the routines EffK_EffectiveForce_Sp() or EffK_EffectiveForce_PS() should be used instead.
 *
 * \pre
 * - All elements of type \c MatrixVector_t must be properly initialised through the MatrixVector_Create()
 *   routine.
 * - The dimensions of the vectors and matrices must be coherent since it will not be checked in the routine.
 * - The mass and viscous damping matrices must be symmetrical and only the upper part will be referenced
 *   (lower part in FORTRAN routines).
 * - The integration constants must be properly initialised.
 *
 * \param[in]     Mass       The Mass matrix \f$\mathcal{M}\f$.
 * \param[in]     Damp       The Viscous Damping matrix \f$\mathcal{C}\f$.
 * \param[in]     DispT      The displacement vector \f$\vec u^t\f$.
 * \param[in]     VelT       The velocity vector \f$\dot{\vec u}^t\f$.
 * \param[in]     AccT       The acceleration vector \f$\ddot{\vec u}^t\f$.
 * \param[in,out] Tempvec    Temporal vector that helps in the calculations. This is included in order to
 *                           avoid allocating and deallocating memory each step. On input only the number of
 *                           rows is used.
 * \param[in]     a0         The integration constant \f$a_0\f$ (\cite Dorka_1998).
 * \param[in]     a1         The integration constant \f$a_1\f$ (\cite Dorka_1998).
 * \param[in]     a2         The integration constant \f$a_2\f$ (\cite Dorka_1998).
 * \param[in]     a3         The integration constant \f$a_3\f$ (\cite Dorka_1998).
 * \param[in]     a4         The integration constant \f$a_4\f$ (\cite Dorka_1998).
 * \param[in]     a5         The integration constant \f$a_5\f$ (\cite Dorka_1998).
 * \param[out]    Eff_ForceT The effective force vector \f$\vec f_{eff}^t\f$.
 *
 * \post
 * - \c Eff_Force is the result of the operation:
 *
 * \f[\vec f_{eff}^t = \mathcal{M}(a_0\vec u^t + a_2\dot{\vec u}^t + a_3\ddot{\vec u}^t) + \mathcal{C}(a_1\vec
 * u^t + a_4\dot{\vec u}^t + a_5\ddot{\vec u}^t)\f]
 *
 * \sa MatrixVector_t.
 */
void EffK_EffectiveForce( const MatrixVector_t *const Mass, const MatrixVector_t *const Damp,
			  const MatrixVector_t *const DispT, const MatrixVector_t *const VelT,
			  const MatrixVector_t *const AccT, MatrixVector_t *const Tempvec, const HYSL_FLOAT a0,
			  const HYSL_FLOAT a1, const HYSL_FLOAT a2, const HYSL_FLOAT a3, const HYSL_FLOAT a4,
			  const HYSL_FLOAT a5, MatrixVector_t *const Eff_ForceT );

/**
 * \brief Calculates the effective force vector according to the formulation using the effective stiffness
 * matrix. Packed storage version.
 * 
 * The effective force vector \f$\vec f_{eff}^t\f$ is calculated according to the formulation using the
 * effective stiffness matrix described in page 53 (\cite Dorka_1998):
 *
 * \f[\vec f_{eff}^t = \mathcal{M}(a_0\vec u^t + a_2\dot{\vec u}^t + a_3\ddot{\vec u}^t) + \mathcal{C}(a_1\vec
 * u^t + a_4\dot{\vec u}^t + a_5\ddot{\vec u}^t)\f]
 *
 * where:
 * - \f$\vec f_{eff}^t\f$ is the effective force vector at time \f$t\f$,
 * - \f$\mathcal{M}\f$ is the mass matrix,
 * - \f$\mathcal{C}\f$ is the proportional viscous damping matrix,
 * - \f$\vec u^t\f$ is the displacement vector at time \f$t\f$,
 * - \f$\dot{\vec u}^t\f$ is the velocity vector at time \f$t\f$,
 * - \f$\ddot{\vec u}^t\f$ is the acceleration vector at time \f$t\f$,
 * - and \f$a_0\f$, \f$a_1\f$, \f$a_2\f$, \f$a_3\f$, \f$a_4\f$ and \f$a_5\f$ are integration constants (see
 *   \cite Dorka_1998).
 *
 * It makes use of BLAS routines to perform the linear algebra operations. For sparse matrices or matrices in
 * general storage the routines EffK_EffectiveForce_Sp() or EffK_EffectiveForce() should be used instead.
 *
 * \pre
 * - \c Mass and \c Damp should be symmetrical matrices in packed storage format with the upper triangular
 *   part referenced (lower part in FORTRAN). They should also be properly initialised through the
 *   MatrixVector_Create_PS() routine.
 * - The rest of the elements of type \c MatrixVector_t must be properly initialised through the
 *   MatrixVector_Create() routine.
 * - The dimensions of the vectors and matrices must be coherent since it will not be checked in the routine.
 * - The integration constants must be properly initialised.
 *
 * \param[in]     Mass       The Mass matrix \f$\mathcal{M}\f$.
 * \param[in]     Damp       The Viscous Damping matrix \f$\mathcal{C}\f$.
 * \param[in]     DispT      The displacement vector \f$\vec u^t\f$.
 * \param[in]     VelT       The velocity vector \f$\dot{\vec u}^t\f$.
 * \param[in]     AccT       The acceleration vector \f$\ddot{\vec u}^t\f$.
 * \param[in,out] Tempvec    Temporal vector that helps in the calculations. This is included in order to
 *                           avoid allocating and deallocating memory each step. On input only the number of
 *                           rows is used.
 * \param[in]     a0         The integration constant \f$a_0\f$ (\cite Dorka_1998).
 * \param[in]     a1         The integration constant \f$a_1\f$ (\cite Dorka_1998).
 * \param[in]     a2         The integration constant \f$a_2\f$ (\cite Dorka_1998).
 * \param[in]     a3         The integration constant \f$a_3\f$ (\cite Dorka_1998).
 * \param[in]     a4         The integration constant \f$a_4\f$ (\cite Dorka_1998).
 * \param[in]     a5         The integration constant \f$a_5\f$ (\cite Dorka_1998).
 * \param[out]    Eff_ForceT The effective force vector \f$\vec f_{eff}^t\f$.
 *
 * \post
 * - \c Eff_Force is the result of the operation:
 *
 * \f[\vec f_{eff}^t = \mathcal{M}(a_0\vec u^t + a_2\dot{\vec u}^t + a_3\ddot{\vec u}^t) + \mathcal{C}(a_1\vec
 * u^t + a_4\dot{\vec u}^t + a_5\ddot{\vec u}^t)\f]
 *
 * \sa MatrixVector_t.
 */
void EffK_EffectiveForce_PS( const MatrixVector_t *const Mass, const MatrixVector_t *const Damp,
			     const MatrixVector_t *const DispT, const MatrixVector_t *const VelT,
			     const MatrixVector_t *const AccT, MatrixVector_t *const Tempvec, const HYSL_FLOAT a0,
			     const HYSL_FLOAT a1, const HYSL_FLOAT a2, const HYSL_FLOAT a3, const HYSL_FLOAT a4, const HYSL_FLOAT a5,
			     MatrixVector_t *const Eff_ForceT );

void EffK_EffectiveForce_HHT( const MatrixVector_t *const Mass, const MatrixVector_t *const Damp,
			      const MatrixVector_t *const Stiff, const MatrixVector_t *const DispT,
			      const MatrixVector_t *const VelT, const MatrixVector_t *const AccT,
			      MatrixVector_t *const Tempvec, const HYSL_FLOAT a0, const HYSL_FLOAT a1,
			      const HYSL_FLOAT a2, const HYSL_FLOAT a3, const HYSL_FLOAT a4,
			      const HYSL_FLOAT a5, const HYSL_FLOAT Alpha_HHT, MatrixVector_t *const Eff_ForceT );
void EffK_EffectiveForce_HHT_PS( const MatrixVector_t *const Mass, const MatrixVector_t *const Damp,
				 const MatrixVector_t *const Stiff, const MatrixVector_t *const DispT,
				 const MatrixVector_t *const VelT, const MatrixVector_t *const AccT,
				 MatrixVector_t *const Tempvec, const HYSL_FLOAT a0, const HYSL_FLOAT a1,
				 const HYSL_FLOAT a2, const HYSL_FLOAT a3, const HYSL_FLOAT a4,
				 const HYSL_FLOAT a5, const HYSL_FLOAT Alpha_HHT, MatrixVector_t *const Eff_ForceT );

/**
 * \brief Calculates the effective force vector according to the formulation using the effective stiffness
 * matrix. Sparse version.
 * 
 * The effective force vector \f$\vec f_{eff}^t\f$ is calculated according to the formulation using the
 * effective stiffness matrix described in page 53 (\cite Dorka_1998):
 *
 * \f[\vec f_{eff}^t = \mathcal{M}(a_0\vec u^t + a_2\dot{\vec u}^t + a_3\ddot{\vec u}^t) + \mathcal{C}(a_1\vec
 * u^t + a_4\dot{\vec u}^t + a_5\ddot{\vec u}^t)\f]
 *
 * where:
 * - \f$\vec f_{eff}^t\f$ is the effective force vector at time \f$t\f$,
 * - \f$\mathcal{M}\f$ is the mass matrix,
 * - \f$\mathcal{C}\f$ is the proportional viscous damping matrix,
 * - \f$\vec u^t\f$ is the displacement vector at time \f$t\f$,
 * - \f$\dot{\vec u}^t\f$ is the velocity vector at time \f$t\f$,
 * - \f$\ddot{\vec u}^t\f$ is the acceleration vector at time \f$t\f$,
 * - and \f$a_0\f$, \f$a_1\f$, \f$a_2\f$, \f$a_3\f$, \f$a_4\f$ and \f$a_5\f$ are integration constants (see
 *   \cite Dorka_1998).
 *
 * It makes use of the BLAS and Sparse BLAS routines from the Intel Math Kernel Library (\cite MKL_2013) to
 * perform the linear algebra operations. For dense matrices, the routine EffK_EffectiveForce() should be used
 * instead.
 *
 * \pre
 * - All elements of type \c MatrixVector_t must be properly initialised through the MatrixVector_Create()
 *   routine.
 * - All elements of type \c MatrixVector_Sp_t must be properly intialised through the
 *   MatrixVector_Create_Sp() routine.
 * - The matrices must in Intel's MKL CSR-\em three \em array \em variation and in one based index.
 * - The matrices must be symmetrical with only the upper part referenced (lower part in FORTRAN routines).
 * - The dimensions of the vectors and matrices must be coherent since it will not be checked in the routine.
 * - The integration constants must be properly initialised.
 *
 * \param[in]     Mass       The Mass matrix \f$\mathcal{M}\f$.
 * \param[in]     Damp       The Viscous Damping matrix \f$\mathcal{C}\f$.
 * \param[in]     DispT      The displacement vector \f$\vec u^t\f$.
 * \param[in]     VelT       The velocity vector \f$\dot{\vec u}^t\f$.
 * \param[in]     AccT       The acceleration vector \f$\ddot{\vec u}^t\f$.
 * \param[in,out] Tempvec    Temporal vector that helps in the calculations. This is included in order to
 *                           avoid allocating and deallocating memory each step. On input only the number of
 *                           rows is used.
 * \param[in]     a0         The integration constant \f$a_0\f$ (\cite Dorka_1998).
 * \param[in]     a1         The integration constant \f$a_1\f$ (\cite Dorka_1998).
 * \param[in]     a2         The integration constant \f$a_2\f$ (\cite Dorka_1998).
 * \param[in]     a3         The integration constant \f$a_3\f$ (\cite Dorka_1998).
 * \param[in]     a4         The integration constant \f$a_4\f$ (\cite Dorka_1998).
 * \param[in]     a5         The integration constant \f$a_5\f$ (\cite Dorka_1998).
 * \param[out]    Eff_ForceT The effective force vector \f$\vec f_{eff}^t\f$.
 *
 * \post
 * - \c Eff_Force is the result of the operation:
 *
 * \f[\vec f_{eff}^t = \mathcal{M}(a_0\vec u^t + a_2\dot{\vec u}^t + a_3\ddot{\vec u}^t) + \mathcal{C}(a_1\vec
 * u^t + a_4\dot{\vec u}^t + a_5\ddot{\vec u}^t)\f]
 *
 * \sa MatrixVector_t and MatrixVector_Sp_t.
 */
void EffK_EffectiveForce_Sp( const MatrixVector_Sp_t *const Mass, const MatrixVector_Sp_t *const Damp,
			     const MatrixVector_t *const DispT, const MatrixVector_t *const VelT,
			     const MatrixVector_t *const AccT, MatrixVector_t *const Tempvec, const HYSL_FLOAT a0,
			     const HYSL_FLOAT a1, const HYSL_FLOAT a2, const HYSL_FLOAT a3, const HYSL_FLOAT a4,
			     const HYSL_FLOAT a5, MatrixVector_t *const Eff_ForceT );

void EffK_EffectiveForce_HHT_Sp( const MatrixVector_Sp_t *const Mass, const MatrixVector_Sp_t *const Damp,
				 const MatrixVector_Sp_t *const Stiff, const MatrixVector_t *const DispT,
				 const MatrixVector_t *const VelT, const MatrixVector_t *const AccT,
				 MatrixVector_t *const Tempvec, const HYSL_FLOAT a0, const HYSL_FLOAT a1,
				 const HYSL_FLOAT a2, const HYSL_FLOAT a3, const HYSL_FLOAT a4, const HYSL_FLOAT a5,
				 const HYSL_FLOAT Alpha_HHT, MatrixVector_t *const Eff_ForceT );

/**
 * \brief Calculates the effective force vector according to the formulation using the effective
 * stiffness matrix. MPI version.
 * 
 * The effective force vector \f$\vec f_{eff}^t\f$ is calculated according to the formulation using the
 * effective stiffness matrix described in page 53 (\cite Dorka_1998):
 *
 * \f[\vec f_{eff}^t = \mathcal{M}(a_0\vec u^t + a_2\dot{\vec u}^t + a_3\ddot{\vec u}^t) + \mathcal{C}(a_1\vec
 * u^t + a_4\dot{\vec u}^t + a_5\ddot{\vec u}^t)\f]
 *
 * where:
 * - \f$\vec f_{eff}^t\f$ is the effective force vector at time \f$t\f$,
 * - \f$\mathcal{M}\f$ is the mass matrix,
 * - \f$\mathcal{C}\f$ is the proportional viscous damping matrix,
 * - \f$\vec u^t\f$ is the displacement vector at time \f$t\f$,
 * - \f$\dot{\vec u}^t\f$ is the velocity vector at time \f$t\f$,
 * - \f$\ddot{\vec u}^t\f$ is the acceleration vector at time \f$t\f$,
 * - and \f$a_0\f$, \f$a_1\f$, \f$a_2\f$, \f$a_3\f$, \f$a_4\f$ and \f$a_5\f$ are integration constants (see
 *   \cite Dorka_1998).
 *
 * It makes use of PBLAS routines to perform the linear algebra operations. For general, packed or sparse
 * formats the routines EffK_EffectiveForce(), EffK_EffectiveForce_PS() or EffK_EffectiveForce_Sp() should be
 * used instead.
 *
 * \pre
 * - All elements of type \c PMatrixVector_t must be properly initialised through the PMatrixVector_Create()
 *   routine.
 * - The dimensions of the vectors and matrices must be coherent since it will not be checked in the routine.
 * - The mass and viscous damping matrices must be symmetrical and only the upper part will be referenced
 *   (lower part in FORTRAN routines).
 * - The integration constants must be properly initialised.
 *
 * \param[in]     Mass       The Mass matrix \f$\mathcal{M}\f$.
 * \param[in]     Damp       The Viscous Damping matrix \f$\mathcal{C}\f$.
 * \param[in]     DispT      The displacement vector \f$\vec u^t\f$.
 * \param[in]     VelT       The velocity vector \f$\dot{\vec u}^t\f$.
 * \param[in]     AccT       The acceleration vector \f$\ddot{\vec u}^t\f$.
 * \param[in,out] Tempvec    Temporal vector that helps in the calculations. This is included in order to
 *                           avoid allocating and deallocating memory each step. On input only the number of
 *                           rows is used.
 * \param[in]     a0         The integration constant \f$a_0\f$ (\cite Dorka_1998).
 * \param[in]     a1         The integration constant \f$a_1\f$ (\cite Dorka_1998).
 * \param[in]     a2         The integration constant \f$a_2\f$ (\cite Dorka_1998).
 * \param[in]     a3         The integration constant \f$a_3\f$ (\cite Dorka_1998).
 * \param[in]     a4         The integration constant \f$a_4\f$ (\cite Dorka_1998).
 * \param[in]     a5         The integration constant \f$a_5\f$ (\cite Dorka_1998).
 * \param[out]    Eff_ForceT The effective force vector \f$\vec f_{eff}^t\f$.
 *
 * \post
 * - \c Eff_Force is the result of the operation:
 *
 * \f[\vec f_{eff}^t = \mathcal{M}(a_0\vec u^t + a_2\dot{\vec u}^t + a_3\ddot{\vec u}^t) + \mathcal{C}(a_1\vec
 * u^t + a_4\dot{\vec u}^t + a_5\ddot{\vec u}^t)\f]
 *
 * \sa PMatrixVector_t.
 */
void EffK_EffectiveForce_MPI( PMatrixVector_t *const Mass, PMatrixVector_t *const Damp,
			      PMatrixVector_t *const DispT, PMatrixVector_t *const VelT,
			      PMatrixVector_t *const AccT, PMatrixVector_t *const Tempvec,
			      const HYSL_FLOAT a0, const HYSL_FLOAT a1, const HYSL_FLOAT a2, const HYSL_FLOAT a3,
			      const HYSL_FLOAT a4, const HYSL_FLOAT a5, PMatrixVector_t *const Eff_ForceT );

/**
 * \brief Computes the new acceleration according to the formulation using the effective stiffness matrix.
 *
 * The acceleration at \f$t+\Delta t\f$ is calculated according to the formulation using the effective
 * stiffness matrix described in page 53 (\cite Dorka_1998) \f$\ddot{\vec{u}}^t\f$:
 *
 * \f[\ddot{\vec u}^{t+\Delta t} = a_0 (\vec u^{t+\Delta t} - \vec u) - a_2\dot{\vec u}^t -a_3\ddot{\vec
 * u}^t\f]
 *
 * where:
 * - \f$\ddot{\vec u}^{t+\Delta t}\f$ is the acceleration vector at time \f$t + \Delta t\f$,
 * - \f$\vec u^t\f$ is the displacement vector at time \f$t\f$,
 * - \f$\vec u^{t+\Delta t}\f$ is the displacement vector at time \f$t + \Delta t\f$,
 * - \f$\dot{\vec u}^t\f$ is the velocity vector at time \f$t\f$,
 * - \f$\ddot{\vec u}^t\f$ is the acceleration vector at time \f$t\f$,
 * - and \f$a_0\f$, \f$a_2\f$ and \f$a_3\f$ are integration constants (see \cite Dorka_1998).
 *
 * It makes use of BLAS routines to perform the linear algebra operations.
 *
 * \pre
 * - All elements of type \c MatrixVector_t must be properly initialised through the MatrixVector_Create()
 *   routine.
 * - The dimensions of the vectors must be identical.
 * - The integration constants must be properly initialised.
 * 
 * \param[in]     DispTdT The displacement vector \f$\vec u^{t+\Delta t}\f$.
 * \param[in]     DispT   The displacement vector \f$\vec u^t\f$.
 * \param[in]     VelT    The velocity vector \f$\dot{\vec u}^t\f$.
 * \param[in]     AccT    The acceleration vector \f$\ddot{\vec u}^t\f$.
 * \param[in]     a0      The integration constant \f$a_0\f$ (\cite Dorka_1998).
 * \param[in]     a2      The integration constant \f$a_2\f$ (\cite Dorka_1998).
 * \param[in]     a3      The integration constant \f$a_3\f$ (\cite Dorka_1998).
 * \param[in,out] AccTdT  The acceleration vector \f$\ddot{\vec u}^{t+\Delta t}\f$. On input only the number
 *                        of rows is used.
 *
 * \post
 * - \c AccTdT is the result of the operation:
 *
 * \f[\ddot{\vec u}^{t+\Delta t} = a_0 (\vec u^{t+\Delta t} - \vec u) - a_2\dot{\vec u}^t -a_3\ddot{\vec
 * u}^t\f]
 *
 * \sa MatrixVector_t.
 */
void EffK_ComputeAcceleration( const MatrixVector_t *const DispTdT, const MatrixVector_t *const DispT,
				const MatrixVector_t *const VelT, const MatrixVector_t *const AccT,
				const HYSL_FLOAT a0, const HYSL_FLOAT a2, const HYSL_FLOAT a3, MatrixVector_t *const AccTdT );

/**
 * \brief Computes the new acceleration according to the formulation using the effective stiffness matrix. MPI
 * version.
 *
 * The acceleration at \f$t+\Delta t\f$ is calculated according to the formulation using the effective
 * stiffness matrix described in page 53 (\cite Dorka_1998) \f$\ddot{\vec{u}}^t\f$:
 *
 * \f[\ddot{\vec u}^{t+\Delta t} = a_0 (\vec u^{t+\Delta t} - \vec u) - a_2\dot{\vec u}^t -a_3\ddot{\vec
 * u}^t\f]
 *
 * where:
 * - \f$\ddot{\vec u}^{t+\Delta t}\f$ is the acceleration vector at time \f$t + \Delta t\f$,
 * - \f$\vec u^t\f$ is the displacement vector at time \f$t\f$,
 * - \f$\vec u^{t+\Delta t}\f$ is the displacement vector at time \f$t + \Delta t\f$,
 * - \f$\dot{\vec u}^t\f$ is the velocity vector at time \f$t\f$,
 * - \f$\ddot{\vec u}^t\f$ is the acceleration vector at time \f$t\f$,
 * - and \f$a_0\f$, \f$a_2\f$ and \f$a_3\f$ are integration constants (see \cite Dorka_1998).
 *
 * It makes use of PBLAS routines to perform the linear algebra operations.
 *
 * \pre
 * - All elements of type \c PMatrixVector_t must be properly initialised through the PMatrixVector_Create()
 *   routine.
 * - The dimensions of the vectors must be identical.
 * - The integration constants must be properly initialised.
 * 
 * \param[in]     DispTdT The displacement vector \f$\vec u^{t+\Delta t}\f$.
 * \param[in]     DispT   The displacement vector \f$\vec u^t\f$.
 * \param[in]     VelT    The velocity vector \f$\dot{\vec u}^t\f$.
 * \param[in]     AccT    The acceleration vector \f$\ddot{\vec u}^t\f$.
 * \param[in]     a0      The integration constant \f$a_0\f$ (\cite Dorka_1998).
 * \param[in]     a2      The integration constant \f$a_2\f$ (\cite Dorka_1998).
 * \param[in]     a3      The integration constant \f$a_3\f$ (\cite Dorka_1998).
 * \param[in,out] AccTdT  The acceleration vector \f$\ddot{\vec u}^{t+\Delta t}\f$. On input only the number
 *                        of rows is used.
 *
 * \post
 * - \c AccTdT is the result of the operation:
 *
 * \f[\ddot{\vec u}^{t+\Delta t} = a_0 (\vec u^{t+\Delta t} - \vec u) - a_2\dot{\vec u}^t -a_3\ddot{\vec
 * u}^t\f]
 *
 * \sa PMatrixVector_t.
 */
void EffK_ComputeAcceleration_MPI( PMatrixVector_t *const DispTdT, PMatrixVector_t *const DispT,
				   PMatrixVector_t *const VelT, PMatrixVector_t *const AccT,
				   const HYSL_FLOAT a0, const HYSL_FLOAT a2, const HYSL_FLOAT a3,
				   PMatrixVector_t *const AccTdT );

/**
 * \brief Computes the new velocity according to the formulation using the effective stiffness matrix.
 *
 * The velocity at \f$t+\Delta t\f$ is calculated according to the formulation using the effective stiffness
 * matrix described in page 53 (\cite Dorka_1998):
 *
 * \f[\dot{\vec u}^{t+\Delta t} = a_1(\vec u^{t+\Delta t} - vec u^t) -a_4\dot{\vec u}^t - a_5\ddot{\vec u}^t\f]
 *
 * where:
 * - \f$\dot{\vec u}^{t+\Delta t}\f$ is the velocity vector at time \f$t\f$,
 * - \f$\dot{\vec u}^t\f$ is the velocity vector at time \f$t\f$,
 * - \f$\ddot{\vec u}^t\f$ is the acceleration vector at time \f$t\f$,
 * - \f$\ddot{\vec u}^{t+\Delta t}\f$ is the acceleration vector at time \f$t + \Delta t\f$,
 * - and \f$a_1\f$, \f$a_4\f$ and \f$a_5\f$ are integration constants (see \cite Dorka_1998).
 *
 * It makes use of BLAS routines to perform the linear algebra operations.
 *
 * \pre
 * - All elements of type \c MatrixVector_t must be properly initialised through the MatrixVector_Create()
 *   routine.
 * - The dimensions of the vectors must be identical.
 * - The integration constants must be properly initialised.
 * 
 * \param[in]     DispTdT The displacement vector \f$\vec u^{t+\Delta t}\f$.
 * \param[in]     DispT   The displacement vector \f$\vec u^t\f$.
 * \param[in]     VelT    The velocity vector \f$\dot{\vec u}^t\f$.
 * \param[in]     AccT    The acceleration vector \f$\ddot{\vec u^t}\f$.
 * \param[in]     a1      The integration constant \f$a_1\f$ (\cite Dorka_1998).
 * \param[in]     a4      The integration constant \f$a_4\f$ (\cite Dorka_1998).
 * \param[in]     a5      The integration constant \f$a_5\f$ (\cite Dorka_1998).
 * \param[in,out] VelTdT  The velocity vector \f$\dot{\vec u}^{t+\Delta t}\f$. On input only the number of rows
 *                        is used.
 *
 * \post
 * - \c VelTdT is the result of the operation:
 *
 * \f[\dot{\vec u}^{t+\Delta t} = a_1(\vec u^{t+\Delta t} - vec u^t) -a_4\dot{\vec u}^t - a_5\ddot{\vec u}^t\f]
 *
 * \sa MatrixVector_t.
 */
void EffK_ComputeVelocity( const MatrixVector_t *const DispTdT, const MatrixVector_t *const DispT,
			   const MatrixVector_t *const VelT, const MatrixVector_t *const AccT,
			   const HYSL_FLOAT a1, const HYSL_FLOAT a4, const HYSL_FLOAT a5, MatrixVector_t *const VelTdT );

/**
 * \brief Computes the new velocity according to the formulation using the effective stiffness matrix. MPI
 * version.
 *
 * The velocity at \f$t+\Delta t\f$ is calculated according to the formulation using the effective stiffness
 * matrix described in page 53 (\cite Dorka_1998):
 *
 * \f[\dot{\vec u}^{t+\Delta t} = \dot{\vec u}^t + a_6\ddot{\vec u}^t + a_7\ddot{\vec u}^{t+\Delta t}\f]
 *
 * where:
 * - \f$\dot{\vec u}^{t+\Delta t}\f$ is the velocity vector at time \f$t\f$,
 * - \f$\dot{\vec u}^t\f$ is the velocity vector at time \f$t\f$,
 * - \f$\ddot{\vec u}^t\f$ is the acceleration vector at time \f$t\f$,
 * - \f$\ddot{\vec u}^{t+\Delta t}\f$ is the acceleration vector at time \f$t + \Delta t\f$,
 * - and \f$a_6\f$ and \f$a_7\f$ are integration constants (see \cite Dorka_1998).
 *
 * It makes use of BLAS routines to perform the linear algebra operations.
 *
 * \pre
 * - All elements of type \c PMatrixVector_t must be properly initialised through the PMatrixVector_Create()
 *   routine.
 * - The dimensions of the vectors must be identical.
 * - The integration constants must be properly initialised.
 * 
 * \param[in]     VelT   The velocity vector \f$\dot{\vec u}^t\f$.
 * \param[in]     AccT   The acceleration vector \f$\ddot{\vec u^t}\f$.
 * \param[in]     AccTdT The acceleration vector \f$\ddot{\vec u}^{t+\Delta t}\f$.
 * \param[in]     a6     The integration constant \f$a_0\f$ (\cite Dorka_1998).
 * \param[in]     a7     The integration constant \f$a_1\f$ (\cite Dorka_1998).
 * \param[in,out] VelTdT The velocity vector \f$\dot{\vec u}^{t+\Delta t}\f$. On input only the number of rows
 *                       is used.
 *
 * \post
 * - \c VelTdT is the result of the operation:
 *
 * \f[\dot{\vec u}^{t+\Delta t} = \dot{\vec u}^t + a_6\ddot{\vec u}^t + a_7\ddot{\vec u}^{t+\Delta t}\f]
 *
 * \sa PMatrixVector_t.
 */
void EffK_ComputeVelocity_MPI( PMatrixVector_t *const VelT, PMatrixVector_t *const AccT,
			       PMatrixVector_t *const AccTdT, const HYSL_FLOAT a6, const HYSL_FLOAT a7,
			       PMatrixVector_t *const VelTdT );

#endif /* EFFK_FORMULATION_H_ */
