/**
 * \file EffK_Formulation_Sp.h
 * \author Ferran Obón Santacana
 * \version 1.0 PARDISO solver for matrix inversion is deprecated.
 * \date 1st of March 2013
 * 
 * \todo Add support for packed storage to reduce memory use.
 *
 * \brief Rayleigh damping routines.
 *
 * Routines for calculating the proportional viscous damping matrix using Rayleigh Damping. The routines
 * make use of the BLAS library to perform the linear algebra operations and they support both single and
 * double precision. Sparse BLAS operations are supported through the Intel MKL library.
 */

#ifndef EFFK_FORMULATION_SP_H_
#define EFFK_FORMULATION_SP_H_

#include "MatrixVector_Sp.h"

/**
 * \brief Calculates the effective force vector according to the formulation using the effective
 * stiffness matrix. Sparse version.
 * 
 * The effective force vector \f$\vec f_{eff}^t\f$ is calculated according to the formulation using the effective stiffness matrix described
 * in page 53 (\cite Dorka_1998):
 *
 * \f[\vec f_{eff}^t = \mathcal{M}(a_0\vec u^t + a_2\dot{\vec u}^t + a_3\ddot{\vec u}^t) + \mathcal{C}(a_1\vec u^t + a_4\dot{\vec u}^t + a_5\ddot{\vec u}^t)\f]
 *
 * where:
 * - \f$\vec f_{eff}^t\f$ is the effective force vector at time \f$t\f$,
 * - \f$\mathcal{M}\f$ is the mass matrix,
 * - \f$\mathcal{C}\f$ is the proportional viscous damping matrix,
 * - \f$\vec u^t\f$ is the displacement vector at time \f$t\f$,
 * - \f$\dot{\vec u}^t\f$ is the velocity vector at time \f$t\f$,
 * - \f$\ddot{\vec u}^t\f$ is the acceleration vector at time \f$t\f$,
 * - and \f$a_0\f$, \f$a_1\f$, \f$a_2\f$, \f$a_3\f$, \f$a_4\f$ and \f$a_5\f$ are integration constants (see \cite Dorka_1998).
 *
 * It makes use of the BLAS and Sparse BLAS routines from the Intel Math Kernel Library (\cite MKL_2013) to
 * perform the linear algebra operations. For dense matrices, the routine EffK_EffectiveForce() should be used
 * instead.
 *
 * \pre
 * - All elements of type \c MatrixVector_t must be properly initialised through the MatrixVector_Create() routine.
 * - All elements of type \c MatrixVector_Sp_t must be properly intialised through the MatrixVector_Create_Sp() routine.
 * - The matrices must in Intel's MKL CSR-\em three \em array \em variation and in one based index.
 * - The matrices must be symmetrical with only the upper part referenced (lower part in FORTRAN routines).
 * - The dimensions of the vectors and matrices must be coherent since it will not be checked in the routine.
 * - The integration constants must be properly initialised.
 *
 * \param[in] Mass is the Mass matrix \f$\mathcal{M}\f$.
 * \param[in] Damp is the Viscous Damping matrix \f$\mathcal{C}\f$.
 * \param[in] DispT The displacement vector \f$\vec u^t\f$.
 * \param[in] VelT is the velocity vector \f$\dot{\vec u}^t\f$.
 * \param[in] AccT is the acceleration vector \f$\ddot{\vec u}^t\f$.
 * \param[in,out] Tempvec is a temporal vector that helps in the calculations. This is included in order to avoid
 *                allocating and deallocating memory each step. On input only the number of rows is used.
 * \param[in] a0 The integration constant \f$a_0\f$ (\cite Dorka_1998).
 * \param[in] a1 The integration constant \f$a_1\f$ (\cite Dorka_1998).
 * \param[in] a2 The integration constant \f$a_2\f$ (\cite Dorka_1998).
 * \param[in] a3 The integration constant \f$a_3\f$ (\cite Dorka_1998).
 * \param[in] a4 The integration constant \f$a_4\f$ (\cite Dorka_1998).
 * \param[in] a5 The integration constant \f$a_5\f$ (\cite Dorka_1998).
 * \param[out] Eff_ForceT is the effective force vector \f$\vec f_{eff}^t\f$.
 *
 * \post
 * - \c Eff_Force is the result of the operation:
 *
 * \f[\vec f_{eff}^t = \mathcal{M}(a_0\vec u^t + a_2\dot{\vec u}^t + a_3\ddot{\vec u}^t) + \mathcal{C}(a_1\vec u^t + a_4\dot{\vec u}^t + a_5\ddot{\vec u}^t)\f]
 *
 * \sa MatrixVector_t and MatrixVector_Sp_t.
 */
void EffK_EffectiveForce_Sp( const MatrixVector_Sp_t *const Mass, const MatrixVector_Sp_t *const Damp, const MatrixVector_t *const DispT,
			     const MatrixVector_t *const VelT, const MatrixVector_t *const AccT, MatrixVector_t *const Tempvec,
			     const double a0, const double a1, const double a2, const double a3, const double a4,
			     const double a5, MatrixVector_t *const Eff_ForceT );


#endif /* EFFK_FORMULATION_SP_H_ */
