/**
 * \file Common_Formulation_Sp.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 26th of February 2013
 * 
 * \todo Add support for packed storage to reduce memory use.
 *
 * \brief Rayleigh damping routines.
 *
 * Routines for calculating the proportional viscous damping matrix using Rayleigh Damping. The routines
 * make use of the BLAS library to perform the linear algebra operations and they support both single and
 * double precision. Sparse BLAS operations are supported through the Intel MKL library.
 */

#ifndef COMMON_FORMULATION_SP_H_
#define COMMON_FORMULATION_SP_H_

#include "MatrixVector.h"
#include "MatrixVector_Sp.h"

/**
 * \brief Calculates the input load as absolute values. Sparse version.
 *
 * This routine calculates the input load as absolute values that is required at the beginning of each step through:
 *
 * \f[\vec l_i^t = \mathcal{C} \dot{\vec u}_g + \mathcal{K} \vec u_g\f]
 *
 * where:
 * - \f$\vec l_i^t\f$ is the input load vector at time \f$t\f$
 * - \f$\mathcal{C}\f$ is the viscous damping matrix,
 * - \f$\dot{\vec u}_g\f$ is a vector with the ground velocity values,
 * - \f$\mathcal{K}\f$ is the stiffness matrix,
 * - \f$\vec u_g\f$ is a vector with the ground displacement values,
 *
 * It makes use of the BLAS and Sparse BLAS routines from the Intel Math Kernel Library (\cite MKL_2013) to
 * perform the linear algebra operations.
 *
 * \pre
 * - All elements of type \c MatrixVector_t must be properly initialised through the MatrixVector_Create() routine.
 * - All elements of type \c MatrixVector_Sp_t must be properly intialised through the MatrixVector_Create_Sp() routine.
 * - The matrices must in Intel's MKL CSR-\em three \em array \em variation and in one based index.
 * - The matrices must be symmetrical and only the upper part will be referenced (lower part in FORTRAN routines).
 * - The dimensions of the matrices must be the identical.
 * - The dimensions of the vectors and matrices must be coherent since it will not be checked in the routine.
 * - The ground displacement and velocity vectors should have already the right values on the right position
 *   since the routine does not check if they are applied to the right degree of freedom.
 *
 * \param[in] Stiff The stiffness matrix \f$\mathcal K\f$.
 * \param[in] Damp The viscous damping matrix \f$\mathcal C\f$.
 * \param[in] GDisp Vector containing the ground displacement of the earthquake at a certain step \f$\vec u_g\f$.
 * \param[in] GVel Vector containing the ground velocity of the earthquake at a certain step \f$\dot{\vec u}_g\f$.
 * \param[in,out] InLoad The input load vector \f$\vec l_i^t\f$ as an absolute value. As an input, only the size
 * of the vector is referenced, not its elements.
 *
 * \post \c InLoad is the input load vector considering absolute values.
 *
 * \f[\vec l_i^t = \mathcal{C} \dot{\vec u}_g + \mathcal{K} \vec u_g\f]
 *
 * \sa MatrixVector_t and MatrixVector_Sp_t.
 */
void InputLoad_AbsValues_Sp( const MatrixVector_Sp_t *const Stiff, const MatrixVector_Sp_t *const Damp,
			     const MatrixVector_t *const GDisp, const MatrixVector_t *const GVel,
			     MatrixVector_t *const InLoad );

/**
 * \brief Calculates the input load as relative value. Sparse version.
 *
 * \warning This routine requires the Intel Math Kernel Library (\cite MKL_2013).
 *
 * This routine calculates the input load as relative values that is required at the beginning of each step through:
 *
 * \f[\vec l_i^t = -\mathcal{M} \ddot{\vec u}_g\f]
 *
 * where:
 * - \f$\vec l_i^t\f$ is the input load vector at time \f$t\f$
 * - \f$\mathcal{M}\f$ is the mass matrix,
 * - \f$\ddot{\vec u}_g\f$ is a vector with the ground acceleration values,
 *
 * It makes use of the BLAS and Sparse BLAS routines from the Intel Math Kernel Library (\cite MKL_2013) to
 * perform the linear algebra operations.
 *
 * \pre
 * - All elements of type \c MatrixVector_t must be properly initialised through the MatrixVector_Create() routine.
 * - \c Mass must be properly intialised through the MatrixVector_Create_Sp() routine.
 * - \c Mass must in Intel's MKL CSR-\em three \em array \em variation and in one based index.
 * - \c Mass must be symmetrical and only the upper part will be referenced (lower part in FORTRAN routines).
 * - The dimensions of the vectors and matrices must be coherent since it will not be checked in the routine.
 * - The ground acceleration vector should have already the right values on the right position
 *   since the routine does not check if they are applied to the right degree of freedom.
 *
 * \param[in] Mass The mass matrix \f$\mathcal M\f$.
 * \param[in] GAcc Vector containing the ground acceleration of the earthquake at a certain step \f$\ddot{\vec u}_g\f$.
 * \param[in,out] InLoad The input load vector \f$\vec l_i^t\f$ as a relative value. As an input, only the size
 * of the vector is referenced, not its elements.
 *
 * \post \c InLoad is the input load vector considering relative values.
 *
 * \f[\vec l_i^t = -\mathcal{M} \ddot{\vec u}_g\f]
 *
 * \sa MatrixVector_t and MatrixVector_Sp_t.
 */
void InputLoad_RelValues_Sp( const MatrixVector_Sp_t *const Mass, const MatrixVector_t *const GAcc, 
			     MatrixVector_t *const InLoad );

#endif /* COMMON_FORMULATION_SP_H_ */
