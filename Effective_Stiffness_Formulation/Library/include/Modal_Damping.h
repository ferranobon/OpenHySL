/**
 * \file Modal_Damping.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 31st of May 2014
 *
 * \brief Modal damping routines.
 *
 *
 * Routines for calculating the modal damping. Sparse BLAS operations are supported through the Intel MKL
 * library.
 */

#ifndef MODAL_DAMPING_H_
#define MODAL_DAMPING_H_

#include "MatrixVector.h"
#include "MatrixVector_Sp.h"
#include "MatrixVector_MPI.h"

#include "Definitions.h"

/**
 * \brief Computes the proportional viscous damping matrix using Rayleigh damping. General storage version.
 *
 * This routine calculates the proportional viscous damping matrix using Rayleigh Damping. The operation
 * performed is (see \cite Clough_1975 p. 234):
 *
 * \f[\mathcal{C} = \alpha\mathcal{M} \cdot \beta \mathcal{K}\f]
 *
 * where:
 * - \f$\mathcal{C}\f$ is the proportional viscous damping matrix,
 * - \f$\mathcal{M}\f$ is the mass matrix,
 * - \f$\mathcal{K}\f$ is the stiffness matrix,
 * - and \f$\alpha\f$ and \f$beta\f$ are the parameters that multiply the mass and stiffness matrices
 *   respectively.
 * 
 * It makes use of BLAS and LAPACK routines to perform the linear algebra operations.
 *
 * \pre 
 * - All elements of type \c MatrixVector_t must be properly intialised through the MatrixVector_Create()
 *   routine.
 * - The matrices must be symmetrical and only the upper part will be referenced (lower part in FORTRAN
 *   routines).
 * - The dimensions of the matrices must be the identical.
 * - The Rayleigh constants must be properly initialised.
 *
 * \param[in]     Mass     The mass matrix \f$\mathcal{M}\f$.
 * \param[in]     Stiff    The stiffness matrix \f$\mathcal{K}\f$.
 * \param[in,out] Damp     The proportional viscous damping matrix \f$\mathcal{C}\f$. As an input, only the
 *                         size of the matrix is referenced.
 * \param[in]     Rayleigh Contains the values of the proportional coefficients \f$\alpha\f$ (\c
 *                         Rayleigh.Alpha) and \f$\beta\f$ (\c Rayleigh.Beta).
 *
 * \post \c Damp is a symmetric matrix in general storage with only the upper part referenced (Lower part in
 *       FORTRAN routines). It contains the result of: \f[\mathcal{C} = \alpha \mathcal{M} \cdot \beta
 *       \mathcal{K}\f]
 *
 * \sa MatrixVector_t and Rayleigh_t.
 */
void Modal_Damping( const MatrixVector_t *const Mass, const MatrixVector_t *const Stiff, MatrixVector_t *const Damp,
		    double DampFactor );

#endif /* MODAL_DAMPING_H_ */
