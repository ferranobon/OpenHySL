/**
 * \file Initiation.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 9th of February 2013
 * 
 * \todo Add support for packed storage to reduce memory use.
 *
 * \brief Rayleigh damping routines.
 *
 * Routines for calculating the proportional viscous damping matrix using Rayleigh Damping. The routines
 * make use of the BLAS library to perform the linear algebra operations and they support both single and
 * double precision. Sparse BLAS operations are supported through the Intel MKL library.
 */

#ifndef INITIATION_H_
#define INITIATION_H_

#include "MatrixVector.h"

/**
 * \brief Stores the Rayleigh coefficients alpha and beta.
 */
typedef struct {
     double Alpha; /*!< \brief Coefficient that multiplies the mass matrix.*/
     double Beta;  /*!< \brief Coefficient that multiplies the stiffness matrix.*/
} Rayleigh_t;

/**
 * \brief Computes the proportional viscous damping matrix using Rayleigh damping.
 *
 * This routine calculates the proportional viscous damping matrix using Rayleigh Damping. The operation
 * performed is (see \cite Clough 1975 p. 234):
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
 * - All elements of type \c MatrixVector_t must be properly intialised through the MatrixVector_Create() routine.
 * - The matrices must be symmetrical and only the upper part will be referenced (lower part in FORTRAN routines).
 * - The dimensions of the matrices must be the identical.
 * - The Rayleigh constants must be properly initialised.
 *
 * \param[in] Mass The mass matrix \f$\mathcal{M}\f$.
 * \param[in] Stiff The stiffness matrix \f$\mathcal{K}\f$.
 * \param[in,out] Damp The proportional viscous damping matrix \f$\mathcal{C}\f$. As an input, only the
                       size of the matrix is referenced.
 * \param[in] Rayleigh Contains the values of the proportional coefficients \f$\alpha\f$ (\c Rayleigh.Alpha)
 *                     and \f$\beta\f$ (\c Rayleigh.Beta).
 *
 * \post \c Damp is a symmetric matrix in general storage with only the upper part referenced
 *       (Lower part in FORTRAN routines). It contains the result of:
 *       \f[\mathcal{C} = \alpha \mathcal{M} \cdot \beta \mathcal{K}\f]
 *
 * \sa MatrixVector_t and Rayleigh_t.
 */
void Rayleigh_Damping( const MatrixVector_t *const Mass, const MatrixVector_t *const Stiff, MatrixVector_t *const Damp,
		       const Rayleigh_t *const Rayleigh );

/**
 * \brief Computes the gain matrix and its inverse.
 *
 * This routine calculates the so called Gain Matrix and its inverse through:
 *
 * \f[\mathcal{G}^{-1} = \lambda(\alpha\mathcal{M} + \beta\mathcal{C} + \gamma\mathcal{C})^{-1}\f]
 *
 * where:
 * - \f$\mathcal{M}\f$ is the mass matrix,
 * - \f$\mathcal{C}\f$ is the proportional viscous damping matrix,
 * - \f$\mathcal{K}\f$ is the stiffness matrix,
 * - \f$\mathcal{G}\f$ is the gain matrix,
 * - and \f$\alpha\f$, \f$beta\f$, \f$\gamma\f$ and \f$\lambda\f$ are the parameters that multiply the mass, damping, stiffness and
 *   matrices and the result of the matrix inversion respectively.
 *
 * It makes use of BLAS and LAPACK routines to perform the linear algebra operations.
 *
 * \pre 
 * - All elements of type \c MatrixVector_t must be properly intialised through the MatrixVector_Create() routine.
 * - The matrices must be symmetrical and only the upper part will be referenced (lower part in FORTRAN routines).
 * - The dimensions of the matrices must be the identical.
 * - The values in scalars must be properly initialised and they depend on the used formulation (see \cite Dorka_1998).
 *
 * \param[in,out] IGain The inverse of the gain matrix. As an input, only the size of the matrix is referenced.
 * \param[in] Mass The mass matrix.
 * \param[in] Damp The proportional viscous damping matrix.
 * \param[in] Stiff The stiffness matrix.
 * \param[in] Const Scalars that multiply the mass (\f$\alpha\f$), \c Const.Aplha), damping (\f$\beta\f$), \c Const.Beta),
 * stiffness (\f$\gamma\f$), \c Const.Gamma) matrices and the inverted matrix (\f$\lambda\f$), \c Const.Lambda) respectively.
 *
 * \post \c IGain is a symmetric matrix in general storage with only the upper part referenced
 *       (Lower part in FORTRAN routines). It contains the result of:
 *       \f[\mathcal{G}^{-1} = \lambda(\alpha\mathcal{M} + \beta\mathcal{C} + \gamma\mathcal{C})^{-1}\f]
 *
 * \sa MatrixVector_t and Scalars_t.
 */
void IGainMatrix( MatrixVector_t *const IGain, const MatrixVector_t *const Mass, const MatrixVector_t *const Damp,
		  const MatrixVector_t *const Stiff, const Scalars_t Const );

#endif /* INITIATION_H_ */

