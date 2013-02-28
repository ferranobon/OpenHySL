/**
 * \file Initiation_Sp.h
 * \author Ferran Obón Santacana
 * \version 1.0 PARDISO solver for matrix inversion is deprecated.
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

#ifndef INITIATION_SP_H_
#define INITIATION_SP_H_

#include "MatrixVector.h"
#include "MatrixVector_Sp.h"

/**
 * \brief Computes the proportional viscous damping matrix using Rayleigh damping. Sparse version.
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
 * It makes use of the BLAS, LAPACK and Sparse BLAS routines from the Intel Math Kernel Library (\cite MKL_2013) to
 * perform the linear algebra operations. This routine requires the MKL library.
 *
 * \pre 
 * - All elements of type \c MatrixVector_Sp_t must be properly intialised through the MatrixVector_Create_Sp() routine.
 * - The matrices must in Intel's MKL CSR-\em three \em array \em variation and in one based index.
 * - The dimensions of the matrices must be identical.
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
 * \sa MatrixVector_Sp_t and Rayleigh_t.
 */
void Rayleigh_Damping_Sp( const MatrixVector_Sp_t *const Mass, const MatrixVector_Sp_t *const Stiff, MatrixVector_Sp_t *const Damp, const Rayleigh_t *const Rayleigh );

/**
 * \brief Computes the gain matrix and its inverse. Sparse version.
 *
 * \warning This routine requires the Intel Math Kernel Library (\cite MKL_2013).
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
 * It makes use of the BLAS, LAPACK and Sparse BLAS routines from the Intel Math Kernel Library (\cite MKL_2013) to
 * perform the linear algebra operations.
 *
 * \pre 
 * - The gain matrix must be properly initialised through the MatrixVector_Create() routine.
 * - All elements of type \c MatrixVector_Sp_t must be properly intialised through the MatrixVector_Create_Sp() routine.
 * - The matrices must in Intel's MKL CSR-\em three \em array \em variation and in one based index.
 * - The matrices must be symmetrical with only the upper part referenced (lower part in FORTRAN routines).
 * - The dimensions of the matrices must be the identical.
 * - The maximum number or non zero elements of the proportional viscous damping matrix must be equal to the
 *   number of non zero elements in the stiffness matrix.
 * - The values in scalars must be properly initialised and they depend on the used formulation (see \cite Dorka_1998).
 *
 * \param[in,out] IGain The gain matrix. As an input, only the size of the matrix is referenced, not its elements.
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
 * \sa MatrixVector_t, MatrixVector_Sp_t and Scalars_t.
 */
void IGainMatrix_Sp( MatrixVector_t *const IGain, const MatrixVector_Sp_t *const Mass, const MatrixVector_Sp_t *const Damp,
		    const MatrixVector_Sp_t *const Stiff, const Scalars_t Const );

/**
 * \brief Computes the gain matrix and its inverse. PARDISO version.
 *
 * \warning This routine requires the Intel Math Kernel Library (\cite MKL_2013).
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
 * It makes use of BLAS libraries to perform the linear algebra operations and the PARDISO solver from the Intel Math Kernel
 * Library (\cite MKL_2013) to compute the matrix inversion.
 *
 * \pre 
 * - All elements of type \c MatrixVector_t must be properly intialised through the MatrixVector_Create() routine.
 * - The matrices must be symmetrical and only the upper part will be referenced (lower part in FORTRAN routines).
 * - The dimensions of the matrices must be the identical.
 * - The values in scalars must be properly initialised and they depend on the used formulation (see \cite Dorka_1998).
 *
 * \param[in,out] IGain The gain matrix. As an input, only the size of the matrix is referenced.
 * \param[in] Mass The mass matrix.
 * \param[in] Damp The proportional viscous damping matrix.
 * \param[in] Stiff The stiffness matrix.
 * \param[in] Const Scalars that multiply the mass (\f$\alpha\f$), \c Const.Aplha), damping (\f$\beta\f$), \c Const.Beta),
 * stiffness (\f$\gamma\f$), \c Const.Gamma) matrices and the inverted matrix (\f$\lambda\f$), \c Const.Lambda) respectively.
 *
 * \post \c IGain is a symmetric matrix in general storage with only the upper part referenced
 *       (Lower part in FORTRAN routines). It contains the result of:
 *       \f[\mathcal{G} = \lambda(\alpha\mathcal{M} + \beta\mathcal{C} + \gamma\mathcal{C})^{-1}\f]
 *
 * \sa MatrixVector_t and Scalars_t.
 */
void IGainMatrix_Pardiso( MatrixVector_t *const IGain, const MatrixVector_t *const Mass, const MatrixVector_t *const Damp,
			 const MatrixVector_t *const Stiff, const Scalars_t Const );

/**
 * \brief Computes the gain matrix and its inverse. Sparse PARDISO version.
 *
 * \warning This routine requires the Intel Math Kernel Library (\cite MKL_2013).
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
 * It makes use of the BLAS and Sparse BLAS routines from the Intel Math Kernel Library (\cite MKL_2013) to
 * perform the linear algebra operations and the PARDISO solver to compute the matrix inversion.
 *
 * \pre
 * - The gain matrix must be properly initialised through the MatrixVector_Create() routine.
 * - All elements of type \c MatrixVector_Sp_t must be properly intialised through the MatrixVector_Create_Sp() routine.
 * - The matrices must in Intel's MKL CSR-\em three \em array \em variation and in one based index.
 * - The matrices must be symmetrical with only the upper part referenced (lower part in FORTRAN routines).
 * - The dimensions of the matrices must be the identical.
 * - The maximum number or non zero elements of the gain matrix must be equal to the number of non zero elements
 *   in the stiffness matrix.
 * - The values in scalars must be properly initialised and they depend on the used formulation (see \cite Dorka_1998).
 *
 * \param[in,out] IGain The gain matrix. As an input, only the size of the matrix is referenced, not its elements.
 * \param[in] Mass The mass matrix.
 * \param[in] Damp The proportional viscous damping matrix.
 * \param[in] Stiff The stiffness matrix.
 * \param[in] Const Scalars that multiply the mass (\f$\alpha\f$), \c Const.Aplha), damping (\f$\beta\f$), \c Const.Beta)
 * and stiffness matrices (\f$\gamma\f$), \c Const.Gamma).
 *
 * \post \c IGain is a symmetric matrix in general storage with only the upper part referenced
 *       (Lower part in FORTRAN routines). It contains the result of:
 *       \f[\mathcal{G} = \lambda(\alpha\mathcal{M} + \beta\mathcal{C} + \gamma\mathcal{C})^{-1}\f]
 *
 * \sa MatrixVector_t, MatrixVector_Sp_t and Scalars_t.
 */
void IGainMatrix_Pardiso_Sp( MatrixVector_t *const IGain, const MatrixVector_Sp_t *const Mass, const MatrixVector_Sp_t *const Damp,
			     const MatrixVector_Sp_t *const Stiff, const Scalars_t Const );
#endif /* INITIATION_SP_H_ */
