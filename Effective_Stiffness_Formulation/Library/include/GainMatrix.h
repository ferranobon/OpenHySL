/**
 * \file GainMatrix.h
 * \author Ferran Obón Santacana
 * \version 1.0 PARDISO solver for matrix inversion is deprecated.
 * \date 15th of March 2013
 *
 * \brief Routines to compute the gain matrix.
 *
 * \todo Add the single precision routines.
 *
 * Routines for calculating the gain matrix. The routines make use of the BLAS and LAPACK libraries to perform
 * the linear algebra operations, including the matrix inversion in single and double precision. Sparse BLAS
 * operations are supported through the Intel MKL library, but since matrix inversion is a dense operations
 * they still rely on LAPACK for this operation. The PARDISO solver is no longer supported since it requires
 * more memory and time than the LAPACK equivalent routines.
 */

#ifndef GAINMATRIX_H_
#define GAINMATRIX_H_

#include "MatrixVector.h"
#include "MatrixVector_PS.h"
#include "MatrixVector_Sp.h"
#include "MatrixVector_MPI.h"

/**
 * \brief Computes the gain matrix and its inverse. General storage version.
 *
 * This routine calculates the so called Gain Matrix and its inverse through:
 *
 * \f[\mathcal{G}^{-1} = \lambda(\alpha\mathcal{M} + \beta\mathcal{C} +
 * \gamma\mathcal{C})^{-1}\f]
 *
 * where:
 * - \f$\mathcal{M}\f$ is the mass matrix,
 * - \f$\mathcal{C}\f$ is the proportional viscous damping matrix,
 * - \f$\mathcal{K}\f$ is the stiffness matrix,
 * - \f$\mathcal{G}\f$ is the gain matrix,
 * - and \f$\alpha\f$, \f$beta\f$, \f$\gamma\f$ and \f$\lambda\f$ are the parameters that multiply the mass,
 *   damping, stiffness and matrices and the result of the matrix inversion respectively.
 *
 * It makes use of BLAS and LAPACK routines to perform the linear algebra operations.
 *
 * \pre 
 * - All elements of type \c MatrixVector_t must be properly intialised through the  MatrixVector_Create()
 *   routine.
 * - The matrices must be symmetrical and only the upper part will be referenced (lower part in FORTRAN
 *   routines).
 * - The dimensions of the matrices must be the identical.
 * - The values in scalars must be properly initialised and they depend on the used formulation
 *   \cite Dorka_1998.
 *
 * \param[in,out] IGain The inverse of the gain matrix. As an input, only the size of the matrix is
 *                      referenced.
 * \param[in]     Mass  The mass matrix.
 * \param[in]     Damp  The proportional viscous damping matrix.
 * \param[in]     Stiff The stiffness matrix.
 * \param[in]     Const Scalars that multiply the mass (\f$\alpha\f$), \c Const.Aplha), damping (\f$\beta\f$),
 *                      \c Const.Beta), stiffness (\f$\gamma\f$), \c Const.Gamma) matrices and the inverted
 *                      matrix (\f$\lambda\f$), \c Const.Lambda) respectively.
 *
 * \post \c IGain is a symmetric matrix in general storage with only the upper part referenced (Lower part in
 *       FORTRAN routines). It contains the result of: \f[\mathcal{G}^{-1} = \lambda(\alpha\mathcal{M} +
 *       \beta\mathcal{C} + \gamma\mathcal{C})^{-1}\f]
 *
 * \sa MatrixVector_t and Scalars_t.
 */
void IGainMatrix( MatrixVector_t *const IGain, const MatrixVector_t *const Mass,
		  const MatrixVector_t *const Damp, const MatrixVector_t *const Stiff, const Scalars_t Const );

/**
 * \brief Computes the gain matrix and its inverse. Packed storage version.
 *
 * This routine calculates the so called Gain Matrix and its inverse using packed storage through:
 *
 * \f[\mathcal{G}^{-1} = \lambda(\alpha\mathcal{M} + \beta\mathcal{C} + \gamma\mathcal{C})^{-1}\f]
 *
 * where:
 * - \f$\mathcal{M}\f$ is the mass matrix,
 * - \f$\mathcal{C}\f$ is the proportional viscous damping matrix,
 * - \f$\mathcal{K}\f$ is the stiffness matrix,
 * - \f$\mathcal{G}\f$ is the gain matrix,
 * - and \f$\alpha\f$, \f$beta\f$, \f$\gamma\f$ and \f$\lambda\f$ are the parameters that multiply the mass,
 *   damping, stiffness and matrices and the result of the matrix inversion respectively.
 *
 *
 * It makes use of BLAS and LAPACK routines to perform the linear algebra operations.
 *
 * \pre 
 * - All elements of type \c MatrixVector_t must be properly intialised through the  MatrixVector_Create_PS()
 *   routine.
 * - The matrices must be symmetrical and in packed storage. The upper triangular part (lower triangular part
 *   in FORTRAN) must be present.
 * - The dimensions of the matrices must be the identical.
 * - The values in scalars must be properly initialised and they depend on the used formulation
 *   \cite Dorka_1998
 *
 * \param[in,out] IGain The inverse of the gain matrix in packed storage. As an input, only the size of the
 *                      matrix is referenced.
 * \param[in]     Mass  The mass matrix.
 * \param[in]     Damp  The proportional viscous damping matrix.
 * \param[in]     Stiff The stiffness matrix.
 * \param[in]     Const Scalars that multiply the mass (\f$\alpha\f$), \c Const.Aplha), damping (\f$\beta\f$),
 *                      \c Const.Beta), stiffness (\f$\gamma\f$), \c Const.Gamma) matrices and the inverted
 *                      matrix (\f$\lambda\f$), \c Const.Lambda) respectively.
 *
 * \post \c IGain is a symmetric matrix in packed storage with only the upper part (lower part in FORTRAN). It
 *       contains the result of: \f[\mathcal{G}^{-1} = \lambda(\alpha\mathcal{M} + \beta\mathcal{C} +
 *       \gamma\mathcal{C})^{-1}\f]
 *
 * \sa MatrixVector_t and Scalars_t.
 */
void IGainMatrix_PS( MatrixVector_t *const IGain, const MatrixVector_t *const Mass,
		     const MatrixVector_t *const Damp, const MatrixVector_t *const Stiff,
		     const Scalars_t Const );

/**
 * \brief Computes the gain matrix and its inverse. Sparse version and general storage.
 *
 * \warning This routine requires the Intel Math Kernel Library \cite MKL_2013.
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
 * - and \f$\alpha\f$, \f$beta\f$, \f$\gamma\f$ and \f$\lambda\f$ are the parameters that multiply the mass,
 *   damping, stiffness and matrices and the result of the matrix inversion respectively.
 *
 * It makes use of the BLAS, LAPACK and Sparse BLAS routines from the Intel Math Kernel Library \cite MKL_2013
 * to perform the linear algebra operations.
 *
 * \pre 
 * - The gain matrix must be properly initialised through the MatrixVector_Create() routine.
 * - All elements of type \c MatrixVector_Sp_t must be properly intialised through the
 *   MatrixVector_Create_Sp() routine.
 * - The sparse matrices must in Intel's MKL CSR-\em three \em array \em variation and in one based index.
 * - The matrices must be symmetrical with only the upper part referenced (lower part in FORTRAN routines).
 * - The dimensions of the matrices must be the identical.
 * - The maximum number or non zero elements of the proportional viscous damping matrix must be equal to the
 *   number of non zero elements in the stiffness matrix.
 * - The values in scalars must be properly initialised and they depend on the used formulation
 *   \cite Dorka_1998.
 *
 * \param[in,out] IGain The gain matrix. As an input, only the size of the matrix is referenced, not its
 *                      elements.
 * \param[in]     Mass  The mass matrix.
 * \param[in]     Damp  The proportional viscous damping matrix.
 * \param[in]     Stiff The stiffness matrix.
 * \param[in]     Const Scalars that multiply the mass (\f$\alpha\f$), \c Const.Aplha), damping (\f$\beta\f$),
 *                      \c Const.Beta), stiffness (\f$\gamma\f$), \c Const.Gamma) matrices and the inverted
 *                      matrix (\f$\lambda\f$), \c Const.Lambda) respectively.

 *
 * \post \c IGain is a symmetric matrix in general storage with only the upper part referenced (Lower part in
 *       FORTRAN routines). It contains the result of: \f[\mathcal{G}^{-1} = \lambda(\alpha\mathcal{M} +
 *       \beta\mathcal{C} + \gamma\mathcal{C})^{-1}\f]
 *
 * \sa MatrixVector_t, MatrixVector_Sp_t and Scalars_t.
 */
void IGainMatrix_Sp( MatrixVector_t *const IGain, const MatrixVector_Sp_t *const Mass,
		     const MatrixVector_Sp_t *const Damp, const MatrixVector_Sp_t *const Stiff,
		     const Scalars_t Const );

/**
 * \brief Computes the gain matrix and its inverse. Sparse version and packed storage.
 *
 * \warning This routine requires the Intel Math Kernel Library \cite MKL_2013.
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
 * - and \f$\alpha\f$, \f$beta\f$, \f$\gamma\f$ and \f$\lambda\f$ are the parameters that multiply the mass,
 *   damping, stiffness and matrices and the result of the matrix inversion respectively.
 *
 * It makes use of the BLAS, LAPACK and Sparse BLAS routines from the Intel Math Kernel Library \cite MKL_2013
 * to perform the linear algebra operations.
 *
 * \pre 
 * - The gain matrix must be in packed storage and properly initialised through the  MatrixVector_Create_PS()
 *   routine.
 * - All elements of type \c MatrixVector_Sp_t must be properly intialised through the
 *   MatrixVector_Create_Sp() routine The sparse matrices must in Intel's MKL CSR-\em three \em array \em
 *   variation and in one based index.
 * - The sparse matrices must be symmetrical with only the upper part referenced (lower part in FORTRAN
 *   routines).
 * - The dimensions of the matrices must be the identical.
 * - The maximum number or non zero elements of the proportional viscous damping matrix must be equal to the
 *   number of non zero elements in the stiffness matrix.
 * - The values in scalars must be properly initialised and they depend on the used formulation 
 *   \cite Dorka_1998.
 *
 * \param[in,out] IGain The gain matrix. As an input, only the size of the matrix is referenced, not its elements.
 * \param[in]     Mass  The mass matrix.
 * \param[in]     Damp  The proportional viscous damping matrix.
 * \param[in]     Stiff The stiffness matrix.
 * \param[in]     Const Scalars that multiply the mass (\f$\alpha\f$), \c Const.Aplha), damping (\f$\beta\f$),
 *                      \c Const.Beta), stiffness (\f$\gamma\f$), \c Const.Gamma) matrices and the inverted
 *                      matrix (\f$\lambda\f$),\c Const.Lambda) respectively.
 *
 * \post \c IGain is a symmetric matrix in packed storage with only the upper part referenced (Lower part in
 *       FORTRAN routines). It contains the result of: \f[\mathcal{G}^{-1} = \lambda(\alpha\mathcal{M} +
 *       \beta\mathcal{C} + \gamma\mathcal{C})^{-1}\f]
 *
 * \sa MatrixVector_t, MatrixVector_Sp_t and Scalars_t.
 */
void IGainMatrix_Sp_PS( MatrixVector_t *const IGain, const MatrixVector_Sp_t *const Mass,
			const MatrixVector_Sp_t *const Damp, const MatrixVector_Sp_t *const Stiff,
			const Scalars_t Const );

/**
 * \brief Computes the gain matrix and its inverse. PARDISO version.
 *
 * \warning This routine requires the Intel Math Kernel Library \cite MKL_2013.
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
 * - and \f$\alpha\f$, \f$beta\f$, \f$\gamma\f$ and \f$\lambda\f$ are the parameters that multiply the mass,
 *   damping, stiffness and matrices and the result of the matrix inversion respectively.
 *
 * It makes use of BLAS libraries to perform the linear algebra operations and the PARDISO solver from the
 * Intel Math Kernel Library \cite MKL_2013 to compute the matrix inversion.
 *
 * \pre 
 * - All elements of type \c MatrixVector_t must be properly intialised through the MatrixVector_Create()
 *   routine.
 * - The matrices must be symmetrical and only the upper part will be referenced (lower part in FORTRAN
 *   routines).
 * - The dimensions of the matrices must be the identical.
 * - The values in scalars must be properly initialised and they depend on the used formulation 
 *   \cite Dorka_1998.
 *
 * \param[in,out] IGain The gain matrix. As an input, only the size of the matrix is referenced.
 * \param[in]     Mass  The mass matrix.
 * \param[in]     Damp  The proportional viscous damping matrix.
 * \param[in]     Stiff The stiffness matrix.
 * \param[in]     Const Scalars that multiply the mass (\f$\alpha\f$), \c Const.Aplha), damping (\f$\beta\f$),
 *                      \c Const.Beta), stiffness (\f$\gamma\f$), \c Const.Gamma) matrices and the inverted
 *                      matrix (\f$\lambda\f$), \c Const.Lambda) respectively.
 *
 * \post \c IGain is a symmetric matrix in general storage with only the upper part referenced (Lower part in
 *       FORTRAN routines). It contains the result of: \f[\mathcal{G} = \lambda(\alpha\mathcal{M} +
 *       \beta\mathcal{C} + \gamma\mathcal{C})^{-1}\f]
 *
 * \sa MatrixVector_t and Scalars_t.
 */
void IGainMatrix_Pardiso( MatrixVector_t *const IGain, const MatrixVector_t *const Mass,
			  const MatrixVector_t *const Damp, const MatrixVector_t *const Stiff,
			  const Scalars_t Const );

/**
 * \brief Computes the gain matrix and its inverse. Sparse PARDISO version.
 *
 * \warning This routine requires the Intel Math Kernel Library \cite MKL_2013.
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
 * - and \f$\alpha\f$, \f$beta\f$, \f$\gamma\f$ and \f$\lambda\f$ are the parameters that multiply the mass,
 *   damping, stiffness and matrices and the result of the matrix inversion respectively.
 *
 * It makes use of the BLAS and Sparse BLAS routines from the Intel Math Kernel Library \cite MKL_2013 to
 * perform the linear algebra operations and the PARDISO solver to compute the matrix inversion.
 *
 * \pre
 * - The gain matrix must be properly initialised through the MatrixVector_Create() routine.
 * - All elements of type \c MatrixVector_Sp_t must be properly intialised through the
 *   MatrixVector_Create_Sp() routine.
 * - The matrices must in Intel's MKL CSR-\em three \em array \em variation and in one based index.
 * - The matrices must be symmetrical with only the upper part referenced (lower part in FORTRAN routines).
 * - The dimensions of the matrices must be the identical.
 * - The maximum number or non zero elements of the gain matrix must be equal to the number of non zero
 *   elements in the stiffness matrix.
 * - The values in scalars must be properly initialised and they depend on the used formulation
 *   \cite Dorka_1998.
 *
 * \param[in,out] IGain The gain matrix. As an input, only the size of the matrix is referenced, not its
 *                      elements.
 * \param[in]     Mass  The mass matrix.
 * \param[in]     Damp  The proportional viscous damping matrix.
 * \param[in]     Stiff The stiffness matrix.
 * \param[in]     Const Scalars that multiply the mass (\f$\alpha\f$), \c Const.Aplha), damping (\f$\beta\f$),
 *                      \c Const.Beta) and stiffness matrices (\f$\gamma\f$), \c Const.Gamma).
 *
 * \post \c IGain is a symmetric matrix in general storage with only the upper part referenced (Lower part in
 *       FORTRAN routines). It contains the result of: \f[\mathcal{G} = \lambda(\alpha\mathcal{M} +
 *       \beta\mathcal{C} + \gamma\mathcal{C})^{-1}\f]
 *
 * \sa MatrixVector_t, MatrixVector_Sp_t and Scalars_t.
 */
void IGainMatrix_Pardiso_Sp( MatrixVector_t *const IGain, const MatrixVector_Sp_t *const Mass,
			     const MatrixVector_Sp_t *const Damp, const MatrixVector_Sp_t *const Stiff,
			     const Scalars_t Const );

/**
 * \brief Computes the gain matrix and its inverse. MPI version.
 *
 * This routine calculates the so called Gain Matrix and its inverse through:
 *
 * \f[\mathcal{G}^{-1} = \lambda(\alpha\mathcal{M} + \beta\mathcal{C} +
 * \gamma\mathcal{C})^{-1}\f]
 *
 * where:
 * - \f$\mathcal{M}\f$ is the mass matrix,
 * - \f$\mathcal{C}\f$ is the proportional viscous damping matrix,
 * - \f$\mathcal{K}\f$ is the stiffness matrix,
 * - \f$\mathcal{G}\f$ is the gain matrix,
 * - and \f$\alpha\f$, \f$beta\f$, \f$\gamma\f$ and \f$\lambda\f$ are the parameters that multiply the mass,
 *   damping, stiffness and matrices and the result of the matrix inversion respectively.
 *
 * It makes use of BLAS and LAPACK routines to perform the linear algebra operations.
 *
 * \pre 
 * - All elements of type \c PMatrixVector_t must be properly intialised through the PMatrixVector_Create()
 *   routine.
 * - The matrices must be symmetrical and only the upper part will be referenced (lower part in FORTRAN
 *   routines).
 * - The dimensions of the matrices must be the identical.
 * - The values in scalars must be properly initialised and they depend on the used formulation
 *   \cite Dorka_1998.
 *
 * \param[in,out] IGain The inverse of the gain matrix. As an input, only the size of the matrix is
 *                      referenced.
 * \param[in]     Mass  The mass matrix.
 * \param[in]     Damp  The proportional viscous damping matrix.
 * \param[in]     Stiff The stiffness matrix.
 * \param[in]     Const Scalars that multiply the mass (\f$\alpha\f$), \c Const.Aplha), damping (\f$\beta\f$),
 *                      \c Const.Beta), stiffness (\f$\gamma\f$), \c Const.Gamma) matrices and the inverted
 *                      matrix (\f$\lambda\f$), \c Const.Lambda) respectively.
 *
 * \post \c IGain is a distributed symmetric matrix in general storage with only the upper part referenced
 *       (Lower part in FORTRAN routines). It contains the result of: \f[\mathcal{G}^{-1} =
 *       \lambda(\alpha\mathcal{M} + \beta\mathcal{C} + \gamma\mathcal{C})^{-1}\f]
 *
 * \sa PMatrixVector_t and Scalars_t.
 */
void IGainMatrix_MPI( PMatrixVector_t *const IGain, PMatrixVector_t *const Mass, PMatrixVector_t *const Damp,
		      PMatrixVector_t *const Stiff, const Scalars_t Const );

void IGainMatrix_Float2Double( MatrixVector_t *const IGain, const MatrixVector_t *const Mass, 
		  const MatrixVector_t *const Damp, const MatrixVector_t *const Stiff,
			       const Scalars_t Const );

void IGainMatrix_Float2Double_PS( MatrixVector_t *const IGain, const MatrixVector_t *const Mass,
		     const MatrixVector_t *const Damp, const MatrixVector_t *const Stiff,
				  const Scalars_t Const );

#endif /* GAINMATRIX_H_ */

