/**
 * \file Rayleigh.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 15th of March 2013
 *
 * \brief Rayleigh damping routines.
 *
 * \todo Add single precision routines.
 *
 * Routines for calculating the proportional viscous damping matrix using Rayleigh Damping. The routines make
 * use of the BLAS library to perform the linear algebra operations and they support both single and double
 * precision. Sparse BLAS operations are supported through the Intel MKL library.
 */

#ifndef RAYLEIGH_H_
#define RAYLEIGH_H_

#include "MatrixVector.h"
#include "MatrixVector_Sp.h"
#include "MatrixVector_MPI.h"

/**
 * \brief Stores the Rayleigh coefficients alpha and beta.
 */
typedef struct {
     double Alpha; /*!< \brief Coefficient that multiplies the mass matrix.*/
     double Beta;  /*!< \brief Coefficient that multiplies the stiffness matrix.*/
} Rayleigh_t;

/**
 * \brief Computes the proportional viscous damping matrix using Rayleigh damping. General storage version.
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
void Rayleigh_Damping( const MatrixVector_t *const Mass, const MatrixVector_t *const Stiff, MatrixVector_t *const Damp,
		       const Rayleigh_t *const Rayleigh );

/**
 * \brief Computes the proportional viscous damping matrix using Rayleigh damping. Packed storage version.
 *
 * This routine calculates the proportional viscous damping matrix using Rayleigh Damping and packed
 * storage. The operation performed is (see \cite Clough 1975 p. 234):
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
 * It makes use of BLAS routines to perform the linear algebra operations.
 *
 * \pre 
 * - All elements of type \c MatrixVector_t must be properly intialised through the MatrixVector_Create_PS()
 *   routine.
 * - The matrices must be symmetrical and in packed storage.
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
 * \post \c Damp is a symmetric matrix in packed storage. It contains the result of: \f[\mathcal{C} = \alpha
 *       \mathcal{M} \cdot \beta \mathcal{K}\f]
 *
 * \sa MatrixVector_t and Rayleigh_t.
 */
void Rayleigh_Damping_PS( const MatrixVector_t *const Mass, const MatrixVector_t *const Stiff, MatrixVector_t *const Damp,
			  const Rayleigh_t *const Rayleigh );

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
 * It makes use of the BLAS, LAPACK and Sparse BLAS routines from the Intel Math Kernel Library (\cite
 * MKL_2013) to perform the linear algebra operations. This routine requires the MKL library.
 *
 * \pre 
 * - All elements of type \c MatrixVector_Sp_t must be properly intialised through the
 *   MatrixVector_Create_Sp() routine.
 * - The matrices must in Intel's MKL CSR-\em three \em array \em variation and in one based index.
 * - The dimensions of the matrices must be identical.
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
 * \sa MatrixVector_Sp_t and Rayleigh_t.
 */
void Rayleigh_Damping_Sp( const MatrixVector_Sp_t *const Mass, const MatrixVector_Sp_t *const Stiff,
			  MatrixVector_Sp_t *const Damp, const Rayleigh_t *const Rayleigh );

/**
 * \brief Computes the proportional viscous damping matrix using Rayleigh damping. MPI version.
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
 * It makes use of PBLAS and ScaLAPACK routines to perform the linear algebra operations.
 *
 * \pre 
 * - All elements of type \c PMatrixVector_t must be properly intialised through the PMatrixVector_Create()
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
 * \post \c Damp is a distributed symmetric matrix in general storage with only the upper part referenced
 *       (Lower part in FORTRAN routines). It contains the result of: \f[\mathcal{C} = \alpha \mathcal{M}
 *       \cdot \beta \mathcal{K}\f]
 *
 * \sa PMatrixVector_t and Rayleigh_t.
 */
void Rayleigh_Damping_MPI( PMatrixVector_t *const Mass, PMatrixVector_t *const Stiff,
			   PMatrixVector_t *const Damp, const Rayleigh_t *const Rayleigh );

#endif /* RAYLEIGH_H_ */
