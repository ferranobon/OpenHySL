/**
 * \file Auxiliary_Math.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 28th of February 2013
 *
 * \brief MatrixVector_t creation and manipulation prototypes.
 *
 * This file contains the prototypes of those functions involved in creating/destroying
 * dense matrices and vectors.
 *
 * \todo{Document norm and signum functions}
 */

#ifndef AUXILIARY_MATH_H_
#define AUXILIARY_MATH_H_

#include "MatrixVector.h"

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0 - EPS)

/**
 * \brief Returns the maximum of two values
 *
 * The maximum of two integer values is returned at the end of the function.
 *
 * \param a First value.
 * \param b Second value.
 * \return max(a,b).
 */
int Max ( const int a, const int b );

/**
 * \brief Returns the minimum of two values
 *
 * The minimum of two integer values is returned at the end of the function.
 *
 * \param a First value.
 * \param b Second value.
 * \return min(a,b).
 */
int Min ( const int a, const int b );

/**
 * \brief Returns the sign of a floating point value.
 *
 * The function returns a 1.0, 0.0 or a -1.0 depending on the sign of the entered floating point value:
 * - \f$1.0~if~num > 0.0\f$.
 * - \f$0.0~if~num = 0.0\f$.
 * - \f$-1.0~if~num < 0.0\f$.

 * \param[in] num Number whose sign needs to be analysed.
 * \returns       Returns the sign of the entered floating point value.
 */
HYSL_FLOAT signum ( const HYSL_FLOAT num );

/**
 * \brief Calculates the norm of a vector.
 *
 * The norm of a vector \f$u\f$ of length \f$n\f$ is calculated:
 *
 * \f[||u|| = \sum_1^n u_i^2\f]
 *
 * \pre
 * - Vector must point to an allocated memory of \f$n \geq Length\f$ elements.
 *
 * \param length Length of the vector.
 * \param Vector Vector whose norm must be computed.
 * \return The norm of the given vector.
 * 
 */
HYSL_FLOAT norm ( const int length, const HYSL_FLOAT *const Vector );

/**
 * \brief Generation of a Identity Matrix.
 *
 * The identity matrix with sizes equal to \c Rows and \c Cols (number of rows and columns respectively) is
 * generated. The output format is in general storage.
 * 
 * \pre The number of rows must be equal to the number of columns.
 *
 * \param[in] Rows The number of rows.
 * \param[in] Cols The number of columns.
 * \return A identity matrix of \f$Size = Rows*Cols\f$ in general storage.
 *
 * \post The deallocation of memory should be performed through MatrixVector_Destroy().
 */
MatrixVector_t Generate_IdentityMatrix( int Rows, int Cols );

/**
 * \brief Returns a 0-based index of a matrix in packed storage containing the upper triangular part and in
 * row major order given its coordinates.
 *
 * Given the coordinates of an element (1-based index), it returns the position (0-based index) where the
 * element is stored assuming that the upper triangular part is packed in rows. Therefore:
 *
 * If \f$i \geq j\f$, \f$\mathcal a_{ij}\f$ is stored in \f$\mathcal A(i + (2n-j)(j-1)/2)\f$ or in \f$\mathcal
 * A(j + (2n-i)(i-1)/2)\f$ if \f$i < j\f$
 *
 * \pre
 * - \f$1 \leq RowIndex \geq n\f$ and \f$1 \leq ColIndex \geq n\f$.
 * - \em n is the number of rows or columns of the matrix in packed storage.
 *
 * \param[in] RowIndex Row coordinate \em i.
 * \param[in] ColIndex Column coordinate \em j.
 * \param[in] n        Number of rows or columns of the matrix in packed storage.
 * \returns   Position (0-based index) where \f$a_{ij}\f$ is stored in the packed matrix \f$\mathcal A\f$.
 */ 
unsigned int MatrixVector_ReturnIndex_UPS( const unsigned int RowIndex, const unsigned int ColIndex, const int n );

/**
 * \brief Returns a 0-based index of a matrix in packed storage containing the lower triangular part and in
 * row major order given its coordinates.
 *
 * Given the coordinates of an element (1-based index), it returns the position (0-based index) where the
 * element is stored assuming that the lower triangular part is packed in rows. Therefore:
 *
 * If \f$j \geq i\f$, \f$\mathcal a_{ij}\f$ is stored in \f$\mathcal A(i + j*(j-1)/2)\f$ or in \f$\mathcal A(j
 * + i*(i-1)/2)\f$ if \f$j < i\f$
 *
 * \pre
 * - \f$1 \leq RowIndex \geq n\f$ and \f$1 \leq ColIndex \geq n\f$.
 * - \em n is the number of rows or columns of the matrix in packed storage.
 *
 * \param[in] RowIndex Row coordinate \em i.
 * \param[in] ColIndex Column coordinate \em j.
 * \returns   Position (0-based index) where \f$a_{ij}\f$ is stored in the packed matrix \f$\mathcal A\f$.
 */
unsigned int MatrixVector_ReturnIndex_LPS( const unsigned int RowIndex, const unsigned int ColIndex );

/**
 * \brief Computes the Eigenvalues and Eigenvectors the problem \f$\mathcal A*\mathcal B*\vec x = \lambda \vec
 * x\f$.
 *
 * The routine computes de eigenvalues and eigenvectors of a real symmetric-definite eigenproblem, of the form
 * \f$\mathcal A*\mathcal B*\vec x = \lambda \vec x\f$. The eigenvalues and the corresponding eigenvectors are
 * sorted in ascendent order.
 *
 * It makes use of BLAS and LAPACK routines to perform the linear algebra operations.
 *
 * \pre
 * - All elements of type \c MatrixVector_t must be properly intialised through the MatrixVector_Create()
 *   routine.
 * - The matrices must be symmetrical and only the upper part will be referenced (lower part in FORTRAN
 *   routines).
 * - The dimensions of the vectors and matrices must be coherent since it will not be checked in the routine.
 *
 * \param[in]  MatrixA       Matrix \f$\mathcal A\f$.
 * \param[in]  MatrixB       Matrix \f$\mathcal B\f$.
 * \param[out] Eigenvalues  Vector containing the eigenvalues of the \f$\mathcal A*\mathcal B*\vec x =
 *                           \lambda \vec x\f$ problem in ascendant order.
 * \param[out] Eigenvectors Matrix containing the eigenvectors of the \f$\mathcal A*\mathcal B*\vec x =
 *                           \lambda \vec x\f$ problem in ascendant order.
 *
 * \post
 * - \c Eigenvalues and \c Eigenvectors contain the eigenvalues and eigenvectors of the symmetric-definite
 *   eigen problem of the form \f$\mathcal A*\mathcal B*\vec x = \lambda \vec x\f$. They are sorted in
 *   increasing order.
 * 
 * \sa MatrixVector_t.
 */
void Compute_Eigenvalues_Eigenvectors ( MatrixVector_t *const MatrixA, MatrixVector_t *const MatrixB, MatrixVector_t *const Eigenvalues, MatrixVector_t *const Eigenvectors );

/**
 * From Numerical receipes in C \cite[Num_Recipes]..
 */
HYSL_FLOAT Gaussian_Deviate( const HYSL_FLOAT *const mu, const HYSL_FLOAT *const sigma, long int *const idum );
/**
 * From Numerical receipes in C \cite[Num_Recipes]..
 */
HYSL_FLOAT RandomNumber( long int *const idum );

/**
 * Extrapolates a value from a given set of data. Routine based on \c polint from \cite[Num_Recipes].
 * 
 * Given arrayx \f$X[1..n]\f$ and \f$Y[1..n]\f$ and given a value \f$x\f$, this routine returns a value
 * \f$y\f$, and an error estimate \f$d_y\f$. If \f$P(x)\f$ is the polynomial of degree \f$N-1\f$ such that
 * \f$P(X_i) \f$ = Y_i\f$, \f$ = 1,...,n\f$, then the returned value \f$y = P(x)\f$.
 *
 * \pre 
 * - All elements of type \c MatrixVector_t must be properly intialised through the MatrixVector_Create()
 *   routine and be of the same length.
 * 
 * \param[in]  X  Array of x values to interpolate.
 * \param[in]  Y  Array of y values to interpolate.
 * \param[in]  x  Value whose result has to be extrapolated.
 * \param[out] y  Extrapolated value \f$P(x) = y\f$.
 * \param[out] dy Error estimate of the operation \f$P(x) = y\f$.
 *
 * \sa MatrixVector_t.
 */
void Interpolate_Extrapolate( const MatrixVector_t *const X, const MatrixVector_t *const Y, const HYSL_FLOAT x, HYSL_FLOAT *const y, HYSL_FLOAT *const dy );

#endif /* AUXILIARY_MATH_H_ */