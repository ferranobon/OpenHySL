/**
 * \file Auxiliary_Math.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 28th of February 2013
 * \todo Add support for packaged storages to decrease the memory use
 *
 * \brief MatrixVector_t creation and manipulation prototypes.
 *
 * This file contains the prototypes of those functions involved in creating/destroying
 * dense matrices and vectors.
 */

#ifndef AUXILIARY_MATH_H_
#define AUXILIARY_MATH_H_

#include "MatrixVector.h"

/**
 * \brief Returns the maximum of two values
 *
 * The maximum of two double values is returned at the end of the function.
 *
 * \param a First value.
 * \param b Second value.
 * \return max(a,b).
 */
int Max ( const int a, const int b );

/**
 * \brief Returns the minimum of two values
 *
 * The minimum of two double values is returned at the end of the function.
 *
 * \param a First value.
 * \param b Second value.
 * \return min(a,b).
 */
int Min ( const int a, const int b );

/**
 * \brief Generation of a Identity Matrix.
 *
 * The identity matrix with sizes equal to \c Rows and \c Cols (number of rows and columns
 * respectively) is generated. The output format is in general storage.
 * 
 * \pre The number of rows must be equal to the number of columns.
 *
 * \param[in] Rows The number of rows.
 * \param[in] Cols The number of columns.
 * \return A identity matrix of \f$Size= Rows*Cols\f$ in general storage.
 *
 * \post The deallocation of memory should be performed through MatrixVector_Destroy().
 */
MatrixVector_t Generate_IdentityMatrix( int Rows, int Cols );

/**
 * \brief Returns a 0-based index of a matrix in packed storage containing the upper
 * triangular part and in row major order given its coordinates.
 *
 * Given the coordinates of an element (1-based index), it returns the position (0-based
 * index) where the element is stored assuming that the upper triangular part is packed in
 * rows. Therefore:
 *
 * If \f$i \geq j\f$, \f$\mathcal a_{ij}\f$ is stored in \f$\mathcal A(i + (2n-j)(j-1)/2)\f$ or in
 * \f$\mathcal A(j + (2n-i)(i-1)/2)\f$ if \f$i < j\f$
 *
 * \pre
 * - \f$1 \leq RowIndex \geq n\f$ and \f$1 \leq ColIndex \geq n\f$.
 * - \em n is the number of rows or columns of the matrix in packed storage.
 *
 * \param[in] RowIndex Row coordinate \em i.
 * \param[in] ColIndex Column coordinate \em j.
 * \param[in] n        Number of rows or columns of the matrix in packed storage.
 * \returns   Position (0-based index) where \f$a_{ij}\f$ is stored in the packed matrix
 *            \f$\mathcal A\f$.
 */ 
unsigned int MatrixVector_ReturnIndex_UPS( const unsigned int RowIndex, const unsigned int ColIndex, const int n );

/**
 * \brief Returns a 0-based index of a matrix in packed storage containing the lower
 * triangular part and in row major order given its coordinates.
 *
 * Given the coordinates of an element (1-based index), it returns the position (0-based
 * index) where the element is stored assuming that the lower triangular part is packed in
 * rows. Therefore:
 *
 * If \f$j \geq i\f$, \f$\mathcal a_{ij}\f$ is stored in \f$\mathcal A(i + j*(j-1)/2)\f$ or in
 * \f$\mathcal A(j + i*(i-1)/2)\f$ if \f$j < i\f$
 *
 * \pre
 * - \f$1 \leq RowIndex \geq n\f$ and \f$1 \leq ColIndex \geq n\f$.
 * - \em n is the number of rows or columns of the matrix in packed storage.
 *
 * \param[in] RowIndex Row coordinate \em i.
 * \param[in] ColIndex Column coordinate \em j.
 * \param[in] n        Number of rows or columns of the matrix in packed storage.
 * \returns   Position (0-based index) where \f$a_{ij}\f$ is stored in the packed matrix
 *            \f$\mathcal A\f$.
 */
unsigned int MatrixVector_ReturnIndex_LPS( const unsigned int RowIndex, const unsigned int ColIndex );

#endif /* AUXILIARY_MATH_H_ */
