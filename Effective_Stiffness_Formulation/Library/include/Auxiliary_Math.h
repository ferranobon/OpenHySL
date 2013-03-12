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

#endif /* AUXILIARY_MATH_H_ */
