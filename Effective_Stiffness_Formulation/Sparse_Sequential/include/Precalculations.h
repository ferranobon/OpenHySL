/**
 * \file Precalculations.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 10th of August 2011
 *
 * \brief Prototypes of the functions used in the precalculations phase.
 *
 * This file contains the prototypes of the functions that are called during the precalculations phase. This includes,
 * reading the earthquake data from a file, and perform some operations like Copying the diagonal values of the Mass matrix
 * and calculate the input load.
 */

#include "MatrixVector.h"

#ifndef PRECALCULATIONS_H_
#define PRECALCULATIONS_H_

/**
 * \brief Reads the data of an earthquake from a file.
 *
 * The data from an earthquake is readen from the specified file and its components are stored in different variables, namely:
 * displacement, velocity and acceleration for further use. If the specified file does not exist, an appropriate error is displayed
 * and the program exists with \c EXIT_FAILURE.
 *
 * \pre
 * - The input text file should be in ASCII format and it must contain no header and four columns in the following order:
 * 	-# It can contain anything, since it will not be stored, it is recommended that it has the step number to detect any possible
 * errors (see the next precondition).
 * 	-# Acceleration in \f$[m/s^2]\f$.
 *	-# Velocity in \f$[m/s]\f$.
 * 	-# Displacement in \f$[m]\f$.
 * - Note also that the number of steps in the file should be equal to the number of steps defined in InitConstants(). No check
 * is performed.
 * - The arrays must be properly initialised and its length should be \f$l = Number~of~Steps\f$.
 *
 * \param[out] Acceleration Array that will store the second column of the data file, that is, the acceleration of the recorded earthquake.
 * \param[out] Velocity Array that will store the third column of the data file, that is, the velocity of the recorded earthquake..
 * \param[out] Displacement Array that will store the fourth column of the data file, that is, the displacement of the recorded earthquake..
 * \param[in] NumSteps The number of steps in which the earthquake is divided.
 * \param[in] Filename The name of the file where the earthquake is recorded.
 *
 * \post
 * - \c Acceleration, \c Velocity and \c Displacement will store the acceleration, velocity and displacement of the recorded earthquake.
 * \sa InitConstants()
 */
void ReadDataEarthquake( float *Acceleration, float *Velocity, float *Displacement, const int NumSteps, const char *Filename );

/**
 * \brief Copies the diagonal values of a matrix into a vector.
 *
 * This routine copies the diagonal values of a square matrix into a vector. It makes use of the BLAS routine dcopy_().
 *
 * \pre
 * - The input matrix must be square and must be stored in general storage.
 * - The vector must be of size \f$l = Order\f$ where \e Order is the number of rows or columns of the input matrix.
 * - The number of rows of the vector must be indicative of its length.
 *
 * \param[in] Mat The matrix that contains the elements to be copied.
 * \param[in,out] Vec The vector where the diagonal elements of \c Mat will be stored. As an input, only the size of the vector
 * is referenced, not its elements.
 *
 * \post The vector \c Vec contains the diagonal elements of the matrix \c Mat.
 *
 * \sa Dense_MatrixVector.
 */
void CopyDiagonalValues( const Dense_MatrixVector *const Mat, Dense_MatrixVector *const Vec );

/**
 * \brief Calculates the input load.
 *
 * This routine calculates the input load that is required at the beginning of each step, that is
 * \f$ \{l_{i+1}\} = \{ld_g\} + \{ld_v\} + \{ld_a\} -[M]\cdot [I]\cdot \{Acceleration_i\}\f$, where:
 * - \f$\{ld_g\} = [M]\cdot{Displacement_i}\f$ is the load caused by the ground motion,
 * - \f$\{ld_v\} = [C]\cdot{Velocity_i}\f$ is the load caused by ground velocity,
 * - \f$\{ld_a\} = [C]\cdot{Acceleration_i}\f$ is the load caused by ground acceleration and
 * - \e i denotes the step.
 * It makes use of the level 2 BLAS routine dsymv_().
 *
 * \pre
 * - All elements of type Dense_MatrixVector must be properly initialised through the Init_Dense_MatrixVector() routine.
 * - In the case of matrices, they are supposed to be symmetrical and they must contain at least the upper elements in
 * general storage format.
 * - The number of rows of the vectors must be indicative of their length.
 * - The size of the elements must be coeherent, since it will not be checked in the routine: the number of rows of the
 * vectors (\c InLoad, \c DiagM, \c D, \c V and \c A) must be \f$NumRows_{Vec} = NumRows_{Mat}\f$.
 *
 * \param[in,out] InLoad Will contain the input load at the exit of the function. As an input, only the size of the vector
 * is referenced, not its elements.
 * \param[in] Stif The Stiffness matrix.
 * \param[in] Damp The Viscous Damping matrix.
 * \param[in] Mass The Mass matrix.
 * \param[in] DiagM Vector that has the diagonal elements of the Mass matrix.
 * \param[in] D Vector containing the ground motion of the earthquake at a certain step.
 * \param[in] V Vector containing the ground velocity of the earthquake at a certain step.
 * \param[in] A Vector containing the ground acceleration of the earthquake at a certain step.
 *
 * \post \c InLoad has the value of \f$ \{l_{i+1}\} = \{ld_g\} + \{ld_v\} + \{ld_a\} -[M]\cdot [I]\cdot\{Acceleration_i\}\f$.
 *
 * \sa Dense_MatrixVector.
 */
void Calc_Input_Load( Dense_MatrixVector *const InLoad, const Dense_MatrixVector *const Stif, const Dense_MatrixVector *const Damp, const Dense_MatrixVector *const Mass, const Dense_MatrixVector *const DiagM, const Dense_MatrixVector *const D, const Dense_MatrixVector *const V, const Dense_MatrixVector *const A );

#endif /* PRECALCULATIONS_H_ */
