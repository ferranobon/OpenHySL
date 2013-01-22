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

#ifndef PRECALCULATIONS_H_
#define PRECALCULATIONS_H_

/**
 * \brief Reads the velocity and displacement of an earthquake from a file.
 *
 * The data from an earthquake is readen from the specified file and its components are stored in different variables, namely:
 * displacement and velocity. If the specified file does not exist, an appropriate error is displayed
 * and the program exists with \c EXIT_FAILURE.
 *
 * \pre
 * - The input text file should be in ASCII format and it must contain no header and four columns in the following order:
 *     -# It can contain anything, since it will not be stored, it is recommended that it has the step number to detect any possible
 * errors (see the next precondition).
 *     -# Acceleration in \f$[m/s^2]\f$.
 *     -# Velocity in \f$[m/s]\f$.
 *     -# Displacement in \f$[m]\f$.
 * - Note also that the number of steps in the file should be equal to the number of steps defined in InitConstants(). No check
 * is performed.
 * - The arrays must be properly initialised and its length should be \f$l = Number~of~Steps\f$.
 *
 * \param[out] Velocity Array that will store the third column of the data file, that is, the velocity of the recorded earthquake.
 * \param[out] Displacement Array that will store the fourth column of the data file, that is, the displacement of the recorded earthquake.
 * \param[in] NumSteps The number of steps in which the earthquake is divided.
 * \param[in] Filename The name of the file where the earthquake is recorded.
 *
 * \post
 * - \c Velocity and \c Displacement will store the velocity and displacement of the recorded earthquake.
 * \sa InitConstants()
 */
void ReadDataEarthquake_AbsValues( double *Velocity, double *Displacement, const unsigned int NumSteps, const char *Filename );

/**
 * \brief Reads the accelerations of an earthquake from a file.
 *
 * The acceleration record from an earthquake is readen from the specified file. If it does not exist, an appropriate
 * error is displayed and the program exists with \c EXIT_FAILURE.
 *
 * \pre
 * - The input text file should be in ASCII format and it must contain no header and four columns in the following order:
 *      -# It can contain anything, since it will not be stored, it is recommended that it has the step number to detect any possible
 * errors (see the next precondition).
 *      -# Acceleration in \f$[m/s^2]\f$.
 *      -# Velocity in \f$[m/s]\f$.
 *      -# Displacement in \f$[m]\f$.
 * - Note also that the number of steps in the file should be equal to the number of steps defined in InitConstants(). No check
 * is performed.
 * - The arrays must be properly initialised and its length should be \f$l = Number~of~Steps\f$.
 *
 * \param[out] Acceleration Array that will store the second column of the data file, that is, the acceleration.

 * \param[in] NumSteps The number of steps in which the earthquake is divided.
 * \param[in] Filename The name of the file where the earthquake is recorded.
 *
 * \post
 * - \c Acceleration will store the velocity and displacement of the recorded earthquake.
 * \sa InitConstants()
 */
void ReadDataEarthquake_RelValues( double *Acceleration, const unsigned int NumSteps, const char *Filename );

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
 * \sa MatrixVector.
 */
void CopyDiagonalValues( const MatrixVector *const Mat, MatrixVector *const Vec );

/**
 * \brief Calculates the input load as an absolute value.
 *
 * This routine calculates the input load that is required at the beginning of each step, that is
 * \f$ \{l_{i+1}\} = \{ld_g\} + \{ld_v\}\f$, where:
 * - \f$\{ld_g\} = [M]\cdot{Displacement_i}\f$ is the load caused by the ground motion,
 * - \f$\{ld_v\} = [C]\cdot{Velocity_i}\f$ is the load caused by ground velocity,
 * - \e i denotes the step.
 * It makes use of the level 2 BLAS routine dsymv_().
 *
 * \pre
 * - All elements of type MatrixVector must be properly initialised through the Init_MatrixVector() routine.
 * - In the case of matrices, they are supposed to be symmetrical and they must contain at least the upper elements in
 * general storage format.
 * - The number of rows of the vectors must be indicative of their length.
 * - The size of the elements must be coeherent, since it will not be checked in the routine: the number of rows of the
 * vectors (\c InLoad, \c D and \c V) must be \f$NumRows_{Vec} = NumRows_{Mat}\f$.
 *
 * \param[in,out] InLoad Will contain the input load as an absolute value at the exit of the function. As an input,
 * only the size of the vector is referenced, not its elements.
 *
 * \param[out] InLoad On output it contains the input load \f$ \{l_{i+1}\} = \{ld_g\} + \{ld_v\}\f$.
 * \param[in] Stif The Stiffness matrix.
 * \param[in] Damp The Viscous Damping matrix.
 * \param[in] D Vector containing the ground motion of the earthquake at a certain step.
 * \param[in] V Vector containing the ground velocity of the earthquake at a certain step.
 *
 * \post \c InLoad has the value of \f$ \{l_{i+1}\} = \{ld_g\} + \{ld_v\}\f$.
 *
 * \sa MatrixVector.
 */
void Calc_Input_Load_AbsValues( MatrixVector *const InLoad, const MatrixVector *const Stif, const MatrixVector *const Damp, const MatrixVector *const D, const MatrixVector *const V );

/**
 * \brief Calculates the input load as a relative value.
 *
 * This routine calculates the input load that is required at the beginning of each step, that is
 * \f$ \{l_{i+1}\} = \{ld_a\}\f$, where:
 * - \f$\{ld_a\} = [M]\cdot{Acceleration_i}\f$ is the load caused by the ground acceleration,
 * - \e i denotes the step.
 * It makes use of the level 2 BLAS routine dsymv_().
 *
 * \pre
 * - All elements of type MatrixVector must be properly initialised through the Init_MatrixVector() routine.
 * - In the case of matrices, they are supposed to be symmetrical and they must contain at least the upper elements in
 * general storage format.
 * - The number of rows of the vectors must be indicative of their length.
 * - The size of the elements must be coeherent, since it will not be checked in the routine: the number of rows of the
 * vectors (\c InLoad and \c A) must be \f$NumRows_{Vec} = NumRows_{Mat}\f$.
 *
 * \param[in,out] InLoad Will contain the input load as a relative value at the exit of the function. As an input,
 * only the size of the vector is referenced, not its elements.
 * \param[out] InLoad On output it contains the input load \f$ \{l_{i+1}\} = \{ld_g\} + \{ld_v\}\f$.
 * \param[in] Mass The Mass matrix.
 * \param[in] A Vector containing the ground acceleration of the earthquake at a certain step.
 *
 * \post \c InLoad has the value of \f$ \{l_{i+1}\} = \{ld_a\}\f$.
 *
 * \sa MatrixVector.
 */
void Calc_Input_Load_RelValues( MatrixVector *const InLoad, const MatrixVector *const Mass, const MatrixVector *const A );

void Apply_LoadVectorForm ( MatrixVector *const Vector, const MatrixVector *const LoadForm, const double Value );

#if _SPARSE_
void Calc_Input_Load_AbsValues_Sparse( MatrixVector *const InLoad, const Sp_MatrixVector *const Stif, const Sp_MatrixVector *const Damp, const MatrixVector *const D, const MatrixVector *const V );

void Calc_Input_Load_RelValues_Sparse( MatrixVector *const InLoad, const Sp_MatrixVector *const Mass, const MatrixVector *const A );
#endif

#endif /* PRECALCULATIONS_H_ */
