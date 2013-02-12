/**
 * \file Common_Formulation.h
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

#ifndef COMMON_FORMULATION_H_
#define COMMON_FORMULATION_H_

#include "MatrixVector.h"

/**
 * \brief Calculates the new state.
 *
 * The displacement is computed according to the formulation using the effective stiffness matrix in page 53 of (\cite Dorka_1998). This is the implicit displacement.
 * 
 * \f[\overrightarrow{u_0} = \mathcal{K}^{-1}_e(\overrightarrow{f_{eff}} + \overrightarrow{l_i} + \overrightarrow{f(e)})\f]
 * 
 * \pre
 * - All elements of type \c MatrixVector must be properly initialised through the Init_MatrixVector() routine.
 * - \c Keinv must be a symmetrical matrix in general storage. Only the upper part will be referenced (lower part in FORTRAN).
 * - The size of the vectors and matrices must be coherent, since it will not be checked in the routine.
 * 
 * \param[in] Eff_ForceT the effective force vector \f$\overrightarrow{f_{eff}}\f$.
 * \param[in] In_Load the input load vector \f$\overrightarrow{l_i}\f$.
 * \param[in] Err_Force the error force vector coming from the PID controller \f$\overrightarrow{f(e)}\f$.
 * \param[in] PID_P the proportional part of the PID controller.
 * \param[in] Keinv the inverted effective stiffness matrix \f$\mathcal{K}^{1}_e\f$.
 * \param[in,out] Tempvec is a temporal vector that helps in the calculations. This is included in order to avoid allocating and deallocating memory each step. On input
 * only the number of rows is used.
 * \param[out] Disp0 is new implicit displacement \f$\overrightarrow{u_0}\f$.
 *
 * \post
 * - \c Disp0 is the result of the operation \f$\overrightarrow{u_0} = \mathcal{K}^{-1}_e(\overrightarrow{f_{eff}} + \overrightarrow{l_i} + \overrightarrow{f(e)})\f$.
 *
 * \sa MatrixVector and EffK_Calc_Effective_Force().
 */
void EffK_ComputeU0( const MatrixVector *const Eff_Force, const MatrixVector *const In_Load,
		     const MatrixVector *const Err_Force, const double PID_P, const MatrixVector *const Keinv,
		     MatrixVector *const Tempvec, MatrixVector *const Disp0 );

/**
 * \brief Joins the non-coupling of a vector.
 *
 * The non-coupling (\f$Order-OrderC\f$) part of the vector is added to the global vector of size (\f$Order\f$).
 * This routine makes use of the level 1 BLAS dcopy_().
 *
 * \pre
 * - The global vector (length \f$Order\f$) must be properly initialised through the Init_MatrixVector() routine.
 * - The vector of length \f$Order- OrderC\f$ must contain the non-coupling part that the global vector requires.
 * - The coupling nodes are assumed to be consecutive (the first node is assumed to be PosCouple).
 * - The number of rows of the vectors must be indicative of their length.
 * - \e Order is the number of rows and columns of the input matrix.
 * - \e OrderC is the number of coupling degrees of freedom.
 *
 * \param[in] VecXm The non-coupling part of the vector.
 * \param[in,out] Vec The global vector. As an input, only the size of the vector is referenced, not its elements.
 * \param[in] CNodes Structure containing the coupling nodes.
 *
 * \post
 * - \c Vec contains the non-coupling part of the vector, leaving the coupling nodes untouched.
 *
 * \sa MatrixVector.
 */
void JoinNonCouplingPart( MatrixVector *const VecXm, const MatrixVector *const Keinv_m, const MatrixVector *const fcprevsub,
			  MatrixVector *const Vec, const Coupling_Node *const CNodes );


#endif /* COMMON_FORMULATION_H_ */
