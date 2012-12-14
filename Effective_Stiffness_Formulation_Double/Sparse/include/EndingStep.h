/**
 * \file EndingStep.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 10th of August 2011
 *
 * \brief Prototypes of the functions used during the Compute U0 phase.
 *
 * This file contains the prototypes of the functions that are called during the Ending Step phase of the substructure
 * algorithm. For the moment, only the joining of the non-coupling and the coupling part routine has been considered.
 *
 * \todo Evaluate the performance impact (overhead) of writing routines to calculate the acceleration, velocity and displacement in this
 * phase as a separate routines.
 */

#ifndef ENDINGSTEP_H_
#define ENDINGSTEP_H_

#include "MatrixVector.h"
#include "Initiation.h"

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
void JoinNonCouplingPart( MatrixVector *const VecXm, const MatrixVector *const Keinv_m, const MatrixVector *const fcprevsub, MatrixVector *const Vec, const Coupling_Node *const CNodes );

/**
 * \brief Computes the new acceleration.
 *
 * The new acceleration is computed according to...........
 *
 * \param[in] DispTdT
 * \param[in] DispT
 * \param[in] VelT
 * \param[in] AccT
 * \param[in] a0
 * \param[in] a2
 * \param[in] a3
 * \param[out] AccTdT
 */
void Compute_Acceleration( const MatrixVector *const DispTdT, const MatrixVector *const DispT, const MatrixVector *const VelT,
			   const MatrixVector *const AccT, const double a0, const double a2, const double a3,
			   MatrixVector *const AccTdT );

void Compute_Velocity( const MatrixVector *const VelT, const MatrixVector *const AccT, const MatrixVector *const AccTdT,
		       const double a6, const double a7, MatrixVector *const VelTdT );

void Compute_Force_Error( const MatrixVector *const Mass, const MatrixVector *const Damp, const MatrixVector *Stiff,
			  const MatrixVector *const AccTdT, const MatrixVector *const VelTdT, const MatrixVector *const DispTdT,
			  const MatrixVector *const fc, const MatrixVector *const LoadTdT, MatrixVector *const fu );

#if _SPARSE_
void JoinNonCouplingPart_Sparse( MatrixVector *const VecXm, const Sp_MatrixVector *const Keinv_m, const MatrixVector *const fcprevsub, MatrixVector *const Vec, const Coupling_Node *const CNodes );

void Compute_Force_Error_Sparse( const Sp_MatrixVector *const Mass, const Sp_MatrixVector *const Damp, const Sp_MatrixVector *Stiff,
			  const MatrixVector *const AccTdT, const MatrixVector *const VelTdT, const MatrixVector *const DispTdT,
				 const MatrixVector *const fc, const MatrixVector *const LoadTdT, MatrixVector *const fu );
#endif

#endif /* ENDINGSTEP_H_ */
