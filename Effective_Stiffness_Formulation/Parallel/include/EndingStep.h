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

#include "PMatrixVector.h"

/**
 * \brief Joins the non-coupling of a vector.
 *
 * The non-coupling (\f$Order-OrderC\f$) part of the vector is added to the global vector of size (\f$Order\f$).
 * This routine makes use of the level 1 PBLAS pscopy_().
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
 * \param[in] PosCouple The position of the first coupling node.
 * \param[in] OrderC The number of coupling nodes. In case that \f$OrderC > 1\f$, the routine assumes that they are consecutive.
 *
 * \post
 * - \c Vec contains the non-coupling part of the vector, leaving the coupling nodes untouched.
 *
 * \sa MatrixVector.
 */
void JoinNonCouplingPart( PMatrixVector *const VecXm, PMatrixVector *const Keinv_m, PMatrixVector *const fcprevsub, PMatrixVector *const Vec, const int PosCouple, const int OrderC );

void Compute_Acceleration( PMatrixVector *const DispTdT, PMatrixVector *const DispT, PMatrixVector *const VelT,
			   PMatrixVector *const AccT, const float a0, const float a2, const float a3,
			   PMatrixVector *const AccTdT );

void Compute_Velocity( PMatrixVector *const VelT, PMatrixVector *const AccT, PMatrixVector *const AccTdT,
		       const float a6, const float a7, PMatrixVector *const VelTdT );

void Compute_Force_Error( PMatrixVector *const Mass, PMatrixVector *const Damp, PMatrixVector *Stiff,
			  PMatrixVector *const AccTdT, PMatrixVector *const VelTdT, PMatrixVector *const DispTdT,
			  PMatrixVector *const fc, PMatrixVector *const LoadTdT, PMatrixVector *const fu );


#endif /* ENDINGSTEP_H_ */
