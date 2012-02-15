/**
 * \file ComputeU0.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 10th of August 2011
 *
 * \brief Prototypes of the functions used during the Compute U0 phase.
 *
 * This file contains the prototypes of the functions that are called during the Compute U0 phase of the substructure
 * algorithm. This includes a function to calculate the input load and functions to decouple vectors into their non-coupling
 * and coupling part.
 *
 */

#ifndef COMPUTEU0_H_
#define COMPUTEU0_H_

#include "MatrixVector.h"

/**
 * \brief Calculate force vector.
 *
 * The force vector is calculated according to \f$f_i = l_i + f_c + \Delta_f\f$, where:
 * - \f$f_i\f$ is the force vector.
 * - \f$l_i\f$ is the input load.
 * - \f$\Delta_f\f$ is the force due to the PID error compensator.
 * - \f$i\f$ denotes the step.
 *
 * \pre
 * - All elements of type MatrixVector must be properly initialised through the Init_MatrixVector() routine.
 * - The size of the vectors must be coeherent, since it will not be checked in the routine.
 * - The number of rows of the vector \f$f_i\f$ must be indicative of its length.
 *
 * \param[in, out] fi The force vector at step \e i. As an input, only the size of the vector is referenced, not its elements.
 * \param[in] fc The coupling force at step \e i.
 * \param[in] li The input load at step \e i.
 * \param[in] Deltaf The force due to the PID error compensator.
 *
 * \post
 * - \f$f_i\f$ is the result of the operation \f$f_i = l_i + f_c + \Delta_f\f$.
 *
 * \sa MatrixVector.
 */
void Calculatefi( MatrixVector *const fi, const MatrixVector *const fc, const MatrixVector *const li, const MatrixVector *const Deltaf );

void EffK_Calc_Effective_Force( const MatrixVector *const Mass, const MatrixVector *const Damp,
				const MatrixVector *const Disp, const MatrixVector *const Vel,
				const MatrixVector *const Acc, MatrixVector *const Tempvec,
				const float a0, const float a1, const float a2,
				const float a3, const float a4, const float a5,
			       MatrixVector *const Eff_Force );
void EffK_ComputeU0( const MatrixVector *const Eff_Force, const MatrixVector *const In_Load,
		     const MatrixVector *const Err_Force, const float PID_P, const MatrixVector *const Keinv,
		     MatrixVector *const Tempvec, MatrixVector *const Disp0 );

/**
 * \brief Copies the non-coupling part a vector.
 *
 * The non-coupling part a vector (\f$Order - OrderC\f$) is copied. This routine makes use of the level 1 BLAS dcopy_().
 *
 * \pre
 * - The global vector (length \f$Order\f$) must be properly initialised through the Init_MatrixVector() routine.
 * - The non-coupling vector must be of length \f$Order- OrderC\f$ and properly initialised through the Init_MatrixVector() routine.
 * - The number of rows of the vectors must be indicative of their length.
 * - The coupling nodes are assumed to be consecutive (the first node is assumed to be PosCouple).
 * - \e Order is the number of rows of the input vector.
 * - \e OrderC is the number of coupling degrees of freedom.
 *
 * \param[in] VectorX The global vector.
 * \param[in,out] VectorXm The vector that will contain the non-coupling elements of \c VectorX. As an input,
 * only the size of the vector is referenced, not its elements.
 * \param[in] PosCouple The position of the first coupling node.
 * \param[in] OrderC The number of coupling nodes. In case that \f$OrderC > 1\f$, the routine assumes that they are consecutive.
 *
 * \post
 * - \c VectorXm contains only the non-coupling nodes of \c Vec.
 *
 * \sa MatrixVector.
 */
void CreateVectorXm( const MatrixVector *const VectorX, MatrixVector *const VectorXm, const int PosCouple, const int OrderC );

/**
 * \brief Copies the coupling part a vector.
 *
 * The coupling part a vector (\f$OrderC\f$) is copied. This routine makes use of the level 1 BLAS dcopy_().
 *
 * \pre
 * - The global vector (length \f$Order\f$) must be properly initialised through the Init_MatrixVector() routine.
 * - The number of rows of the vectors must be indicative of their length.
 * - The coupling vector must be of length \f$OrderC\f$ and properly initialised.
 * - The coupling nodes are assumed to be consecutive (the first node is assumed to be PosCouple).
 * - \e Order is the number of rows of the input vector.
 * - \e OrderC is the number of coupling degrees of freedom.
 *
 * \param[in] VecX The global vector.
 * \param[out] VecXc The vector that will contain the coupling elements of \c VectorX.
 * \param[in] PosCouple The position of the first coupling node.
 * \param[in] OrderC The number of coupling nodes. In case that \f$OrderC > 1\f$, the routine assumes that they are consecutive.
 *
 * \post
 * - \c VectorXc contains only the coupling nodes of \c Vec.
 *
 * \sa MatrixVector.
 */
void CreateVectorXc( const MatrixVector *const VecX, float *VecXc, const int PosCouple, int OrderC );

#endif /* COMPUTEU0_H_ */
