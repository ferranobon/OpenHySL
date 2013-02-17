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
 * The explicit part of the new state is computed. This routine does not depend on a specific formulation and therefore it can
 * be used in displacement, velocity and acceleration control formulations.
 * 
 * \f[\vec n_0^{t + \Delta t} = \mathcal{G}^{-1}(\vec f_{eff}^t + \vec l_i^t - \vec f_{err}^t)\f]
 *
 * where:
 * - \f$\vec n_0^{t + \Delta t}\f$ is the explicit part of the new state vector at time \f$t + \Delta t\f$,
 * - \f$\mathcal{G}\f$ is the gain matrix,
 * - \f$\vec f_{eff}^t\f$ is the effective force vector at time \f$t\f$,
 * - \f$\vec l_i^t\f$ is the input load vector at time \f$t\f$,
 * - and \f$\vec f_{err}^t\f$ is the error compensation force at time \f$t\f$.
 *
 * It makes use of BLAS routines to perform the lineal algebra operations.
 * 
 * \pre
 * - All elements of type \c MatrixVector must be properly initialised through the Init_MatrixVector() routine.
 * - \c Gain must be a symmetrical matrix in general storage. Only the upper part will be referenced (lower part in FORTRAN).
 * - The size of the vectors and matrices must be coherent since it will not be checked in the routine.
 * 
 * \param[in] IGain the inverted gain matrix \f$\mathcal{G}^{1}\f$.
 * \param[in] Eff_ForceT the effective force vector \f$\vec f_{eff}^t\f$.
 * \param[in] In_LoadT the input load vector \f$\vec l_i^t\f$.
 * \param[in] Err_ForceT the error force vector \f$\vec f_{err}^t\f$.
 * \param[in,out] Tempvec is a temporal vector that helps in the calculations. This is included in order to avoid
 * allocating and deallocating memory each step. On input only the number of rows is used.
 * \param[out] VecTdT_0 is new explicit displacement \f$\vec n_0^{t + \Delta t}\f$.
 *
 * \post
 * - \c VecTdT_0 is the result of the operation:
 *
 * \f[\vec n_0^{t + \Delta t} = \mathcal{G}^{-1}(\vec f_{eff}^t + \vec l_i^t + \vec f_{err}^t)\f]
 *
 * \sa MatrixVector and EffK_Calc_Effective_Force().
 */
void Compute_NewState( const MatrixVector *const IGain, const MatrixVector *const Eff_ForceT, const MatrixVector *const In_LoadT,
		       const MatrixVector *const Err_ForceT, MatrixVector *const Tempvec, MatrixVector *const VecTdT_0 );

/**
 * \brief Joins the non-coupling of a vector.
 *
 * The non-coupling (\f$Order-OrderC\f$) part of the vector is added to the global vector of size (\f$Order\f$). It performs the
 * operation:
 *
 * \f[\hat{\vec n}^{t + \Delta t} = \mathcal{\hat G}^{-1} (\hat{\vec f}_r^{n_{sub}-1} + \hat{\vec f}_s^{n_{sub}-1}) = \mathcal{\hat G}^{-1} \hat{\vec f}_c^{n_{sub}-1}\f]
 * \f[\hat{\vec n}^{t + \Delta t} \longrightarrow \vec n^{t + \Delta t}\f]
 *
 * where:
 * - \f$\hat{\vec n}^{t + \Delta t}\f$ is the new displacement of the non-coupling nodes after applying the restoring forces,
 * - \f$\mathcal{G}\f$ is a part of the gain matrix with the non-coupling elements of a column with a node marked as coupling node.
 * - \f$f_r^{i-1}\f$ and \f$f_s^{i-1}\f$ are the calculated and measured force vectors at sub-step \f$n_{sub}-1\f$, with $n_{sub} being the number of sub-steps,
 * - and \f$n^{t + \Delta t}\f$ is the new displacement vector including both parts: the coupling and non-coupling part.
 *
 * It makes use of BLAS routines to perform the lineal algebra operations.
 *
 * \pre
 * - All elements of type \c MatrixVector must be properly initialised through the Init_MatrixVector() routine.
 * - \c Gain_m must be a matrix of size \f$Order\f$x\f$OrderC\f$ and have the output of Build_MatrixXm().
 * - \c VecTdT_m must be of length \f$Order - OrderC\f$.
 * - \c VecTdT must be of length \f$Order\f$.
 * - \c fcprevsub must be of length \f$OrderC\f$.
 * - The coupling nodes are assumed to be consecutive and in increasing equation number.
 * - \e Order is the number of rows and columns of the global matrix.
 * - \e OrderC is the number of coupling degrees of freedom.
 *
 * It makes use of BLAS routines to perform the linear algebra operations.
 *
 * \param[in] VecTdT_m The non-coupling part of the vector.
 * \param[in] Gain_m The non-coupling part of the vector.
 * \param[in] fcprevsub The non-coupling part of the vector.
 * \param[in] CNodes Structure containing the coupling nodes.
 * \param[in,out] VecTdT The non-coupling part of the vector. Only the non-coupling positions are modified.
 *
 * \post
 * - \c VecTdT contains the non-coupling part of the vector, leaving the coupling nodes untouched. This is accomplished through: 
 *
 * \f[\hat{\vec n}^{t + \Delta t} = \mathcal{\hat G}^{-1} (\hat{\vec f}_r^{n_{sub}-1} + \hat{\vec f}_s^{n_{sub}-1}) = \mathcal{\hat G}^{-1} \hat{\vec f}_c^{n_{sub}-1}\f]
 * \f[\hat{\vec n}^{t + \Delta t} \longrightarrow \vec n^{t + \Delta t}\f]
 *
 *
 * \sa MatrixVector.
 */
void Join_NonCouplingPart( MatrixVector *const VecTdT_m, const MatrixVector *const Gain_m,
			   const MatrixVector *const fcprevsub, const Coupling_Node *const CNodes,
			   MatrixVector *const VecTdT );


#endif /* COMMON_FORMULATION_H_ */
