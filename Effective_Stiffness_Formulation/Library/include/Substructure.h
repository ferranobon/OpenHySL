#ifndef SUBSTRUCTURE_H_
#define SUBSTRUCTURE_H_

#include "MatrixVector.h"
#include "Substructure_CouplingNodes.h"

#define NUM_TYPE_SUB  8  /*!< Number of recognized sub-structure types */

enum Substructure_Id { SIM_EXACT,     /*!< Simulate the substructure using the exact solution. */
		       SIM_UHYDE,     /*!< Simulate the substructure using the UHYDE-fbr device. */
		       SIM_MEASURED,  /*!< Simulate the substructure using measured values. */
		       EXP_ADWIN,     /*!< Run using ADwin */
		       REMOTE_TCP,    /*!< Remote node using TCP/IP connection. */
		       REMOTE_UDP,    /*!< Remote node using UDP connection. */
		       REMOTE_NSEP,   /*!< Remote node using NSEP protocol. */
		       REMOTE_OF      /*!< Remote node using OpenFresco. */
};

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
 * - All elements of type \c MatrixVector_t must be properly initialised through the MatrixVector_Create() routine.
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
 * \param[in] VecTdT_m The non-coupling part of the vector (displacement, velocity or acceleration depending on the formulation).
 * \param[in] Gain_m The non-coupling part of the gain matrix \f$\mathcal G\f$.
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
 * \sa MatrixVector_t.
 */
void Join_NonCouplingPart( MatrixVector_t *const VecTdT_m, const MatrixVector_t *const Gain_m,
			   const MatrixVector_t *const fcprevsub, const CouplingNode_t *const CNodes,
			   MatrixVector_t *const VecTdT );

void Substructure_Simulate( CouplingNode_t *const CNodes, double *Keinv, double *const u0c, double *const uc, double *const fcprev, double *const fc, const unsigned int NSubstep, const double DeltaT_Sub );

#endif /* SUBSTRUCTURE_H_ */
