/**
 * \file Substructure.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 13th of September 2013
 * 
 * \brief Routines for dealing with substructures and the sub-stepping process. Definition of supported
 * substructure types.
 *
 * This file contains routines that deal with the sub-stepping process or that involve substructures,
 * including sending the gain matrix to remote facilities or to the control hardware. It defines the
 * recognised types of substructures.
 *
 * \todo Implement sending the gain matrix when using OpenFresco in Substructure_SendGainMatrix(). What is
 * coded is an ugly hack.
 */
#ifndef SUBSTRUCTURE_H_
#define SUBSTRUCTURE_H_

#include "Substructure_CouplingNodes.h"

#if _MPI_
#include <mpi.h>
#include "MatrixVector_MPI.h"
#define MPI_NEW_STATE 1   /*!< \brief To send/receive MPI messages that contain the new value. */
#define MPI_PREV_FORCE 2  /*!< \brief To send/receive MPI messages that contain the coupling force at the
			   * previous substep. */
#define MPI_FORCE 3       /*!< \brief To send/receive MPI messages that contain the coupling force at the last
			   * substep. */
#endif

/**
 * \brief Sends the Gain matrix to the corresponding facility/controller.

 * \pre
 * - The Gain matrix \f$\mathcal G_c\f$ (order \f$n_c x n_c\f$ where \f$n_c\f$ is the number of coupling DOFs)
 *   must have been created through Substructure_MatrixXc(), Substructure_MatrixXc_PS(),
 *   Substructure_MatrixXc_MPI() routines if the complete gain matrix \f$\mathcal G\f$ (of order \f$n\f$ where
 *   \f$n > n_c\f$) is stored in general, packed storage or is a distributed matrix respectively. It must be a
 *   symmetric matrix with all its elements referenced.
 * - \c Order is the number of substructures/coupling DOFs.
 * - \c Substructure must be properly initialised through Substructure_ReadCouplingNodes().
 *
 * The portion of the gain matrix \f$\mathcal G_c\f$ that is affected by substructures is sent to the remote
 * facility or controller depending on the substructure type. Those DOFs that have a numerical substructure
 * assigned remain untouched. 
 * 
 * \param[in] Gain           Symmetric matrix of order \c Order. This is usually the The part of the gain matrix
 *                           \f$\mathcal G\f$ that affects the coupling degrees of freedom \f$\mathcal G_c\f$.
 * \param[in] Order          Order of the matrix \f$G_c\f$.
 * \param[in] Substructure   Array of substructures of order \c Order. Used to identify its type and act
 *                           accordingly.
 *
 * \post Each remote facility or controller has received the portion of the gain matrix that they require.
 *
 * \sa Substructure_t
 */
void Substructure_SendGainMatrix( const double *const Gain, const unsigned int Order, const Substructure_t *const Substructure );

/**
 * \brief Performs the sub-stepping process.
 * 
 * This routine arranges the sub-stepping process depending on the type of substructure and will perform a
 * different for each type. The degrees of freedom assigned to a particular substructure may be sent to a
 * remote site (OpenFresco, NSEP, Celestina, ...), to the controller (Experimental) or treated in a separate
 * way if they are numerical substructures.
 *
 * \pre
 * - \c CNodes must be properly initialised through Substructure_ReadCouplingNodes().
 * - The Gain matrix \f$\mathcal G_c\f$ (order \f$n_c x n_c\f$ where \f$n_c\f$ is the number of coupling DOFs)
 *   must have been created through Substructure_MatrixXc() or Substructure_MatrixXc_PS() routines if the
 *   complete gain matrix \f$\mathcal G\f$ (of order \f$n\f$ where \f$n > n_c\f$) is stored in general or
 *   packed storage. It must be a symmetric matrix with all its elements referenced.
 * - \c VecTdT0_c should be generated through the Substructure_VectorXc() routine and be of length \f$n_c\f$.
 * - \c VecTdT, \c CoupForcePrev and \c CoupForce are of length \f$n > n_c\f$.
 *
 * \param[in]  CNodes        Information regarding the coupling nodes.
 * \param[in]  IGain         Symmetric matrix of order \c CNodes.Order. The part of the gain matrix
 *                           \f$\mathcal G\f$ that affects the coupling degrees of freedom \f$\mathcal G_c\f$.
 * \param[in]  VecTdT0_c     Vector with the initial value problem
 * \param[in]  Time          Only used in case of remote substructures using the NSEP protocol and it
 *                           indicates the ellapsed time from since the start of the simulation. It does not
 *                           indicate actual computation time but it is usually computed as \f$\Delta t*i\f$,
 *                           where \f$i\f$ is the step count and \f$\Delta t\f$ the time increment.
 * \param[in]  GAcc          Ground acceleration
 * \param[in]  NSubstep      The number of sub-steps \f$n_{Sub}\f$.
 * \param[in]  DeltaT_Sub    Time increment of a sub-step \f$\Delta t_{Sub}\f$.
 * \param[out] VecTdT        Vector with the updated displacement, velocity or acceleration.
 * \param[out] CoupForcePrev Vector containing the coupling force at the sub-step \f$i - 1\f$, with \f$i =
 *                           n_{Sub}\f$.
 * \param[out] CoupForce     Vector containing the coupling force at the sub-step \f$n_{Sub}\f$.
 *
 * \post The vectors \c VecTdT \c CoupForcePrev and \c CoupForce have the updated displacements at the
 * equations with substructures acting on them.
 *
 * \sa CouplingNote_t.
 */
void Substructure_Substepping( const CouplingNode_t *const CNodes, const double *const IGain,
			       const double *const VecTdT0_c, const double Time, const double GAcc,
			       const unsigned int NSubstep, const double DeltaT_Sub, double *const VecTdT,
			       double *const CoupForcePrev, double *const CoupForce );

#if _MPI_
/**
 * \brief Performs the sub-stepping process. MPI version.
 * 
 * This routine arranges the sub-stepping process depending on the type of substructure and will perform a
 * different for each type. The degrees of freedom assigned to a particular substructure may be sent to a
 * remote site (OpenFresco, NSEP, Celestina, ...), to the controller (Experimental) or treated in a separate
 * way if they are numerical substructures. This is the version that is called when using MPI.
 *
 * \pre
 * - \c CNodes must be properly initialised through Substructure_ReadCouplingNodes().
 * - The Gain matrix \f$\mathcal G_c\f$ (order \f$n_c x n_c\f$ where \f$n_c\f$ is the number of coupling DOFs)
 *   must have been created through Substructure_MatrixXc_MPI() routine. It must be a symmetric matrix with
 *   all its elements referenced.
 * - \c VecTdT0_c should be generated through the Substructure_VectorXc_MPI() routine and be of length
 *   \f$n_c\f$.
 * - \c VecTdT, \c CoupForcePrev and \c CoupForce are of length \f$n > n_c\f$.
 *
 * \param[in]  CNodes        Information regarding the coupling nodes.
 * \param[in]  IGain         Symmetric matrix of order \c CNodes.Order. The part of the gain matrix
 *                           \f$\mathcal G\f$ that affects the coupling degrees of freedom \f$\mathcal G_c\f$.
 * \param[in]  VecTdT0_c     Vector with the initial value problem
 * \param[in]  Time          Only used in case of remote substructures using the NSEP protocol and it
 *                           indicates the ellapsed time from since the start of the simulation. It does not
 *                           indicate actual computation time but it is usually computed as \f$\Delta t*i\f$,
 *                           where \f$i\f$ is the step count and \f$\Delta t\f$ the time increment.
 * \param[in]  GAcc          Ground acceleration
 * \param[in]  NSubstep      The number of sub-steps \f$n_{Sub}\f$.
 * \param[in]  DeltaT_Sub    Time increment of a sub-step \f$\Delta t_{Sub}\f$.
 * \param[out] VecTdT        Vector with the updated displacement, velocity or acceleration.
 * \param[out] CoupForcePrev Vector containing the coupling force at the sub-step \f$i - 1\f$, with \f$i =
 *                           n_{Sub}\f$.
 * \param[out] CoupForce     Vector containing the coupling force at the sub-step \f$n_{Sub}\f$.
 *
 * \post The vectors \c VecTdT \c CoupForcePrev and \c CoupForce have the updated displacements at the
 * equations with substructures acting on them.
 *
 * \sa CouplingNote_t, PMatrixVector_t.
 */
void Substructure_Substepping_MPI( const CouplingNode_t *const CNodes, const double *const IGain,
				   const double *const VecTdT0_c, const double Time, const double GAcc,
				   const unsigned int NSubstep, const double DeltaT_Sub,
				   const MPI_Comm Comm, const int *const LRowIndex_Coupling,
				   const int *const RowProcess_Coupling, const int *const ColProcess_Coupling,
				   PMatrixVector_t *const VecTdT, PMatrixVector_t *const CoupForcePrev,
				   PMatrixVector_t *const CoupForce );
#endif

/**
 * \brief Performs the sub-stepping process for those DOFs that have numerical substructures acting on them.
 * 
 * The sub-stepping process is performed on those substructures that involve numerical simulation. The initial
 * value vector (displacement, velocity or acceleration) is applied in small increments as a ramp function in
 * order to calculate the influence of the substructure in it. The process follows the next equation:
 *
 * \f[\vec x_c^{t + \Delta t} = \vec x_{0c}^t\biggl(1 - \frac{i}{k_{sub}}\biggr) + \vec x_{0c}^{t+\Delta
 * t}\biggl(\frac{i}{k_{sub}}\biggr) + \mathcal G_c\vec f_{s}\f]
 *
 * where:
 * - \f$\vec x_0\f$ is the explicit part of the new state vector;
 * - \f$\mathcal G_c\f$ is the portion of the gain matrix \f$\mathcal G\f$ that affects the coupling degrees
 * of freedom with numerical substructures;
 * - \f$\vec x_c\f$ is the new state vector with the influence of the applied numerical substructures;
 * - \f$k_{sub}\f$ is the number of sub-steps;
 * - \f$i\f$ is the current sub-step \f$1\leq i\leq k_{sub}\f$.
 *
 * \pre 
 * - \c CNodes must be properly initialised through Substructure_ReadCouplingNodes().
 * - The Gain matrix \f$\mathcal G_c\f$ (order \f$n_c x n_c\f$ where \f$n_c\f$ is the number of coupling DOFs)
 *   must have been created through Substructure_MatrixXc(), Substructure_MatrixXc_PS(),
 *   Substructure_MatrixXc_MPI() routines if the complete gain matrix \f$\mathcal G\f$ (of order \f$n\f$ where
 *   \f$n > n_c\f$) is stored in general, packed storage or is a distributed matrix respectively. It must be a
 *   symmetric matrix and only the upper part will be referenced (lower part in FORTRAN routines).
 * - \c VecTdT0_c should be generated through the Substructure_VectorXc() or Substructure_VectorXc_MPI()
 *   routines and be of length \f$n_c\f$. It must contain only those elements of the original vector that are
 *   being affected by a numerical substructure.
 * - \c VecTdT_c, \c CoupForcePrev_c and \c CoupForce_c are of length \f$n_c\f$ and must contain only those
 *   elements of the original vector that are being affected by a numerical substructure.
 * 
 * \param[in]  CNodes          Information regarding the coupling nodes.
 * \param[in]  IGain           Symmetric matrix of order \c CNodes.Order. The part of the gain matrix 
 *                             \f$\mathcal G\f$ that affects the coupling degrees of freedom \f$\mathcal
 *                             G_c\f$.
 * \param[in]  VecTdT0_c       Vector of length \c CNodes.Order. The initial value problem.
 * \param[in]  GAcc            Ground acceleration.
 * \param[in]  NSubstep        The number of sub-steps \f$n_{Sub}\f$.
 * \param[in]  DeltaT_Sub      Time increment of a sub-step \f$\Delta t_{Sub}\f$.
 * \param[out] VecTdT_c        Vector of length \c CNodes.Order with the updated displacement, velocity or
 *                             acceleration.
 * \param[out] CoupForcePrev_c Vector of length \c CNodes.Order containing the coupling force at the sub-step
 *                             \f$i - 1\f$, with \f$i = n_{Sub}\f$.
 * \param[out] CoupForce_c     Vector of length \c CNodes.Order containing the coupling force at the sub-step
 *                             \f$n_{Sub}\f$.
 *
 * \post The vectors \c VecTdT_c \c CoupForcePrev_c and \c CoupForce_c have the updated displacements at the
 * equations with numerical substructures acting on them as a result of the following equation:
 * 
 * \f[\vec x_c^{t + \Delta t} = \vec x_{0c}^t\biggl(1 - \frac{i}{k_{sub}}\biggr) + \vec x_{0c}^{t+\Delta
 * t}\biggl(\frac{i}{k_{sub}}\biggr) + \mathcal G_c\vec f_{s}\f]
 *
 * \sa CouplingNode_t.
 */
void Substructure_Simulate( const CouplingNode_t *const CNodes, const double *const IGain,
			    const double *const VecTdT0_c, const double GAcc, const unsigned int NSubstep,
			    const double DeltaT_Sub, double *const VecTdT_c, double *const CoupForcePrev_c,
			    double *const CoupForce_c );

#endif /* SUBSTRUCTURE_H_ */
