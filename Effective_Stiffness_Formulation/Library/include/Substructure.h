#ifndef SUBSTRUCTURE_H_
#define SUBSTRUCTURE_H_

#include "MatrixVector.h"

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

//extern const char *Substructure_Type[]; /*!< Has all the substructure types. The order of the substructures
//					 * must be the same as in Substructure_Id type.
//					 * \sa Substructure_Id.*/

#define MAX_LINE 200          /*!< Maximum input line length. */
#define MAX_SUBTYPE 20        /*!< Maximum length of the substructure type */
#define MAX_DESCRIPTION 80    /*!< Maximum description length for a substructure. */
#define MAX_FILENAME 20       /*!< Maximum file name length. */

typedef struct Substructure{
     void *SimStruct;
     int Type;
} Substructure_t;

/**
 * \brief Structure to store the coupling nodes.
 * 
 * This structure is used in order to store the coupling nodes that will be used
 * during a test. The nodes are stored sequentially and in increasing order.
 */
typedef struct CouplingNode {
     int *Array;  /*!< \brief Array containing the coupling nodes */
     int Order;   /*!< \brief Number of coupling nodes */
     double *u0c0;
     Substructure_t *Sub;
} CouplingNode_t;

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

/**
 * \brief Deallocates the memory in the \c CouplingNode_t.
 * 
 * \pre \c CNodes must be properly initialised through Substructure_ReadCouplingNodes().
 *
 * \param[in,out] CNodes On entry only the number of substructures is referenced.
 * 
 * \post The dynamically allocated memory inside \c CNodes is freed and <tt>CNodes.Order = 0</tt>.
 */
void Substructure_DeleteCouplingNodes( CouplingNode_t *CNodes );

/**
 * \brief Identifies the type of substructure.
 *
 * \param[in] Type Substructure type to be validated/identified.
 * \param[out] Identity_Num Substructure identifier.
 *
 * \post If the input substructure type in \c Type is not recognised. The program will terminate.
 */
void Substructure_Identify( char *const Type, int *const Identity_Num );

void Substructure_ReadCouplingNodes( CouplingNode_t *const CNodes, const unsigned int NSteps, const unsigned int NSubsteps,
				     const int OrderSub, const double DeltaTSub, const char *Filename );

/**
 * \brief Perform */
void Substructure_SimMeasuredValues( const char *FileName, const double *const Keinv, const double *const u0c, double *const uc, double *const fcprev, double *const fc, const unsigned int OrderC, const unsigned int NSub );
void Substructure_Simulate( CouplingNode_t *const CNodes, double *Keinv, double *const u0c, double *const uc, double *const fcprev, double *const fc, const unsigned int NSubstep, const double DeltaT_Sub );

#endif /* SUBSTRUCTURE_H_ */
