/**
 * \file Substructure_CouplingNodes.h
 * \author Ferran Obón Santacana
 * \version 1.0 
 * \date 21st of June 2013
 *
 * \brief Routines to deal with coupling nodes.
 *
 * These routines deal with the creation and destriction of experimental substructures. They do not depend on
 * specific hardware (eg. ADwin).
 */
#ifndef SUBSTRUCTURE_COUPLINGNODES_H_
#define SUBSTRUCTURE_COUPLINGNODES_H_

#include <stdio.h>
#include "Algorithm_Aux.h"
#include "Substructure.h"

#define MAX_LINE 200          /*!< \brief Maximum input line length. */
#define MAX_SUBTYPE 20        /*!< \brief Maximum length of the substructure type */
#define MAX_DESCRIPTION 80    /*!< \brief Maximum description length for a substructure. */
#define MAX_FILENAME 20       /*!< \brief Maximum file name length. */

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
     int *Array;             /*!< \brief Array containing the coupling nodes of length \c Order.*/
     int Order;              /*!< \brief Number of coupling nodes. */
     double *VecTdT0_c0;     /*!< \brief Previous values of the explicit displacement. They are used in the
			      * sub-stepping process (ramp function). */
     Substructure_t *Sub;    /*!< \brief Array of substructures of length \c Order.*/
} CouplingNode_t;


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
 * \param[in]  Type         Substructure type to be validated/identified.
 * \param[out] Identity_Num Substructure identifier.
 *
 * \post If the input substructure type in \c Type is not recognised. The program will terminate.
 */
void Substructure_Identify( char *const Type, int *const Identity_Num );

void Substructure_GetDescription( FILE *const InFile, const int LineNum, char *const Description );

void Substructure_ReadCouplingNodes( const AlgConst_t *const InitCnt, CouplingNode_t *const CNodes );

/**
 * \brief Sort the coupling nodes in ascending order.
 *
 * \param[in,out] CNodes Structure containing information regarding the coupling nodes. On output, the
 *                       coupling nodes are sorted in ascending order.
 *
 * \post \f$\forall i \in [1..N]~A(i-1) \leq A(i)\f$ where:
 * - \f$N\f$ is the number of coupling nodes.
 * - \f$A\f$ is the array containing the coupling nodes in <tt>struct CNodes (CNodes.Array)</tt>.
 */
 void Substructure_SortCouplingNodes( CouplingNode_t *const CNodes );

/**
 * \brief Perform */
void Substructure_SimMeasuredValues( const char *FileName, const double *const IGain, const double *const VecTdT0_c, double *const VecTdT_c, double *const fcprev, double *const fc, const unsigned int OrderC, const unsigned int NSub );


void Substructure_BroadCastCouplingNodes( CouplingNode_t *const CNodes );

#endif /* SUBSTRUCTURE_COUPLINGNODES_H_ */

