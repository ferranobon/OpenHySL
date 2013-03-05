#ifndef SUBSTRUCTURE_COUPLINGNODES_H_
#define SUBSTRUCTURE_COUPLINGNODES_H_

#include "Substructure.h"

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
     double *VecTdT0_c0;
     Substructure_t *Sub;
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
void Substructure_SimMeasuredValues( const char *FileName, const double *const IGain, const double *const VecTdT0_c, double *const VecTdT_c, double *const fcprev, double *const fc, const unsigned int OrderC, const unsigned int NSub );


#endif /* SUBSTRUCTURE_COUPLINGNODES_H_ */

