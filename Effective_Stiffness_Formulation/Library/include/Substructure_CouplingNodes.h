/**
 * \file Substructure_CouplingNodes.h
 * \author Ferran Obón Santacana
 * \version 1.0 
 * \date 21st of June 2013
 *
 * \brief Routines to deal with coupling nodes.
 *
 * \todo Strip Substructure_ReadCouplingNodes() from depending on AlgConst_t.
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

#define NUM_TYPE_SUB  7  /*!< \brief Number of recognized sub-structure types */

/**
 * \brief Supported sub-structure types.
 */
enum Substructure_Id { SIM_EXACT_MDOF, /*!< \brief Simulate the substructure using the exact solution. MDOF. */
		       SIM_EXACT_SDOF, /*!< \brief Simulate the substructure using the exact solution. SDOF. */
		       SIM_EXACT_ESP,  /*!< \brief Simulate the sub-structure using an exact integration
					* method \cite Exact */
		       SIM_UHYDE,      /*!< \brief Simulate the substructure using the UHYDE-\em fbr device
					* (1-Dimension).*/
		       SIM_MEASURED,   /*!< \brief Simulate the substructure using measured values. */
		       EXP_ADWIN,      /*!< \brief Run using ADwin. */
		       REMOTE,         /*!< \brief Remote substructure. */
};

/**
 * \brief Structure to accomodate different types of substructures.
 */
typedef struct Substructure{
     void *SimStruct;   /*!< \brief Pointer to a recognised substructure type */
     int Type;          /*!< \brief Type of Substructure. \sa Substructure_Id */
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
 * \brief Deallocates the memory in the \c CouplingNode_t by calling the appropiate routines.
 * 
 * \pre \c CNodes must be properly initialised through Substructure_ReadCouplingNodes().
 *
 * \param[in,out] CNodes On entry only the number of substructures is referenced.
 * 
 * \post The dynamically allocated memory inside \c CNodes is freed and <tt>CNodes.Order = 0</tt>.
 *
 * \sa CouplingNode_t, Substructure_ExactSolutionMDOF_Destroy(), Substructure_ExactSolutionSDOF_Destroy(),
 * Substructure_ExactSolutionESP_Destroy(), Substructure_SimUHYDE_1D_Destroy(),
 * Substructure_SimMeasured_Init(), SubstruSubstructure_Experimental_Destroy(), Substructure_Remote_Destroy().
 */
void Substructure_DeleteCouplingNodes( CouplingNode_t *CNodes );

/**
 * \brief Identifies the type of substructure.
 *
 * \param[in]  Type         Substructure type to be validated/identified.
 * \param[out] Identity_Num Substructure identifier.
 *
 * \post If the input substructure type in \c Type is not recognised. The program will terminate.
 *
 * \sa Substructure_Id.
 */
void Substructure_Identify( char *const Type, int *const Identity_Num );

/**
 * \brief Gets a description from a substructure declaration.
 * 
 * \pre
 * - \c InFile must point to an opened file.
 * - \c LineNum should be greater than 0 and should be 0-based index.
 * - \c Description must point to an allocated memory of length \f$l \geq MAX\_DESCRIPTION\f$.
 * 
 * \param[in]  InFile      File identifier.
 * \param[in]  LineNum     Line number where the description should be read from. Only used for error handling
 *                         purposes in the input file.
 * \param[out] Description Description of the substructure.
 *
 * \post
 * - \c Description contains a string not longer than MAX_DESCRIPTION.
 *
 * \sa MAX_DESCRIPTION.
 */
void Substructure_GetDescription( FILE *const InFile, const int LineNum, char *const Description );

/**
 * \brief Reads the declaration of substructures and the corresponding coupling nodes from a file.
 * 
 * \warning In the case of the MPI version, this routine should only be called by a single MPI process
 * (usually process 0) and then broadcasted to the rest through the Substructure_BroadCastCouplingNodes()
 * routine.

 * The specified file in \c InitCnt is parsed in order to obtain the information regarding the
 * substructures. The routine initialises each substructure depending on its type (allocating necessary
 * memory, opening sockets in the case of remote substructures etc) through calling the appropiate routines.
 *
 * Each type of substructure conditions the way a line in the file has to be formated and the first line in
 * the specified file must be equal to the number of substructures to initialise (or number of lines to
 * read). Each line should terminate with a ";".
 *
 * - <tt>SIM_EXACT_MDOF: Sim_Exact_MDOF, \<ndof\> \<dof1\> ... \<dofn\>, 3 \<num_dof\> \<mass_file\>
 *   \<stiffness_file\> \<Rayleigh_alpha\> \<Rayleigh_beta\>,Optional Description;</tt>
 * - <tt>SIM_EXACT_SDOF: Sim_Exact_SDOF, \<ndof\> \<dof1\> ... \<dofn\>, 3 \<mass\> \<damp\>
 *   \<stiffness\>,Optional Description;</tt>
 * - <tt>SIM_EXACT_ESP: Sim_Exact_ESP, \<ndof\> \<dof1\> ... \<dofn\>, 3 \<mass\> \<damp\>
 *   \<stiffness\>,Optional Description;</tt>
 * - <tt>SIM_UHYDE: Sim_UHYDEfbr, \<ndof\> \<dof1\> ... \<dofn\>, 3 \<yield_displacement\> \<yield_factor\>
 *   \<friction_force\>,Optional Description;</tt>
 * - <tt>SIM_MEASURED: Sim_Measured, \<ndof\> \<dof1\> ... \<dofn\>, \<Input_File\>,Optional Description;</tt>
 * - <tt>EXP_ADWIN: Exp_ADwin, \<ndof\> \<dof1\> ... \<dofn\>,Optional Description;</tt>
 * - <tt>REMOTE: Remote, \<ndof\> \<dof1\> ... \<dofn\>, \<TCP|UDP|NSEP|OpenFresco|Celestina\> \<IP\>
 *   \<Port\>,Optional Description;</tt>
 *
 * \pre 
 * - \c InitCnt must be properly initialised through Algorithm_Init().
 *
 * \param[in]  InitCnt A structure that comprises of several constants.
 * \param[out] CNodes  Structure containing information regarding the coupling nodes sorted in ascending
 *                     order.
 *
 * \post \f$\forall i \in [1..N]~A(i-1) \leq A(i)\f$ where:
 * - \f$N\f$ is the number of coupling nodes.
 * - \f$A\f$ is the array containing the coupling nodes in <tt>struct CNodes (CNodes.Array)</tt>.
 * - The memory should be deallocated through Substructure_DeleteCouplingNodes().
 *
 * \sa AlgConst_t, Coupling_Node_t, Substructure_ExactSolutionMDOF_Init(),
 * Substructure_ExactSolutionSDOF_Init(), Substructure_ExactSolutionESP_Init(),
 * Substructure_SimUHYDE_1D_Init(), Substructure_SimMeasured_Init(), Substructure_Experimental_Init(),
 * Substructure_Remote_Init().
 */
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
 *
 * \sa CouplingNode_t.
 */
void Substructure_SortCouplingNodes( CouplingNode_t *const CNodes );

/**
 * \brief Broadcasts the CouplingNode_t struct to the rest of the MPI processes.
 * 
 * \warning This routine requires MPI.
 * 
 * \pre \c CNodes must be properly initialised through Substructure_ReadCouplingNodes() in MPI process 0.
 *
 * \sa CouplingNode_t.
 */
void Substructure_BroadCastCouplingNodes( CouplingNode_t *const CNodes );

#endif /* SUBSTRUCTURE_COUPLINGNODES_H_ */

