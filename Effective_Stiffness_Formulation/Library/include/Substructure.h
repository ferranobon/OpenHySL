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
 */
#ifndef SUBSTRUCTURE_H_
#define SUBSTRUCTURE_H_

#include "Substructure_CouplingNodes.h"

#define NUM_TYPE_SUB  7  /*!< \brief Number of recognized sub-structure types */

/**
 * \brief Supported sub-structure types.
 */
enum Substructure_Id { SIM_EXACT_MDOF, /*!< \brief Simulate the substructure using the exact solution. MDOF. */
		       SIM_EXACT_SDOF, /*!< \brief Simulate the substructure using the exact solution. SDOF. */
		       SIM_EXACT_ESP,  /*!< \brief Simulate the sub-structure using an exact integration method \cite{Exact} */
		       SIM_UHYDE,      /*!< \brief Simulate the substructure using the UHYDE-fbr device. */
		       SIM_MEASURED,   /*!< \brief Simulate the substructure using measured values. */
		       EXP_ADWIN,      /*!< \brief Run using ADwin */
		       REMOTE,         /*!< \brief Remote substructure. */
};

void Substructure_SendGainMatrix( const double *const Gain, const unsigned int Order, const Substructure_t *const Substructure );

void Substructure_Substepping( const double *const IGain, const double *const DispTdT0_c, const double Time, const double GAcc, const unsigned int NSubstep,
			       const double DeltaT_Sub, const CouplingNode_t *const CNodes, double *const DispTdT,
			       double *const fcprevsub, double *const fc );

void Substructure_Simulate( const CouplingNode_t *const CNodes, const double *const IGain, const double *const VecTdT0_c, const double GAcc, const unsigned int NSubstep, const double DeltaT_Sub, double *const VecTdT_c, double *const fcprev, double *const fc );

#endif /* SUBSTRUCTURE_H_ */
