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

#define NUM_TYPE_SUB  7  /*!< Number of recognized sub-structure types */

/**
 * \brief Supported sub-structure types.
 */
enum Substructure_Id { SIM_EXACT_MDOF, /*!< Simulate the substructure using the exact solution. MDOF. */
		       SIM_EXACT_SDOF, /*!< Simulate the substructure using the exact solution. SDOF. */
		       SIM_EXACT_ESP,  /*!< Simulate the sub-structure using an exact integration method \cite */
		       SIM_UHYDE,      /*!< Simulate the substructure using the UHYDE-fbr device. */
		       SIM_MEASURED,   /*!< Simulate the substructure using measured values. */
		       EXP_ADWIN,      /*!< Run using ADwin */
		       REMOTE,         /*!< Remote substructure. */
};


void Substructure_SendGainMatrix( double *Gain, unsigned int Order, const Substructure_t *const Substructure );

void Substructure_Substepping( double *const IGain, double *const DispTdT0_c, const double Time, const double GAcc, const unsigned int NSubstep,
			       const double DeltaT_Sub, CouplingNode_t *const CNodes, double *const DispTdT,
			       double *const fcprevsub, double *const fc );

void Substructure_Simulate( CouplingNode_t *const CNodes, double *IGain, double *const VecTdT0_c, const double GAcc, const unsigned int NSubstep, const double DeltaT_Sub, double *const VecTdT_c, double *const fcprev, double *const fc );

#endif /* SUBSTRUCTURE_H_ */
