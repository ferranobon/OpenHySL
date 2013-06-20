#ifndef SUBSTRUCTURE_H_
#define SUBSTRUCTURE_H_

#include "Substructure_CouplingNodes.h"

#define NUM_TYPE_SUB  7  /*!< Number of recognized sub-structure types */

enum Substructure_Id { SIM_EXACT_MDOF, /*!< Simulate the substructure using the exact solution. MDOF. */
		       SIM_EXACT_SDOF, /*!< Simulate the substructure using the exact solution. SDOF. */
		       SIM_EXACT_ESP,
		       SIM_UHYDE,      /*!< Simulate the substructure using the UHYDE-fbr device. */
		       SIM_MEASURED,   /*!< Simulate the substructure using measured values. */
		       EXP_ADWIN,      /*!< Run using ADwin */
		       REMOTE,         /*!< Remote node using TCP/IP connection. */
};

void Substructure_SendGainMatrix( double *Gain, unsigned int Order, const Substructure_t *const Substructure );

void Substructure_Substepping( double *const IGain, double *const DispTdT0_c, const double Time, const double GAcc, const unsigned int NSubstep,
			       const double DeltaT_Sub, CouplingNode_t *const CNodes, double *const DispTdT,
			       double *const fcprevsub, double *const fc );

void Substructure_Simulate( CouplingNode_t *const CNodes, double *IGain, double *const VecTdT0_c, const double GAcc, const unsigned int NSubstep, const double DeltaT_Sub, double *const VecTdT_c, double *const fcprev, double *const fc );

#endif /* SUBSTRUCTURE_H_ */
