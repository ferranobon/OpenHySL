#ifndef SUBSTRUCTURE_H_
#define SUBSTRUCTURE_H_

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

void Substructure_Substepping( double *const IGain, double *const DispTdT0_c, const double Time, const unsigned int NSubstep,
			       const double DeltaT_Sub, CouplingNode_t *const CNodes, double *const DispTdT,
			       double *const fcprevsub, double *const fc );

void Substructure_Simulate( CouplingNode_t *const CNodes, double *IGain, double *const VecTdT0_c, const unsigned int NSubstep, const double DeltaT_Sub, double *const VecTdT_c, double *const fcprev, double *const fc );

#endif /* SUBSTRUCTURE_H_ */
