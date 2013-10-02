#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>

#include "Substructure.h"
#include "Substructure_Exact.h"
#include "Substructure_Remote.h"
#include "Substructure_Remote_OpenFresco.h"
#include "Substructure_Remote_NSEP.h"
#include "Substructure_UHYDEfbr.h"
#include "Substructure_SimMeasured.h"
#include "Substructure_Experimental.h"
#include "Substructure_CouplingNodes.h"

#if _ADWIN_
#include "ADwin_Routines.h"
#endif

void Substructure_Substepping_MPI( const CouplingNode_t *const CNodes, const double *const IGain,
				   const double *const VecTdT0_c, const double Time, const double GAcc,
				   const unsigned int NSubstep, const double DeltaT_Sub, double *const VecTdT,
				   double *const CoupForcePrev, double *const CoupForce )
{

     int i;
     bool Called_Sub = false;
     double *Recv = NULL;
     double *Send = NULL;

     Remote_t *Remote;

     Recv = (double *) calloc( (size_t) 3*(size_t)CNodes->Order, sizeof(double) );

     /* Copy the older coupling force. This is necessary for simulations */
     for( i = 0; i < CNodes->Order; i++ ){
	  Recv[2*CNodes->Order + i] = CoupForce[CNodes->Array[i] - 1];
     }

     for( i = 0; i < CNodes->Order; i++ ){
	  switch ( CNodes->Sub[i].Type ){
	  case SIM_EXACT_MDOF:
	       /* This is the same case as SIM_MEASURED. All the simulated substructures are treated together
		* in the same routine.*/
	  case SIM_EXACT_SDOF:
	       /* This is the same case as SIM_MEASURED. All the simulated substructures are treated together
		* in the same routine.*/
	  case SIM_EXACT_ESP:
	  case SIM_UHYDE:
	       /* This is the same case as SIM_MEASURED. All the simulated substructures are treated together
		* in the same routine.*/
	  case SIM_MEASURED:
	       /* Call the Simulate_Substructures() function only once. All the simulated substructures are
		* handled together in this routine */
	       if( !Called_Sub ){
		    Substructure_Simulate( CNodes, IGain, VecTdT0_c, GAcc, NSubstep, DeltaT_Sub, &Recv[0], &Recv[CNodes->Order], &Recv[2*CNodes->Order] );
		    Called_Sub = true;
	       }
	       break;
#if _ADWIN_
	  case EXP_ADWIN:
	       /* Tell ADwin to perform the substepping process */
	       ADwin_Substep( VecTdT0_c, (unsigned int) CNodes->Order, 0.75, &Recv[0], &Recv[1], &Recv[2] );
	       break;
#endif
	  case REMOTE:
	       Remote = (Remote_t *) CNodes->Sub[i].SimStruct;

	       if( Remote->Type == REMOTE_TCP || Remote->Type == REMOTE_UDP || Remote->Type == REMOTE_CELESTINA ){
		    Send = (double *) calloc( (size_t) 1+(size_t)CNodes->Order, sizeof(double) );
		    for( i = 0; i < CNodes->Order; i++ ){
			 Send[i] = VecTdT0_c[i];
		    }
		    Send[CNodes->Order] = GAcc;

		    Substructure_Remote_Send( Remote->Socket, (unsigned int) CNodes->Order + 1, sizeof(double), (char *const) Send );

		    Substructure_Remote_Receive( Remote->Socket, 3*(unsigned int) CNodes->Order, sizeof(double), (char *const) Recv );
		    free( Send );
	       } else if( Remote->Type == REMOTE_NSEP ){
		    /* Using NSEP Protocol */
		    Substructure_Remote_NSEP( Remote, NSEP_SEND_CMD, Time, Remote->NSub, Send, Recv );
		    /* Receive the force from the PNSE server. WhatToDo = 2 */
		    Substructure_Remote_NSEP( Remote, NSEP_REQUEST_CSIG, Time, 3*Remote->NSub, Send, Recv );
	       } else if( Remote->Type == REMOTE_OF ){
		    /* Using OpenFresco */
		    /* Send the trial response */
		    Substructure_Remote_OpenFresco( Remote->Socket, OF_REMOTE_SET_TRIAL_RESPONSE, Remote->NSub, Send, NULL );
		    /* Get the DAQ response */
		    Substructure_Remote_OpenFresco( Remote->Socket, OF_REMOTE_GET_DAQ_RESPONSE, Remote->NSub, NULL, Recv );
	       } else assert( Remote->Type >= 0 || Remote->Type < NUM_REMOTE_TYPE );
	       break;
	  }
     }

#pragma omp parallel for
     for ( i = 0; i < CNodes->Order; i++ ){
	  VecTdT[i] = Recv[i];
	  CoupForcePrev[i] = Recv[CNodes->Order + i];
	  CoupForce[i] = Recv[2*CNodes->Order + i];
     }

     free( Recv );
}
