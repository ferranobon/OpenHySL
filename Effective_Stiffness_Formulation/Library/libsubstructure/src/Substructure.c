#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>

#include "Print_Messages.h"

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

#if _MKL_
#include <mkl_blas.h>
#else
#include "Netlib.h"
#endif

void Substructure_SendGainMatrix( const double *const Gain, unsigned int Order, const Substructure_t *const Substructure )
{
     Remote_t *Remote = NULL;
     unsigned int i, j;       /* Counters */
     double *Send, *Recv;     /* Only used in case of of substructures of remote type REMOTE_NSEP */

     if( Substructure->Type == EXP_ADWIN ){
#if _ADWIN_

	  /* Send matrix Gc to ADwin. */
	  /* In the code of ADWIN, the variable G is stored in DATA_1 */
	  ADwin_SendArray( 1, Gain, Order*Order );
	  
	  Print_Header( SUCCESS );
	  printf("Gain Matrix successfully sent to ADwin system.\n" );
#else
	  Print_Header( ERROR );
	  fprintf( stderr, "The support for ADwin was disabled at runtime.\n");
	  exit( EXIT_FAILURE );
#endif
     } else if ( Substructure->Type == REMOTE ){
	  Remote = (Remote_t *) Substructure->SimStruct;

	  if( Remote->Type == REMOTE_TCP || Remote->Type == REMOTE_UDP || Remote->Type == REMOTE_CELESTINA ){
	       Substructure_Remote_Send( Remote->Socket, Order*Order, sizeof(double), (const char *const) Gain );
	       Print_Header( SUCCESS );
	       printf("Gain Matrix successfully sent to Remote site %s:%s (%s protocol).\n", Remote->IP, Remote->Port,
		      Substructure_Remote_Type[Remote->Type] );
	  } else if( Remote->Type == REMOTE_NSEP ){
	  /* Send the matrix Gc to the PNSE server in order to reach the FCM */
	  /*
	   * Note that the CGM can only send a NSEP_CMD message with the size of order, therefore it is needed
	   * to send the matrix Gc per rows to fullfill this requisite. This also means, that the FCM must
	   * send as many NSEP_CSIG packets as the number of rows to keep everything synchronised
	   */
	       Send = (double *) calloc( (size_t) Order, sizeof(double) );
	       Recv = (double *) calloc( (size_t) Order, sizeof(double) );
	       for ( i = 0; i < Order; i++ ){
		    for ( j = 0; j < Order; j++ ){
			 Send[j] = Gain[i*Order + j];
		    }
		    /* Send the matrix Keinv_c to PNSE Server */
		    Substructure_Remote_NSEP( Remote, NSEP_SEND_CMD, 0.0, Remote->NSub, Send, Recv );
		    /* This is done so that PNSE do not overtake the first step */
		    Substructure_Remote_NSEP( Remote, NSEP_REQUEST_CSIG, 0.0, Remote->NSub, Send, Recv );
	       }
	       
	       free( Send );
	       free( Recv );
	  } else if( Remote->Type == REMOTE_OF ){
	       /* TODO Implement Send the Matrix G in OpenFresco. Wait for the answer from Andreas. What
		* follows is an ugly hack. */
	       Substructure_Remote_OpenFresco( Remote->Socket, OF_REMOTE_SET_TRIAL_RESPONSE,
					       Remote->NSub*Remote->NSub, Gain, NULL );
	  } else assert( Remote->Type >= 0 || Remote->Type < NUM_REMOTE_TYPE );
     } else assert( Substructure->Type == EXP_ADWIN || Substructure->Type == REMOTE );
}

void Substructure_Substepping( const CouplingNode_t *const CNodes, const double *const IGain,
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
#if _MPI_
	  VecTdT[i] = Recv[i];
	  CoupForcePrev[i] = Recv[CNodes->Order + i];
	  CoupForce[i] = Recv[2*CNodes->Order + i];
#else
	  VecTdT[CNodes->Array[i] - 1] = Recv[i];
	  CoupForcePrev[i] = Recv[CNodes->Order + i];
	  CoupForce[CNodes->Array[i] - 1] = Recv[2*CNodes->Order + i];
#endif
     }

     free( Recv );
}

void Substructure_Simulate( const CouplingNode_t *const CNodes, const double *IGain, const double *const VecTdT0_c, const double GAcc, 
			    const unsigned int NSubstep, const double DeltaT_Sub, double *const VecTdT_c,
			    double *const CoupForcePrev_c, double *const CoupForce_c )
{

     unsigned int i, Substep;
     double ramp0, ramp;
     int incx = 1, incy = 1;
     double One;
     int Length;
     char uplo = 'L';
     ExactSim_t *Exact;
     ExactSimESP_t *ExactEsp;
     UHYDEfbrSim_t *UHYDE;  
     MeasuredSim_t *Measured;

     Length = (int) CNodes->Order;
     One = 1.0;

     for ( Substep = 1; Substep <= NSubstep; Substep++ ){

	  /* Backup data so that CoupForcePrev_c contains always the last coupling force */
	  dcopy( &Length, CoupForce_c, &incx, CoupForcePrev_c, &incy );
	       
	  ramp = (double) Substep / (double) NSubstep;

	  ramp0 = 1.0 - ramp;   

	  if ( CNodes->Order > 1 ){
	       dcopy( &Length, CNodes->VecTdT0_c0, &incx, VecTdT_c, &incy );
	       dscal( &Length, &ramp0, VecTdT_c, &incx );
	       daxpy( &Length, &ramp, VecTdT0_c, &incx, VecTdT_c, &incy );
	       dsymv( &uplo, &Length, &One, IGain, &Length, CoupForce_c, &incx, &One, VecTdT_c, &incy ); 
	  } else {
	       VecTdT_c[0] = ramp0*CNodes->VecTdT0_c0[0] + ramp*VecTdT0_c[0] + IGain[0]*CoupForce_c[0];
	  }
	  
	  /* Compute the new CoupForce_c */
	  for( i = 0; i < (unsigned int) CNodes->Order; i ++ ){
	       switch( CNodes->Sub[i].Type ){
	       case SIM_EXACT_MDOF:
		    Exact = (ExactSim_t *) CNodes->Sub[i].SimStruct;
		    Substructure_ExactSolutionMDOF( VecTdT_c[i], ramp, GAcc, DeltaT_Sub, Exact, &CoupForce_c[i] );
		    break;
	       case SIM_EXACT_SDOF:
		    Exact = (ExactSim_t *) CNodes->Sub[i].SimStruct;
		    Substructure_ExactSolutionSDOF( VecTdT_c[i], ramp, GAcc, DeltaT_Sub, Exact, &CoupForce_c[i] );
		    break;
	       case SIM_EXACT_ESP:
		    ExactEsp = (ExactSimESP_t *) CNodes->Sub[i].SimStruct;
		    Substructure_ExactSolutionESP_SDOF( VecTdT_c[i], DeltaT_Sub, ExactEsp, &CoupForce_c[i] );
		    break;
	       case SIM_UHYDE:
		    UHYDE = (UHYDEfbrSim_t *) CNodes->Sub[i].SimStruct;
		    Substructure_SimUHYDE_1D( VecTdT_c[i], DeltaT_Sub, UHYDE, &CoupForce_c[i] );
		    break;
	       case SIM_MEASURED:
		    Measured = (MeasuredSim_t *) CNodes->Sub[i].SimStruct;
		    Substructure_SimMeasured( Measured, &CoupForce_c[i] );
		    break;
	       }
	  }	  
     }

     /* Backup VecTdT0_c */
     dcopy( &Length, VecTdT0_c, &incx, CNodes->VecTdT0_c0, &incy );
     for( i = 0; i < (unsigned int) CNodes->Order; i ++ ){
	  switch( CNodes->Sub[i].Type ){
	  case SIM_EXACT_MDOF:
	       /* Same as SIM_EXACT_SDOF */
	  case SIM_EXACT_SDOF:
	       Exact->Acc0 = Exact->AccT; Exact->AccT = Exact->AccTdT;
	       Exact->Vel0 = Exact->VelT; Exact->VelT = Exact->VelTdT;
	       Exact->Disp0 = Exact->DispT; Exact->DispT = VecTdT_c[0];
	       break;
	  case SIM_EXACT_ESP:
	       break;
	  case SIM_UHYDE:
	       break;
	  case SIM_MEASURED:
	       break;
	  }
     }
}
