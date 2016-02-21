#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>

#include "Print_Messages.h"

#include "Substructure.h"
#include "Substructure_Exact.h"
#include "Substructure_Newmark.h"
#include "Substructure_Remote.h"
#include "Substructure_Remote_OpenFresco.h"
#include "Substructure_Remote_NSEP.h"
#include "Substructure_UHYDEfbr.h"
#include "Substructure_SimMeasured.h"
#include "Substructure_Experimental.h"
#include "Substructure_CouplingNodes.h"
#include "Substructure_StoneDrums.h"

#include "Definitions.h"

#if _ADWIN_
#include "ADwin_Routines.h"
#endif

#if _MKL_
#include <mkl_blas.h>
#else
#include "Netlib.h"
#endif

void Substructure_SendGainMatrix( const HYSL_FLOAT *const Gain, unsigned int Order, const Substructure_t *const Substructure )
{
     Remote_t *Remote = NULL;
     unsigned int i, j;       /* Counters */
     HYSL_FLOAT *Send, *Recv;     /* Only used in case of of substructures of remote type REMOTE_NSEP */

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
	       Substructure_Remote_Send( Remote->Socket, Order*Order, sizeof(HYSL_FLOAT), (const char *const) Gain );
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
	       Send = (HYSL_FLOAT *) calloc( (size_t) Order, sizeof(HYSL_FLOAT) );
	       Recv = (HYSL_FLOAT *) calloc( (size_t) Order, sizeof(HYSL_FLOAT) );
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

void Substructure_Substepping( const HYSL_FLOAT *const IGain, const HYSL_FLOAT *const VecTdT0_c, const HYSL_FLOAT Time,
			       const HYSL_FLOAT GAcc, const unsigned int NSubstep, const HYSL_FLOAT DeltaT_Sub,
			       CouplingNode_t *const CNodes, HYSL_FLOAT *const VecTdT, HYSL_FLOAT *const CoupForcePrev,
			       HYSL_FLOAT *const CoupForce )
{

     int i;
     bool Called_Sub = false;
     HYSL_FLOAT *Recv = NULL;
     HYSL_FLOAT *Send = NULL;

     Remote_t *Remote;

     Recv = (HYSL_FLOAT *) calloc( (size_t) 3*(size_t)CNodes->Order, sizeof(HYSL_FLOAT) );

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
	  case SIM_NEWMARK:
	       /* This is the same case as SIM_MEASURED. All the simulated substructures are treated together
		* in the same routine.*/
	  case SIM_UHYDE:
	       /* This is the same case as SIM_MEASURED. All the simulated substructures are treated together
		* in the same routine.*/
	  case SIM_STONEDRUMS:
	       /* This is the same case as SIM_MEASURED. All the simulated substructures are treated together
		* in the same routine.*/
	  case SIM_MEASURED:
	       /* Call the Simulate_Substructures() function only once. All the simulated substructures are
		* handled together in this routine */
	       if( !Called_Sub ){
		    Substructure_Simulate( IGain, VecTdT0_c, GAcc, NSubstep, DeltaT_Sub, CNodes, &Recv[0],
					   &Recv[CNodes->Order], &Recv[2*CNodes->Order] );
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
		    Send = (HYSL_FLOAT *) calloc( (size_t) 1+(size_t)CNodes->Order, sizeof(HYSL_FLOAT) );
		    for( i = 0; i < CNodes->Order; i++ ){
			 Send[i] = VecTdT0_c[i];
		    }
		    Send[CNodes->Order] = GAcc;

		    Substructure_Remote_Send( Remote->Socket, (unsigned int) CNodes->Order + 1, sizeof(HYSL_FLOAT), (char *const) Send );

		    Substructure_Remote_Receive( Remote->Socket, 3*(unsigned int) CNodes->Order, sizeof(HYSL_FLOAT), (char *const) Recv );
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
	  VecTdT[CNodes->Array[i] - 1] = Recv[i];
	  CoupForcePrev[i] = Recv[CNodes->Order + i];
	  CoupForce[CNodes->Array[i] - 1] = Recv[2*CNodes->Order + i];
     }

     free( Recv );
}

void Substructure_Simulate( const HYSL_FLOAT *IGain, const HYSL_FLOAT *const VecTdT0_c, const HYSL_FLOAT GAcc, 
			    const unsigned int NSubstep, const HYSL_FLOAT DeltaT_Sub, CouplingNode_t *const CNodes,
			    HYSL_FLOAT *const VecTdT_c, HYSL_FLOAT *const CoupForcePrev_c, HYSL_FLOAT *const CoupForce_c )
{

     unsigned int i, Substep, temp;
     HYSL_FLOAT ramp0, ramp;
     int incx = 1, incy = 1;
     HYSL_FLOAT One;
     int Length;
     char uplo = 'L';
     ExactSim_t *Exact;
     ExactSimESP_t *ExactEsp;
     NewmarkSim_t *Newmark;
     UHYDEfbrSim_t *UHYDE;  
     MeasuredSim_t *Measured;
     StoneDrums_t *StoneDrums, *StoneDrums_Prev;

     Length = (int) CNodes->Order;
     One = 1.0;

     for ( Substep = 1; Substep <= NSubstep; Substep++ ){

	  /* Backup data so that CoupForcePrev_c contains always the last coupling force */
	  hysl_copy( &Length, CoupForce_c, &incx, CoupForcePrev_c, &incy );
	       
	  ramp = (HYSL_FLOAT) Substep / (HYSL_FLOAT) NSubstep;

#if _FLOAT_
	  ramp0 = 1.0f - ramp;   
#else
	  ramp0 = 1.0 - ramp;   
#endif

	  if ( CNodes->Order > 1 ){
	       hysl_copy( &Length, CNodes->VecTdT0_c0, &incx, VecTdT_c, &incy );
	       hysl_scal( &Length, &ramp0, VecTdT_c, &incx );
	       hysl_axpy( &Length, &ramp, VecTdT0_c, &incx, VecTdT_c, &incy );
	       hysl_symv( &uplo, &Length, &One, IGain, &Length, CoupForce_c, &incx, &One, VecTdT_c, &incy ); 
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
		    Substructure_ExactSolutionESP_SDOF( VecTdT_c[i], ramp, GAcc, DeltaT_Sub, ExactEsp, &CoupForce_c[i] );
		    break;
	       case SIM_NEWMARK:
		    Newmark = (NewmarkSim_t *) CNodes->Sub[i].SimStruct;
		    Substructure_Newmark_SDOF( VecTdT_c[i], ramp, GAcc, Newmark, &CoupForce_c[i] );
		    break;
	       case SIM_UHYDE:
		    UHYDE = (UHYDEfbrSim_t *) CNodes->Sub[i].SimStruct;
		    Substructure_SimUHYDE_1D( VecTdT_c[i], DeltaT_Sub, UHYDE, &CoupForce_c[i] );
		    break;
	       case SIM_MEASURED:
		    Measured = (MeasuredSim_t *) CNodes->Sub[i].SimStruct;
		    Substructure_SimMeasured( Measured, &CoupForce_c[i] );
		    break;
	       case SIM_STONEDRUMS:
		    StoneDrums = (StoneDrums_t *) CNodes->Sub[i].SimStruct;
		    
		    if (StoneDrums->PrevDOF >= 0){
			 // Substructure_StoneDrums( VecTdT_c[i] - VecTdT_c[Substructure_FindPosition(
			 // StoneDrums->PrevDOF, CNodes ) - 1], ramp, StoneDrums, &CoupForce_c[i] );
			 Substructure_StoneDrums( VecTdT_c[i], ramp, StoneDrums, &CoupForce_c[i] );
		    } else {
			 Substructure_StoneDrums( VecTdT_c[i], ramp, StoneDrums, &CoupForce_c[i] );
		    }
		    printf("%d %lE %lE %lE %lE\n", Substep, StoneDrums->AccTdT, StoneDrums->VelTdT, VecTdT_c[i], CoupForce_c[i]);
		    break;
	       }
	  }

	  for( i = 0; i < (unsigned int) CNodes->Order; i ++ ){
	       if ( CNodes->Sub[i].Type == SIM_STONEDRUMS ) {
		    if (StoneDrums->PrevDOF >= 0){
			 temp = Substructure_FindPosition( StoneDrums->PrevDOF, CNodes ) - 1;
			 StoneDrums = (StoneDrums_t *) CNodes->Sub[i].SimStruct;
			 StoneDrums_Prev = (StoneDrums_t *) CNodes->Sub[temp].SimStruct;
			 CoupForce_c[temp] = StoneDrums_Prev->Result_Fc - CoupForce_c[i];
		    }
	       }
	  }
     }

       
     /* Backup VecTdT0_c */
     hysl_copy( &Length, VecTdT0_c, &incx, CNodes->VecTdT0_c0, &incy );
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
	       ExactEsp->Acc0 = ExactEsp->AccT; ExactEsp->AccT = ExactEsp->AccTdT;
	       ExactEsp->Vel0 = ExactEsp->VelT; ExactEsp->VelT = ExactEsp->VelTdT;
	       ExactEsp->Disp0 = ExactEsp->DispT; ExactEsp->DispT = VecTdT_c[0];
	       break;
	  case SIM_NEWMARK:
	       Newmark->Acc0 = Newmark->AccT; Newmark->AccT = Newmark->AccTdT;
	       Newmark->Vel0 = Newmark->VelT; Newmark->VelT = Newmark->VelTdT;
	       Newmark->Disp0 = Newmark->DispT; Newmark->DispT = VecTdT_c[0];
	       break;
	  case SIM_UHYDE:
	       break;
	  case SIM_MEASURED:
	       break;
	  case SIM_STONEDRUMS:
	       StoneDrums->Acc0 = StoneDrums->AccT; StoneDrums->AccT = StoneDrums->AccTdT;
	       StoneDrums->Vel0 = StoneDrums->VelT; StoneDrums->VelT = StoneDrums->VelTdT;
	       StoneDrums->Disp0 = StoneDrums->DispT; StoneDrums->DispT = VecTdT_c[0];
	       break;
	  }
     }
}
