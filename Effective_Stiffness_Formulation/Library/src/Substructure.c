#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>

#include "Substructure.h"
#include "Substructure_Exact.h"
#include "Substructure_UHYDEfbr.h"
#include "Substructure_SimMeasured.h"
#include "Substructure_Experimental.h"
#include "Substructure_CouplingNodes.h"

#if _ADWIN_
#include "ADwinRoutines.h"
#endif

#if _MKL_
#include <mkl_blas.h>
#else
#include "Netlib.h"
#endif

void Substructure_Substepping( double *const IGain, double *const DispTdT0_c, const double Time, const unsigned int NSubstep,
			       const double DeltaT_Sub, CouplingNode_t *const CNodes, double *const DispTdT,
			       double *const fcprevsub, double *const fc )
{

     int i;
     bool Called_Sub = false;
     double *Recv = NULL;

     Recv = (double *) calloc( (size_t) 3*(size_t)CNodes->Order, sizeof(double) );

     for( i = 0; i < CNodes->Order; i++ ){
	  switch ( CNodes->Sub[i].Type ){
	  case SIM_EXACT:
	       /* This is the same case as SIM_MEASURED. All the simulated substructures are treated together in
		* the same routine.*/
	  case SIM_UHYDE:
	       /* This is the same case as SIM_MEASURED. All the simulated substructures are treated together in
		* the same routine.*/
	  case SIM_MEASURED:
	       /* Call the Simulate_Substructures() function only once. All the simulated substructures are handled
		* together in this routine */
	       if( !Called_Sub ){
		    Substructure_Simulate( CNodes, IGain, DispTdT0_c, NSubstep, DeltaT_Sub, &Recv[0], &Recv[CNodes->Order], &Recv[2*CNodes->Order] );
		    Called_Sub = true;
	       }
	       break;
#if _ADWIN_
	  case EXP_ADWIN:
	       /* Tell ADwin to perform the substepping process */
	       ADWIN_Substep( DispTdT0_c, &Recv[0], &Recv[1], &Recv[2], CNodes->Order );
	       break;
#endif
	  case REMOTE_TCP:
	       /* Using TCP communication protocol */
//	       Send_Data( DispTdT0_c, CNodes->Order, Socket );

//	       Receive_Data( Recv, 3*CNodes->Order, Socket );
	       break;
	  case REMOTE_UDP:
	       /* Using UDP communication protocol */

//	       Send_Data( DispTdT0_c, CNodes->Order, Socket );
//	       if ( recv( Socket, Recv, sizeof(double)*3*(size_t) CNodes->Order,0) != (int) sizeof(double)*3*CNodes->Order ){    /* sizeof returns an unsigned integer ? */
//		    PrintErrorAndExit( "recv() failed in connected UDP mode" );
//	       }
	       break;
	  case REMOTE_NSEP:
	       /* Using NSEP Protocol */
//	       Communicate_With_PNSE( 1, Time, DispTdT0_c, Recv, CNodes->Order );
	       /* Receive the force from the PNSE server. WhatToDo = 2 */
//	       Communicate_With_PNSE( 2, Time, DispTdT0_c, Recv, 3*CNodes->Order );
	       break;
	  case REMOTE_OF:
	       /* Using OpenFresco */
//	       Communicate_With_OpenFresco( DispTdT0_c, Recv, CNodes->Order, 3 ); 
	       break;
	  }
     }

#pragma omp parallel for
     for ( i = 0; i < CNodes->Order; i++ ){
#if _MPI_

	  DispTdT[i] = Recv[i];
	  fcprevsub[i] = Recv[CNodes->Order + i];
	  fc[i] = Recv[2*CNodes->Order + i];
#else
	  DispTdT[CNodes->Array[i] - 1] = Recv[i];
	  fcprevsub[i] = Recv[CNodes->Order + i];
	  fc[CNodes->Array[i] - 1] = Recv[2*CNodes->Order + i];
#endif
     }

     free( Recv );
}

void Substructure_Simulate( CouplingNode_t *const CNodes, double *IGain, double *const VecTdT0_c, const unsigned int NSubstep, const double DeltaT_Sub, double *const VecTdT_c, double *const fcprev, double *const fc )
{

     unsigned int i, Substep;
     double ramp0, ramp;
     int incx = 1, incy = 1;
     double One;
     int Length;
     char uplo = 'L';
     ExactSim_t *Exact;
     UHYDEfbrSim_t *UHYDE;  
     MeasuredSim_t *Measured;

     Length = (int) CNodes->Order;
     One = 1.0;

     for ( Substep = 1; Substep <= NSubstep; Substep++ ){

	  /* Backup data so that fcprev contains always the last coupling force */
	  dcopy( &Length, fc, &incx, fcprev, &incy );
	       
	  ramp = (double) Substep / (double) NSubstep;

	  ramp0 = 1.0 - ramp;   

	  if ( CNodes->Order > 1 ){
	       dcopy( &Length, CNodes->VecTdT0_c0, &incx, VecTdT_c, &incy );
	       dscal( &Length, &ramp0, VecTdT_c, &incx );
	       daxpy( &Length, &ramp, VecTdT0_c, &incx, VecTdT_c, &incy );
	       dsymv( &uplo, &Length, &One, IGain, &Length, fc, &incx, &One, VecTdT_c, &incy ); 
	  } else {
	       VecTdT_c[0] = ramp0*CNodes->VecTdT0_c0[0] + ramp*VecTdT0_c[0] + IGain[0]*fc[0];
	  }
	  
	  /* Compute the new fc */
	  for( i = 0; i < (unsigned int) CNodes->Order; i ++ ){
	       switch( CNodes->Sub[i].Type ){
	       case SIM_EXACT:
		    Exact = (ExactSim_t *) CNodes->Sub[i].SimStruct;
		    Substructure_ExactSolution_SDOF( VecTdT0_c[i], DeltaT_Sub, Exact, &fc[i] );
		    break;
	       case SIM_UHYDE:
		    UHYDE = (UHYDEfbrSim_t *) CNodes->Sub[i].SimStruct;
		    Substructure_SimUHYDE_1D( VecTdT0_c[i], DeltaT_Sub, UHYDE, &fc[i] );
		    break;
	       case SIM_MEASURED:
		    Measured = (MeasuredSim_t *) CNodes->Sub[i].SimStruct;
		    Substructure_SimMeasured( Measured, &fc[i] );
		    break;
	       }
	  }
	  
     }

     /* Backup VecTdT0_c */
     dcopy( &Length, VecTdT0_c, &incx, CNodes->VecTdT0_c0, &incy );
}
