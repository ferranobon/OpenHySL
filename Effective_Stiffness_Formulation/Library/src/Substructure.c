#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>

#include "Substructure.h"
#include "Substructure_Exact.h"
#include "Substructure_UHYDEfbr.h"
#include "Substructure_SimMeasured.h"
#include "Substructure_Experimental.h"

#include "Auxiliary_Math.h"

#if _ADWIN_
#include "ADwinRoutines.h"
#endif

#if _MKL_
#include <mkl_blas.h>
#else
#include "Netlib.h"
#endif

void Join_NonCouplingPart( MatrixVector_t *const VecTdT_m, const MatrixVector_t *const Gain_m,
			   const MatrixVector_t *const fcprevsub, const CouplingNode_t *const CNodes,
			   MatrixVector_t *const VecTdT )			  
{
     static int icoup;                 /* Counter for the coupling nodes */
     static int incx, incy;            /* Stride in the vectors */
     static double Alpha, Beta;        /* Constants for the BLAS routines */
     static char trans;                /* Use or not the transpose */
     static int Rows, Cols;            /* Number of Rows and columns */
     static int lda;                   /* Leading dimension */
     static int Length, PosX, PosXm;   /* Length and position counters */
     
     incx = 1; incy = 1;
     trans = 'N';
     Alpha = 1.0; Beta = 1.0;
     Rows = Gain_m->Rows;
     Cols = Gain_m->Cols;
     lda = Max( 1, Gain_m->Rows);

     /* Update the VecTdT_m displacments to include the effects of the coupling force */
     /* BLAS: VecTdT_m = Gain_m*fcprevsub */
     dgemv( &trans, &Rows, &Cols, &Alpha, Gain_m->Array, &lda,
	     fcprevsub->Array, &incx, &Beta, VecTdT_m->Array, &incy );

     /* Copy the updated values into the complete displacement vector */
     PosX = 0; PosXm = 0;
     for ( icoup = 0; icoup < CNodes->Order; icoup++ ){
	  Length = CNodes->Array[icoup] - PosX -1;
	  dcopy( &Length, &VecTdT_m->Array[PosXm], &incx, &VecTdT->Array[PosX], &incy );
	  PosX = CNodes->Array[icoup];
	  PosXm = PosXm + Length;
     }

     /* Add the elements between the final coupling node and the final element
      * of the complete displacement vector */
     Length = VecTdT->Rows - CNodes->Array[CNodes->Order -1];
     dcopy( &Length, &VecTdT_m->Array[PosXm], &incx, &VecTdT->Array[PosX], &incy );	
}

void Substructure_Substepping( double *const IGain, double *const DispTdT0_c, const double Time, const int NSubstep,
			       const double DeltaT_Sub, CouplingNode_t *const CNodes, double *const DispTdT,
			       double *const fcprevsub, double *const fc )
{

     unsigned int i;
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
	  for( i = 0; i < CNodes->Order; i ++ ){
	       switch( CNodes->Sub[i].Type )
	       if( CNodes->Sub[i].Type == SIM_EXACT ){
		    Exact = (ExactSim_t *) CNodes->Sub[i].SimStruct;
		    Substructure_ExactSolution_SDOF( VecTdT0_c[i], DeltaT_Sub, Exact, &fc[i] );
	       } else if ( CNodes->Sub[i].Type == SIM_UHYDE ){
		    UHYDE = (UHYDEfbrSim_t *) CNodes->Sub[i].SimStruct;
		    Substructure_SimUHYDE_1D( VecTdT0_c[i], DeltaT_Sub, UHYDE, &fc[i] );
	       } else if ( CNodes->Sub[i].Type == SIM_MEASURED ){
		    Measured = (MeasuredSim_t *) CNodes->Sub[i].SimStruct;
		    Substructure_SimMeasured( Measured, &fc[i] );
	       }
	  }
	  
     }

     /* Backup VecTdT0_c */
     dcopy( &Length, VecTdT0_c, &incx, CNodes->VecTdT0_c0, &incy );
}
