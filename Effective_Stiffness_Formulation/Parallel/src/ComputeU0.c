/*
 * ComputeU0.c
 *
 *  Created on: 22/07/2011
 *      Author: ferran
 */

#include <mpi.h>

#include "ComputeU0.h"
#include "Netlib.h"
#include "PMatrixVector.h"
#include "Initiation.h"

void EffK_Calc_Effective_Force( PMatrixVector *const Mass, PMatrixVector *const Damp,
				PMatrixVector *const Disp, PMatrixVector *const Vel,
				PMatrixVector *const Acc, PMatrixVector *const Tempvec,
				const float a0, const float a1, const float a2,
				const float a3, const float a4, const float a5,
				PMatrixVector *const Eff_Force )
{
     static int ione = 1;
     static int incx = 1, incy = 1;
     static char uplo = 'L';
     static float Alpha, Beta;

     /* PBLAS: tempvec = Disp */
     pscopy_( &Tempvec->GlobalSize.Row, Disp->Array, &ione, &ione, Disp->Desc, &ione, Tempvec->Array, &ione, &ione, Tempvec->Desc, &ione );
     /* PBLAS: tempvec = a0*Disp = a0*tempvec */
     Alpha = a0;
     psscal_( &Tempvec->GlobalSize.Row, &Alpha, Tempvec->Array, &ione, &ione, Tempvec->Desc, &incx );
     /* BLAS: tempvec = a0*Disp + a2*Vel = tempvec + a2*Vel */
     Alpha = a2;
     psaxpy_( &Tempvec->GlobalSize.Row, &Alpha, Vel->Array, &ione, &ione, Vel->Desc, &incx, Tempvec->Array, &ione, &ione, Tempvec->Desc, &incy );
     /* PBLAS: tempvec = a0*Disp + a2*Vel + a3*Acc = tempvec + a3*Acc */
     Alpha = a3;
     psaxpy_( &Tempvec->GlobalSize.Row, &Alpha, Acc->Array, &ione, &ione, Acc->Desc, &incx, Tempvec->Array, &ione, &ione, Tempvec->Desc, &incy );
     /* PBLAS: Eff_Force = Mass*(a0*Disp + a2*Vel + a3*Acc) = Mass*tempvec */
     Alpha = 1.0; Beta = 0.0;
     pssymv_( &uplo, &Tempvec->GlobalSize.Row, &Alpha, Mass->Array, &ione, &ione, Mass->Desc, Tempvec->Array, &ione, &ione, Tempvec->Desc, &incx, &Beta, Eff_Force->Array, &ione, &ione, Eff_Force->Desc, &incy );

     /* PBLAS: tempvec = Disp */
     pscopy_( &Tempvec->GlobalSize.Row, Disp->Array, &ione, &ione, Disp->Desc, &ione, Tempvec->Array, &ione, &ione, Tempvec->Desc, &ione );
     /* PBLAS: tempvec = a1*Disp = a0*tempvec */
     Alpha = a1;
     psscal_( &Tempvec->GlobalSize.Row, &Alpha, Tempvec->Array, &ione, &ione, Tempvec->Desc, &incx );
     /* PBLAS: tempvec = a1*Disp + a4*Vel = tempvec + a4*Vel */
     Alpha = a4;
     psaxpy_( &Tempvec->GlobalSize.Row, &Alpha, Vel->Array, &ione, &ione, Vel->Desc, &incx, Tempvec->Array, &ione, &ione, Tempvec->Desc, &incy );
     /* PBLAS: tempvec = a1*Disp + a4*Vel + a5*Acc = tempvec + a5*Acc */
     Alpha = a5;
     psaxpy_( &Tempvec->GlobalSize.Row, &Alpha, Acc->Array, &ione, &ione, Acc->Desc, &incx, Tempvec->Array, &ione, &ione, Tempvec->Desc, &incy );
     /* PBLAS: Eff_Force = Mass*(a0*Disp + a2*Vel + a3*Acc) + Damp*(a1*Disp + a4*Vel + a5*Acc) = Eff_Force + Damp*tempvec */
     Alpha = 1.0; Beta = 1.0;
     pssymv_( &uplo, &Tempvec->GlobalSize.Row, &Alpha, Damp->Array, &ione, &ione, Damp->Desc, Tempvec->Array, &ione, &ione, Tempvec->Desc, &incx, &Beta, Eff_Force->Array, &ione, &ione, Eff_Force->Desc, &incy );
}

void EffK_ComputeU0( PMatrixVector *const Eff_Force, PMatrixVector *const In_Load,
		     PMatrixVector *const Err_Force, const float PID_P, PMatrixVector *const Keinv,
		     PMatrixVector *const Tempvec, PMatrixVector *const Disp0 )
{
     int ione;
     int incx, incy;
     float Alpha, Beta;
     char uplo;

     ione = 1;
     incx = 1; incy = 1;
     Alpha = 1.0; Beta = 0.0;
     uplo = 'L';

     /* PBLAS: tempvec = Eff_Force */
     pscopy_( &Tempvec->GlobalSize.Row, Eff_Force->Array, &ione, &ione, Eff_Force->Desc, &incx,
	      Tempvec->Array, &ione, &ione, Tempvec->Desc, &incy );

     /* PBLAS: tempvec = Eff_Force + LoadTdT = tempvec + LoadTdT */
     psaxpy_( &Tempvec->GlobalSize.Row, &Alpha, In_Load->Array, &ione, &ione, In_Load->Desc, &incx,
	      Tempvec->Array, &ione, &ione, Tempvec->Desc, &incy );

     /* BLAS: tempvec = Eff_Force + LoadTdT + Err_Force = tempvec + Err_Force. The sign of Err_Force was already applied when calculating it. */
     Alpha = PID_P;
     psaxpy_( &Tempvec->GlobalSize.Row, &Alpha, Err_Force->Array, &ione, &ione, Err_Force->Desc, &incx,
	      Tempvec->Array,  &ione, &ione, Tempvec->Desc, &incy );

     /* BLAS: Disp0 = Keinv*(Eff_Force + LoadTdT + Err_Force) = Keinv*Tempvec */
     Alpha = 1.0;
     pssymv_( &uplo, &Tempvec->GlobalSize.Row, &Alpha, Keinv->Array, &ione, &ione, Keinv->Desc, 
	      Tempvec->Array, &ione, &ione, Tempvec->Desc, &incx,
	      &Beta, Disp0->Array, &ione, &ione, Disp0->Desc, &incy );
}

void CreateVectorXm( MPI_Comm Comm, PMatrixVector *const VectorX, PMatrixVector *const VectorXm, const Coupling_Node *const CNodes )
{

	static int incx, incy;
	static int Length;
	static int PosX_Row, PosX_Col, PosXm_Row, PosXm_Col;
	static int icoup;    /* Counter for the coupling nodes */

	incx = 1; incy = 1;
	PosX_Col = 1;
	PosXm_Col = 1;

	PosX_Row = 1;
	PosXm_Row = 1;
	for( icoup = 0; icoup < CNodes->Order; icoup++ ){
	     Length = CNodes->Array[icoup] - PosX_Row;

	     /* Copy the part of the vector between two positions */
	     pscopy_( &Length, VectorX->Array, &PosX_Row, &PosX_Col, VectorX->Desc, &incx,
		      VectorXm->Array, &PosXm_Row, &PosXm_Col, VectorXm->Desc, &incy );

	     /* Update the values of the position in the vectors */
	     PosX_Row = CNodes->Array[icoup] + 1;
	     PosXm_Row = PosXm_Row + Length;
	}

	/* Copy the elements from the last position until the end of the vector */
	Length = VectorX->GlobalSize.Row - CNodes->Array[CNodes->Order-1];
	pscopy_( &Length, VectorX->Array, &PosX_Row, &PosX_Col, VectorX->Desc, &incx,
		 VectorXm->Array, &PosXm_Row, &PosXm_Col, VectorXm->Desc, &incy );
}

void CreateVectorXc( MPI_Comm Comm, PMatrixVector *const VecX, float *VecXc, const Coupling_Node *const CNodes )
{
     int icoup;
     int nprow, npcol, myrow, mycol;
     int rank;
     int GRowIndex, GColIndex;
     int LRowIndex, LColIndex;
     int RowProcess, ColProcess;

     MPI_Status status;

     MPI_Comm_rank( Comm, &rank );

     /* Get grid info */
     Cblacs_gridinfo( (*VecX).Desc[1], &nprow, &npcol, &myrow, &mycol );
     
     GColIndex = 1;

     for( icoup = 0; icoup < CNodes->Order; icoup++ ){

	  GRowIndex = CNodes->Array[icoup];   /* 1 based */
	  /* Given the global index of an element (GRowIndex, GColIndex) returns the local index of the element
	     (LRowIndex, LColIndex) and the coordinates of the process (Row Process, ColProcess) */
	  infog2l_( &GRowIndex, &GColIndex, (*VecX).Desc, &nprow, &npcol, &myrow, &mycol, &LRowIndex, &LColIndex,
		    &RowProcess, &ColProcess );

	  if ( Cblacs_pnum( (*VecX).Desc[1], RowProcess, ColProcess ) == 0 && myrow == RowProcess && mycol == ColProcess ){ 
	       VecXc[icoup] = (*VecX).Array[(LRowIndex - 1)*(*VecX).LocalSize.Col + (LColIndex - 1)];
		    
	  } else {
	       if ( myrow == RowProcess && mycol == ColProcess ){
			 
		    MPI_Send( &(*VecX).Array[(LRowIndex - 1)*(*VecX).LocalSize.Col + (LColIndex - 1)], 1, MPI_FLOAT, 0, 1, Comm );
	       }
	       
	       if ( rank == 0 ){
		    /* Store the diagonal elements */
		    MPI_Recv( &VecXc[icoup], 1, MPI_FLOAT, Cblacs_pnum( (*VecX).Desc[1], RowProcess, ColProcess ), 1, Comm, &status );
	       }
	  }	  
     }     
}
