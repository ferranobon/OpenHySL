/**
 * \file EndingStep.c
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 10th of August 2011
 *
 * \brief Source code of the functions used during the Compute U0 phase.
 *
 * This file contains the source code of the functions that are called during the Ending Step phase of the substructure
 * algorithm. For the moment, only the joining of the non-coupling and the coupling part routine has been considered.
 *
 * \todo Evaluate the performance impact (overhead) of writing routines to calculate the acceleration, velocity and displacement in this
 * phase as a separate routines.
 */

#include "EndingStep.h"
#include "PMatrixVector.h"
#include "Initiation.h"
#include "Netlib.h"

void JoinNonCouplingPart( PMatrixVector *const VecXm, PMatrixVector *const Keinv_m, PMatrixVector *const fcprevsub,
			  PMatrixVector *const Vec, const Coupling_Node *const CNodes )
{

     static int incx, incy, ione;
     int icoup, Length;
     float Alpha, Beta;
     static char trans = 'N';
     static int Rows, Cols;
     static int PosX_Row, PosX_Col, PosXm_Row, PosXm_Col;

     incx = 1; incy = 1;
     ione = 1;
     Rows = Keinv_m->GlobalSize.Row;
     Cols = Keinv_m->GlobalSize.Col;
     Alpha = 1.0f; Beta = 1.0f;
     trans = 'N';

     psgemv_( &trans, &Rows, &Cols, &Alpha, Keinv_m->Array, &ione, &ione, Keinv_m->Desc, fcprevsub->Array, &ione, &ione, fcprevsub->Desc, &incx, &Beta, VecXm->Array, &ione, &ione, VecXm->Desc, &incy );

     PosX_Col = 1;
     PosXm_Col = 1;
     
     PosX_Row = 1;
     PosXm_Row = 1;
     for( icoup = 0; icoup < CNodes->Order; icoup++ ){
	  Length = CNodes->Array[icoup] - PosX_Row;
	  
	  /* Copy the part of the vector between two positions */
	  pscopy_( &Length, VecXm->Array, &PosXm_Row, &PosXm_Col, VecXm->Desc, &incx,
		   Vec->Array, &PosX_Row, &PosX_Col, Vec->Desc, &incy );
	  
	  /* Update the values of the position in the vectors */
	  PosX_Row = CNodes->Array[icoup] + 1; /* 1 based index */
	  PosXm_Row = PosXm_Row + Length;
	}

	/* Copy the elements from the last position until the end of the vector */
	Length = Vec->GlobalSize.Row - CNodes->Array[CNodes->Order-1];
	pscopy_( &Length, VecXm->Array, &PosXm_Row, &PosXm_Col, VecXm->Desc, &incx,
		 Vec->Array, &PosX_Row, &PosX_Col, Vec->Desc, &incy );

}

void Compute_Acceleration( PMatrixVector *const DispTdT, PMatrixVector *const DispT, PMatrixVector *const VelT,
			   PMatrixVector *const AccT, const float a0, const float a2, const float a3,
			   PMatrixVector *const AccTdT )
{

     static int ione = 1;
     static int incx = 1, incy = 1;
     static float Alpha;

     /* PBLAS: AccTdT = DispTdT */
     pscopy_( &AccTdT->GlobalSize.Row, DispTdT->Array, &ione, &ione, DispTdT->Desc, &ione, AccTdT->Array, &ione, &ione, AccTdT->Desc, &ione );
     /* PBLAS: AccTdT = DispTdT - DispT = AccTdT - DispT */
     Alpha = -1.0;
     psaxpy_( &AccTdT->GlobalSize.Row, &Alpha, DispT->Array, &ione, &ione, DispT->Desc, &incx, AccTdT->Array, &ione, &ione, AccTdT->Desc, &incy );
     /* PBLAS: AccTdT = a0*(DispTdT - DispT) = a0*AccTdT */
     Alpha = a0;
     psscal_( &AccTdT->GlobalSize.Row, &Alpha, AccTdT->Array, &ione, &ione, AccTdT->Desc, &incx );
     /* PBLAS: AccTdT = a0*(DispTdT - DispT) - a2*VelT = AccTdT - a2*VelT */
     Alpha = -a2;
     psaxpy_( &AccTdT->GlobalSize.Row, &Alpha, VelT->Array, &ione, &ione, VelT->Desc, &incx, AccTdT->Array, &ione, &ione, AccTdT->Desc, &incy );
     /* OBLAS: AccTdT = a0*(DispTdT - DispT) - a2*VelT - a3*AccT = AccTdT - a3*AccT */
     Alpha = -a3;
     psaxpy_( &AccTdT->GlobalSize.Row, &Alpha, AccT->Array, &ione, &ione, AccT->Desc, &incx, AccTdT->Array, &ione, &ione, AccTdT->Desc, &incy );
}

void Compute_Velocity( PMatrixVector *const VelT, PMatrixVector *const AccT, PMatrixVector *const AccTdT,
		       const float a6, const float a7, PMatrixVector *const VelTdT )
{
     static int ione = 1;
     static int incx = 1, incy= 1;
     static float Alpha;

     /* BLAS: VelTdT = VelT */
     pscopy_( &VelTdT->GlobalSize.Row, VelT->Array, &ione, &ione, VelT->Desc, &incx, VelTdT->Array, &ione, &ione, VelTdT->Desc, &incy );
     /* BLAS: VelTdT = VelT + a6*AccT = VelTdT + a6*AccT */
     Alpha = a6;
     psaxpy_( &VelTdT->GlobalSize.Row, &Alpha, AccT->Array, &ione, &ione, AccT->Desc, &incx, VelTdT->Array, &ione, &ione, VelTdT->Desc, &incy );
     /* BLAS: VelTdT = VelT + a6*AccT + a7*AccTdT = VelTdT + a7*AccTdT */
     Alpha = a7;
     psaxpy_( &VelTdT->GlobalSize.Row, &Alpha, AccTdT->Array, &ione, &ione, AccTdT->Desc, &incx, VelTdT->Array, &ione, &ione, VelTdT->Desc, &incy );
}

void Compute_Force_Error( PMatrixVector *const Mass, PMatrixVector *const Damp, PMatrixVector *Stiff,
			  PMatrixVector *const AccTdT, PMatrixVector *const VelTdT, PMatrixVector *const DispTdT,
			  PMatrixVector *const fc, PMatrixVector *const LoadTdT, PMatrixVector *const fu )
{

     static int ione = 1;
     static int incx = 1, incy = 1;
     static float Alpha, Beta;
     static char uplo = 'L';

     /* PBLAS: fu = Mass*AccTdT */
     Alpha = 1.0; Beta = 0.0;
     pssymv_( &uplo, &fu->GlobalSize.Row, &Alpha, Mass->Array, &ione, &ione, Mass->Desc, AccTdT->Array, &ione, &ione, AccTdT->Desc, &incx, &Beta, fu->Array, &ione, &ione, fu->Desc, &incy );
     /* PBLAS: fu = Mass*AccTdT + Damp*VelTdT = fu + Damp*VelTdT */
     Beta = 1.0;
     pssymv_( &uplo, &fu->GlobalSize.Row, &Alpha, Damp->Array, &ione, &ione, Damp->Desc, VelTdT->Array, &ione, &ione, VelTdT->Desc, &incx, &Beta, fu->Array, &ione, &ione, fu->Desc, &incy );
     /* BLAS: fu = Mass*AccTdT + Damp*VelTdT + Stiff*DispTdT = fu + Stiff*DispTdT */
     pssymv_( &uplo, &fu->GlobalSize.Row, &Alpha, Stiff->Array, &ione, &ione, Stiff->Desc, DispTdT->Array, &ione, &ione, DispTdT->Desc, &incx, &Beta, fu->Array, &ione, &ione, fu->Desc, &incy );
     /* BLAS: fu = -(Mass*AccTdT + Damp*VelTdT + Stiff*DispTdT) = -fu */
     Alpha = -1.0;
     psscal_( &fu->GlobalSize.Row, &Alpha, fu->Array, &ione, &ione, fu->Desc, &incx );
     /* BLAS: fu = fc -(Mass*AccTdT + Damp*VelTdT + Stiff*DispTdT) = fc -fu */
     Alpha = 1.0;
     psaxpy_( &fu->GlobalSize.Row, &Alpha, fc->Array, &ione, &ione, fc->Desc, &incx, fu->Array, &ione, &ione, fu->Desc, &incy );
     /* BLAS: fu = LoatTdT + fc -(Mass*AccTdT + Damp*VelTdT + Stiff*DispTdT) = LoadTdT -fu */
     psaxpy_( &fu->GlobalSize.Row, &Alpha, LoadTdT->Array, &ione, &ione, LoadTdT->Desc, &incx, fu->Array, &ione, &ione, fu->Desc, &incy );

}
