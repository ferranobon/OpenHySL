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
#include "MatrixVector.h"
#include "Netlib.h"

#include <mkl_spblas.h>

void JoinNonCouplingPart_Dense( Dense_MatrixVector *const VecXm, const Dense_MatrixVector *const Keinv_m, const Dense_MatrixVector *const fcprevsub,
			  Dense_MatrixVector *const Vec, const int PosCouple, const int OrderC )
{

	static int incx, incy;
	static float Alpha, Beta;
	static char trans;
	static int Rows, Cols, TempSize;
	static int lda;

	incx = 1; incy = 1;
	trans = 'N';
	Alpha = 1.0; Beta = 0.0;
	Rows = Keinv_m->Rows;
	Cols = Keinv_m->Cols;
	lda = Max( 1, Keinv_m->Rows);

	sgemv_( &trans, &Rows, &Cols, &Alpha, Keinv_m->Array, &lda,
		&fcprevsub->Array[PosCouple - 1], &incx, &Beta, VecXm->Array, &incy );

	/* Copy the first elements */
	TempSize = PosCouple - 1;
	scopy_( &TempSize, (*VecXm).Array, &incx, (*Vec).Array, &incy );

	/* Join the part after the coupling position */
	TempSize = (*Vec).Rows - ( PosCouple + OrderC - 1 );
	scopy_( &TempSize, &(*VecXm).Array[PosCouple - 1], &incx, &(*Vec).Array[PosCouple + OrderC - 1], &incy );

}

void JoinNonCouplingPart_Sparse( Dense_MatrixVector *const VecXm, const Sp_MatrixVector *const Keinv_m, const Dense_MatrixVector *const fcprevsub,
			  Dense_MatrixVector *const Vec, const int PosCouple, const int OrderC )
{

	static int incx, incy;
	static float Alpha, Beta;
	static int Rows,Cols, TempSize;
	static char trans = 'N';
	static char matdescra[6] = {'G', 'U', 'N', 'C'};

	Alpha = 1.0; Beta = 0.0;
	incx = 1; incy = 1;
	Rows = Keinv_m->Rows;
	Cols = Keinv_m->Cols;

	mkl_scsrmv( &trans, &Rows, &Cols, &Alpha, matdescra, Keinv_m->Values, Keinv_m->Columns, Keinv_m->RowIndex, &Keinv_m->RowIndex[1], &fcprevsub->Array[PosCouple -1],&Beta, VecXm->Array );

	/* Copy the first elements */
	TempSize = PosCouple - 1;
	scopy_( &TempSize, (*VecXm).Array, &incx, (*Vec).Array, &incy );

	/* Join the part after the coupling position */
	TempSize = (*Vec).Rows - ( PosCouple + OrderC - 1 );
	scopy_( &TempSize, &(*VecXm).Array[PosCouple - 1], &incx, &(*Vec).Array[PosCouple + OrderC - 1], &incy );

}

void Compute_Acceleration( const Dense_MatrixVector *const DispTdT, const Dense_MatrixVector *const DispT, const Dense_MatrixVector *const VelT,
			   const Dense_MatrixVector *const AccT, const float a0, const float a2, const float a3,
			   Dense_MatrixVector *const AccTdT )
{
     static int incx = 1, incy = 1;
     static float Alpha;

     /* BLAS: AccTdT = DispTdT */
     scopy_( &AccTdT->Rows, DispTdT->Array, &incx, AccTdT->Array, &incy ); 
     /* BLAS: AccTdT = DispTdT - DispT = AccTdT - DispT */
     Alpha = -1.0;
     saxpy_( &AccTdT->Rows, &Alpha, DispT->Array, &incx, AccTdT->Array, &incy );
     /* BLAS: AccTdT = a0*(DispTdT - DispT) = a0*AccTdT */
     Alpha = a0;
     sscal_( &AccTdT->Rows, &Alpha, AccTdT->Array, &incx );
     /* BLAS: AccTdT = a0*(DispTdT - DispT) - a2*VelT = AccTdT - a2*VelT */
     Alpha = -a2;
     saxpy_( &AccTdT->Rows, &Alpha, VelT->Array, &incx, AccTdT->Array, &incy );
     /* BLAS: AccTdT = a0*(DispTdT - DispT) - a2*VelT - a3*AccT = AccTdT - a3*AccT */
     Alpha = -a3;
     saxpy_( &AccTdT->Rows, &Alpha, AccT->Array, &incx, AccTdT->Array, &incy );
}

void Compute_Velocity( const Dense_MatrixVector *const VelT, const Dense_MatrixVector *const AccT, const Dense_MatrixVector *const AccTdT,
		       const float a6, const float a7, Dense_MatrixVector *const VelTdT )
{
     static int incx = 1, incy= 1;
     static float Alpha;

     /* BLAS: VelTdT = VelT */
     scopy_( &VelTdT->Rows, VelT->Array, &incx, VelTdT->Array, &incy );
     /* BLAS: VelTdT = VelT + a6*AccT = VelTdT + a6*AccT */
     Alpha = a6;
     saxpy_( &VelTdT->Rows, &Alpha, AccT->Array, &incx, VelTdT->Array, &incy );
     /* BLAS: VelTdT = VelT + a6*AccT + a7*AccTdT = VelTdT + a7*AccTdT */
     Alpha = a7;
     saxpy_( &VelTdT->Rows, &Alpha, AccTdT->Array, &incx, VelTdT->Array, &incy );
}

void Compute_Force_Error( const Sp_MatrixVector *const Mass, const Sp_MatrixVector *const Damp, const Sp_MatrixVector *Stiff,
			  const Dense_MatrixVector *const AccTdT, const Dense_MatrixVector *const VelTdT, const Dense_MatrixVector *const DispTdT,
			  const Dense_MatrixVector *const fc, const Dense_MatrixVector *const LoadTdT, Dense_MatrixVector *const fu, Dense_MatrixVector *const Tempvec )
{

     static int incx = 1, incy = 1;
     static float Alpha, Beta;
     static char uplo = 'L';
     static char trans = 'N';
     static char matdescra[6] = {'S', 'U', 'N', 'C'};



     /* BLAS: fu = Mass*AccTdT */
     Alpha = 1.0; Beta = 0.0;
     mkl_scsrmv( &trans, &fu->Rows, &fu->Rows, &Alpha, matdescra, Mass->Values, Mass->Columns, Mass->RowIndex, &Mass->RowIndex[1], AccTdT->Array, &Beta, fu->Array );
     /* BLAS: fu = Mass*AccTdT + Damp*VelTdT = fu + Damp*VelTdT */
     Beta = 1.0;
     mkl_scsrmv( &trans, &fu->Rows, &fu->Rows, &Alpha, matdescra, Damp->Values, Damp->Columns, Damp->RowIndex, &Damp->RowIndex[1], VelTdT->Array, &Beta, fu->Array );
     /* BLAS: fu = Mass*AccTdT + Damp*VelTdT + Stiff*DispTdT = fu + Stiff*DispTdT */
     mkl_scsrmv( &trans, &fu->Rows, &fu->Rows, &Alpha, matdescra, Stiff->Values, Stiff->Columns, Stiff->RowIndex, &Stiff->RowIndex[1], DispTdT->Array, &Beta, fu->Array );
     /* BLAS: fu = -(Mass*AccTdT + Damp*VelTdT + Stiff*DispTdT) = -fu */
     Alpha = -1.0;
     sscal_( &fu->Rows, &Alpha, fu->Array, &incx );
     /* BLAS: fu = fc -(Mass*AccTdT + Damp*VelTdT + Stiff*DispTdT) = fc -fu */
     Alpha = 1.0;
     saxpy_( &fu->Rows, &Alpha, fc->Array, &incx, fu->Array, &incy );
     /* BLAS: fu = LoatTdT + fc -(Mass*AccTdT + Damp*VelTdT + Stiff*DispTdT) = LoadTdT -fu */
     saxpy_( &fu->Rows, &Alpha, LoadTdT->Array, &incx, fu->Array, &incy );
}
