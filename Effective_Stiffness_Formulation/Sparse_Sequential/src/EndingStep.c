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

void JoinNonCouplingPart( MatrixVector *const VecXm, const MatrixVector *const Keinv_m, const MatrixVector *const fcprevsub,
			  MatrixVector *const Vec, const int PosCouple, const int OrderC )
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
	lda = Max( 1, Keinv_m->Rows - OrderC );

	sgemv_( &trans, &Rows, &Cols, &Alpha, Keinv_m->Array, &lda,
		  &fcprevsub->Array[PosCouple - 1], &incx, &Beta, VecXm->Array, &incy );

	/* Copy the first elements */
	TempSize = PosCouple - 1;
	scopy_( &TempSize, (*VecXm).Array, &incx, (*Vec).Array, &incy );

	/* Join the part after the coupling position */
	TempSize = (*Vec).Rows - ( PosCouple + OrderC - 1 );
	scopy_( &TempSize, &(*VecXm).Array[PosCouple - 1], &incx, &(*Vec).Array[PosCouple + OrderC - 1], &incy );

}

void Compute_Acceleration( const MatrixVector *const DispTdT, const MatrixVector *const DispT, const MatrixVector *const VelT,
			   const MatrixVector *const AccT, const float a0, const float a2, const float a3,
			   MatrixVector *const AccTdT )
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

void Compute_Velocity( const MatrixVector *const VelT, const MatrixVector *const AccT, const MatrixVector *const AccTdT,
		       const float a6, const float a7, MatrixVector *const VelTdT )
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

void Compute_Force_Error( const MatrixVector *const Mass, const MatrixVector *const Damp, const MatrixVector *Stiff,
			  const MatrixVector *const AccTdT, const MatrixVector *const VelTdT, const MatrixVector *const DispTdT,
			  const MatrixVector *const fc, const MatrixVector *const LoadTdT, MatrixVector *const fu )
{

     static int incx = 1, incy = 1;
     static float Alpha, Beta;
     static char uplo = 'L';

     /* BLAS: fu = Mass*AccTdT */
     Alpha = 1.0; Beta = 0.0;
     ssymv_( &uplo, &fu->Rows, &Alpha, Mass->Array, &fu->Rows, AccTdT->Array, &incx, &Beta, fu->Array, &incy );
     /* BLAS: fu = Mass*AccTdT + Damp*VelTdT = fu + Damp*VelTdT */
     Beta = 1.0;
     ssymv_( &uplo, &fu->Rows, &Alpha, Damp->Array, &fu->Rows, VelTdT->Array, &incx, &Beta, fu->Array, &incy );
     /* BLAS: fu = Mass*AccTdT + Damp*VelTdT + Stiff*DispTdT = fu + Stiff*DispTdT */
     ssymv_( &uplo, &fu->Rows, &Alpha, Stiff->Array, &fu->Rows, DispTdT->Array, &incx, &Beta, fu->Array, &incy );
     /* BLAS: fu = -(Mass*AccTdT + Damp*VelTdT + Stiff*DispTdT) = -fu */
     Alpha = -1.0;
     sscal_( &fu->Rows, &Alpha, fu->Array, &incx );
     /* BLAS: fu = fc -(Mass*AccTdT + Damp*VelTdT + Stiff*DispTdT) = fc -fu */
     Alpha = 1.0;
     saxpy_( &fu->Rows, &Alpha, fc->Array, &incx, fu->Array, &incy );
     /* BLAS: fu = LoatTdT + fc -(Mass*AccTdT + Damp*VelTdT + Stiff*DispTdT) = LoadTdT -fu */
     saxpy_( &fu->Rows, &Alpha, LoadTdT->Array, &incx, fu->Array, &incy );

}
