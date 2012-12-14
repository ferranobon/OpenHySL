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

#if _SPARSE_
#include <mkl_spblas.h>
#endif

#include "EndingStep.h"
#include "MatrixVector.h"
#include "Netlib.h"
#include "Initiation.h"

void JoinNonCouplingPart( MatrixVector *const VecXm, const MatrixVector *const Keinv_m, const MatrixVector *const fcprevsub,
			  MatrixVector *const Vec, const Coupling_Node *const CNodes )
{
	static int icoup;    /* Counter for the coupling nodes */
	static int incx, incy;
	static double Alpha, Beta;
	static char trans;
	static int Rows, Cols;
	static int lda;
	static int Length, PosX, PosXm;

	incx = 1; incy = 1;
	trans = 'N';
	Alpha = 1.0; Beta = 1.0;
	Rows = Keinv_m->Rows;
	Cols = Keinv_m->Cols;
	lda = Max( 1, Keinv_m->Rows);

	dgemv_( &trans, &Rows, &Cols, &Alpha, Keinv_m->Array, &lda,
		fcprevsub->Array, &incx, &Beta, VecXm->Array, &incy );

	PosX = 0; PosXm = 0;
	for ( icoup = 0; icoup < CNodes->Order; icoup++ ){
	     Length = CNodes->Array[icoup] - PosX -1;
	     dcopy_( &Length, &VecXm->Array[PosXm], &incx, &Vec->Array[PosX], &incy );
	     PosX = CNodes->Array[icoup];
	     PosXm = PosXm + Length;
	}

	/* Add the elements between the final coupling node and the final element
	 * in the vector */
	Length = Vec->Rows - CNodes->Array[CNodes->Order -1];
	dcopy_( &Length, &VecXm->Array[PosXm], &incx, &Vec->Array[PosX], &incy );	
}

void Compute_Acceleration( const MatrixVector *const DispTdT, const MatrixVector *const DispT, const MatrixVector *const VelT,
			   const MatrixVector *const AccT, const double a0, const double a2, const double a3,
			   MatrixVector *const AccTdT )
{
     static int incx = 1, incy = 1;
     static double Alpha;

     /* BLAS: AccTdT = DispTdT */
     dcopy_( &AccTdT->Rows, DispTdT->Array, &incx, AccTdT->Array, &incy ); 
     /* BLAS: AccTdT = DispTdT - DispT = AccTdT - DispT */
     Alpha = -1.0;
     daxpy_( &AccTdT->Rows, &Alpha, DispT->Array, &incx, AccTdT->Array, &incy );
     /* BLAS: AccTdT = a0*(DispTdT - DispT) = a0*AccTdT */
     Alpha = a0;
     dscal_( &AccTdT->Rows, &Alpha, AccTdT->Array, &incx );
     /* BLAS: AccTdT = a0*(DispTdT - DispT) - a2*VelT = AccTdT - a2*VelT */
     Alpha = -a2;
     daxpy_( &AccTdT->Rows, &Alpha, VelT->Array, &incx, AccTdT->Array, &incy );
     /* BLAS: AccTdT = a0*(DispTdT - DispT) - a2*VelT - a3*AccT = AccTdT - a3*AccT */
     Alpha = -a3;
     daxpy_( &AccTdT->Rows, &Alpha, AccT->Array, &incx, AccTdT->Array, &incy );
}

void Compute_Velocity( const MatrixVector *const VelT, const MatrixVector *const AccT, const MatrixVector *const AccTdT,
		       const double a6, const double a7, MatrixVector *const VelTdT )
{
     static int incx = 1, incy= 1;
     static double Alpha;

     /* BLAS: VelTdT = VelT */
     dcopy_( &VelTdT->Rows, VelT->Array, &incx, VelTdT->Array, &incy );
     /* BLAS: VelTdT = VelT + a6*AccT = VelTdT + a6*AccT */
     Alpha = a6;
     daxpy_( &VelTdT->Rows, &Alpha, AccT->Array, &incx, VelTdT->Array, &incy );
     /* BLAS: VelTdT = VelT + a6*AccT + a7*AccTdT = VelTdT + a7*AccTdT */
     Alpha = a7;
     daxpy_( &VelTdT->Rows, &Alpha, AccTdT->Array, &incx, VelTdT->Array, &incy );
}

void Compute_Force_Error( const MatrixVector *const Mass, const MatrixVector *const Damp, const MatrixVector *Stiff,
			  const MatrixVector *const AccTdT, const MatrixVector *const VelTdT, const MatrixVector *const DispTdT,
			  const MatrixVector *const fc, const MatrixVector *const LoadTdT, MatrixVector *const fu )
{

     static int incx = 1, incy = 1;
     static double Alpha, Beta;
     static char uplo = 'L';

     /* BLAS: fu = Mass*AccTdT */
     Alpha = 1.0; Beta = 0.0;
     dsymv_( &uplo, &fu->Rows, &Alpha, Mass->Array, &fu->Rows, AccTdT->Array, &incx, &Beta, fu->Array, &incy );
     /* BLAS: fu = Mass*AccTdT + Damp*VelTdT = fu + Damp*VelTdT */
     Beta = 1.0;
     dsymv_( &uplo, &fu->Rows, &Alpha, Damp->Array, &fu->Rows, VelTdT->Array, &incx, &Beta, fu->Array, &incy );
     /* BLAS: fu = Mass*AccTdT + Damp*VelTdT + Stiff*DispTdT = fu + Stiff*DispTdT */
     dsymv_( &uplo, &fu->Rows, &Alpha, Stiff->Array, &fu->Rows, DispTdT->Array, &incx, &Beta, fu->Array, &incy );
     /* BLAS: fu = -(Mass*AccTdT + Damp*VelTdT + Stiff*DispTdT) = -fu */
     Alpha = -1.0;
     dscal_( &fu->Rows, &Alpha, fu->Array, &incx );
     /* BLAS: fu = fc -(Mass*AccTdT + Damp*VelTdT + Stiff*DispTdT) = fc -fu */
     Alpha = 1.0;
     daxpy_( &fu->Rows, &Alpha, fc->Array, &incx, fu->Array, &incy );
     /* BLAS: fu = LoatTdT + fc -(Mass*AccTdT + Damp*VelTdT + Stiff*DispTdT) = LoadTdT -fu */
     daxpy_( &fu->Rows, &Alpha, LoadTdT->Array, &incx, fu->Array, &incy );

}

#if _SPARSE_

void JoinNonCouplingPart_Sparse( MatrixVector *const VecXm, const Sp_MatrixVector *const Keinv_m, const MatrixVector *const fcprevsub,
			  MatrixVector *const Vec, const Coupling_Node *const CNodes )
{
	static int icoup;    /* Counter for the coupling nodes */
	static int incx, incy;
	static double Alpha, Beta;
	static char trans;
	static int Rows, Cols;
	static int lda;
	static int Length, PosX, PosXm;
	static char matdescra[6] = {'G', 'U', 'N', 'F'};

	incx = 1; incy = 1;
	trans = 'N';
	Alpha = 1.0; Beta = 1.0;
	Rows = Keinv_m->Rows;
	Cols = Keinv_m->Cols;
	lda = Max( 1, Keinv_m->Rows);

	mkl_dcsrmv( &trans, &Rows, &Cols, &Alpha, matdescra, Keinv_m->Values, Keinv_m->Columns, Keinv_m->RowIndex, &Keinv_m->RowIndex[1], &fcprevsub->Array,&Beta, VecXm->Array );

	PosX = 0; PosXm = 0;
	for ( icoup = 0; icoup < CNodes->Order; icoup++ ){
	     Length = CNodes->Array[icoup] - PosX -1;
	     dcopy_( &Length, &VecXm->Array[PosXm], &incx, &Vec->Array[PosX], &incy );
	     PosX = CNodes->Array[icoup];
	     PosXm = PosXm + Length;
	}

	/* Add the elements between the final coupling node and the final element
	 * in the vector */
	Length = Vec->Rows - CNodes->Array[CNodes->Order -1];
	dcopy_( &Length, &VecXm->Array[PosXm], &incx, &Vec->Array[PosX], &incy );	
}

void Compute_Force_Error_Sparse( const Sp_MatrixVector *const Mass, const Sp_MatrixVector *const Damp, const Sp_MatrixVector *Stiff,
			  const MatrixVector *const AccTdT, const MatrixVector *const VelTdT, const MatrixVector *const DispTdT,
			  const MatrixVector *const fc, const MatrixVector *const LoadTdT, MatrixVector *const fu )
{

     static int incx = 1, incy = 1;
     static double Alpha, Beta;
     static char trans = 'N';
     static char matdescra[6] = {'S', 'U', 'N', 'F'};



     /* BLAS: fu = Mass*AccTdT */
     Alpha = 1.0; Beta = 0.0;
     mkl_dcsrmv( &trans, &fu->Rows, &fu->Rows, &Alpha, matdescra, Mass->Values, Mass->Columns, Mass->RowIndex, &Mass->RowIndex[1], AccTdT->Array, &Beta, fu->Array );
     /* BLAS: fu = Mass*AccTdT + Damp*VelTdT = fu + Damp*VelTdT */
     Beta = 1.0;
     mkl_dcsrmv( &trans, &fu->Rows, &fu->Rows, &Alpha, matdescra, Damp->Values, Damp->Columns, Damp->RowIndex, &Damp->RowIndex[1], VelTdT->Array, &Beta, fu->Array );
     /* BLAS: fu = Mass*AccTdT + Damp*VelTdT + Stiff*DispTdT = fu + Stiff*DispTdT */
     mkl_dcsrmv( &trans, &fu->Rows, &fu->Rows, &Alpha, matdescra, Stiff->Values, Stiff->Columns, Stiff->RowIndex, &Stiff->RowIndex[1], DispTdT->Array, &Beta, fu->Array );
     /* BLAS: fu = -(Mass*AccTdT + Damp*VelTdT + Stiff*DispTdT) = -fu */
     Alpha = -1.0;
     dscal_( &fu->Rows, &Alpha, fu->Array, &incx );
     /* BLAS: fu = fc -(Mass*AccTdT + Damp*VelTdT + Stiff*DispTdT) = fc -fu */
     Alpha = 1.0;
     daxpy_( &fu->Rows, &Alpha, fc->Array, &incx, fu->Array, &incy );
     /* BLAS: fu = LoatTdT + fc -(Mass*AccTdT + Damp*VelTdT + Stiff*DispTdT) = LoadTdT -fu */
     daxpy_( &fu->Rows, &Alpha, LoadTdT->Array, &incx, fu->Array, &incy );
}

#endif
