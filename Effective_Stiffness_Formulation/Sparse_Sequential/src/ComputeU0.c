/**
 * \file ComputeU0.c
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 10th of August 2011
 *
 * \brief Source code of the functions used during the Compute U0 phase.
 *
 * This file contains the source code of the functions that are called during the Compute U0 phase of the substructure
 * algorithm. This includes a function to calculate the input load and functions to decouple vectors into their non-coupling
 * and coupling part.
 *
 */
#include <stdio.h>

#include "ComputeU0.h"
#include "MatrixVector.h"
#include "Netlib.h"

#include <mkl_spblas.h>

void Calculatefi( Dense_MatrixVector *const fi, const Dense_MatrixVector *const fc, const Dense_MatrixVector *const li, const Dense_MatrixVector *const Deltaf )
{

	static int incx, incy;
	static float Alpha;

	incx = 1; incy = 1;
	Alpha = 1.0;

	/* BLAS: fi = Deltaf */
	scopy_( &(*fi).Rows, (*Deltaf).Array, &incx, (*fi).Array, &incy );

	/* BLAS fi = fc + fi = fc + Deltaf */
	saxpy_( &(*fi).Rows, &Alpha, (*fc).Array, &incx, (*fi).Array, &incy );
	/* BLAS fi = li + fi = li + fc + Deltaf */
	saxpy_( &(*fi).Rows, &Alpha, (*li).Array, &incx, (*fi).Array, &incy );
}

void EffK_Calc_Effective_Force( const Sp_MatrixVector *const Mass, const Sp_MatrixVector *const Damp,
				const Dense_MatrixVector *const Disp, const Dense_MatrixVector *const Vel,
				const Dense_MatrixVector *const Acc, Dense_MatrixVector *const Tempvec,
				const float a0, const float a1, const float a2,
				const float a3, const float a4, const float a5,
				Dense_MatrixVector *const Eff_Force )
{

     static int incx = 1, incy = 1;
     static char trans = 'N';
     static char matdescra[6] = {'S', 'U', 'N', 'C'};
     static float Alpha, Beta;

     /* BLAS: tempvec = Disp */
     scopy_( &Tempvec->Rows, Disp->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = a0*Disp = a0*tempvec */
     Alpha = a0;
     sscal_( &Tempvec->Rows, &Alpha, Tempvec->Array, &incx );
     /* BLAS: tempvec = a0*Disp + a2*Vel = tempvec + a2*Vel */
     Alpha = a2;
     saxpy_( &Tempvec->Rows, &Alpha, Vel->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = a0*Disp + a2*Vel + a3*Acc = tempvec + a3*Acc */
     Alpha = a3;
     saxpy_( &Tempvec->Rows, &Alpha, Acc->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: Eff_Force = Mass*(a0*Disp + a2*Vel + a3*Acc) = Mass*tempvec */
     Alpha = 1.0; Beta = 0.0;
     mkl_scsrmv( &trans, &Tempvec->Rows, &Tempvec->Rows, &Alpha, matdescra, Mass->Values, Mass->Columns, Mass->RowIndex, &Mass->RowIndex[1], Tempvec->Array, &Beta, Eff_Force->Array );

     /* BLAS: tempvec = Disp */
     scopy_( &Tempvec->Rows, Disp->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = a1*Disp = a0*tempvec */
     Alpha = a1;
     sscal_( &Tempvec->Rows, &Alpha, Tempvec->Array, &incx );
     /* BLAS: tempvec = a1*Disp + a4*Vel = tempvec + a4*Vel */
     Alpha = a4;
     saxpy_( &Tempvec->Rows, &Alpha, Vel->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = a1*Disp + a4*Vel + a5*Acc = tempvec + a5*Acc */
     Alpha = a5;
     saxpy_( &Tempvec->Rows, &Alpha, Acc->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: Eff_Force = Mass*(a0*Disp + a2*Vel + a3*Acc) + Damp*(a1*Disp + a4*Vel + a5*Acc) = Eff_Force + Damp*tempvec */
     Alpha = 1.0; Beta = 1.0;
     mkl_scsrmv( &trans, &Tempvec->Rows, &Tempvec->Rows, &Alpha, matdescra, Damp->Values, Damp->Columns, Damp->RowIndex, &Damp->RowIndex[1], Tempvec->Array, &Beta, Eff_Force->Array );
     
}

void EffK_ComputeU0( const Dense_MatrixVector *const Eff_Force, const Dense_MatrixVector *const In_Load,
		     const Dense_MatrixVector *const Err_Force, const float PID_P, const Sp_MatrixVector *const Keinv, Dense_MatrixVector *const Tempvec, Dense_MatrixVector *const Disp0 )
{
     static int incx = 1, incy = 1;
     static float Alpha = 1.0, Beta = 0.0;
     static char trans = 'N';
     static char matdescra[6] = {'S', 'U', 'N', 'C'};

     /* BLAS: tempvec = Eff_Force */
     scopy_( &Tempvec->Rows, Eff_Force->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = Eff_Force + LoadTdT = tempvec + LoadTdT */
     saxpy_( &Tempvec->Rows, &Alpha, In_Load->Array, &incx, Tempvec->Array, &incy );
     
     /* BLAS: tempvec = Eff_Force + LoadTdT + Err_Force = tempvec + Err_Force. The sign of Err_Force was already applied when calculating it. */
     Alpha = PID_P;
     saxpy_( &Tempvec->Rows, &Alpha, Err_Force->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: Disp0 = Keinv*(Eff_Force + LoadTdT + Err_Force) = Keinv*Tempvec */
     Alpha = 1.0;
     mkl_scsrmv( &trans, &Tempvec->Rows, &Tempvec->Rows, &Alpha, matdescra, Keinv->Values, Keinv->Columns, Keinv->RowIndex, &Keinv->RowIndex[1], Tempvec->Array, &Beta, Disp0->Array );
}

void CreateVectorXm( const Dense_MatrixVector *const VectorX, Dense_MatrixVector *const VectorXm, const int PosCouple, const int OrderC )
{

	static int incx, incy;
	static int length;
	static int posX, posXm;

	incx = 1; incy = 1;

	/* Copy the first part of the vector (from 1 to position PosCoupl -1) */
	length = PosCouple - 1;
	scopy_( &length, (*VectorX).Array, &incx, (*VectorXm).Array, &incy );

	/* Copy the second part of the vector (from position PosCoupl + 1 until the end) */
	length = (*VectorX).Rows - (PosCouple + OrderC - 1);
	posX = PosCouple + OrderC - 1;
	posXm = PosCouple - 1;

	scopy_( &length, &(*VectorX).Array[posX], &incx, &(*VectorXm).Array[posXm], &incy );

}

void CreateVectorXc( const Dense_MatrixVector *const VecX, float *VecXc, const int PosCouple, int OrderC )
{


	static int incx, incy;

	incx = 1; incy = 1;

	scopy_( &OrderC, &(*VecX).Array[PosCouple - 1], &incx, VecXc, &incy );
}

