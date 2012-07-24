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

#include "ComputeU0.h"
#include "MatrixVector.h"
#include "Netlib.h"

void Calculatefi( MatrixVector *const fi, const MatrixVector *const fc, const MatrixVector *const li, const MatrixVector *const Deltaf )
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

void EffK_Calc_Effective_Force( const MatrixVector *const Mass, const MatrixVector *const Damp,
				const MatrixVector *const Disp, const MatrixVector *const Vel,
				const MatrixVector *const Acc, MatrixVector *const Tempvec,
				const float a0, const float a1, const float a2,
				const float a3, const float a4, const float a5,
				MatrixVector *const Eff_Force )
{

     static int incx = 1, incy = 1;
     static char uplo = 'L';
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
     ssymv_( &uplo, &Tempvec->Rows, &Alpha, Mass->Array, &Tempvec->Rows, Tempvec->Array, &incx, &Beta, Eff_Force->Array, &incy );

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
     ssymv_( &uplo, &Tempvec->Rows, &Alpha, Damp->Array, &Tempvec->Rows, Tempvec->Array, &incx, &Beta, Eff_Force->Array, &incy );
     
}

void EffK_ComputeU0( const MatrixVector *const Eff_Force, const MatrixVector *const In_Load,
		     const MatrixVector *const Err_Force, const float PID_P, const MatrixVector *const Keinv,
		     MatrixVector *const Tempvec, MatrixVector *const Disp0 )
{
     static int incx = 1, incy = 1;
     static float Alpha = 1.0, Beta = 0.0;
     static char uplo = 'L';

     /* BLAS: tempvec = Eff_Force */
     scopy_( &Tempvec->Rows, Eff_Force->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = Eff_Force + LoadTdT = tempvec + LoadTdT */
     saxpy_( &Tempvec->Rows, &Alpha, In_Load->Array, &incx, Tempvec->Array, &incy );
     
     /* BLAS: tempvec = Eff_Force + LoadTdT + Err_Force = tempvec + Err_Force. The sign of Err_Force was already applied when calculating it. */
     Alpha = PID_P;
     saxpy_( &Tempvec->Rows, &Alpha, Err_Force->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: Disp0 = Keinv*(Eff_Force + LoadTdT + Err_Force) = Keinv*Tempvec */
     Alpha = 1.0;
     ssymv_( &uplo, &Tempvec->Rows, &Alpha, Keinv->Array, &Tempvec->Rows, Tempvec->Array, &incx, &Beta, Disp0->Array, &incy );
}

void CreateVectorXm( const MatrixVector *const VectorX, MatrixVector *const VectorXm, const Coupling_Node *const CNodes )
{

	static int incx, incy;
	static int Length;
	static int PosX, PosXm;
	static int icoup;    /* Counter for the coupling nodes */

	incx = 1; incy = 1;
	PosX = 0; PosXm = 0;
	for( icoup = 0; icoup < CNodes->Order; icoup++ ){

	     /* Copy the part of the vector between twwo positions */
	     Length = CNodes->Array[icoup] - PosX - 1;
	     scopy_( &Length, &VectorX->Array[PosX], &incx, &VectorXm->Array[PosXm], &incy );
	     /* Update the values of the position in the vectors */
	     PosX = CNodes->Array[icoup];
	     PosXm = PosXm + Length;
	}

	/* Copy the elements from the last position until the end of the vector */
	Length = VectorX->Rows - CNodes->Array[CNodes->Order-1];
	scopy_( &Length, &VectorX->Array[PosX], &incx, &VectorXm->Array[PosXm], &incy );

}

void CreateVectorXc( const MatrixVector *const VecX, float *VecXc, const Coupling_Node *const CNodes )
{
	static int icoup;    /* Counter for the coupling nodes */

#pragma omp parallel for
	for( icoup = 0; icoup < CNodes->Order; icoup++ ){
	     VecXc[icoup] = VecX->Array[CNodes->Array[icoup] - 1];
	}
}

