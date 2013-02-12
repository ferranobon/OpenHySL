#include "MatrixVector.h"
#include "Common_Formulation.h"

void Compute_NewState( const MatrixVector *const Eff_Force, const MatrixVector *const In_Load,
		       const MatrixVector *const Err_Force, const double PID_P, const MatrixVector *const Gain,
		       MatrixVector *const Tempvec, MatrixVector *const New_State )
{
     static int incx = 1, incy = 1;           /* Stride in the vectors */
     static double Alpha = 1.0, Beta = 0.0;   /* Constants for the BLAS routines */
     static char uplo = 'L';                  /* The lower part (upper part in C) will be used and the upper part (lower part
					       * in C) will strictly not be referenced */

     /* BLAS: tempvec = Eff_Force */
     dcopy_( &Tempvec->Rows, Eff_Force->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = Eff_Force + LoadTdT = tempvec + LoadTdT */
     daxpy_( &Tempvec->Rows, &Alpha, In_Load->Array, &incx, Tempvec->Array, &incy );
     
     /* BLAS: tempvec = Eff_Force + LoadTdT + Err_Force = tempvec + Err_Force. The sign of Err_Force was already applied
      * when calculating it. */
     Alpha = PID_P;
     daxpy_( &Tempvec->Rows, &Alpha, Err_Force->Array, &incx, Tempvec->Array, &incy );

     /* BLAS: Disp0 = Keinv*(Eff_Force + LoadTdT + Err_Force) = Keinv*Tempvec */
     Alpha = 1.0;
     dsymv_( &uplo, &Tempvec->Rows, &Alpha, Gain->Array, &Tempvec->Rows, Tempvec->Array, &incx, &Beta,
	     New_State->Array, &incy );
}

void JoinNonCouplingPart( MatrixVector *const DispTdT_m, const MatrixVector *const Keinv_m,
			  const MatrixVector *const fcprevsub, MatrixVector *const DispTdT,
			  const Coupling_Node *const CNodes )
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
	Rows = Keinv_m->Rows;
	Cols = Keinv_m->Cols;
	lda = Max( 1, Keinv_m->Rows);

	/* Update the DispTdT_m displacments to include the effects of the coupling force */
	/* BLAS: DispTdT_m = Keinv_m*fcprevsub */
	dgemv_( &trans, &Rows, &Cols, &Alpha, Keinv_m->Array, &lda,
		fcprevsub->Array, &incx, &Beta, DispTdT_m->Array, &incy );

	/* Copy the updated values into the complete displacement vector */
	PosX = 0; PosXm = 0;
	for ( icoup = 0; icoup < CNodes->Order; icoup++ ){
	     Length = CNodes->Array[icoup] - PosX -1;
	     dcopy_( &Length, &DispTdT_m->Array[PosXm], &incx, &DispTdT->Array[PosX], &incy );
	     PosX = CNodes->Array[icoup];
	     PosXm = PosXm + Length;
	}

	/* Add the elements between the final coupling node and the final element
	 * of the complete displacement vector */
	Length = DispTdT->Rows - CNodes->Array[CNodes->Order -1];
	dcopy_( &Length, &DispTdT_m->Array[PosXm], &incx, &DispTdT->Array[PosX], &incy );	
}
