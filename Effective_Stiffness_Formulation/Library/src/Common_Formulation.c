#include "MatrixVector.h"
#include "Common_Formulation.h"

void Compute_NewState( const MatrixVector *const IGain, const MatrixVector *const Eff_ForceT, const MatrixVector *const In_LoadT,
		       const MatrixVector *const Err_ForceT, MatrixVector *const Tempvec, MatrixVector *const VecTdT_0 )
{
     static int incx = 1, incy = 1;           /* Stride in the vectors */
     static double Alpha = 1.0, Beta = 0.0;   /* Constants for the BLAS routines */
     static char uplo = 'L';                  /* The lower part (upper part in C) will be used and the upper part (lower part
					       * in C) will strictly not be referenced */

     /* BLAS: tempvec = Eff_Force */
     dcopy_( &Tempvec->Rows, Eff_ForceT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = Eff_Force + LoadTdT = tempvec + LoadTdT */
     daxpy_( &Tempvec->Rows, &Alpha, In_LoadT->Array, &incx, Tempvec->Array, &incy );
     
     /* BLAS: tempvec = Eff_Force + LoadTdT - Err_Force = tempvec - Err_Force. */
     Alpha = -1.0;
     daxpy_( &Tempvec->Rows, &Alpha, Err_ForceT->Array, &incx, Tempvec->Array, &incy );

     /* BLAS: Disp0 = Keinv*(Eff_Force + LoadTdT + Err_Force) = Keinv*Tempvec */
     Alpha = 1.0;
     dsymv_( &uplo, &Tempvec->Rows, &Alpha, Gain->Array, &Tempvec->Rows, Tempvec->Array, &incx, &Beta,
	     VecTdT_0->Array, &incy );
}

void Join_NonCouplingPart( MatrixVector *const VecTdT_m, const MatrixVector *const Gain_m,
			   const MatrixVector *const fcprevsub, const Coupling_Node *const CNodes,
			   MatrixVector *const VecTdT )			  
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
     dgemv_( &trans, &Rows, &Cols, &Alpha, Gain_m->Array, &lda,
	     fcprevsub->Array, &incx, &Beta, VecTdT_m->Array, &incy );

     /* Copy the updated values into the complete displacement vector */
     PosX = 0; PosXm = 0;
     for ( icoup = 0; icoup < CNodes->Order; icoup++ ){
	  Length = CNodes->Array[icoup] - PosX -1;
	  dcopy_( &Length, &VecTdT_m->Array[PosXm], &incx, &VecTdT->Array[PosX], &incy );
	  PosX = CNodes->Array[icoup];
	  PosXm = PosXm + Length;
     }

     /* Add the elements between the final coupling node and the final element
      * of the complete displacement vector */
     Length = VecTdT->Rows - CNodes->Array[CNodes->Order -1];
     dcopy_( &Length, &VecTdT_m->Array[PosXm], &incx, &VecTdT->Array[PosX], &incy );	
}
