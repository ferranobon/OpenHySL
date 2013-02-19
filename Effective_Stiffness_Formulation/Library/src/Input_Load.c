#include "Common_Formulation.h"
#include "MatrixVector.h"

void InputLoad_AbsValues( MatrixVector *const InLoad, const MatrixVector *const Stiff,
			  const MatrixVector *const Damp, const MatrixVector *const GDisp,
			  const MatrixVector *const GVel )
{

     static int incx, incy;       /* Stride in the vectors for BLAS library */
     static double Alpha, Beta;   /* Constants to use in the BLAS library */
     static char uplo;            /* Character to use in the BLAS library */

     incx = 1; incy = 1;
     Alpha = 1.0; Beta = 0.0;
     uplo = 'L';                  /* The lower part (FORTRAN) of the symmetric matrix is
				   * referenced (upper part in C) */


     /* BLAS: li = K*ug */
     dsymv_( &uplo, &InLoad->Rows, &Alpha, Stiff->Array, &InLoad->Rows, GDisp->Array, &incx, &Beta,
	     InLoad->.Array, &incy );

     /* BLAS: li = K*ug + C*vg = li + C*vg */
     Beta = 1.0;
     dsymv_( &uplo, &InLoad->Rows, &Alpha, Damp->Array, &InLoad->Rows, GVel->Array, &incx, &Beta,
	     InLoad->Array, &incy );
}

void InputLoad_RelValues( MatrixVector *const InLoad, const MatrixVector *const Mass,
			  const MatrixVector *const GAcc )
{

     static int incx = 1, incy = 1;           /* Stride in the vectors for BLAS library */
     static double Alpha = -1.0, Beta = 0.0;  /* Constants to use in the BLAS library */
     static char uplo = 'L';                  /* Character defining that the lower part (FORTRAN) of the
					       * symmetric matrix is referenced (upper part in C) */

     /* BLAS: li = -M*ag */
     dsymv_( &uplo, &InLoad->Rows, &Alpha, Mass->Array, &InLoad->Rows, GAcc->Array, &incx, &Beta,
	     InLoad->Array, &incy );
}

#if _SPARSE_
void InputLoad_AbsValues_Sp( MatrixVector *const InLoad, const MatrixVector_Sp *const Stiff,
				   const MatrixVector_Sp *const Damp, const MatrixVector *const GDisp,
				   const MatrixVector *const GVel )
{

     static double Alpha, Beta;        /* Constants to use in the Sparse BLAS routines */
     static char trans = 'N';          /* No transpose operation */
     static char matdescra[6] = {'S',  /* The matrix is symmetric */
				 'U',  /* The upper part is referenced */
				 'N',  /* Non-unit values in the diagonal */
				 'F'}; /* One based index */

     Alpha = 1.0; Beta = 0.0;

     /* Sparse BLAS: li = K*ug */
     mkl_dcsrmv( &trans, &InLoad->Rows, &InLoad->Rows, &Alpha, matdescra, Stiff->Values, Stiff->Columns,
		 Stiff->RowIndex, &Stiff->RowIndex[1], GDisp->Array, &Beta, InLoad->Array );

     /* Sparse BLAS: li = K*ug + C*vg = li + C*vg */
     Beta = 1.0;
     mkl_dcsrmv( &trans, &InLoad->Rows, &InLoad->Rows, &Alpha, matdescra, Damp->Values, Damp->Columns,
		 Damp->RowIndex, &Damp->RowIndex[1], GVel->Array, &Beta, InLoad->Array );
}

void InputLoad_RelValues_Sp( MatrixVector *const InLoad, const MatrixVector_Sp *const Mass,
			     const MatrixVector *const GAcc )
{
     
     static double Alpha, Beta;        /* Constants to use in the Sparse BLAS routines */
     static char trans = 'N';          /* No transpose operation */
     static char matdescra[6] = {'S',  /* The matrix is symmetric */
				 'U',  /* The upper part is referenced */
				 'N',  /* Non-unit values in the diagonal */
				 'F'}; /* One based index */

     Alpha = -1.0; Beta = 0.0;

     /* Sparse BLAS: li = -M*ag */
     mkl_dcsrmv( &trans, &InLoad->Rows, &InLoad->Rows, &Alpha, matdescra, Mass->Values, Mass->Columns,
		 Mass->RowIndex, &Mass->RowIndex[1], GAcc->Array, &Beta, InLoad->Array );
}

#endif
