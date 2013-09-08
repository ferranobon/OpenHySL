#include "PCMethods.h"
#include "MatrixVector.h"
#include "MatrixVector_Sp.h"

#include <mkl_blas.h>
#include <mkl_spblas.h>

void PC_Correct_Acceleration_Sp( const MatrixVector_Sp_t *const MInv, const MatrixVector_t *const In_LoadT,
				 const MatrixVector_t *const fc, MatrixVector_t *const Tempvec,
				 MatrixVector_t *const AccTdT_Corr )
{
     int incx = 1, incy = 1;    /* Stride in the vectors */
     double Alpha, Beta;        /* Constants to use in the Sparse BLAS routines */
     char trans = 'N';          /* No transpose operation */
     char matdescra[6] = {'S',  /* The matrix is symmetric */
			  'U',  /* The upper part is referenced */
			  'N',  /* Non-unit values in the diagonal */
			  'F'}; /* One based index */

     Alpha = 1.0; Beta = 0.0;

     /* BLAS: tempvec = In_LoadT */
     dcopy( &Tempvec->Rows, In_LoadT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = In_LoadT + fc = tempvec + fc */
     daxpy( &Tempvec->Rows, &Alpha, fc->Array, &incx, Tempvec->Array, &incy );

     /* Sparse BLAS: AccTdT_Corr = MInv*(In_LoadT + fc) = MInv*Tempvec */
     mkl_dcsrmv( &trans, &Tempvec->Rows, &Tempvec->Rows, &Alpha, matdescra, MInv->Values, MInv->Columns,
		 MInv->RowIndex, &MInv->RowIndex[1], Tempvec->Array, &Beta, AccTdT_Corr->Array );
}
