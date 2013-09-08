#include "MatrixVector.h"
#include "PCMethods.h"

#if _MKL_
#include <mkl_blas.h>
#else
#include "Netlib.h"
#endif

void PC_Correct_Acceleration( const MatrixVector_t *const MInv, const MatrixVector_t *const In_LoadT,
			      const MatrixVector_t *const fc, MatrixVector_t *const Tempvec,
			      MatrixVector_t *const AccTdT_Corr )
{

     int incx = 1, incy = 1;           /* Stride in the vectors */
     double Alpha = 1.0, Beta = 0.0;   /* Constants for the BLAS routines */
     char uplo = 'L';                  /* The lower part (upper part in C) will be used and the upper part
					* (lower part in C) will strictly not be referenced */

     /* BLAS: tempvec = In_LoadT */
     dcopy( &Tempvec->Rows, In_LoadT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = In_LoadT + fc = tempvec + fc */
     daxpy( &Tempvec->Rows, &Alpha, fc->Array, &incx, Tempvec->Array, &incy );

     /* BLAS: AccTdT_Corr = MInv*(In_LoadT + fc) = MInv*Tempvec */
     dsymv( &uplo, &Tempvec->Rows, &Alpha, MInv->Array, &Tempvec->Rows, Tempvec->Array, &incx, &Beta,
	    AccTdT_Corr->Array, &incy );
}

void PC_Correct_Acceleration_PS( const MatrixVector_t *const MInv, const MatrixVector_t *const In_LoadT,
				 const MatrixVector_t *const fc, MatrixVector_t *const Tempvec,
				 MatrixVector_t *const AccTdT_Corr )
{

     int incx = 1, incy = 1;           /* Stride in the vectors */
     double Alpha = 1.0, Beta = 0.0;   /* Constants for the BLAS routines */
     char uplo = 'L';                  /* The lower part (upper part in C) will be used and the upper part
					* (lower part in C) will strictly not be referenced */

     /* BLAS: tempvec = In_LoadT */
     dcopy( &Tempvec->Rows, In_LoadT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = In_LoadT + fc = tempvec + fc */
     daxpy( &Tempvec->Rows, &Alpha, fc->Array, &incx, Tempvec->Array, &incy );

     /* BLAS: AccTdT_Corr = MInv*(In_LoadT + fc) = MInv*Tempvec */
     dspmv( &uplo, &Tempvec->Rows, &Alpha, MInv->Array, Tempvec->Array, &incx, &Beta,
	    AccTdT_Corr->Array, &incy );
}
