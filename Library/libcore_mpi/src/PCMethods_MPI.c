#include "PCMethods.h"
#include "MatrixVector_MPI.h"
#include "Definitions.h"

#if _MKL_
#include <mkl_pblas.h>
#else
#include "Netlib.h"
#endif

void PC_Correct_Acceleration_MPI( const PMatrixVector_t *const MInv, const PMatrixVector_t *const In_LoadT,
				  const PMatrixVector_t *const fc, PMatrixVector_t *const Tempvec,
				  PMatrixVector_t *const AccTdT_Corr )
{

     int ione = 1;
     int incx = 1, incy = 1;           /* Stride in the vectors */
     hysl_float_t Alpha = 1.0, Beta = 0.0;   /* Constants for the BLAS routines */
     char uplo = 'L';                  /* The lower part (upper part in C) will be used and the upper part
					* (lower part in C) will strictly not be referenced */

     /* PBLAS: tempvec = In_LoadT */
     hysl_pcopy( &Tempvec->GlobalSize.Row, In_LoadT->Array, &ione, &ione, In_LoadT->Desc, &incx,
	      Tempvec->Array, &ione, &ione, Tempvec->Desc, &incy );
     /* PBLAS: tempvec = In_LoadT + fc = tempvec + fc */
     hysl_pcopy( &Tempvec->GlobalSize.Row, fc->Array, &ione, &ione, fc->Desc, &incx,
	      Tempvec->Array, &ione, &ione, Tempvec->Desc, &incy );

     /* PBLAS: AccTdT_Corr = MInv*(In_LoadT + fc) = MInv*Tempvec */
     hysl_psymv( &uplo, &Tempvec->GlobalSize.Row, &Alpha, MInv->Array, &ione, &ione, MInv->Desc, 
	      Tempvec->Array, &ione, &ione, Tempvec->Desc, &incx,
	      &Beta, AccTdT_Corr->Array, &ione, &ione, AccTdT_Corr->Desc, &incy );
}
