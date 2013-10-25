#include "MatrixVector.h"
#include "New_State.h"
#include "Definitions.h"

#if _MKL_
#include <mkl_blas.h>
#else
#include "Netlib.h"
#endif

void Compute_NewState( const MatrixVector_t *const IGain, const MatrixVector_t *const Eff_ForceT,
		       const MatrixVector_t *const In_LoadT, const MatrixVector_t *const Err_ForceT,
		       MatrixVector_t *const Tempvec, MatrixVector_t *const VecTdT_0 )
{
     int incx = 1, incy = 1;           /* Stride in the vectors */
     HYSL_FLOAT Alpha = 1.0, Beta = 0.0;   /* Constants for the BLAS routines */
     char uplo = 'L';                  /* The lower part (upper part in C) will be used and the upper part
					* (lower part in C) will strictly not be referenced */

     /* BLAS: tempvec = Eff_Force */
     hysl_copy( &Tempvec->Rows, Eff_ForceT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = Eff_Force + LoadTdT = tempvec + LoadTdT */
     hysl_axpy( &Tempvec->Rows, &Alpha, In_LoadT->Array, &incx, Tempvec->Array, &incy );
     
     /* BLAS: tempvec = Eff_Force + LoadTdT - Err_Force = tempvec - Err_Force. */
     Alpha = -1.0;
     hysl_axpy( &Tempvec->Rows, &Alpha, Err_ForceT->Array, &incx, Tempvec->Array, &incy );

     /* BLAS: VecTdT_0 = Keinv*(Eff_Force + LoadTdT - Err_Force) = Keinv*Tempvec */
     Alpha = 1.0;
     hysl_symv( &uplo, &Tempvec->Rows, &Alpha, IGain->Array, &Tempvec->Rows, Tempvec->Array, &incx, &Beta,
	    VecTdT_0->Array, &incy );
}

void Compute_NewState_PS( const MatrixVector_t *const IGain, const MatrixVector_t *const Eff_ForceT,
			  const MatrixVector_t *const In_LoadT, const MatrixVector_t *const Err_ForceT,
			  MatrixVector_t *const Tempvec, MatrixVector_t *const VecTdT_0 )
{
     int incx = 1, incy = 1;           /* Stride in the vectors */
     HYSL_FLOAT Alpha = 1.0, Beta = 0.0;   /* Constants for the BLAS routines */
     char uplo = 'L';                  /* The lower part (upper part in C) will be used and the upper part
					* (lower part in C) will strictly not be referenced */

     /* BLAS: tempvec = Eff_Force */
     hysl_copy( &Tempvec->Rows, Eff_ForceT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = Eff_Force + LoadTdT = tempvec + LoadTdT */
     hysl_axpy( &Tempvec->Rows, &Alpha, In_LoadT->Array, &incx, Tempvec->Array, &incy );
     
     /* BLAS: tempvec = Eff_Force + LoadTdT - Err_Force = tempvec - Err_Force. */
     Alpha = -1.0;
     hysl_axpy( &Tempvec->Rows, &Alpha, Err_ForceT->Array, &incx, Tempvec->Array, &incy );

     /* BLAS: Disp0 = Keinv*(Eff_Force + LoadTdT + Err_Force) = Keinv*Tempvec */
     Alpha = 1.0;
     hysl_spmv( &uplo, &Tempvec->Rows, &Alpha, IGain->Array, Tempvec->Array, &incx, &Beta,
	    VecTdT_0->Array, &incy );
}
