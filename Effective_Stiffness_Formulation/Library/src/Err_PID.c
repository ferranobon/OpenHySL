#include "Error_Compensation.h"
#include "MatrixVector.h"

#if _MKL_
#include <mkl_blas.h>
#else
#include "Netlib.h"
#endif

void ErrorForce_PID( const MatrixVector_t *const Mass, const MatrixVector_t *const Damp, const MatrixVector_t *Stiff,
		     const MatrixVector_t *const AccTdT, const MatrixVector_t *const VelTdT, const MatrixVector_t *const DispTdT,
		     const MatrixVector_t *const fc, const MatrixVector_t *const LoadTdT, const PID_t *const PID, MatrixVector_t *const fe )
{

     int incx = 1, incy = 1;  /* Stride in the vectors for BLAS routines */
     double Alpha, Beta;      /* Constants to use in the Sparse BLAS routines */
     char uplo = 'L';

     /* BLAS: fe = fc */
     dcopy( &fe->Rows, fc->Array, &incx, fe->Array, &incy );
     /* BLAS: fe = fc + LoadTdT */
     Alpha = 1.0;
     daxpy( &fe->Rows, &Alpha, LoadTdT->Array, &incx, fe->Array, &incy );

     /* BLAS: fe =  fc + LoadTdT - Mass*AccTdT = fe - Mass*AccTdT */
     Alpha = -1.0; Beta = 1.0;
     dsymv( &uplo, &fe->Rows, &Alpha, Mass->Array, &fe->Rows, AccTdT->Array, &incx, &Beta, fe->Array, &incy );
     /* BLAS: fe =  fc + LoadTdT - Mass*AccTdT - Damp*VelTdT = fe - Damp*VelTdT */
     dsymv( &uplo, &fe->Rows, &Alpha, Damp->Array, &fe->Rows, VelTdT->Array, &incx, &Beta, fe->Array, &incy );
     /* BLAS: fe =  fc + LoadTdT - Mass*AccTdT - Damp*VelTdT - Stiff*DispTdT = fe - Stiff*DispTdT */
     dsymv( &uplo, &fe->Rows, &Alpha, Stiff->Array, &fe->Rows, DispTdT->Array, &incx, &Beta, fe->Array, &incy );

     /* Apply the PID */
     /* BLAS: fe = P*fe */
     Alpha = -PID->P;
     dscal( &fe->Rows, &Alpha, fe->Array, &incx );
}

void ErrorForce_PID_PS( const MatrixVector_t *const Mass, const MatrixVector_t *const Damp, const MatrixVector_t *Stiff,
		     const MatrixVector_t *const AccTdT, const MatrixVector_t *const VelTdT, const MatrixVector_t *const DispTdT,
		     const MatrixVector_t *const fc, const MatrixVector_t *const LoadTdT, const PID_t *const PID, MatrixVector_t *const fe )
{

     int incx = 1, incy = 1;  /* Stride in the vectors for BLAS routines */
     double Alpha, Beta;      /* Constants to use in the Sparse BLAS routines */
     char uplo = 'L';

     /* BLAS: fe = fc */
     Alpha = 1.0;
     daxpy( &fe->Rows, &Alpha, fc->Array, &incx, fe->Array, &incy );
     /* BLAS: fe = fc + LoadTdT */
     daxpy( &fe->Rows, &Alpha, LoadTdT->Array, &incx, fe->Array, &incy );

     /* BLAS: fe =  fc + LoadTdT - Mass*AccTdT = fe - Mass*AccTdT */
     Alpha = -1.0; Beta = 1.0;
     dspmv( &uplo, &fe->Rows, &Alpha, Mass->Array, AccTdT->Array, &incx, &Beta, fe->Array, &incy );
     /* BLAS: fe =  fc + LoadTdT - Mass*AccTdT - Damp*VelTdT = fe - Damp*VelTdT */
     dspmv( &uplo, &fe->Rows, &Alpha, Damp->Array, VelTdT->Array, &incx, &Beta, fe->Array, &incy );
     /* BLAS: fe =  fc + LoadTdT - Mass*AccTdT - Damp*VelTdT - Stiff*DispTdT = fe - Stiff*DispTdT */
     dspmv( &uplo, &fe->Rows, &Alpha, Stiff->Array, DispTdT->Array, &incx, &Beta, fe->Array, &incy );

     /* Apply the PID */
     /* BLAS: fe = P*fe */
     Alpha = -PID->P;
     dscal( &fe->Rows, &Alpha, fe->Array, &incx );
}
