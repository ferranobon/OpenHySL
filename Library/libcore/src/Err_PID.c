#include "Error_Compensation.h"
#include "MatrixVector.h"
#include "Definitions.h"

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
     HYSL_FLOAT Alpha, Beta;      /* Constants to use in the Sparse BLAS routines */
     char uplo = 'L';

     /* BLAS: fe = fc */
     hysl_copy( &fe->Rows, fc->Array, &incx, fe->Array, &incy );
     /* BLAS: fe = fc + LoadTdT */
     Alpha = 1.0;
     hysl_axpy( &fe->Rows, &Alpha, LoadTdT->Array, &incx, fe->Array, &incy );

     /* BLAS: fe =  fc + LoadTdT - Mass*AccTdT = fe - Mass*AccTdT */
     Alpha = -1.0; Beta = 1.0;
     hysl_symv( &uplo, &fe->Rows, &Alpha, Mass->Array, &fe->Rows, AccTdT->Array, &incx, &Beta, fe->Array, &incy );
     /* BLAS: fe =  fc + LoadTdT - Mass*AccTdT - Damp*VelTdT = fe - Damp*VelTdT */
     hysl_symv( &uplo, &fe->Rows, &Alpha, Damp->Array, &fe->Rows, VelTdT->Array, &incx, &Beta, fe->Array, &incy );
     /* BLAS: fe =  fc + LoadTdT - Mass*AccTdT - Damp*VelTdT - Stiff*DispTdT = fe - Stiff*DispTdT */
     hysl_symv( &uplo, &fe->Rows, &Alpha, Stiff->Array, &fe->Rows, DispTdT->Array, &incx, &Beta, fe->Array, &incy );

     /* Apply the PID */
     /* BLAS: fe = P*fe */
     Alpha = -PID->P;
     hysl_scal( &fe->Rows, &Alpha, fe->Array, &incx );
}

void ErrorForce_PID_PS( const MatrixVector_t *const Mass, const MatrixVector_t *const Damp, const MatrixVector_t *Stiff,
		     const MatrixVector_t *const AccTdT, const MatrixVector_t *const VelTdT, const MatrixVector_t *const DispTdT,
		     const MatrixVector_t *const fc, const MatrixVector_t *const LoadTdT, const PID_t *const PID, MatrixVector_t *const fe )
{

     int incx = 1, incy = 1;  /* Stride in the vectors for BLAS routines */
     HYSL_FLOAT Alpha, Beta;      /* Constants to use in the Sparse BLAS routines */
     char uplo = 'L';

     /* BLAS: fe = fc */
     hysl_copy( &fe->Rows, fc->Array, &incx, fe->Array, &incy );
     /* BLAS: fe = fc + LoadTdT */
     Alpha = 1.0;
     hysl_axpy( &fe->Rows, &Alpha, LoadTdT->Array, &incx, fe->Array, &incy );

     /* BLAS: fe =  fc + LoadTdT - Mass*AccTdT = fe - Mass*AccTdT */
     Alpha = -1.0; Beta = 1.0;
     hysl_spmv( &uplo, &fe->Rows, &Alpha, Mass->Array, AccTdT->Array, &incx, &Beta, fe->Array, &incy );
     /* BLAS: fe =  fc + LoadTdT - Mass*AccTdT - Damp*VelTdT = fe - Damp*VelTdT */
     hysl_spmv( &uplo, &fe->Rows, &Alpha, Damp->Array, VelTdT->Array, &incx, &Beta, fe->Array, &incy );
     /* BLAS: fe =  fc + LoadTdT - Mass*AccTdT - Damp*VelTdT - Stiff*DispTdT = fe - Stiff*DispTdT */
     hysl_spmv( &uplo, &fe->Rows, &Alpha, Stiff->Array, DispTdT->Array, &incx, &Beta, fe->Array, &incy );

     /* Apply the PID */
     /* BLAS: fe = P*fe */
     Alpha = -PID->P;
     hysl_scal( &fe->Rows, &Alpha, fe->Array, &incx );
}

void ErrorForce_PID_HHT( const MatrixVector_t *const Mass, const MatrixVector_t *const Damp, const MatrixVector_t *Stiff,
			 const MatrixVector_t *const VelT, const MatrixVector_t *const DispT,
			 const MatrixVector_t *const AccTdT, const MatrixVector_t *const VelTdT, const MatrixVector_t *const DispTdT,
			 const MatrixVector_t *const fc, const MatrixVector_t *const LoadT, const MatrixVector_t *const LoadTdT,
			 const HYSL_FLOAT Alpha_HHT, const PID_t *const PID, MatrixVector_t *const fe )
{

     int incx = 1, incy = 1;  /* Stride in the vectors for BLAS routines */
     HYSL_FLOAT Alpha, Beta;      /* Constants to use in the Sparse BLAS routines */
     char uplo = 'L';

     /* BLAS: fe = fc */
     hysl_copy( &fe->Rows, fc->Array, &incx, fe->Array, &incy );
     /* BLAS: fe = fc + (1 + Alpha_HHT)*LoadTdT = fe + (1 + Alpha_HHT)*LoadTdT */
     Alpha = (1.0 + Alpha_HHT);
     hysl_axpy( &fe->Rows, &Alpha, LoadTdT->Array, &incx, fe->Array, &incy );
     /* BLAS: fe = fc + (1 + Alpha_HHT)*LoadTdT - Alpha_HHT*LoadT = fe - Alpha_HHT*LoadT */
     Alpha = -Alpha_HHT;
     hysl_axpy( &fe->Rows, &Alpha, LoadTdT->Array, &incx, fe->Array, &incy );
     
     /* BLAS: fe =  fc + (1 + Alpha_HHT)*LoadTdT - Alpha_HHT*LoadT - Mass*AccTdT = fe - Mass*AccTdT */
     Alpha = -1.0; Beta = 1.0;
     hysl_symv( &uplo, &fe->Rows, &Alpha, Mass->Array, &fe->Rows, AccTdT->Array, &incx, &Beta, fe->Array, &incy );
     /* BLAS: fe =  fc + (1 + Alpha_HHT)*LoadTdT - Alpha_HHT*LoadT - Mass*AccTdT - (1 + Alpha_HHT) Damp*VelTdT =
      * fe - (1 + Alpha_HHT)*Damp*VelTdT */
     Alpha = -(1.0 + Alpha_HHT);
     hysl_symv( &uplo, &fe->Rows, &Alpha, Damp->Array, &fe->Rows, VelTdT->Array, &incx, &Beta, fe->Array, &incy );
     /* BLAS: fe =  fc + (1 + Alpha_HHT)*LoadTdT - Alpha_HHT*LoadT - Mass*AccTdT - (1 + Alpha_HHT)*Damp*VelTdT
      * - (1 + Alpha_HHT)*Stiff*DispTdT = fe - (1 + Alpha_HHT)*Stiff*DispTdT */
     hysl_symv( &uplo, &fe->Rows, &Alpha, Stiff->Array, &fe->Rows, DispTdT->Array, &incx, &Beta, fe->Array, &incy );

     /* BLAS: fe =  fc + (1 + Alpha_HHT)*LoadTdT - Alpha_HHT*LoadT - Mass*AccTdT - (1 + Alpha_HHT)*Damp*VelTdT
      * - (1 + Alpha_HHT)*Stiff*DispTdT + Alpha_HHT*Damp*VelT = fe + Alpha_HHT*Damp*VelTdT */
     Alpha = Alpha_HHT;
     hysl_symv( &uplo, &fe->Rows, &Alpha, Damp->Array, &fe->Rows, VelT->Array, &incx, &Beta, fe->Array, &incy );
     /* BLAS: fe =  fc + (1 + Alpha_HHT)*LoadTdT - Alpha_HHT*LoadT - Mass*AccTdT - (1 + Alpha_HHT)*Damp*VelTdT
      * - (1 + Alpha_HHT)*Stiff*DispTdT + Alpha_HHT*Damp*VelT + Alpha_HHT*Damp*DispTdT = fe + Alpha_HHT*Damp*DispTdT */
     hysl_symv( &uplo, &fe->Rows, &Alpha, Stiff->Array, &fe->Rows, DispT->Array, &incx, &Beta, fe->Array, &incy );


     /* Apply the PID */
     /* BLAS: fe = P*fe */
     Alpha = -PID->P;
     hysl_scal( &fe->Rows, &Alpha, fe->Array, &incx );
}
