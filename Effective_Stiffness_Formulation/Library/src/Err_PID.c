#include "Error_Compensation.h"

#if _SPARSE_
#include <mkl.h>
#else
#include "Netlib.h"
#endif

void ErrorForce_PID( const MatrixVector_t *const Mass, const MatrixVector_t *const Damp, const MatrixVector *Stiff,
		     const MatrixVector_t *const AccTdT, const MatrixVector_t *const VelTdT, const MatrixVector_t *const DispTdT,
		     const MatrixVector_t *const fc, const MatrixVector_t *const LoadTdT, const PID_t *const PID, MatrixVector_t *const fe )
{

     static int incx = 1, incy = 1;  /* Stride in the vectors for BLAS routines */
     static double Alpha, Beta;      /* Constants to use in the Sparse BLAS routines */
     static char uplo = 'L';

     /* BLAS: fe = fc */
     daxpy_( &fe->Rows, &Alpha, fc->Array, &incx, fe->Array, &incy );
     /* BLAS: fe = fc + LoadTdT */
     daxpy_( &fe->Rows, &Alpha, LoadTdT->Array, &incx, fe->Array, &incy );

     /* BLAS: fe =  fc + LoadTdT - Mass*AccTdT = fe - Mass*AccTdT */
     Alpha = -1.0; Beta = 1.0;
     dsymv_( &uplo, &fe->Rows, &Alpha, Mass->Array, &fe->Rows, AccTdT->Array, &incx, &Beta, fe->Array, &incy );
     /* BLAS: fe =  fc + LoadTdT - Mass*AccTdT - Damp*VelTdT = fe - Damp*VelTdT */
     dsymv_( &uplo, &fe->Rows, &Alpha, Damp->Array, &fe->Rows, VelTdT->Array, &incx, &Beta, fe->Array, &incy );
     /* BLAS: fe =  fc + LoadTdT - Mass*AccTdT - Damp*VelTdT - Stiff*DispTdT = fe - Stiff*DispTdT */
     dsymv_( &uplo, &fe->Rows, &Alpha, Stiff->Array, &fe->Rows, DispTdT->Array, &incx, &Beta, fe->Array, &incy );

     /* Apply the PID */
     /* BLAS: fe = P*fe */
     Alpha = -PID->P;
     dscal( &fe->Rows, &Alpha, fe->Array, &incx )
}

#if _SPARSE_
void ErrorForce_PID_Sp( const MatrixVector_Sp_t *const Mass, const MatrixVector_Sp_t *const Damp, const MatrixVector_Sp *Stiff,
			const MatrixVector_t *const AccTdT, const MatrixVector_t *const VelTdT, const MatrixVector_t *const DispTdT,
			const MatrixVector_t *const fc, const MatrixVector_t *const LoadTdT, const PID_t *const PID, MatrixVector_t *const fe )
{

     static int incx = 1, incy = 1;    /* Stride in the vectors for BLAS routines */
     static double Alpha, Beta;        /* Constants to use in the Sparse BLAS routines */
     static char trans = 'N';          /* No transpose operation */
     static char matdescra[6] = {'S',  /* The matrix is symmetric */
				 'U',  /* The upper part is referenced */
				 'N',  /* Non-unit values in the diagonal */
				 'F'}; /* One based index */

     /* BLAS: fe = Mass*AccTdT */
     Alpha = 1.0; Beta = 0.0;
     mkl_dcsrmv( &trans, &fe->Rows, &fe->Rows, &Alpha, matdescra, Mass->Values, Mass->Columns, Mass->RowIndex, &Mass->RowIndex[1], AccTdT->Array, &Beta, fe->Array );
     /* BLAS: fe = Mass*AccTdT + Damp*VelTdT = fe + Damp*VelTdT */
     Beta = 1.0;
     mkl_dcsrmv( &trans, &fe->Rows, &fe->Rows, &Alpha, matdescra, Damp->Values, Damp->Columns, Damp->RowIndex, &Damp->RowIndex[1], VelTdT->Array, &Beta, fe->Array );
     /* BLAS: fe = Mass*AccTdT + Damp*VelTdT + Stiff*DispTdT = fe + Stiff*DispTdT */
     mkl_dcsrmv( &trans, &fe->Rows, &fe->Rows, &Alpha, matdescra, Stiff->Values, Stiff->Columns, Stiff->RowIndex, &Stiff->RowIndex[1], DispTdT->Array, &Beta, fe->Array );
     /* BLAS: fe = -(Mass*AccTdT + Damp*VelTdT + Stiff*DispTdT) = -fe */
     Alpha = -1.0;
     dscal_( &fe->Rows, &Alpha, fe->Array, &incx );
     /* BLAS: fe = fc -(Mass*AccTdT + Damp*VelTdT + Stiff*DispTdT) = fc -fe */
     Alpha = 1.0;
     daxpy_( &fe->Rows, &Alpha, fc->Array, &incx, fe->Array, &incy );
     /* BLAS: fe = LoatTdT + fc -(Mass*AccTdT + Damp*VelTdT + Stiff*DispTdT) = LoadTdT -fe */
     daxpy_( &fe->Rows, &Alpha, LoadTdT->Array, &incx, fe->Array, &incy );
}

#endif
