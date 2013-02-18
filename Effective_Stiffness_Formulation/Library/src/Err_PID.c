#include "Error_Compensation.h"

#if _SPARSE_
#include <mkl.h>
#else
#include "Netlib.h"
#endif

void Compute_Force_Error( const MatrixVector *const Mass, const MatrixVector *const Damp, const MatrixVector *Stiff,
			  const MatrixVector *const AccTdT, const MatrixVector *const VelTdT, const MatrixVector *const DispTdT,
			  const MatrixVector *const fc, const MatrixVector *const LoadTdT, const PID_t *const PID, MatrixVector *const fe )
{

     static int incx = 1, incy = 1;
     static double Alpha, Beta;
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
void Compute_Force_Error_Sparse( const Sp_MatrixVector *const Mass, const Sp_MatrixVector *const Damp, const Sp_MatrixVector *Stiff,
				 const MatrixVector *const AccTdT, const MatrixVector *const VelTdT, const MatrixVector *const DispTdT,
				 const MatrixVector *const fc, const MatrixVector *const LoadTdT, const PID_t *const PID, MatrixVector *const fe )
{

     static int incx = 1, incy = 1;
     static double Alpha, Beta;
     static char trans = 'N';
     static char matdescra[6] = {'S', 'U', 'N', 'F'};



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
