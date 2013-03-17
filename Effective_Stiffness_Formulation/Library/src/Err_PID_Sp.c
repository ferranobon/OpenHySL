#include "Error_Compensation.h"

#include "MatrixVector.h"
#include "MatrixVector_Sp.h"

#include <mkl_blas.h>
#include <mkl_spblas.h>

void ErrorForce_PID_Sp( const MatrixVector_Sp_t *const Mass, const MatrixVector_Sp_t *const Damp, const MatrixVector_Sp_t *const Stiff,
			const MatrixVector_t *const AccTdT, const MatrixVector_t *const VelTdT, const MatrixVector_t *const DispTdT,
			const MatrixVector_t *const fc, const MatrixVector_t *const LoadTdT, const PID_t *const PID, MatrixVector_t *const fe )
{

     int incx = 1, incy = 1;    /* Stride in the vectors for BLAS routines */
     double Alpha, Beta;        /* Constants to use in the Sparse BLAS routines */
     char trans = 'N';          /* No transpose operation */
     char matdescra[6] = {'S',  /* The matrix is symmetric */
			  'U',  /* The upper part is referenced */
			  'N',  /* Non-unit values in the diagonal */
			  'F'}; /* One based index */

     /* BLAS: fe = fc */
     dcopy( &fe->Rows, fc->Array, &incx, fe->Array, &incy );
     /* BLAS: fe = fc + LoadTdT */
     Alpha = 1.0;
     daxpy( &fe->Rows, &Alpha, LoadTdT->Array, &incx, fe->Array, &incy );

     /* Sparse BLAS: fe =  fc + LoadTdT - Mass*AccTdT = fe - Mass*AccTdT */
     Alpha = -1.0; Beta = 1.0;
     mkl_dcsrmv( &trans, &fe->Rows, &fe->Rows, &Alpha, matdescra, Mass->Values, Mass->Columns, Mass->RowIndex, &Mass->RowIndex[1], AccTdT->Array, &Beta, fe->Array );
      /* Sparse BLAS: fe =  fc + LoadTdT - Mass*AccTdT - Damp*VelTdT = fe - Damp*VelTdT */
     mkl_dcsrmv( &trans, &fe->Rows, &fe->Rows, &Alpha, matdescra, Damp->Values, Damp->Columns, Damp->RowIndex, &Damp->RowIndex[1], VelTdT->Array, &Beta, fe->Array );
     /* Sparse BLAS: fe =  fc + LoadTdT - Mass*AccTdT - Damp*VelTdT - Stiff*DispTdT = fe - Stiff*DispTdT */
     mkl_dcsrmv( &trans, &fe->Rows, &fe->Rows, &Alpha, matdescra, Stiff->Values, Stiff->Columns, Stiff->RowIndex, &Stiff->RowIndex[1], DispTdT->Array, &Beta, fe->Array );

     /* Apply the PID */
     /* BLAS: fe = P*fe */
     Alpha = -PID->P;
     dscal( &fe->Rows, &Alpha, fe->Array, &incx );
}
