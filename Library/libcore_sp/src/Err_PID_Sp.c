#include "Error_Compensation.h"

#include "MatrixVector.h"
#include "MatrixVector_Sp.h"

#include "Definitions.h"

#include <mkl_blas.h>
#include <mkl_spblas.h>

void ErrorForce_PID_Sp( const MatrixVector_Sp_t *const Mass, const MatrixVector_Sp_t *const Damp, const MatrixVector_Sp_t *const Stiff,
			const MatrixVector_t *const AccTdT, const MatrixVector_t *const VelTdT, const MatrixVector_t *const DispTdT,
			const MatrixVector_t *const fc, const MatrixVector_t *const LoadTdT, const PID_t *const PID, MatrixVector_t *const fe )
{

     int incx = 1, incy = 1;    /* Stride in the vectors for BLAS routines */
     hysl_float_t Alpha, Beta;        /* Constants to use in the Sparse BLAS routines */
     char trans = 'N';          /* No transpose operation */
     char matdescra[6] = {'S',  /* The matrix is symmetric */
			  'U',  /* The upper part is referenced */
			  'N',  /* Non-unit values in the diagonal */
			  'F'}; /* One based index */

     /* BLAS: fe = fc */
     hysl_copy( &fe->Rows, fc->Array, &incx, fe->Array, &incy );
     /* BLAS: fe = fc + LoadTdT */
     Alpha = 1.0;
     hysl_axpy( &fe->Rows, &Alpha, LoadTdT->Array, &incx, fe->Array, &incy );

     /* Sparse BLAS: fe =  fc + LoadTdT - Mass*AccTdT = fe - Mass*AccTdT */
     Alpha = -1.0; Beta = 1.0;
     hysl_mkl_csrmv( &trans, &fe->Rows, &fe->Rows, &Alpha, matdescra, Mass->Values, Mass->Columns, Mass->RowIndex, &Mass->RowIndex[1], AccTdT->Array, &Beta, fe->Array );
      /* Sparse BLAS: fe =  fc + LoadTdT - Mass*AccTdT - Damp*VelTdT = fe - Damp*VelTdT */
     hysl_mkl_csrmv( &trans, &fe->Rows, &fe->Rows, &Alpha, matdescra, Damp->Values, Damp->Columns, Damp->RowIndex, &Damp->RowIndex[1], VelTdT->Array, &Beta, fe->Array );
     /* Sparse BLAS: fe =  fc + LoadTdT - Mass*AccTdT - Damp*VelTdT - Stiff*DispTdT = fe - Stiff*DispTdT */
     hysl_mkl_csrmv( &trans, &fe->Rows, &fe->Rows, &Alpha, matdescra, Stiff->Values, Stiff->Columns, Stiff->RowIndex, &Stiff->RowIndex[1], DispTdT->Array, &Beta, fe->Array );

     /* Apply the PID */
     /* BLAS: fe = P*fe */
     Alpha = -PID->P;
     hysl_scal( &fe->Rows, &Alpha, fe->Array, &incx );
}

void ErrorForce_PID_HHT_Sp ( const MatrixVector_Sp_t *const Mass, const MatrixVector_Sp_t *const Damp, const MatrixVector_Sp_t *const Stiff,
			     const MatrixVector_t *const VelT, const MatrixVector_t *const DispT,
			     const MatrixVector_t *const AccTdT, const MatrixVector_t *const VelTdT, const MatrixVector_t *const DispTdT,
			     const MatrixVector_t *const fc, const MatrixVector_t *const LoadT, const MatrixVector_t *const LoadTdT,
			     const hysl_float_t Alpha_HHT, const PID_t *const PID, MatrixVector_t *const fe )
{
     int incx = 1, incy = 1;    /* Stride in the vectors for BLAS routines */
     hysl_float_t Alpha, Beta;        /* Constants to use in the Sparse BLAS routines */
     char trans = 'N';          /* No transpose operation */
     char matdescra[6] = {'S',  /* The matrix is symmetric */
			  'U',  /* The upper part is referenced */
			  'N',  /* Non-unit values in the diagonal */
			  'F'}; /* One based index */
     
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
     hysl_mkl_csrmv( &trans, &fe->Rows, &fe->Rows, &Alpha, matdescra, Mass->Values, Mass->Columns, Mass->RowIndex, &Mass->RowIndex[1], AccTdT->Array, &Beta, fe->Array );
     /* BLAS: fe =  fc + (1 + Alpha_HHT)*LoadTdT - Alpha_HHT*LoadT - Mass*AccTdT - (1 + Alpha_HHT) Damp*VelTdT =
      * fe - (1 + Alpha_HHT)*Damp*VelTdT */
     Alpha = -(1.0 + Alpha_HHT);
     hysl_mkl_csrmv( &trans, &fe->Rows, &fe->Rows, &Alpha, matdescra, Damp->Values, Damp->Columns, Damp->RowIndex, &Damp->RowIndex[1], VelTdT->Array, &Beta, fe->Array );
     /* BLAS: fe =  fc + (1 + Alpha_HHT)*LoadTdT - Alpha_HHT*LoadT - Mass*AccTdT - (1 + Alpha_HHT)*Damp*VelTdT
      * - (1 + Alpha_HHT)*Stiff*DispTdT = fe - (1 + Alpha_HHT)*Stiff*DispTdT */
     hysl_mkl_csrmv( &trans, &fe->Rows, &fe->Rows, &Alpha, matdescra, Stiff->Values, Stiff->Columns, Stiff->RowIndex, &Stiff->RowIndex[1], DispTdT->Array, &Beta, fe->Array );

     /* BLAS: fe =  fc + (1 + Alpha_HHT)*LoadTdT - Alpha_HHT*LoadT - Mass*AccTdT - (1 + Alpha_HHT)*Damp*VelTdT
      * - (1 + Alpha_HHT)*Stiff*DispTdT + Alpha_HHT*Damp*VelT = fe + Alpha_HHT*Damp*VelTdT */
     Alpha = Alpha_HHT;
     hysl_mkl_csrmv( &trans, &fe->Rows, &fe->Rows, &Alpha, matdescra, Damp->Values, Damp->Columns, Damp->RowIndex, &Damp->RowIndex[1], VelT->Array, &Beta, fe->Array );
     /* BLAS: fe =  fc + (1 + Alpha_HHT)*LoadTdT - Alpha_HHT*LoadT - Mass*AccTdT - (1 + Alpha_HHT)*Damp*VelTdT
      * - (1 + Alpha_HHT)*Stiff*DispTdT + Alpha_HHT*Damp*VelT + Alpha_HHT*Damp*DispTdT = fe + Alpha_HHT*Damp*DispTdT */
     hysl_mkl_csrmv( &trans, &fe->Rows, &fe->Rows, &Alpha, matdescra, Stiff->Values, Stiff->Columns, Stiff->RowIndex, &Stiff->RowIndex[1], DispT->Array, &Beta, fe->Array );

     /* Apply the PID */
     /* BLAS: fe = P*fe */
     Alpha = -PID->P;
     hysl_scal( &fe->Rows, &Alpha, fe->Array, &incx );
}
