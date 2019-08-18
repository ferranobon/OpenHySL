#include "MatrixVector_Sp.h"
#include "New_State.h"
#include "Definitions.h"

#include <mkl_blas.h>
#include <mkl_spblas.h>

void Compute_NewState_Zienkiewicz_Sp( const MatrixVector_t *const Meff, const MatrixVector_Sp_t *const MatA,
				      const MatrixVector_Sp_t *const MatB, const MatrixVector_t *const DispT,
				      const MatrixVector_t *const DispT0, const MatrixVector_t *const In_LoadT,
				      const MatrixVector_t *const Err_ForceT, const MatrixVector_t *const ForceT,
				      const MatrixVector_t *const ForceT0, const hysl_float_t a8, const hysl_float_t a17,
				      const hysl_float_t a18, MatrixVector_t *const Tempvec, MatrixVector_t *const VecTdT_0 )
{

     int incx = 1, incy = 1;              /* Stride in the vectors */
     char trans = 'N';                    /* No transpose operation */
     char matdescra[6] = {'S',            /* The matrix is symmetric */
			  'U',            /* The upper part is referenced */
			  'N',            /* Non-unit values in the diagonal */
			  'F'};           /* One based index */
     hysl_float_t Alpha = 1.0, Beta = 1.0;  /* Constants for the BLAS routines */
     char uplo = 'L';                     /* The lower part (upper part in C) will be used and the upper part
					   * (lower part in C) will strictly not be referenced */


    /* BLAS: Tempvec = ForceT */
     hysl_copy( &Tempvec->Rows, ForceT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: Tempvec = (1/2 - 2*beta + gamma)*ForceT = a17*ForceT */
     Alpha = a17;
     hysl_scal( &Tempvec->Rows, &Alpha, Tempvec->Array, &incx );
     /* BLAS: Tempvec = (1/2 - 2*beta + gamma)*Deltat^2*ForceT + (1/2 + beta - gamma)*ForceT0 = Tempvec + a18*ForceT0 */
     Alpha = a18;
     hysl_axpy( &Tempvec->Rows, &Alpha, ForceT0->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: Tempvec = 1/2*beta*Deltat^2*In_LoadT +(1/2 - 2*beta + gamma)*ForceT + (1/2 + beta - gamma)*ForceT0 = Tempvec + a8*In_LoadT */
     Alpha = a8;
     hysl_axpy( &Tempvec->Rows, &Alpha, In_LoadT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: Tempvec = MatA*DispT + 1/2*beta*Deltat^2*In_LoadT +(1/2 - 2*beta + gamma)*ForceT + (1/2 + beta - gamma)*ForceT0 = Tempvec + MatA*DispT */
     Alpha = 1.0;
     hysl_mkl_csrmv( &trans, &Tempvec->Rows, &Tempvec->Rows, &Alpha, matdescra, MatA->Values, MatA->Columns, MatA->RowIndex,
		     &MatA->RowIndex[1], DispT->Array, &Beta, Tempvec->Array );
     /* BLAS: Tempvec = MatA*DispT - MatB*DispT0 + 1/2*beta*Deltat^2*In_LoadT +(1/2 - 2*beta + gamma)*ForceT
      * + (1/2 + beta - gamma)*ForceT0 = Tempvec - MatB*DispT0 */
     /* Note: MatB has already the minus sign applied */
     hysl_mkl_csrmv( &trans, &Tempvec->Rows, &Tempvec->Rows, &Alpha, matdescra, MatB->Values, MatB->Columns, MatB->RowIndex,
		     &MatB->RowIndex[1], DispT0->Array, &Beta, Tempvec->Array );
     /* BLAS: Tempvec = MatA*DispT - MatB*DispT0 + 1/2*beta*Deltat^2*In_LoadT +(1/2 - 2*beta + gamma)*ForceT
      * + (1/2 + beta - gamma)*ForceT0 = Tempvec - MatB*DispT0 + 1/2*beta*Deltat^2*Err_ForceT = Tempvec +
      * a8*Err_ForceT */
     Alpha = a8;
     hysl_axpy( &Tempvec->Rows, &Alpha, Err_ForceT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: VecTdT_0 = Meff*Tempvec */
     Alpha = 1.0;
     Beta = 0.0;
     hysl_symv( &uplo, &VecTdT_0->Rows, &Alpha, Meff->Array, &Tempvec->Rows, Tempvec->Array, &incx, &Beta,
		VecTdT_0->Array, &incy );
}
