#include "EffK_Formulation.h"
#include "MatrixVector_Sp.h"

#include "Definitions.h"

#include <mkl_blas.h>
#include <mkl_spblas.h>

void EffK_EffectiveForce_Sp( const MatrixVector_Sp_t *const Mass, const MatrixVector_Sp_t *const Damp, const MatrixVector_t *const DispT,
			     const MatrixVector_t *const VelT, const MatrixVector_t *const AccT, MatrixVector_t *const Tempvec,
			     const HYSL_FLOAT a0, const HYSL_FLOAT a1, const HYSL_FLOAT a2, const HYSL_FLOAT a3, const HYSL_FLOAT a4,
			     const HYSL_FLOAT a5, MatrixVector_t *const Eff_ForceT )
{

     int incx = 1, incy = 1;    /* Stride in the vectors */
     char trans = 'N';          /* No transpose operation */
     char matdescra[6] = {'S',  /* The matrix is symmetric */
			  'U',  /* The upper part is referenced */
			  'N',  /* Non-unit values in the diagonal */
			  'F'}; /* One based index */
     HYSL_FLOAT Alpha, Beta;    /* Constants for the BLAS routines */

     /* BLAS: tempvec = Disp */
     hysl_copy( &Tempvec->Rows, DispT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = a0*Disp = a0*tempvec */
     Alpha = a0;
     hysl_scal( &Tempvec->Rows, &Alpha, Tempvec->Array, &incx );
     /* BLAS: tempvec = a0*Disp + a2*Vel = tempvec + a2*Vel */
     Alpha = a2;
     hysl_axpy( &Tempvec->Rows, &Alpha, VelT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = a0*Disp + a2*Vel + a3*Acc = tempvec + a3*Acc */
     Alpha = a3;
     hysl_axpy( &Tempvec->Rows, &Alpha, AccT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: Eff_Force = Mass*(a0*Disp + a2*Vel + a3*Acc) = Mass*tempvec */
     Alpha = 1.0; Beta = 0.0;
     hysl_mkl_csrmv( &trans, &Tempvec->Rows, &Tempvec->Rows, &Alpha, matdescra, Mass->Values, Mass->Columns, Mass->RowIndex,
		     &Mass->RowIndex[1], Tempvec->Array, &Beta, Eff_ForceT->Array );

     /* BLAS: tempvec = Disp */
     hysl_copy( &Tempvec->Rows, DispT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = a1*Disp = a0*tempvec */
     Alpha = a1;
     hysl_scal( &Tempvec->Rows, &Alpha, Tempvec->Array, &incx );
     /* BLAS: tempvec = a1*Disp + a4*Vel = tempvec + a4*Vel */
     Alpha = a4;
     hysl_axpy( &Tempvec->Rows, &Alpha, VelT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = a1*Disp + a4*Vel + a5*Acc = tempvec + a5*Acc */
     Alpha = a5;
     hysl_axpy( &Tempvec->Rows, &Alpha, AccT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: Eff_Force = Mass*(a0*Disp + a2*Vel + a3*Acc) + Damp*(a1*Disp + a4*Vel + a5*Acc) = Eff_Force + Damp*tempvec */
     Alpha = 1.0; Beta = 1.0;
     hysl_mkl_csrmv( &trans, &Tempvec->Rows, &Tempvec->Rows, &Alpha, matdescra, Damp->Values, Damp->Columns, Damp->RowIndex,
		     &Damp->RowIndex[1], Tempvec->Array, &Beta, Eff_ForceT->Array );
     
}
