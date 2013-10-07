#include "EffM_Formulation.h"
#include "MatrixVector_Sp.h"

#if _MKL_
#include <mkl_blas.h>
#include <mkl_spblas.h>
#else
#include "Netlib.h"
#endif

void EffM_EffectiveForce_Sp( const MatrixVector_Sp_t *const Stiff, const MatrixVector_Sp_t *const Damp, const MatrixVector_t *const DispT,
			     const MatrixVector_t *const VelT, const MatrixVector_t *const AccT, MatrixVector_t *const Tempvec,
			     const double a6, const double a9, const double a10,
			     MatrixVector_t *const Eff_ForceT )
{

     int incx = 1, incy = 1;    /* Stride in the vectors */
     char trans = 'N';          /* No transpose operation */
     char matdescra[6] = {'S',  /* The matrix is symmetric */
			  'U',  /* The upper part is referenced */
			  'N',  /* Non-unit values in the diagonal */
			  'F'}; /* One based index */
     double Alpha, Beta;        /* Constants for the BLAS routines */

     /* BLAS: tempvec = Disp */
     dcopy( &Tempvec->Rows, DispT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = Disp + a9*Vel = tempvec + a9*Vel */
     Alpha = a9;
     daxpy( &Tempvec->Rows, &Alpha, VelT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = Disp + a9*Vel + a10*Acc = tempvec + a10*Acc */
     Alpha = a10;
     daxpy( &Tempvec->Rows, &Alpha, AccT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: Eff_Force = -Stiff*(Disp + a2*Vel + a10*Acc) = -Stiff*tempvec */
     Alpha = -1.0; Beta = 0.0;
     mkl_dcsrmv( &trans, &Tempvec->Rows, &Tempvec->Rows, &Alpha, matdescra, Stiff->Values, Stiff->Columns, Stiff->RowIndex, &Stiff->RowIndex[1], Tempvec->Array, &Beta, Eff_ForceT->Array );

     /* BLAS: tempvec = Vel */
     dcopy( &Tempvec->Rows, VelT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = Vel + a6*Acc = tempvec + a6*Acc */
     Alpha = a6;
     daxpy( &Tempvec->Rows, &Alpha, AccT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: Eff_Force = -Stiff*(Disp + a9*Vel + a10*Acc) - Damp*(Vel + 6*Acc) = Eff_Force - Damp*tempvec */
     Alpha = -1.0; Beta = 1.0;
     mkl_dcsrmv( &trans, &Tempvec->Rows, &Tempvec->Rows, &Alpha, matdescra, Damp->Values, Damp->Columns, Damp->RowIndex, &Damp->RowIndex[1], Tempvec->Array, &Beta, Eff_ForceT->Array );
     
}
