#include "EffM_Formulation.h"
#include "MatrixVector_Sp.h"

#include "Definitions.h"

#include <mkl_blas.h>
#include <mkl_spblas.h>

void EffM_EffectiveForce_Sp (const MatrixVector_Sp_t *const Stiff, const MatrixVector_Sp_t *const Damp, const MatrixVector_t *const DispT, const MatrixVector_t *const VelT,
        const MatrixVector_t *const AccT, MatrixVector_t *const Tempvec, const hysl_float_t a6, const hysl_float_t a9, const hysl_float_t a10, MatrixVector_t *const Eff_ForceT) {

    int incx = 1, incy = 1; /* Stride in the vectors */
    char trans = 'N'; /* No transpose operation */
    char matdescra[6] = { 'S', /* The matrix is symmetric */
    'U', /* The upper part is referenced */
    'N', /* Non-unit values in the diagonal */
    'F' }; /* One based index */
    hysl_float_t Alpha, Beta; /* Constants for the BLAS routines */

    /* BLAS: tempvec = Disp */
    hysl_copy(&Tempvec->Rows, DispT->Array, &incx, Tempvec->Array, &incy);
    /* BLAS: tempvec = Disp + a9*Vel = tempvec + a9*Vel */
    Alpha = a9;
    hysl_axpy(&Tempvec->Rows, &Alpha, VelT->Array, &incx, Tempvec->Array, &incy);
    /* BLAS: tempvec = Disp + a9*Vel + a10*Acc = tempvec + a10*Acc */
    Alpha = a10;
    hysl_axpy(&Tempvec->Rows, &Alpha, AccT->Array, &incx, Tempvec->Array, &incy);
    /* BLAS: Eff_Force = -Stiff*(Disp + a2*Vel + a10*Acc) = -Stiff*tempvec */
    Alpha = -1.0;
    Beta = 0.0;
    hysl_mkl_csrmv(&trans, &Tempvec->Rows, &Tempvec->Rows, &Alpha, matdescra, Stiff->Values, Stiff->Columns, Stiff->RowIndex, &Stiff->RowIndex[1], Tempvec->Array, &Beta, Eff_ForceT->Array);

    /* BLAS: tempvec = Vel */
    hysl_copy(&Tempvec->Rows, VelT->Array, &incx, Tempvec->Array, &incy);
    /* BLAS: tempvec = Vel + a6*Acc = tempvec + a6*Acc */
    Alpha = a6;
    hysl_axpy(&Tempvec->Rows, &Alpha, AccT->Array, &incx, Tempvec->Array, &incy);
    /* BLAS: Eff_Force = -Stiff*(Disp + a9*Vel + a10*Acc) - Damp*(Vel + 6*Acc) = Eff_Force - Damp*tempvec */
    Alpha = -1.0;
    Beta = 1.0;
    hysl_mkl_csrmv(&trans, &Tempvec->Rows, &Tempvec->Rows, &Alpha, matdescra, Damp->Values, Damp->Columns, Damp->RowIndex, &Damp->RowIndex[1], Tempvec->Array, &Beta, Eff_ForceT->Array);

}
