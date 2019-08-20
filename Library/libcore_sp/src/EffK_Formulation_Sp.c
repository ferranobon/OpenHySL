#include "EffK_Formulation.h"
#include "MatrixVector_Sp.h"

#include "Definitions.h"

#include <mkl_blas.h>
#include <mkl_spblas.h>

void EffK_EffectiveForce_Sp (const MatrixVector_Sp_t *const Mass, const MatrixVector_Sp_t *const Damp, const MatrixVector_t *const DispT, const MatrixVector_t *const VelT,
        const MatrixVector_t *const AccT, MatrixVector_t *const Tempvec, const hysl_float_t a0, const hysl_float_t a1, const hysl_float_t a2, const hysl_float_t a3, const hysl_float_t a4,
        const hysl_float_t a5, MatrixVector_t *const Eff_ForceT) {

    int32_t incx = 1;
    int32_t incy = 1; /* Stride in the vectors */
    char trans = 'N'; /* No transpose operation */
    char matdescra[6] = { 'S', /* The matrix is symmetric */
    'U', /* The upper part is referenced */
    'N', /* Non-unit values in the diagonal */
    'F' }; /* One based index */
    hysl_float_t Alpha, Beta; /* Constants for the BLAS routines */

    /* BLAS: tempvec = Disp */
    hysl_copy(&Tempvec->Rows, DispT->Array, &incx, Tempvec->Array, &incy);
    /* BLAS: tempvec = a0*Disp = a0*tempvec */
    Alpha = a0;
    hysl_scal(&Tempvec->Rows, &Alpha, Tempvec->Array, &incx);
    /* BLAS: tempvec = a0*Disp + a2*Vel = tempvec + a2*Vel */
    Alpha = a2;
    hysl_axpy(&Tempvec->Rows, &Alpha, VelT->Array, &incx, Tempvec->Array, &incy);
    /* BLAS: tempvec = a0*Disp + a2*Vel + a3*Acc = tempvec + a3*Acc */
    Alpha = a3;
    hysl_axpy(&Tempvec->Rows, &Alpha, AccT->Array, &incx, Tempvec->Array, &incy);
    /* BLAS: Eff_Force = Mass*(a0*Disp + a2*Vel + a3*Acc) = Mass*tempvec */
    Alpha = 1.0;
    Beta = 0.0;
    hysl_mkl_csrmv(&trans, &Tempvec->Rows, &Tempvec->Rows, &Alpha, matdescra, Mass->Values, Mass->Columns, Mass->RowIndex, &Mass->RowIndex[1], Tempvec->Array, &Beta, Eff_ForceT->Array);

    /* BLAS: tempvec = Disp */
    hysl_copy(&Tempvec->Rows, DispT->Array, &incx, Tempvec->Array, &incy);
    /* BLAS: tempvec = a1*Disp = a0*tempvec */
    Alpha = a1;
    hysl_scal(&Tempvec->Rows, &Alpha, Tempvec->Array, &incx);
    /* BLAS: tempvec = a1*Disp + a4*Vel = tempvec + a4*Vel */
    Alpha = a4;
    hysl_axpy(&Tempvec->Rows, &Alpha, VelT->Array, &incx, Tempvec->Array, &incy);
    /* BLAS: tempvec = a1*Disp + a4*Vel + a5*Acc = tempvec + a5*Acc */
    Alpha = a5;
    hysl_axpy(&Tempvec->Rows, &Alpha, AccT->Array, &incx, Tempvec->Array, &incy);
    /* BLAS: Eff_Force = Mass*(a0*Disp + a2*Vel + a3*Acc) + Damp*(a1*Disp + a4*Vel + a5*Acc) = Eff_Force + Damp*tempvec */
    Alpha = 1.0;
    Beta = 1.0;
    hysl_mkl_csrmv(&trans, &Tempvec->Rows, &Tempvec->Rows, &Alpha, matdescra, Damp->Values, Damp->Columns, Damp->RowIndex, &Damp->RowIndex[1], Tempvec->Array, &Beta, Eff_ForceT->Array);

}

void EffK_EffectiveForce_HHT_Sp (const MatrixVector_Sp_t *const Mass, const MatrixVector_Sp_t *const Damp, const MatrixVector_Sp_t *const Stiff, const MatrixVector_t *const DispT,
        const MatrixVector_t *const VelT, const MatrixVector_t *const AccT, MatrixVector_t *const Tempvec, const hysl_float_t a0, const hysl_float_t a1, const hysl_float_t a2, const hysl_float_t a3,
        const hysl_float_t a4, const hysl_float_t a5, const hysl_float_t Alpha_HHT, MatrixVector_t *const Eff_ForceT) {

    int incx = 1, incy = 1; /* Stride in the vectors */
    char trans = 'N'; /* No transpose operation */
    char matdescra[6] = { 'S', /* The matrix is symmetric */
    'U', /* The upper part is referenced */
    'N', /* Non-unit values in the diagonal */
    'F' }; /* One based index */
    hysl_float_t Alpha, Beta; /* Constants for the BLAS routines */

    /* BLAS: tempvec = Disp */
    hysl_copy(&Tempvec->Rows, DispT->Array, &incx, Tempvec->Array, &incy);
    /* BLAS: tempvec = a0*Disp = a0*tempvec */
    Alpha = a0;
    hysl_scal(&Tempvec->Rows, &Alpha, Tempvec->Array, &incx);
    /* BLAS: tempvec = a0*Disp + a2*Vel = tempvec + a2*Vel */
    Alpha = a2;
    hysl_axpy(&Tempvec->Rows, &Alpha, VelT->Array, &incx, Tempvec->Array, &incy);
    /* BLAS: tempvec = a0*Disp + a2*Vel + a3*Acc = tempvec + a3*Acc */
    Alpha = a3;
    hysl_axpy(&Tempvec->Rows, &Alpha, AccT->Array, &incx, Tempvec->Array, &incy);
    /* BLAS: Eff_Force = Mass*(a0*Disp + a2*Vel + a3*Acc) = Mass*tempvec */
    Alpha = 1.0;
    Beta = 0.0;
    hysl_mkl_csrmv(&trans, &Tempvec->Rows, &Tempvec->Rows, &Alpha, matdescra, Mass->Values, Mass->Columns, Mass->RowIndex, &Mass->RowIndex[1], Tempvec->Array, &Beta, Eff_ForceT->Array);

    /* BLAS: tempvec = Disp */
    hysl_copy(&Tempvec->Rows, DispT->Array, &incx, Tempvec->Array, &incy);
    /* BLAS: tempvec = a1*Disp = a1*tempvec */
    Alpha = a1;
    hysl_scal(&Tempvec->Rows, &Alpha, Tempvec->Array, &incx);
    /* BLAS: tempvec = a1*Disp + a4*Vel = tempvec + a4*Vel */
    Alpha = a4;
    hysl_axpy(&Tempvec->Rows, &Alpha, VelT->Array, &incx, Tempvec->Array, &incy);
    /* BLAS: tempvec = a1*Disp + a4*Vel + a5*Acc = tempvec + a5*Acc */
    Alpha = a5;
    hysl_axpy(&Tempvec->Rows, &Alpha, AccT->Array, &incx, Tempvec->Array, &incy);
    /* BLAS: Eff_Force = Mass*(a0*Disp + a2*Vel + a3*Acc) + (1 + Alpha_HHT)*Damp*(a1*Disp + a4*Vel + a5*Acc)
     * = Eff_Force + (1 + Alpha_HHT)*Damp*tempvec */
    Alpha = 1.0 + Alpha_HHT;
    Beta = 1.0;
    hysl_mkl_csrmv(&trans, &Tempvec->Rows, &Tempvec->Rows, &Alpha, matdescra, Damp->Values, Damp->Columns, Damp->RowIndex, &Damp->RowIndex[1], Tempvec->Array, &Beta, Eff_ForceT->Array);

    /* BLAS: Eff_Force = Mass*(a0*Disp + a2*Vel + a3*Acc) + (1 + Alpha_HHT)*Damp*(a1*Disp + a4*Vel + a5*Acc)
     * + Alpha_HHT*Damp*Vel = Eff_Force + Alpha_HHT*Damp*Vel */
    Alpha = Alpha_HHT;
    hysl_mkl_csrmv(&trans, &Tempvec->Rows, &Tempvec->Rows, &Alpha, matdescra, Damp->Values, Damp->Columns, Damp->RowIndex, &Damp->RowIndex[1], VelT->Array, &Beta, Eff_ForceT->Array);

    /* BLAS: Eff_Force = Mass*(a0*Disp + a2*Vel + a3*Acc) + (1 + Alpha_HHT)*Damp*(a1*Disp + a4*Vel + a5*Acc)
     * + Alpha_HHT*Damp*Vel + Alpha_HHT*Stiff*Disp= Eff_Force + Alpha_HHT*Stiff*Disp */
    hysl_mkl_csrmv(&trans, &Tempvec->Rows, &Tempvec->Rows, &Alpha, matdescra, Stiff->Values, Stiff->Columns, Stiff->RowIndex, &Stiff->RowIndex[1], DispT->Array, &Beta, Eff_ForceT->Array);
}
