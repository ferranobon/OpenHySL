#include "EffM_Formulation.h"
#include "MatrixVector.h"
#include "Definitions.h"

#if _MKL_
#include <mkl_blas.h>
#else
#include "Netlib.h"
#endif

void EffM_EffectiveForce(const MatrixVector_t *const Damp, const MatrixVector_t *const Stiff, const MatrixVector_t *const DispT,
        const MatrixVector_t *const VelT, const MatrixVector_t *const AccT, MatrixVector_t *const Tempvec, const hysl_float_t a6, const hysl_float_t a9,
        const hysl_float_t a10, MatrixVector_t *const Eff_ForceT) {
    int32_t incx = 1;
    int32_t incy = 1; /* Stride in the vectors */

    /* BLAS: tempvec = Disp */
    hysl_copy(&Tempvec->Rows, DispT->Array, &incx, Tempvec->Array, &incy);
    /* BLAS: tempvec = Disp + a9*Vel = tempvec + a9*Vel */
    hysl_float_t Alpha = a9;
    hysl_axpy(&Tempvec->Rows, &Alpha, VelT->Array, &incx, Tempvec->Array, &incy);
    /* BLAS: tempvec = Disp + a9*Vel + a10*Acc = tempvec + a10*Acc */
    Alpha = a10;
    hysl_axpy(&Tempvec->Rows, &Alpha, AccT->Array, &incx, Tempvec->Array, &incy);
    /* BLAS: Eff_Force = -Stiff*(Disp + a2*Vel + a10*Acc) = -Stiff*tempvec */
    Alpha = -1.0;
    hysl_float_t Beta = 0.0;
    char uplo = 'L'; /* The lower part (upper part in C) will be used and the upper part (lower part in C) will strictly not be referenced */
    hysl_symv(&uplo, &Tempvec->Rows, &Alpha, Stiff->Array, &Tempvec->Rows, Tempvec->Array, &incx, &Beta, Eff_ForceT->Array, &incy);

    /* BLAS: tempvec = Vel */
    hysl_copy(&Tempvec->Rows, VelT->Array, &incx, Tempvec->Array, &incy);
    /* BLAS: tempvec = Vel + a6*Acc = tempvec + a6*Acc */
    Alpha = a6;
    hysl_axpy(&Tempvec->Rows, &Alpha, AccT->Array, &incx, Tempvec->Array, &incy);
    /* BLAS: Eff_Force = -Stiff*(Disp + a9*Vel + a10*Acc) - Damp*(Vel + 6*Acc) = Eff_Force - Damp*tempvec */
    Alpha = -1.0;
    Beta = 1.0;
    hysl_symv(&uplo, &Tempvec->Rows, &Alpha, Damp->Array, &Tempvec->Rows, Tempvec->Array, &incx, &Beta, Eff_ForceT->Array, &incy);
}

void EffM_EffectiveForce_PS(const MatrixVector_t *const Stiff, const MatrixVector_t *const Damp, const MatrixVector_t *const DispT,
        const MatrixVector_t *const VelT, const MatrixVector_t *const AccT, MatrixVector_t *const Tempvec, const hysl_float_t a6, const hysl_float_t a9,
        const hysl_float_t a10, MatrixVector_t *const Eff_ForceT) {
    int32_t incx = 1;
    int32_t incy = 1; /* Stride in the vectors */

    /* BLAS: tempvec = Disp */
    hysl_copy(&Tempvec->Rows, DispT->Array, &incx, Tempvec->Array, &incy);
    /* BLAS: tempvec = Disp + a9*Vel = tempvec + a9*Vel */
    hysl_float_t Alpha = a9;
    hysl_axpy(&Tempvec->Rows, &Alpha, VelT->Array, &incx, Tempvec->Array, &incy);
    /* BLAS: tempvec = Disp + a9*Vel + a10*Acc = tempvec + a10*Acc */
    Alpha = a10;
    hysl_axpy(&Tempvec->Rows, &Alpha, AccT->Array, &incx, Tempvec->Array, &incy);
    /* BLAS: Eff_Force = -Stiff*(Disp + a2*Vel + a10*Acc) = -Stiff*tempvec */
    Alpha = -1.0;
    hysl_float_t Beta = 0.0;
    char uplo = 'L'; /* The lower part (upper part in C) will be used and the upper part (lower part in C) will strictly not be referenced */
    hysl_spmv(&uplo, &Tempvec->Rows, &Alpha, Stiff->Array, Tempvec->Array, &incx, &Beta, Eff_ForceT->Array, &incy);

    /* BLAS: tempvec = Vel */
    hysl_copy(&Tempvec->Rows, VelT->Array, &incx, Tempvec->Array, &incy);
    /* BLAS: tempvec = Vel + a6*Acc = tempvec + a6*Acc */
    Alpha = a6;
    hysl_axpy(&Tempvec->Rows, &Alpha, AccT->Array, &incx, Tempvec->Array, &incy);
    /* BLAS: Eff_Force = -Stiff*(Disp + a9*Vel + a10*Acc) - Damp*(Vel + 6*Acc) = Eff_Force - Damp*tempvec */
    Alpha = -1.0;
    Beta = 1.0;
    hysl_spmv(&uplo, &Tempvec->Rows, &Alpha, Damp->Array, Tempvec->Array, &incx, &Beta, Eff_ForceT->Array, &incy);
}

void EffM_ComputeDisplacement(const MatrixVector_t *const DispT, const MatrixVector_t *const VelT, const MatrixVector_t *const AccT,
        const MatrixVector_t *const AccTdT, const hysl_float_t a8, const hysl_float_t a9, const hysl_float_t a10, MatrixVector_t *const DispTdT) {
    int32_t incx = 1;
    int32_t incy = 1; /* Stride in the vectors */

    /* BLAS: DispTdT = DispT */
    hysl_copy(&DispTdT->Rows, DispT->Array, &incx, DispTdT->Array, &incy);
    /* BLAS: DispTdT = DispT + a9*VelT = DispTdT + a9*VelT */
    hysl_float_t Alpha = a9;
    hysl_axpy(&DispTdT->Rows, &Alpha, VelT->Array, &incx, DispTdT->Array, &incy);
    /* BLAS: DispTdT = DispT + a9*VelT + a10*AccT = DispTdT + a10*AccT */
    Alpha = a10;
    hysl_axpy(&DispTdT->Rows, &Alpha, AccT->Array, &incx, DispTdT->Array, &incy);
    /* BLAS: DispTdT = DispT + a9*VelT + a10*AccT + a8*AccTdT = DispTdT + a8*AccTdT */
    Alpha = a8;
    hysl_axpy(&DispTdT->Rows, &Alpha, AccTdT->Array, &incx, DispTdT->Array, &incy);
}

void EffM_ComputeVelocity(const MatrixVector_t *const VelT, const MatrixVector_t *const AccT, const MatrixVector_t *const AccTdT, const hysl_float_t a6,
        const hysl_float_t a7, MatrixVector_t *const VelTdT) {
    int32_t incx = 1;
    int32_t incy = 1; /* Stride in the vectors */

    /* BLAS: VelTdT = VelT */
    hysl_copy(&VelTdT->Rows, VelT->Array, &incx, VelTdT->Array, &incy);
    /* BLAS: VelTdT = VelT + a6*AccT = VelTdT + a6*AccT */
    hysl_float_t Alpha = a6;
    hysl_axpy(&VelTdT->Rows, &Alpha, AccT->Array, &incx, VelTdT->Array, &incy);
    /* BLAS: VelTdT = VelT + a6*AccT + a7*AccTdT = VelTdT + a7*AccTdT */
    Alpha = a7;
    hysl_axpy(&VelTdT->Rows, &Alpha, AccTdT->Array, &incx, VelTdT->Array, &incy);
}

