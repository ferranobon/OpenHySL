#include "EffK_Formulation.h"
#include "MatrixVector.h"

#include "Definitions.h"

#if _MKL_
#include <mkl_blas.h>
#else
#include "Netlib.h"
#endif

void EffK_EffectiveForce(const MatrixVector_t *const Mass, const MatrixVector_t *const Damp, const MatrixVector_t *const DispT,
        const MatrixVector_t *const VelT, const MatrixVector_t *const AccT, MatrixVector_t *const Tempvec, const hysl_float_t a0, const hysl_float_t a1,
        const hysl_float_t a2, const hysl_float_t a3, const hysl_float_t a4, const hysl_float_t a5, MatrixVector_t *const Eff_ForceT) {

    int32_t incx = 1;
    int32_t incy = 1; /* Stride in the vectors */

    /* BLAS: tempvec = Disp */
    hysl_copy(&Tempvec->Rows, DispT->Array, &incx, Tempvec->Array, &incy);
    /* BLAS: tempvec = a0*Disp = a0*tempvec */
    hysl_float_t Alpha = a0;
    hysl_scal(&Tempvec->Rows, &Alpha, Tempvec->Array, &incx);
    /* BLAS: tempvec = a0*Disp + a2*Vel = tempvec + a2*Vel */
    Alpha = a2;
    hysl_axpy(&Tempvec->Rows, &Alpha, VelT->Array, &incx, Tempvec->Array, &incy);
    /* BLAS: tempvec = a0*Disp + a2*Vel + a3*Acc = tempvec + a3*Acc */
    Alpha = a3;
    hysl_axpy(&Tempvec->Rows, &Alpha, AccT->Array, &incx, Tempvec->Array, &incy);
    /* BLAS: Eff_Force = Mass*(a0*Disp + a2*Vel + a3*Acc) = Mass*tempvec */
    Alpha = 1.0;
    hysl_float_t Beta = 0.0;
    char uplo = 'L';  /* The lower part (upper part in C) will be used and the upper part (lower part  in C) will strictly not be referenced */
    hysl_symv(&uplo, &Tempvec->Rows, &Alpha, Mass->Array, &Tempvec->Rows, Tempvec->Array, &incx, &Beta, Eff_ForceT->Array, &incy);

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
    /* BLAS: Eff_Force = Mass*(a0*Disp + a2*Vel + a3*Acc) + Damp*(a1*Disp + a4*Vel + a5*Acc) = Eff_Force +
     * Damp*tempvec */
    Alpha = 1.0;
    Beta = 1.0;
    hysl_symv(&uplo, &Tempvec->Rows, &Alpha, Damp->Array, &Tempvec->Rows, Tempvec->Array, &incx, &Beta, Eff_ForceT->Array, &incy);
}

void EffK_EffectiveForce_PS(const MatrixVector_t *const Mass, const MatrixVector_t *const Damp, const MatrixVector_t *const DispT,
        const MatrixVector_t *const VelT, const MatrixVector_t *const AccT, MatrixVector_t *const Tempvec, const hysl_float_t a0, const hysl_float_t a1,
        const hysl_float_t a2, const hysl_float_t a3, const hysl_float_t a4, const hysl_float_t a5, MatrixVector_t *const Eff_ForceT) {

    int32_t incx = 1;
    int32_t incy = 1; /* Stride in the vectors */

    /* BLAS: tempvec = Disp */
    hysl_copy(&Tempvec->Rows, DispT->Array, &incx, Tempvec->Array, &incy);
    /* BLAS: tempvec = a0*Disp = a0*tempvec */
    hysl_float_t Alpha = a0;
    hysl_scal(&Tempvec->Rows, &Alpha, Tempvec->Array, &incx);
    /* BLAS: tempvec = a0*Disp + a2*Vel = tempvec + a2*Vel */
    Alpha = a2;
    hysl_axpy(&Tempvec->Rows, &Alpha, VelT->Array, &incx, Tempvec->Array, &incy);
    /* BLAS: tempvec = a0*Disp + a2*Vel + a3*Acc = tempvec + a3*Acc */
    Alpha = a3;
    hysl_axpy(&Tempvec->Rows, &Alpha, AccT->Array, &incx, Tempvec->Array, &incy);
    /* BLAS: Eff_Force = Mass*(a0*Disp + a2*Vel + a3*Acc) = Mass*tempvec */
    Alpha = 1.0;
    hysl_float_t Beta = 0.0;
    char uplo = 'L'; /* The lower part (upper part in C) will be used and the upper part (lower part in C) will strictly not be referenced */
    hysl_spmv(&uplo, &Tempvec->Rows, &Alpha, Mass->Array, Tempvec->Array, &incx, &Beta, Eff_ForceT->Array, &incy);

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
    /* BLAS: Eff_Force = Mass*(a0*Disp + a2*Vel + a3*Acc) + Damp*(a1*Disp + a4*Vel + a5*Acc) = Eff_Force +
     * Damp*tempvec */
    Alpha = 1.0;
    Beta = 1.0;
    hysl_spmv(&uplo, &Tempvec->Rows, &Alpha, Damp->Array, Tempvec->Array, &incx, &Beta, Eff_ForceT->Array, &incy);

}

void EffK_EffectiveForce_HHT(const MatrixVector_t *const Mass, const MatrixVector_t *const Damp, const MatrixVector_t *const Stiff,
        const MatrixVector_t *const DispT, const MatrixVector_t *const VelT, const MatrixVector_t *const AccT, MatrixVector_t *const Tempvec,
        const hysl_float_t a0, const hysl_float_t a1, const hysl_float_t a2, const hysl_float_t a3, const hysl_float_t a4, const hysl_float_t a5,
        const hysl_float_t Alpha_HHT, MatrixVector_t *const Eff_ForceT) {

    int32_t incx = 1;
    int32_t incy = 1; /* Stride in the vectors */

    /* BLAS: tempvec = Disp */
    hysl_copy(&Tempvec->Rows, DispT->Array, &incx, Tempvec->Array, &incy);
    /* BLAS: tempvec = a0*Disp = a0*tempvec */
    hysl_float_t Alpha = a0;
    hysl_scal(&Tempvec->Rows, &Alpha, Tempvec->Array, &incx);
    /* BLAS: tempvec = a0*Disp + a2*Vel = tempvec + a2*Vel */
    Alpha = a2;
    hysl_axpy(&Tempvec->Rows, &Alpha, VelT->Array, &incx, Tempvec->Array, &incy);
    /* BLAS: tempvec = a0*Disp + a2*Vel + a3*Acc = tempvec + a3*Acc */
    Alpha = a3;
    hysl_axpy(&Tempvec->Rows, &Alpha, AccT->Array, &incx, Tempvec->Array, &incy);
    /* BLAS: Eff_Force = Mass*(a0*Disp + a2*Vel + a3*Acc) = Mass*tempvec */
    Alpha = 1.0;
    hysl_float_t Beta = 0.0;
    char uplo = 'L';  /* The lower part (upper part in C) will be used and the upper part (lower part in C) will strictly not be referenced */
    hysl_symv(&uplo, &Tempvec->Rows, &Alpha, Mass->Array, &Tempvec->Rows, Tempvec->Array, &incx, &Beta, Eff_ForceT->Array, &incy);

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
    hysl_symv(&uplo, &Tempvec->Rows, &Alpha, Damp->Array, &Tempvec->Rows, Tempvec->Array, &incx, &Beta, Eff_ForceT->Array, &incy);

    /* BLAS: Eff_Force = Mass*(a0*Disp + a2*Vel + a3*Acc) + (1 + Alpha_HHT)*Damp*(a1*Disp + a4*Vel + a5*Acc)
     * + Alpha_HHT*Damp*Vel = Eff_Force + Alpha_HHT*Damp*Vel */
    Alpha = Alpha_HHT;
    hysl_symv(&uplo, &Tempvec->Rows, &Alpha, Damp->Array, &Tempvec->Rows, VelT->Array, &incx, &Beta, Eff_ForceT->Array, &incy);

    /* BLAS: Eff_Force = Mass*(a0*Disp + a2*Vel + a3*Acc) + (1 + Alpha_HHT)*Damp*(a1*Disp + a4*Vel + a5*Acc)
     * + Alpha_HHT*Damp*Vel + Alpha_HHT*Stiff*Disp= Eff_Force + Alpha_HHT*Stiff*Disp */
    hysl_symv(&uplo, &Tempvec->Rows, &Alpha, Stiff->Array, &Tempvec->Rows, DispT->Array, &incx, &Beta, Eff_ForceT->Array, &incy);
}

void EffK_EffectiveForce_HHT_PS(const MatrixVector_t *const Mass, const MatrixVector_t *const Damp, const MatrixVector_t *const Stiff,
        const MatrixVector_t *const DispT, const MatrixVector_t *const VelT, const MatrixVector_t *const AccT, MatrixVector_t *const Tempvec,
        const hysl_float_t a0, const hysl_float_t a1, const hysl_float_t a2, const hysl_float_t a3, const hysl_float_t a4, const hysl_float_t a5,
        const hysl_float_t Alpha_HHT, MatrixVector_t *const Eff_ForceT) {

    int32_t incx = 1;
    int32_t incy = 1; /* Stride in the vectors */

    /* BLAS: tempvec = Disp */
    hysl_copy(&Tempvec->Rows, DispT->Array, &incx, Tempvec->Array, &incy);
    /* BLAS: tempvec = a0*Disp = a0*tempvec */
    hysl_float_t Alpha = a0;
    hysl_scal(&Tempvec->Rows, &Alpha, Tempvec->Array, &incx);
    /* BLAS: tempvec = a0*Disp + a2*Vel = tempvec + a2*Vel */
    Alpha = a2;
    hysl_axpy(&Tempvec->Rows, &Alpha, VelT->Array, &incx, Tempvec->Array, &incy);
    /* BLAS: tempvec = a0*Disp + a2*Vel + a3*Acc = tempvec + a3*Acc */
    Alpha = a3;
    hysl_axpy(&Tempvec->Rows, &Alpha, AccT->Array, &incx, Tempvec->Array, &incy);
    /* BLAS: Eff_Force = Mass*(a0*Disp + a2*Vel + a3*Acc) = Mass*tempvec */
    Alpha = 1.0;
    hysl_float_t Beta = 0.0;
    char uplo = 'L'; /* The lower part (upper part in C) will be used and the upper part (lower part n C) will strictly not be referenced */
    hysl_spmv(&uplo, &Tempvec->Rows, &Alpha, Mass->Array, Tempvec->Array, &incx, &Beta, Eff_ForceT->Array, &incy);

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
    /* BLAS: Eff_Force = Mass*(a0*Disp + a2*Vel + a3*Acc) + (1 - Alpha_HHT)*Damp*(a1*Disp + a4*Vel + a5*Acc)
     * = Eff_Force + (1 - Alpha_HHT)*Damp*tempvec */
    Alpha = 1.0 + Alpha_HHT;
    Beta = 1.0;
    hysl_spmv(&uplo, &Tempvec->Rows, &Alpha, Damp->Array, Tempvec->Array, &incx, &Beta, Eff_ForceT->Array, &incy);

    /* BLAS: Eff_Force = Mass*(a0*Disp + a2*Vel + a3*Acc) + (1 - Alpha_HHT)*Damp*(a1*Disp + a4*Vel + a5*Acc)
     * + Alpha_HHT*Damp*Vel = Eff_Force + Alpha_HHT*Damp*Vel */
    Alpha = Alpha_HHT;
    hysl_spmv(&uplo, &Tempvec->Rows, &Alpha, Damp->Array, VelT->Array, &incx, &Beta, Eff_ForceT->Array, &incy);

    /* BLAS: Eff_Force = Mass*(a0*Disp + a2*Vel + a3*Acc) + (1 - Alpha_HHT)*Damp*(a1*Disp + a4*Vel + a5*Acc)
     * + Alpha_HHT*Damp*Vel + Alpha_HHT*Stiff*Disp= Eff_Force + Alpha_HHT*Stiff*Disp */
    hysl_spmv(&uplo, &Tempvec->Rows, &Alpha, Stiff->Array, DispT->Array, &incx, &Beta, Eff_ForceT->Array, &incy);

}

void EffK_ComputeAcceleration(const MatrixVector_t *const DispTdT, const MatrixVector_t *const DispT, const MatrixVector_t *const VelT,
        const MatrixVector_t *const AccT, const hysl_float_t a0, const hysl_float_t a2, const hysl_float_t a3, MatrixVector_t *const AccTdT) {
    int32_t incx = 1;
    int32_t incy = 1; /* Stride in the vectors */

    /* BLAS: AccTdT = DispTdT */
    hysl_copy(&AccTdT->Rows, DispTdT->Array, &incx, AccTdT->Array, &incy);
    /* BLAS: AccTdT = DispTdT - DispT = AccTdT - DispT */
    hysl_float_t Alpha = -1.0;
    hysl_axpy(&AccTdT->Rows, &Alpha, DispT->Array, &incx, AccTdT->Array, &incy);
    /* BLAS: AccTdT = a0*(DispTdT - DispT) = a0*AccTdT */
    Alpha = a0;
    hysl_scal(&AccTdT->Rows, &Alpha, AccTdT->Array, &incx);
    /* BLAS: AccTdT = a0*(DispTdT - DispT) - a2*VelT = AccTdT - a2*VelT */
    Alpha = -a2;
    hysl_axpy(&AccTdT->Rows, &Alpha, VelT->Array, &incx, AccTdT->Array, &incy);
    /* BLAS: AccTdT = a0*(DispTdT - DispT) - a2*VelT - a3*AccT = AccTdT - a3*AccT */
    Alpha = -a3;
    hysl_axpy(&AccTdT->Rows, &Alpha, AccT->Array, &incx, AccTdT->Array, &incy);
}

void EffK_ComputeVelocity(const MatrixVector_t *const DispTdT, const MatrixVector_t *const DispT, const MatrixVector_t *const VelT,
        const MatrixVector_t *const AccT, const hysl_float_t a1, const hysl_float_t a4, const hysl_float_t a5, MatrixVector_t *const VelTdT) {
    int32_t incx = 1;
    int32_t incy = 1; /* Stride in the vectors */

    /* BLAS: VelTdT = DispTdT */
    hysl_copy(&VelTdT->Rows, DispTdT->Array, &incx, VelTdT->Array, &incy);
    /* BLAS: VelTdT = DispTdT - DispT = VeldT - DispT */
    hysl_float_t Alpha = -1.0;
    hysl_axpy(&VelTdT->Rows, &Alpha, DispT->Array, &incx, VelTdT->Array, &incy);
    /* BLAS: AccTdT = a1*(DispTdT - DispT) = a1*AccTdT */
    Alpha = a1;
    hysl_scal(&VelTdT->Rows, &Alpha, VelTdT->Array, &incx);
    /* BLAS: VelTdT = a1*(DispTdT - DispT) - a4*VelT = VelTdT - a4*VelT */
    Alpha = -a4;
    hysl_axpy(&VelTdT->Rows, &Alpha, VelT->Array, &incx, VelTdT->Array, &incy);
    /* BLAS: VelTdT = a1*(DispTdT - DispT) - a4*VelT - a5*AccT = VelTdT - a5*AccT */
    Alpha = -a5;
    hysl_axpy(&VelTdT->Rows, &Alpha, AccT->Array, &incx, VelTdT->Array, &incy);
}
