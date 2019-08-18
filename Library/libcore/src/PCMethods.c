#include "MatrixVector.h"
#include "PCMethods.h"
#include "Definitions.h"

#if _MKL_
#include <mkl_blas.h>
#else
#include "Netlib.h"
#endif

void PC_PredictorStep_Displacement(const MatrixVector_t *const DispT, const MatrixVector_t *const VelT, const MatrixVector_t *const AccT, const hysl_float_t a9,
        const hysl_float_t a10, MatrixVector_t *const DispTdT_Pred) {

    int incx = 1, incy = 1; /* Stride in the vectors */
    hysl_float_t Alpha; /* Constant for the BLAS routines */

    /* BLAS: DispTdT_Pred = DispT */
    hysl_copy(&DispTdT_Pred->Rows, DispT->Array, &incx, DispTdT_Pred->Array, &incy);
    /* BLAS: DispTdT_Pred = DispT + a9*VelT */
    Alpha = a9;
    hysl_axpy(&DispTdT_Pred->Rows, &Alpha, VelT->Array, &incx, DispTdT_Pred->Array, &incy);
    /* BLAS: DispTdT_Pred = DispT + a9*VelT + a10*AccT */
    Alpha = a10;
    hysl_axpy(&DispTdT_Pred->Rows, &Alpha, AccT->Array, &incx, DispTdT_Pred->Array, &incy);
}

void PC_PredictorStep_Velocity(const MatrixVector_t *const VelT, const MatrixVector_t *const AccT, const hysl_float_t a6, MatrixVector_t *const VelTdT_Pred) {
    int incx = 1, incy = 1; /* Stride in the vectors */
    hysl_float_t Alpha; /* Constant for the BLAS routines */

    /* BLAS: VelTdT_Pred = VelT */
    hysl_copy(&VelTdT_Pred->Rows, VelT->Array, &incx, VelTdT_Pred->Array, &incy);
    /* BLAS: VelTdT_Pred = VelT + a6*AccT */
    Alpha = a6;
    hysl_axpy(&VelTdT_Pred->Rows, &Alpha, AccT->Array, &incx, VelTdT_Pred->Array, &incy);

}

void PC_CorrectorStep_Displacement(const MatrixVector_t *const DispTdT_Pred, const MatrixVector_t *const AccTdT, const hysl_float_t a8,
        MatrixVector_t *const DispTdT) {
    int incx = 1, incy = 1; /* Stride in the vectors */
    hysl_float_t Alpha; /* Constant for the BLAS routines */

    /* BLAS: DispTdT = DispTdT_Pred */
    hysl_copy(&DispTdT->Rows, DispTdT_Pred->Array, &incx, DispTdT->Array, &incy);
    /* BLAS: DispTdT = DispTdT_Pred + a8*AccTdT */
    Alpha = a8;
    hysl_axpy(&DispTdT->Rows, &Alpha, AccTdT->Array, &incx, DispTdT->Array, &incy);
}

void PC_CorrectorStep_Velocity(const MatrixVector_t *const VelTdT_Pred, const MatrixVector_t *const AccTdT, const hysl_float_t a7, MatrixVector_t *const VelTdT) {
    int incx = 1, incy = 1; /* Stride in the vectors */
    hysl_float_t Alpha; /* Constant for the BLAS routines */

    /* BLAS: VelTdT = VelTdT_Pred */
    hysl_copy(&VelTdT->Rows, VelTdT_Pred->Array, &incx, VelTdT->Array, &incy);
    /* BLAS: VelTdT = VelTdT_Pred + a7*AccTdT */
    Alpha = a7;
    hysl_axpy(&VelTdT->Rows, &Alpha, AccTdT->Array, &incx, VelTdT->Array, &incy);
}

void PC_ReactionForces_Numerical(const MatrixVector_t *const DispTdT_Pred, MatrixVector_t *const K, MatrixVector_t *const RForceTdT) {
    int incx = 1, incy = 1; /* Stride in the vectors */
    hysl_float_t Alpha = 1.0, Beta = 0.0; /* Constant for the BLAS routines */
    char uplo = 'L'; /* The lower part (upper part in C) will be used and the upper part
     * (lower part in C) will strictly not be referenced */

    /* BLAS: RForceTdT = K*DispTdT_Pred */
    hysl_symv(&uplo, &RForceTdT->Rows, &Alpha, K->Array, &RForceTdT->Rows, DispTdT_Pred->Array, &incx, &Beta, RForceTdT->Array, &incy);
}

void PC_ReactionForces_Numerical_PS(const MatrixVector_t *const DispTdT_Pred, MatrixVector_t *const K, MatrixVector_t *const RForceTdT) {
    int incx = 1, incy = 1; /* Stride in the vectors */
    hysl_float_t Alpha = 1.0, Beta = 0.0; /* Constant for the BLAS routines */
    char uplo = 'L'; /* The lower part (upper part in C) will be used and the upper part
     * (lower part in C) will strictly not be referenced */

    /* BLAS: RForceTdT = K*DispTdT_Pred */
    hysl_spmv(&uplo, &RForceTdT->Rows, &Alpha, K->Array, DispTdT_Pred->Array, &incx, &Beta, RForceTdT->Array, &incy);

}

void PC_Calculate_Acceleration(const MatrixVector_t *const LoadTdT, const MatrixVector_t *const LoadT, const MatrixVector_t *const RForceTdT,
        const MatrixVector_t *const RForceT, const MatrixVector_t *const VelTdT_Pred, const MatrixVector_t *const VelT_Pred, const MatrixVector_t *const AccT,
        const MatrixVector_t *const K, const MatrixVector_t *const C, const MatrixVector_t *const Meinv, const hysl_float_t alpha_H, const hysl_float_t a7,
        const hysl_float_t a8, MatrixVector_t *const AccTdT) {
    int incx = 1, incy = 1; /* Stride in the vectors */
    hysl_float_t Alpha, Beta; /* Constant for the BLAS routines */
    char uplo = 'L';        /* The lower part (upper part in C) will be used and the upper part
                             * (lower part in C) will strictly not be referenced */

    /* BLAS: AccTdT = LoadTdT */
    hysl_copy(&AccTdT->Rows, LoadTdT->Array, &incx, AccTdT->Array, &incy);
    /* BLAS: AccTdT = (1 + alpha_H)*LoadTdT */
    Alpha = 1.0 + alpha_H;
    hysl_scal(&AccTdT->Rows, &Alpha, AccTdT->Array, &incx);
    /* BLAS: AccTdT = (1 + alpha_H)*LoadTdT - alpha_H*LoadT */
    Alpha = -alpha_H;
    hysl_axpy(&AccTdT->Rows, &Alpha, LoadT->Array, &incx, AccTdT->Array, &incy);
    /* BLAS: AccTdT = (1 + alpha_H)*LoadTdT - alpha_H*LoadT + alpha_H*RForceT */
    Alpha = alpha_H;
    hysl_axpy(&AccTdT->Rows, &Alpha, RForceT->Array, &incx, AccTdT->Array, &incy);
    /* BLAS: AccTdT = (1 + alpha_H)*LoadTdT - alpha_H*LoadT + alpha_H*RForceT - (1 + alpha_H)*RForceTdT*/
    Alpha = -(1.0 + alpha_H);
    hysl_axpy(&AccTdT->Rows, &Alpha, RForceTdT->Array, &incx, AccTdT->Array, &incy);
    /* BLAS: AccTdT = (1 + alpha_H)*LoadTdT - alpha_H*LoadT + alpha_H*RForceT - (1 + alpha_H)*RForceTdT +
     * alpha_H*C*VelT_Pred */
    Alpha = alpha_H;
    Beta = 1.0;
    hysl_symv(&uplo, &AccTdT->Rows, &Alpha, C->Array, &AccTdT->Rows, VelT_Pred->Array, &incx, &Beta, AccTdT->Array, &incy);
    /* BLAS: AccTdT = (1 + alpha_H)*LoadTdT - alpha_H*LoadT + alpha_H*RForceT - (1 + alpha_H)*RForceTdT +
     * alpha_H*C*VelT_Pred - (1 + alpha_H)*c*VelTdT_Pred */
    Alpha = -(1.0 + alpha_H);
    hysl_symv(&uplo, &AccTdT->Rows, &Alpha, C->Array, &AccTdT->Rows, VelTdT_Pred->Array, &incx, &Beta, AccTdT->Array, &incy);
    /* BLAS: AccTdT = (1 + alpha_H)*LoadTdT - alpha_H*LoadT + alpha_H*RForceT - (1 + alpha_H)*RForceTdT +
     * alpha_H*C*VelT_Pred - (1 + alpha_H)*c*VelTdT_Pred + alpha_H*a7*C*AccT. */
    Alpha = alpha_H * a7;
    hysl_symv(&uplo, &AccTdT->Rows, &Alpha, C->Array, &AccTdT->Rows, AccT->Array, &incx, &Beta, AccTdT->Array, &incy);
    /* BLAS: AccTdT = (1 + alpha_H)*LoadTdT - alpha_H*LoadT + alpha_H*RForceT - (1 + alpha_H)*RForceTdT +
     * alpha_H*C*VelT_Pred - (1 + alpha_H)*c*VelTdT_Pred + alpha_H*a7*C*AccT + alpha_H*a8*K*AccT */
    Alpha = alpha_H * a8;
    hysl_symv(&uplo, &AccTdT->Rows, &Alpha, K->Array, &AccTdT->Rows, AccT->Array, &incx, &Beta, AccTdT->Array, &incy);
    /* BLAS: AccTdT = Meinv*((1 + alpha_H)*LoadTdT - alpha_H*LoadT + alpha_H*RForceT - (1 + alpha_H)*RForceTdT +
     * alpha_H*C*VelT_Pred + (1 + alpha_H)*c*VelTdT_Pred + alpha_H*a7*C*AccT + alpha_H*a8*K*AccT) */
    Alpha = 1.0;
    Beta = 0.0;

    MatrixVector_t Tempvec;
    MatrixVector_Create(AccTdT->Rows, AccTdT->Cols, &Tempvec);

    hysl_copy(&AccTdT->Rows, AccTdT->Array, &incx, Tempvec.Array, &incy);
    hysl_symv(&uplo, &AccTdT->Rows, &Alpha, Meinv->Array, &AccTdT->Rows, Tempvec.Array, &incx, &Beta, AccTdT->Array, &incy);

    MatrixVector_Destroy(&Tempvec);

}

void PC_Calculate_Acceleration_PS(const MatrixVector_t *const LoadTdT, const MatrixVector_t *const LoadT, const MatrixVector_t *const RForceTdT,
        const MatrixVector_t *const RForceT, const MatrixVector_t *const VelTdT_Pred, const MatrixVector_t *const VelT_Pred, const MatrixVector_t *const AccT,
        const MatrixVector_t *const K, const MatrixVector_t *const C, const MatrixVector_t *const Meinv, const hysl_float_t alpha_H, const hysl_float_t a7,
        const hysl_float_t a8, MatrixVector_t *const AccTdT) {
    int incx = 1, incy = 1; /* Stride in the vectors */
    hysl_float_t Alpha, Beta; /* Constant for the BLAS routines */
    char uplo = 'L'; /* The lower part (upper part in C) will be used and the upper part
     * (lower part in C) will strictly not be referenced */

    /* BLAS: AccTdT = LoadTdT */
    hysl_copy(&AccTdT->Rows, LoadTdT->Array, &incx, AccTdT->Array, &incy);
    /* BLAS: AccTdT = (1 + alpha_H)*LoadTdT */
    Alpha = 1.0 + alpha_H;
    hysl_scal(&AccTdT->Rows, &Alpha, AccTdT->Array, &incx);
    /* BLAS: AccTdT = (1 + alpha_H)*LoadTdT - alpha_H*LoadT */
    Alpha = -alpha_H;
    hysl_axpy(&AccTdT->Rows, &Alpha, LoadT->Array, &incx, AccTdT->Array, &incy);
    /* BLAS: AccTdT = (1 + alpha_H)*LoadTdT - alpha_H*LoadT + alpha_H*RForceT */
    Alpha = alpha_H;
    hysl_axpy(&AccTdT->Rows, &Alpha, RForceT->Array, &incx, AccTdT->Array, &incy);
    /* BLAS: AccTdT = (1 + alpha_H)*LoadTdT - alpha_H*LoadT + alpha_H*RForceT - (1 + alpha_H)*RForceTdT*/
    Alpha = -(1.0 + alpha_H);
    hysl_axpy(&AccTdT->Rows, &Alpha, RForceTdT->Array, &incx, AccTdT->Array, &incy);
    /* BLAS: AccTdT = (1 + alpha_H)*LoadTdT - alpha_H*LoadT + alpha_H*RForceT - (1 + alpha_H)*RForceTdT +
     * alpha_H*C*VelT_Pred */
    Alpha = alpha_H;
    Beta = 1.0;
    hysl_spmv(&uplo, &AccTdT->Rows, &Alpha, C->Array, VelT_Pred->Array, &incx, &Beta, AccTdT->Array, &incy);
    /* BLAS: AccTdT = (1 + alpha_H)*LoadTdT - alpha_H*LoadT + alpha_H*RForceT - (1 + alpha_H)*RForceTdT +
     * alpha_H*C*VelT_Pred - (1 + alpha_H)*c*VelTdT_Pred */
    Alpha = -(1.0 + alpha_H);
    hysl_spmv(&uplo, &AccTdT->Rows, &Alpha, C->Array, VelTdT_Pred->Array, &incx, &Beta, AccTdT->Array, &incy);
    /* BLAS: AccTdT = (1 + alpha_H)*LoadTdT - alpha_H*LoadT + alpha_H*RForceT - (1 + alpha_H)*RForceTdT +
     * alpha_H*C*VelT_Pred - (1 + alpha_H)*c*VelTdT_Pred + alpha_H*a7*C*AccT. */
    Alpha = alpha_H * a7;
    hysl_spmv(&uplo, &AccTdT->Rows, &Alpha, C->Array, AccT->Array, &incx, &Beta, AccTdT->Array, &incy);
    /* BLAS: AccTdT = (1 + alpha_H)*LoadTdT - alpha_H*LoadT + alpha_H*RForceT - (1 + alpha_H)*RForceTdT +
     * alpha_H*C*VelT_Pred - (1 + alpha_H)*c*VelTdT_Pred + alpha_H*a7*C*AccT + alpha_H*a8*K*AccT */
    Alpha = alpha_H * a8;
    hysl_spmv(&uplo, &AccTdT->Rows, &Alpha, K->Array, AccT->Array, &incx, &Beta, AccTdT->Array, &incy);

    /* BLAS: AccTdT = Meinv*((1 + alpha_H)*LoadTdT - alpha_H*LoadT + alpha_H*RForceT - (1 + alpha_H)*RForceTdT +
     * alpha_H*C*VelT_Pred + (1 + alpha_H)*c*VelTdT_Pred + alpha_H*a7*C*AccT + alpha_H*a8*K*AccT) */
    Alpha = 1.0;
    Beta = 0.0;

    MatrixVector_t Tempvec;
    MatrixVector_Create(AccTdT->Rows, AccTdT->Cols, &Tempvec);

    hysl_copy(&AccTdT->Rows, AccTdT->Array, &incx, Tempvec.Array, &incy);
    hysl_spmv(&uplo, &AccTdT->Rows, &Alpha, Meinv->Array, Tempvec.Array, &incx, &Beta, AccTdT->Array, &incy);

    MatrixVector_Destroy(&Tempvec);
}
