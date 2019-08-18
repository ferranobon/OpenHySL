#include <stdio.h>          /* For printf(), fprintf() */
#include <stdlib.h>         /* For exit() */

#include "Auxiliary_Math.h" /* For Max() */
#include "MatrixVector.h"   /* MatrixVector definition */
#include "Print_Messages.h" /* For Print_Header() */
#include "Rayleigh.h"       /* Rayleigh damping routines */
#include "Definitions.h"

#if _MKL_
#include <mkl_blas.h>
#include <mkl_lapack.h>
#else
#include "Netlib.h"
#endif

void Rayleigh_Damping(const MatrixVector_t *const Mass, const MatrixVector_t *const Stiff, MatrixVector_t *const Damp, const Rayleigh_t *const Rayleigh) {
    int32_t ione = 1;
    hysl_float_t done = 1.0;
    int32_t incx = 1;
    int32_t incy = 1;
    char uplo = 'L'; /* The lower part of the matrix will be used and the upper part will strictly
                 * not be referenced */
    int32_t Rows = Mass->Rows;
    int32_t Cols = Mass->Cols;
    hysl_float_t alpha = Rayleigh->Alpha;
    hysl_float_t beta = Rayleigh->Beta;

    int32_t lda = Max(1, Rows);
    int32_t ldb = Max(1, Damp->Rows);

    /* LAPACK: C = M */
    hysl_lacpy(&uplo, &Rows, &Cols, Mass->Array, &lda, Damp->Array, &ldb);

    int32_t info = -1;
    /* LAPACK: C = Rayleigh.alpha*M = Rayleigh.alpha*C */
    hysl_lascl(&uplo, &ione, &ione, &done, &alpha, &Damp->Rows, &Damp->Cols, Damp->Array, &lda, &info);

    if (info < 0) {
        Print_Header( ERROR);
        fprintf( stderr, "dlascl: The %d-th argument had an illegal value.\n", -info);
        exit( EXIT_FAILURE);
    }

    /* BLAS: C = alpha*M + beta*K = C + beta*K. Only half of the matrix is calculated */
    for (int32_t idx = 0; idx < Damp->Rows; idx++) {
        int32_t Length = Damp->Rows - idx;
        hysl_axpy(&Length, &beta, &Stiff->Array[idx * Stiff->Rows + idx], &incx, &Damp->Array[idx * Damp->Rows + idx], &incy);
    }

    Print_Header( SUCCESS);
    printf("Damping matrix successfully calculated.\n");
}

void Rayleigh_Damping_PS(const MatrixVector_t *const Mass, const MatrixVector_t *const Stiff, MatrixVector_t *const Damp, const Rayleigh_t *const Rayleigh) {
    int32_t incx = 1;
    int32_t incy = 1;

    hysl_float_t alpha = Rayleigh->Alpha;
    hysl_float_t beta = Rayleigh->Beta;

    int32_t Length = (Damp->Rows * Damp->Cols + Damp->Rows) / 2;

    /* BLAS: C = M */
    hysl_copy(&Length, Mass->Array, &incx, Damp->Array, &incy);

    /* BLAS: C = Rayleigh.alpha*M = Rayleigh.alpha*C */
    hysl_scal(&Length, &alpha, Damp->Array, &incx);

    /* BLAS: C = alpha*M + beta*K = C + beta*K. Only half of the matrix is calculated */
    hysl_axpy(&Length, &beta, Stiff->Array, &incx, Damp->Array, &incy);

    Print_Header( SUCCESS);
    printf("Damping matrix successfully calculated.\n");
}
