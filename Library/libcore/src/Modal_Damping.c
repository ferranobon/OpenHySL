#include <stdio.h>          /* For printf(), fprintf() */
#include <stdlib.h>         /* For exit() */
#include <math.h>

#include "Auxiliary_Math.h" /* For Max() */
#include "MatrixVector.h"   /* MatrixVector definition */
#include "Print_Messages.h" /* For Print_Header() */
#include "Modal_Damping.h"       /* Rayleigh damping routines */
#include "Definitions.h"

#if _MKL_
#include <mkl_blas.h>
#include <mkl_lapack.h>
#else
#include "Netlib.h"
#endif

void Modal_Damping(const MatrixVector_t *const Mass, const MatrixVector_t *const Stiff, MatrixVector_t *const Damp, double DampFactor) {
    MatrixVector_t EVectors;
    MatrixVector_t EValues;
    MatrixVector_t temp, temp1;

    Print_Header(WARNING);
    fprintf(stderr, "Modal_Damping(): Untested routine.");

    int32_t ione = 1;
    HYSL_FLOAT done = 1.0;
    int32_t incx = 1;
    int32_t incy = 1;
    char uplo = 'L'; /* The lower part of the matrix will be used and the upper part will strictly
                      * not be referenced */
    int32_t Rows = Mass->Rows;
    int32_t Cols = Mass->Cols;

    MatrixVector_Create(Rows, Cols, &EVectors);
    MatrixVector_Create(Rows, 1,    &EValues);
    MatrixVector_Create(Rows, Cols, &temp);
    MatrixVector_Create(Rows, Cols, &temp1);

    for (int32_t idx = 0; idx < Rows; idx++) {
        if (Mass->Array[idx * Rows + idx] == 0.0) {
            temp1.Array[idx * Rows + idx] = 1E-12;
        } else {
            temp1.Array[idx * Rows + idx] = Mass->Array[idx * Rows + idx];
        }
        printf("%lE %lE\n", Mass->Array[idx * Rows + idx], temp1.Array[idx * Rows + idx]);
    }

    Compute_Eigenvalues_Eigenvectors(Stiff, &temp1, &EValues, &EVectors);
    for (int32_t idx = 0; idx < Rows; idx++) {
        Damp->Array[idx * Rows + idx] = 2.0 * sqrt(EValues.Array[idx]) * DampFactor;
    }

    HYSL_FLOAT alpha = 1.0;
    HYSL_FLOAT beta = 0.0;
    char transa = 'T';
    char transb = 'N';
    dgemm_(&transa, &transb, &Rows, &Cols, &Rows, &alpha, EVectors.Array, &Rows, Damp->Array, &Rows, &beta, temp.Array, &Rows);
    transa = 'N';
    transb = 'T';
    dgemm_(&transa, &transb, &Rows, &Cols, &Rows, &alpha, temp.Array, &Rows, EVectors.Array, &Rows, &beta, Damp->Array, &Rows);

    MatrixVector_ToFile(Damp, "Damp.txt");

    int32_t lda = Max(1, Rows);
    int32_t ldb = Max(1, Damp->Rows);

    MatrixVector_Destroy(&EVectors);
    MatrixVector_Destroy(&EValues);
    MatrixVector_Destroy(&temp);
    MatrixVector_Destroy(&temp1);

    Print_Header( SUCCESS);
    printf("Damping matrix successfully calculated.\n");
}
