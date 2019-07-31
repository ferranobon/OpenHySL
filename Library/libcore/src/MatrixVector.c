#include <stdio.h>
#include <stdlib.h>
#include <string.h>  /* For strcmp() */

#include "MatrixVector.h"
#include "Print_Messages.h"
#include "Auxiliary_Math.h"  /* For max() */
#include "Definitions.h"

#if _MKL_
#include <mkl_blas.h>
#include <mkl_lapack.h>
#else
#include "Netlib.h"
#endif

void MatrixVector_Create(const int32_t Rows, const int32_t Cols, MatrixVector_t *const MatVec) {

    /* Check the input data */
    if ((Rows > 0) && (Cols > 0)) {
        MatVec->Rows = Rows;
        MatVec->Cols = Cols;
    } else {
        Print_Header( ERROR);
        fprintf( stderr, "MatrixVector_Create(): The number of rows or columns must be greater than zero.\n");
        exit( EXIT_FAILURE);
    }

    /* Allocate the memory */
    MatVec->Array = NULL;
    MatVec->Array = (HYSL_FLOAT*) calloc((size_t) MatVec->Rows * (size_t) MatVec->Cols, sizeof(HYSL_FLOAT));
    if (MatVec->Array == NULL) {
        Print_Header( ERROR);
        fprintf( stderr, "MatrixVector_Create(): Out of memory.\n");
        exit( EXIT_FAILURE);
    }
}

void MatrixVector_Set2Value(const HYSL_FLOAT Value, MatrixVector_t *const MatVec) {
    int32_t incx = 0; /* No stride in the vector */
    int32_t incy = 1; /* Stride of one */
    int32_t Length;
    HYSL_FLOAT Val;

    Length = MatVec->Rows * MatVec->Cols;
    Val = Value;

    /* BLAS: All elements of MatVec are equal to Value */
    hysl_copy(&Length, &Val, &incx, MatVec->Array, &incy);
}

void MatrixVector_ModifyElement(const int32_t RowIndex, const int32_t ColIndex, const HYSL_FLOAT Alpha, const char *Operation, MatrixVector_t *const MatVec) {

    const char *OpSet = "Set";
    const char *OpAdd = "Add";
    const char *OpMult = "Multiply";
    const char *OpDiv = "Divide";

    if (strcmp(Operation, OpSet) == 0) {
        MatVec->Array[(RowIndex - 1) * MatVec->Cols + (ColIndex - 1)] = Alpha;
    } else if (strcmp(Operation, OpAdd) == 0) {
        MatVec->Array[(RowIndex - 1) * MatVec->Cols + (ColIndex - 1)] = MatVec->Array[(RowIndex - 1) * MatVec->Cols + (ColIndex - 1)] + Alpha;
    } else if (strcmp(Operation, OpMult) == 0) {
        MatVec->Array[(RowIndex - 1) * MatVec->Cols + (ColIndex - 1)] = MatVec->Array[(RowIndex - 1) * MatVec->Cols + (ColIndex - 1)] * Alpha;
    } else if (strcmp(Operation, OpDiv) == 0) {
        MatVec->Array[(RowIndex - 1) * MatVec->Cols + (ColIndex - 1)] = MatVec->Array[(RowIndex - 1) * MatVec->Cols + (ColIndex - 1)] / Alpha;
    } else {
        Print_Header( ERROR);
        fprintf( stderr, "MatrixVector_ModifyElement: Operation '%s' not identified. Valid operations are:\n", Operation);
        fprintf( stderr, "[......] 1) %s.\n", OpSet);
        fprintf( stderr, "[......] 2) %s.\n", OpAdd);
        fprintf( stderr, "[......] 3) %s.\n", OpMult);
        fprintf( stderr, "[......] 4) %s.\n", OpDiv);
        exit( EXIT_FAILURE);
    }
}

void MatrixVector_Add3Mat(const MatrixVector_t *const MatA, const MatrixVector_t *const MatB, const MatrixVector_t *const MatC, const Scalars_t Const,
        MatrixVector_t *const MatY) {

    int32_t incx, incy;      /* Stride in the vectors for BLAS library */
    int32_t lda, ldy;
    int32_t ione;
    HYSL_FLOAT done, Scalar; /* Constant to use in the BLAS library */
    char uplo;               /* BLAS & LAPACK: Character to specify which part of the matrix has been referenced. */
    int32_t info;            /* LAPACK: Variable to inform if the operations of Cholesky factorization and inverse were successful or not. */

    ione = 1;
    incx = 1;
    incy = 1;
    uplo = 'L'; /* Character defining that the lower part of the symmetric matrix is referenced
                 * (see man dpotrf) */
    done = 1.0;

    lda = Max(1, MatA->Rows);
    ldy = Max(1, MatY->Rows);

    /* LAPACK: Y = A */
    hysl_lacpy(&uplo, &MatY->Rows, &MatY->Cols, MatA->Array, &lda, MatY->Array, &ldy);

    /* LAPACK: Calculates Y = A  */
    Scalar = Const.Alpha;
    hysl_lascl(&uplo, &ione, &ione, &done, &Scalar, &MatY->Rows, &MatY->Cols, MatY->Array, &ldy, &info);

    if (info < 0) {
        Print_Header( ERROR);
        fprintf( stderr, "dlascl: The %d-th argument had an illegal value.\n", -info);
        exit( EXIT_FAILURE);
    }

    /* BLAS: Calculates Y = Beta*B + Y. Only computes half of the matrix */
    Scalar = Const.Beta;
    for (int32_t idx = 0; idx < MatY->Rows; idx++) {
        int32_t Length = MatY->Rows - idx;
        hysl_axpy(&Length, &Scalar, &MatB->Array[(idx * MatB->Rows) + idx], &incx, &MatY->Array[(idx * MatY->Rows) + idx], &incy);
    }

    /* BLAS: Calculates Y = Gamma*C + Y. Only computes half of the matrix */
    Scalar = Const.Gamma;
    for (int32_t idx = 0; idx < MatY->Rows; idx++) {
        int32_t Length = MatY->Rows - idx;
        hysl_axpy(&Length, &Scalar, &MatC->Array[(idx * MatC->Rows) + idx], &incx, &MatY->Array[(idx * MatY->Rows) + idx], &incy);
    }
}

void MatrixVector_Destroy(MatrixVector_t *const MatVec) {
    /* Set the number of rows and columns to 0 */
    MatVec->Rows = 0;
    MatVec->Cols = 0;
    free(MatVec->Array);
}
