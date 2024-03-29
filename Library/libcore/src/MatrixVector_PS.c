#include <stdio.h>
#include <stdlib.h>
#include <string.h>          /* For strcmp() */

#include "MatrixVector.h"
#include "MatrixVector_PS.h"

#include "Print_Messages.h"
#include "Auxiliary_Math.h"  /* For max() */
#include "Definitions.h"

#if _MKL_
#include <mkl_blas.h>
#include <mkl_lapack.h>
#else
#include "Netlib.h"
#endif

void MatrixVector_Create_PS (const int32_t Rows, const int32_t Cols, MatrixVector_t *const Matrix) {

    /* Check the input data */
    if ((Rows > 0) && (Cols > 0)) {
        Matrix->Rows = Rows;
        Matrix->Cols = Cols;
    } else {
        Print_Header( ERROR);
        fprintf( stderr, "MatrixVector_Create_PS(): The number of rows or columns must be greater than zero.\n");
        exit( EXIT_FAILURE);
    }

    /* Allocate the memory */
    Matrix->Array = NULL;
    Matrix->Array = (hysl_float_t*) calloc(((size_t) Matrix->Rows * (size_t) Matrix->Cols + (size_t) Matrix->Rows) / 2, sizeof(hysl_float_t));
    if (Matrix->Array == NULL) {
        Print_Header( ERROR);
        fprintf( stderr, "MatrixVector_Create_PS(): Out of memory.\n");
        exit( EXIT_FAILURE);
    }
}

void MatrixVector_Set2Value_PS (const hysl_float_t Value, MatrixVector_t *const Matrix) {
    int32_t incx = 0; /* No stride in the vector */
    int32_t incy = 1; /* Stride of one */
    int32_t Length = (Matrix->Rows * Matrix->Cols + Matrix->Rows) / 2;
    hysl_float_t Val = Value;

    /* BLAS: X = Val*X */
    hysl_copy(&Length, &Val, &incx, Matrix->Array, &incy);
}

void MatrixVector_ModifyElement_PS (const int RowIndex, const int ColIndex, const hysl_float_t Alpha, const char *Operation, MatrixVector_t *const Matrix) {

    const char *OpSet = "Set";
    const char *OpAdd = "Add";
    const char *OpMult = "Multiply";
    const char *OpDiv = "Divide";

    const int Position = RowIndex + (2 * Matrix->Cols - ColIndex) * (ColIndex - 1) / 2 - 1;

    if (RowIndex < ColIndex) {
        Print_Header( ERROR);
        fprintf( stderr, "MatrixVector_ModifyElement_PS(): Only the upper part of the matrix is stored and therefore the row index should be greater or equal than the column index.\n");
    }

    if (strcmp(Operation, OpSet) == 0) {
        Matrix->Array[Position] = Alpha;
    } else if (strcmp(Operation, OpAdd) == 0) {
        Matrix->Array[Position] = Matrix->Array[Position] + Alpha;
    } else if (strcmp(Operation, OpMult) == 0) {
        Matrix->Array[Position] = Matrix->Array[Position] * Alpha;
    } else if (strcmp(Operation, OpDiv) == 0) {
        Matrix->Array[Position] = Matrix->Array[Position] / Alpha;
    } else {
        Print_Header( ERROR);
        fprintf( stderr, "MatrixVector_ModifyElement_PS(): Operation '%s' not identified. Valid operations are:\n", Operation);
        fprintf( stderr, "[......] 1) %s.\n", OpSet);
        fprintf( stderr, "[......] 2) %s.\n", OpAdd);
        fprintf( stderr, "[......] 3) %s.\n", OpMult);
        fprintf( stderr, "[......] 4) %s.\n", OpDiv);
        exit( EXIT_FAILURE);
    }
}

void MatrixVector_Add3Mat_PS (const MatrixVector_t *const MatA, const MatrixVector_t *const MatB, const MatrixVector_t *const MatC, const Scalars_t Const, MatrixVector_t *const MatY) {
    int32_t incx = 1;
    int32_t incy = 1;

    /* Length of the matrix in packed storage */
    int32_t Length = (MatY->Rows * MatY->Cols + MatY->Rows) / 2;

    /* BLAS: Y = A */
    hysl_copy(&Length, MatA->Array, &incx, MatY->Array, &incy);

    /* BLAS: Calculates Y = A  */
    hysl_float_t Scalar = Const.Alpha;
    hysl_scal(&Length, &Scalar, MatY->Array, &incx);

    /* BLAS: Calculates Y = Beta*B + Y. Only computes half of the matrix */
    Scalar = Const.Beta;
    hysl_axpy(&Length, &Scalar, MatB->Array, &incx, MatY->Array, &incy);

    /* BLAS: Calculates Y = Gamma*C + Y. Only computes half of the matrix */
    Scalar = Const.Gamma;
    hysl_axpy(&Length, &Scalar, MatC->Array, &incx, MatY->Array, &incy);

}

void MatrixVector_Destroy_PS (MatrixVector_t *const Matrix) {
    /* Set the number of rows and columns to 0 */
    Matrix->Rows = 0;
    Matrix->Cols = 0;
    free(Matrix->Array);
}
