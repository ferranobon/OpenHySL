#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "MatrixVector.h"
#include "MatrixVector_MPI.h"
#include "Print_Messages.h"
#include "Auxiliary_Math.h"
#include "Definitions.h"

#if _MKL_
#include "mkl_blas.h"
#include "mkl_pblas.h"
#include "mkl_scalapack.h"
#include "Cblacs.h"
#include "Scalapack_Aux.h"
#else
#include "Netlib.h"
#endif

void PMatrixVector_Add3Mat (PMatrixVector_t *const MatA, PMatrixVector_t *const MatB, PMatrixVector_t *const MatC, const Scalars_t Const, PMatrixVector_t *const MatY) {

    char uplo, trans;
    int ione;
    hysl_float_t ScalarA, ScalarB;

    ione = 1;
    trans = 'N'; /* The operation will not use the transpose matrix */
    uplo = 'L'; /* The lower part of the matrix will be used; the upper part will strictly not be
     * referenced */
    /* ScaLAPACK: Perform Y = A (locally. There is no communication) */
    hysl_placpy(&uplo, &MatY->GlobalSize.Row, &MatY->GlobalSize.Col, MatA->Array, &ione, &ione, MatA->Desc, MatY->Array, &ione, &ione, MatY->Desc);

    /* ScaLAPACK: Perform Y = beta*B + alpha*A = beta*B + alpha*Y */
    ScalarA = Const.Alpha;
    ScalarB = Const.Beta;
    hysl_ptradd(&uplo, &trans, &MatY->GlobalSize.Row, &MatY->GlobalSize.Col, &ScalarB, MatB->Array, &ione, &ione, MatB->Desc, &ScalarA, MatY->Array, &ione, &ione, MatY->Desc);

    /* ScaLAPACK: Perform Y = gamma*C + beta*B = gamma*C + 1.0*Y */
    ScalarA = 1.0;
    ScalarB = Const.Gamma;
    hysl_ptradd(&uplo, &trans, &MatY->GlobalSize.Row, &MatY->GlobalSize.Col, &ScalarB, MatC->Array, &ione, &ione, MatC->Desc, &ScalarA, MatY->Array, &ione, &ione, MatY->Desc);
}

void PMatrixVector_Create (int icntxt, const int Rows, const int Cols, const int BlRows, int const BlCols, PMatrixVector_t *const MatVec) {

    int myrow, mycol; /* Variables to store row and column in the process grid */
    int nprow, npcol;

    int izero = 0; /* Zero of type integer. Used by numroc_( ) */
    int lld, info;

    MatVec->GlobalSize.Row = Rows; /* Number of rows in the global array */
    MatVec->GlobalSize.Col = Cols; /* Number of columns in the global array */

    MatVec->BlockSize.Row = BlRows; /* Block size in the vertical dimension */
    MatVec->BlockSize.Col = BlCols; /* Block size in the horizontal dimension */

    Cblacs_gridinfo(icntxt, &nprow, &npcol, &myrow, &mycol); /* Get information about the grid */

    /* Compute the size of the local matrices */
    MatVec->LocalSize.Row = numroc_(&MatVec->GlobalSize.Row, &MatVec->BlockSize.Row, &myrow, &izero, &nprow);
    MatVec->LocalSize.Col = numroc_(&MatVec->GlobalSize.Col, &MatVec->BlockSize.Col, &mycol, &izero, &npcol);
    lld = Max(1, MatVec->LocalSize.Row);

    descinit_(MatVec->Desc, &MatVec->GlobalSize.Row, &MatVec->GlobalSize.Col, &MatVec->BlockSize.Row, &MatVec->BlockSize.Col, &izero, &izero, &icntxt, &lld, &info);
    if (info < 0) {
        Print_Header( ERROR);
        fprintf( stderr, "descinit: The %d-th argument had an illegal value.\n", -info);
        exit( EXIT_FAILURE);
    }

    /* Allocate memory for the local array */
    MatVec->Array = (hysl_float_t*) calloc((size_t) MatVec->LocalSize.Row * (size_t) MatVec->LocalSize.Col, sizeof(hysl_float_t));
    if (MatVec->Array == NULL) {
        Print_Header( ERROR);
        fprintf( stderr, "PMatrixVector_Create: Out of memory.\n");
        exit( EXIT_FAILURE);
    }
}

void PMatrixVector_Destroy (PMatrixVector_t *const MatVec) {

    /* Set Global and local sizes to 0 */
    MatVec->GlobalSize.Row = 0;
    MatVec->GlobalSize.Col = 0;

    MatVec->LocalSize.Row = 0;
    MatVec->LocalSize.Col = 0;

    /* Deallocate memory */
    free(MatVec->Array);
}

void PMatrixVector_ModifyElement (int GRowIndex, int GColIndex, const hysl_float_t Alpha, const char *Operation, PMatrixVector_t *const MatVec) {

    int myrow, mycol, nprow, npcol;

    const char *OpSet = "Set";
    const char *OpAdd = "Add";
    const char *OpMult = "Multiply";
    const char *OpDiv = "Divide";

    int LRowIndex, LColIndex;
    int RowProcess, ColProcess;

    /* Get grid info */
    Cblacs_gridinfo(MatVec->Desc[1], &nprow, &npcol, &myrow, &mycol);

    /* Given the global index of an element (GRowIndex, GColIndex) returns the local index of the element
     * (LRowIndex, LColIndex) and the coordinates of the process (Row Process, ColProcess) */
    infog2l_(&GRowIndex, &GColIndex, MatVec->Desc, &nprow, &npcol, &myrow, &mycol, &LRowIndex, &LColIndex, &RowProcess, &ColProcess);

    /* Modify the value in the local matrix. -1 is substracted because C starts the array indexes at 0, while
     * FORTRAN starts them at 1 (infog2l is a FORTRAN routine) */
    if (myrow == RowProcess && mycol == ColProcess) {

        if (strcmp(Operation, OpSet) == 0) {
            MatVec->Array[(LRowIndex - 1) * MatVec->LocalSize.Col + (LColIndex - 1)] = Alpha;
        } else if (strcmp(Operation, OpAdd) == 0) {
            MatVec->Array[(LRowIndex - 1) * MatVec->LocalSize.Col + (LColIndex - 1)] = MatVec->Array[(LRowIndex - 1) * MatVec->LocalSize.Col + (LColIndex - 1)] + Alpha;
        } else if (strcmp(Operation, OpMult) == 0) {
            MatVec->Array[(LRowIndex - 1) * MatVec->LocalSize.Col + (LColIndex - 1)] = MatVec->Array[(LRowIndex - 1) * MatVec->LocalSize.Col + (LColIndex - 1)] * Alpha;
        } else if (strcmp(Operation, OpDiv) == 0) {
            MatVec->Array[(LRowIndex - 1) * MatVec->LocalSize.Col + (LColIndex - 1)] = MatVec->Array[(LRowIndex - 1) * MatVec->LocalSize.Col + (LColIndex - 1)] / Alpha;
        } else {
            Print_Header( ERROR);
            fprintf( stderr, "PMatrixVector_ModifyElement: Operation '%s' not identified. Valid operations are:\n", Operation);
            fprintf( stderr, "[......] 1) %s.\n", OpSet);
            fprintf( stderr, "[......] 2) %s.\n", OpAdd);
            fprintf( stderr, "[......] 3) %s.\n", OpMult);
            fprintf( stderr, "[......] 4) %s.\n", OpDiv);
            exit( EXIT_FAILURE);
        }
    }

}

void PMatrixVector_Set2Value (const hysl_float_t Value, PMatrixVector_t *const MatVec) {

    int incx, incy;
    int Length;
    hysl_float_t Val;

    incx = 0;
    incy = 1;

    Length = MatVec->LocalSize.Row * MatVec->LocalSize.Col;
    Val = Value;

    hysl_copy(&Length, &Val, &incx, MatVec->Array, &incy);
}
