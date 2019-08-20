#include "PCMethods.h"
#include "MatrixVector.h"
#include "MatrixVector_Sp.h"

#include "Definitions.h"

#include <mkl_blas.h>
#include <mkl_spblas.h>

void PC_Correct_Acceleration_Sp (const MatrixVector_Sp_t *const MInv, const MatrixVector_t *const In_LoadT, const MatrixVector_t *const fc, MatrixVector_t *const Tempvec,
        MatrixVector_t *const AccTdT_Corr) {
    int32_t incx = 1;
    int32_t incy = 1; /* Stride in the vectors */

    /* BLAS: tempvec = In_LoadT */
    hysl_copy(&Tempvec->Rows, In_LoadT->Array, &incx, Tempvec->Array, &incy);

    /* BLAS: tempvec = In_LoadT + fc = tempvec + fc */
    hysl_float_t Alpha = 1.0;
    hysl_axpy(&Tempvec->Rows, &Alpha, fc->Array, &incx, Tempvec->Array, &incy);

    /* Sparse BLAS: AccTdT_Corr = MInv*(In_LoadT + fc) = MInv*Tempvec */
    char trans = 'N'; /* No transpose operation */
    hysl_float_t Beta = 0.0;
    char matdescra[6] = { 'S', /* The matrix is symmetric */
        'U',                   /* The upper part is referenced */
        'N',                   /* Non-unit values in the diagonal */
        'F'                    /* One based index */
    };
    hysl_mkl_csrmv(&trans, &Tempvec->Rows, &Tempvec->Rows, &Alpha, matdescra, MInv->Values, MInv->Columns, MInv->RowIndex, &MInv->RowIndex[1], Tempvec->Array, &Beta, AccTdT_Corr->Array);
}
