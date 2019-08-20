#include "Input_Load.h"
#include "MatrixVector.h"
#include "Definitions.h"

#if _MKL_
#include <mkl_blas.h>
#else
#include "Netlib.h"
#endif

void InputLoad_AbsValues (const MatrixVector_t *const Stiff, const MatrixVector_t *const Damp, const MatrixVector_t *const GDisp, const MatrixVector_t *const GVel, MatrixVector_t *const InLoad) {
    int32_t incx = 1;
    int32_t incy = 1;
    hysl_float_t Alpha = 1.0;
    hysl_float_t Beta = 0.0;
    char uplo = 'L'; /* The lower part (FORTRAN) of the symmetric matrix is referenced (upper part in C) */

    /* BLAS: li = K*ug */
    hysl_symv(&uplo, &InLoad->Rows, &Alpha, Stiff->Array, &InLoad->Rows, GDisp->Array, &incx, &Beta, InLoad->Array, &incy);

    /* BLAS: li = K*ug + C*vg = li + C*vg */
    Beta = 1.0;
    hysl_symv(&uplo, &InLoad->Rows, &Alpha, Damp->Array, &InLoad->Rows, GVel->Array, &incx, &Beta, InLoad->Array, &incy);
}

void InputLoad_AbsValues_PS (const MatrixVector_t *const Stiff, const MatrixVector_t *const Damp, const MatrixVector_t *const GDisp, const MatrixVector_t *const GVel, MatrixVector_t *const InLoad) {
    int32_t incx = 1;
    int32_t incy = 1;
    hysl_float_t Alpha = 1.0;
    hysl_float_t Beta = 0.0;
    char uplo = 'L'; /* The lower part (FORTRAN) of the symmetric matrix is referenced (upper part in C) */

    /* BLAS: li = K*ug */
    hysl_spmv(&uplo, &InLoad->Rows, &Alpha, Stiff->Array, GDisp->Array, &incx, &Beta, InLoad->Array, &incy);

    /* BLAS: li = K*ug + C*vg = li + C*vg */
    Beta = 1.0;
    hysl_spmv(&uplo, &InLoad->Rows, &Alpha, Damp->Array, GVel->Array, &incx, &Beta, InLoad->Array, &incy);
}

void InputLoad_RelValues (const MatrixVector_t *const Mass, const MatrixVector_t *const GAcc, MatrixVector_t *const InLoad) {

    int32_t incx = 1;
    int32_t incy = 1; /* Stride in the vectors for BLAS library */
    hysl_float_t Alpha = -1.0;
    hysl_float_t Beta = 0.0; /* Constants to use in the BLAS library */
    char uplo = 'L'; /* Character defining that the lower part (FORTRAN) of the symmetric matrix is referenced (upper part in C) */

    /* BLAS: li = -M*ag */
    hysl_symv(&uplo, &InLoad->Rows, &Alpha, Mass->Array, &InLoad->Rows, GAcc->Array, &incx, &Beta, InLoad->Array, &incy);
}

void InputLoad_RelValues_PS (const MatrixVector_t *const Mass, const MatrixVector_t *const GAcc, MatrixVector_t *const InLoad) {

    int32_t incx = 1;
    int32_t incy = 1; /* Stride in the vectors for BLAS library */
    hysl_float_t Alpha = -1.0;
    hysl_float_t Beta = 0.0; /* Constants to use in the BLAS library */
    char uplo = 'L'; /* Character defining that the lower part (FORTRAN) of the symmetric matrix is referenced (upper part in C) */

    /* BLAS: li = -M*ag */
    hysl_spmv(&uplo, &InLoad->Rows, &Alpha, Mass->Array, GAcc->Array, &incx, &Beta, InLoad->Array, &incy);
}

void InputLoad_Generate_LoadVectorForm (int32_t *DOF, MatrixVector_t *const LoadVectorForm) {
    int32_t loadIdx = 0;

    while (loadIdx < LoadVectorForm->Rows) {
        for (int32_t jdx = 1; jdx < DOF[0]; jdx++) {
            LoadVectorForm->Array[loadIdx] = (hysl_float_t) DOF[jdx];
            loadIdx = loadIdx + 1;
        }
    }
}

void InputLoad_Apply_LoadVectorForm (const MatrixVector_t *const LoadForm1, const MatrixVector_t *const LoadForm2, const MatrixVector_t *const LoadForm3, const hysl_float_t Value1,
        const hysl_float_t Value2, const hysl_float_t Value3, MatrixVector_t *const LoadVector) {
    int32_t incx = 1;
    int32_t incy = 1;

    hysl_float_t Scalar = Value1;
    hysl_copy(&LoadVector->Rows, LoadForm1->Array, &incx, LoadVector->Array, &incy);
    hysl_scal(&LoadVector->Rows, &Scalar, LoadVector->Array, &incx);

    Scalar = Value2;
    hysl_axpy(&LoadVector->Rows, &Scalar, LoadForm2->Array, &incx, LoadVector->Array, &incy);

    Scalar = Value3;
    hysl_axpy(&LoadVector->Rows, &Scalar, LoadForm3->Array, &incx, LoadVector->Array, &incy);
}
