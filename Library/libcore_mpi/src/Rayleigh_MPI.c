#include <stdio.h>          /* For printf(), fprintf() */
#include <stdlib.h>         /* For exit() */

#include "Auxiliary_Math.h" /* For Max() */
#include "MatrixVector_MPI.h"   /* MatrixVector definition */
#include "Print_Messages.h" /* For Print_Header() */
#include "Rayleigh.h"       /* Rayleigh damping routines */
#include "Definitions.h"

#if _MKL_
#include <mkl_pblas.h>
#include <mkl_scalapack.h>
#else
#include "Netlib.h"
#endif

void Rayleigh_Damping_MPI (PMatrixVector_t *const Mass, PMatrixVector_t *const Stiff, PMatrixVector_t *const Damp, const Rayleigh_t *const Rayleigh) {
    int32_t ione = 1;

    /* ScaLAPACK: Perform C = M (locally. There is no communication) */
    char uplo = 'L';  /* The lower part of the matrix will be used; the upper part will strictly not be referenced */
    hysl_placpy(&uplo, &Damp->GlobalSize.Row, &Damp->GlobalSize.Col, Mass->Array, &ione, &ione, Mass->Desc, Damp->Array, &ione, &ione, Damp->Desc);

    /* ScaLAPACK: Perform C = alpha*M + beta*K = alpha*C + beta*K */
    hysl_float_t alpha = Rayleigh->Alpha;
    hysl_float_t beta = Rayleigh->Beta;
    char trans = 'N'; /* The operation will not use the transpose matrix */
    hysl_ptradd(&uplo, &trans, &Damp->GlobalSize.Row, &Damp->GlobalSize.Col, &beta, Stiff->Array, &ione, &ione, Stiff->Desc, &alpha, Damp->Array, &ione, &ione, Damp->Desc);
}
