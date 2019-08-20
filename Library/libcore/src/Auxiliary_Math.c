#include <stdio.h>              /* For printf(), fprintf() */
#include <stdlib.h>             /* For exit() */
#include <stdbool.h>            /* For bool, true and false */
#include <math.h>               /* For sqrt(), log() */

#include "Auxiliary_Math.h"
#include "Print_Messages.h"

#include "Definitions.h"

#if _MKL_
#include <mkl_blas.h>
#include <mkl_lapack.h>
#else
#include "Netlib.h"
#endif

int32_t Max (const int32_t a, const int32_t b) {
    if (a >= b) {
        return a;
    } else {
        return b;
    }
}

int32_t Min (const int32_t a, const int32_t b) {
    if (a <= b) {
        return a;
    } else
        return b;
}

hysl_float_t signum (const hysl_float_t num) {
# if _FLOAT_
     if (num > 0.0f ){
	  return 1.0f;
     } else if (num < 0.0f){
	  return -1.0f;
     } else return 0.0f;
#else
    if (num > 0.0) {
        return 1.0;
    } else if (num < 0.0) {
        return -1.0;
    } else
        return 0.0;
#endif
}

hysl_float_t norm (const int32_t length, const hysl_float_t *const Vector) {
    hysl_float_t temp = 0.0;

    for (int32_t idx = 0; idx < length; idx++) {
#if _FLOAT_
        temp = temp + powf(Vector[idx], 2.0f);
#else
        temp = temp + pow(Vector[idx], 2.0);
#endif
    }

#if _FLOAT_
     return sqrtf(temp);
#else
    return sqrt(temp);
#endif
}

MatrixVector_t Generate_IdentityMatrix (const int32_t Rows, const int32_t Cols) {
    MatrixVector_t Identity;

    if (Rows != Cols) {
        Print_Header( ERROR);
        fprintf( stderr, "Generate_IdentityMatrix: The number of rows and columns must be the same.\n");
        exit( EXIT_FAILURE);
    }

    MatrixVector_Create(Rows, Cols, &Identity);

    for (int32_t idx = 0u; idx < Rows; idx++) {
        Identity.Array[idx + (Rows * idx)] = 1.0;
    }

    return Identity;
}

int32_t MatrixVector_ReturnIndex_UPS (const int32_t RowIdx, const int32_t ColIdx, const int32_t numRows) {
    int32_t idx = 0u;

    if (RowIdx >= ColIdx) {
        idx = RowIdx + (((2 * numRows) - ColIdx) * (ColIdx - 1) / 2) - 1;
    } else {
        idx = ColIdx + (((2 * numRows) - RowIdx) * (RowIdx - 1) / 2) - 1;
    }
    return idx;
}

int32_t MatrixVector_ReturnIndex_LPS (const int32_t RowIdx, const int32_t ColIdx) {
    int32_t idx = 0;

    if (ColIdx >= RowIdx) {
        idx = RowIdx + (ColIdx * (ColIdx - 1) / 2) - 1;
    } else {
        idx = ColIdx + (RowIdx * (RowIdx - 1) / 2) - 1;
    }
    return idx;
}

void Compute_Eigenvalues_Eigenvectors (MatrixVector_t *const MatrixA, MatrixVector_t *const MatrixB, MatrixVector_t *const EigenValues, MatrixVector_t *const EigenVectors) {
    if ((MatrixA->Rows != MatrixB->Rows) || (MatrixA->Cols != MatrixB->Cols)) {
        Print_Header( ERROR);
        fprintf( stderr, "Compute_Eigenvalues_Eigenvectors: The matrices must be identical.\n");
        exit(EXIT_FAILURE);
    }

    int32_t lwork = (3 * MatrixA->Rows) - 1;  /* Dimension of the array work */
    int32_t lda = Max(1, MatrixA->Rows);
    int32_t ldb = lda;

    int32_t Length = MatrixA->Rows * MatrixA->Cols;
    hysl_float_t *TempMat = (hysl_float_t*) calloc((size_t) Length, sizeof(hysl_float_t));
    hysl_float_t *work = (hysl_float_t*) calloc((size_t) lwork, sizeof(hysl_float_t));

    if ((TempMat == NULL) || (work == NULL)) {
        Print_Header(ERROR);
        fprintf( stderr, "Compute_Eigenvalues_Eigenvectors(): Out of memory.\n");
        exit(EXIT_FAILURE);
    }

    /* DSYGV_:On Entry EigenVectors must contain the Matrix A */
    int32_t one = 1;
    hysl_copy(&Length, MatrixA->Array, &one, EigenVectors->Array, &one);
    hysl_copy(&Length, MatrixB->Array, &one, TempMat, &one);

    Length = MatrixA->Rows;
    int32_t info = 0;
    hysl_sygv(&one, "V", "L", &Length, EigenVectors->Array, &lda, TempMat, &ldb, EigenValues->Array, work, &lwork, &info);

    if (info == 0) {
        Print_Header( SUCCESS);
        printf("Successfully calculated the eigenvalues and eigenvectors.\n");
    } else if (info < 0) {
        Print_Header( ERROR);
        fprintf( stderr, "Compute_Eigenvalues_Eigenvectors: the %d-th argument of the function hysl_sygv() had an illegal value", info);
        exit( EXIT_FAILURE);
    } else if (info > 0) {
        if (info <= EigenVectors->Rows) {
            Print_Header( ERROR);
            fprintf( stderr, "Compute_Eigenvalues_Eigenvectors: %d off-diagonal elements of an intermediate tridiagonal form did not converge to zero.\n", info);
            exit( EXIT_FAILURE);
        } else {
            Print_Header( ERROR);
            fprintf( stderr,
                    "Compute_Eigenvalues_Eigenvectors: the leading minor of order %d of MatrixB (TempMat) is not positive definite. The factorization of MatrixB could not be completed and no eigenvalues or eigenvectors were computed.\n",
                    info - MatrixB->Rows);
            exit( EXIT_FAILURE);
        }
    }

    for (int32_t idx = 0; idx < Length - 1; idx++) {
        /* Order the Eigenvalues and eigenvectors in ascendent order */
        if (EigenValues->Array[idx] > EigenValues->Array[idx + 1u]) {
            /* Swap Eigenvalues */
            hysl_float_t temp = EigenValues->Array[idx];
            EigenValues->Array[idx] = EigenValues->Array[idx + 1u];
            EigenValues->Array[idx + 1] = temp;
            /* Now Swap Eigenvectors */
            for (int32_t jdx = 0; jdx < Length; jdx++) {
                temp = EigenVectors->Array[(Length * (idx + 1)) + jdx];
                EigenVectors->Array[(Length * idx) + jdx] = EigenVectors->Array[(Length * (idx + 1)) + jdx];
                EigenVectors->Array[(Length * (idx + 1)) + jdx] = temp;
            }
        }
    }

    /* Free the dynamically allocated memory */
    free(TempMat);
    free(work);
}

/* Routine based on ran1 from Numerical Receipes in C */
hysl_float_t RandomNumber (int64_t *const idum) {
    int32_t j;
    int64_t k;
    static int64_t iy = 0;
    static long iv[NTAB];
    hysl_float_t temp;

    if (*idum <= 0 || !iy) {
        if (-(*idum) < 1) {
            *idum = 1;
        } else {
            *idum = -(*idum);
        }

        for (j = NTAB + 7; j >= 0; j--) {
            k = (*idum) / IQ;
            *idum = IA * (*idum - k * IQ) - IR * k;
            if (*idum < 0) {
                *idum += IM;
            }
            if (j < NTAB) {
                iv[j] = *idum;
            }
        }
        iy = iv[0];
    }

    k = (*idum) / (long int) IQ;
    *idum = IA * (*idum - k * IQ) - IR * k;
    if (*idum < 0) {
        *idum += IM;
    }

    j = iy / (long int) NDIV;
    iy = iv[j];
    iv[j] = *idum;
    if ((temp = (hysl_float_t) AM * (hysl_float_t) iy) > (hysl_float_t) RNMX) {
        return (hysl_float_t) RNMX;
    } else {
        return temp;
    }

}

/* Routine based on Gasdev() from Numerical Receipes in C */
hysl_float_t Gaussian_Deviate (const hysl_float_t *const mu, const hysl_float_t *const sigma, int64_t *const idum) {
    static bool iset;
    static hysl_float_t GD_value;
    hysl_float_t fac, rsq, v1, v2;

    if (*idum < 0) { /* Reinitialise */
        iset = false;
    }

    if (iset == false) {
        do {
#if _FLOAT_
	       v1 = 2.0*RandomNumber( idum ) - 1.0;
	       v2 = 2.0*RandomNumber( idum ) - 1.0;
#else
            v1 = 2.0 * RandomNumber(idum) - 1.0;
            v2 = 2.0 * RandomNumber(idum) - 1.0;
#endif
            rsq = v1 * v1 + v2 * v2;
        } while (rsq >= 1.0 || rsq == 0.0);

#if _FLOAT_
	  fac = sqrtf(-2.0*logf(rsq)/rsq);
#else
        fac = sqrt(-2.0 * log(rsq) / rsq);
#endif

        GD_value = (*mu) + (*sigma) * v1 * fac;
        iset = true;
        return (*mu) + (*sigma) * v2 * fac;
    } else {
        iset = false;
        return GD_value;
    }
}

/* Routine based on polint32_t from Numerical Receipes in C */
void Interpolate_Extrapolate (const MatrixVector_t *const X, const MatrixVector_t *const Y, const hysl_float_t x, hysl_float_t *const y,
hysl_float_t *const dy) {

    int32_t i, m, ns = 0;
    hysl_float_t den, dif, dift, ho, hp, w;

    hysl_float_t *c, *d;

    /* Allocate space */
    c = (hysl_float_t*) calloc((size_t) X->Rows, sizeof(hysl_float_t));
    d = (hysl_float_t*) calloc((size_t) X->Rows, sizeof(hysl_float_t));
    if (c == NULL || d == NULL) {
        Print_Header( ERROR);
        fprintf( stderr, "Interpolate_Extrapolate(): Out of memory.\n");
        exit( EXIT_FAILURE);
    }

    /* Find the closest entry in the table and initialise c and d.*/
    dif = hysl_abs(x - X->Array[0]);
    for (i = 0; i < X->Rows; i++) {
        dift = hysl_abs(x - X->Array[i]);
        if (dift < dif) {
            ns = i;
            dif = dift;
        }

        c[i] = Y->Array[i];
        d[i] = Y->Array[i];
    }

    /* Initial approximation to y. */
    *y = Y->Array[ns - 1];

    for (m = 1; m < X->Rows; m++) {
        for (i = 0; i < X->Rows - m; i++) {
            ho = X->Array[i] - x;
            hp = X->Array[i + m] - x;
            den = ho - hp;
            w = c[i + 1] - d[i];
            if (den == 0.0) {
                Print_Header( ERROR);
                fprintf( stderr, "Interpolate_Extrapolate(): Error in the routine.\n");
                exit( EXIT_FAILURE);
            }
            den = w / den;
            d[i] = hp * den;
            c[i] = ho * den;
        }

        *dy = (2 * ns < (X->Rows - m) ? c[ns] : d[ns - 1]);
        *y = *y + *dy;
        /* After each column in the tableau is completed, we decide which correction, c or d, we want to add
         * to our accumulating value of y, i.e., which path to take through the tableau—forking up or
         * down. We do this in such a way as to take the most “straight through the tableau to its apex,
         * updating ns accordingly to keep track of line” routewhere we are. This route keeps the partial
         * approximations centered (insofar as possible) on the target x. The last dy added is thus the
         * error indication.*/

    }

    free(c);
    free(d);
}
