/**
 * \file Netlib.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 10th of August 2011
 *
 * \brief Prototypes of the function for the BLAS and LAPACK routines.
 *
 * This file contains the prototypes of the functions for the BLAS and LAPACK routines. They are only used in
 * the case of the reference libraries or the GOTOBLAS implementation.
 */

#ifndef NETLIB_H_
#define NETLIB_H_

#if _MPI_
#include <mpi.h>
#include "Cblacs.h"
#include "Scalapack_Aux.h"
#endif

#if _MKL_
#else
#define dscal dscal_
#define dcopy dcopy_
#define daxpy daxpy_
#define dgemv dgemv_
#define dsymv dsymv_
#define dspmv dspmv_

#define pdcopy pdcopy_
#define pdaxpy pdaxpy_
#define pdlacpy pdlacpy_
#define pdlascl pdlascl_
#define pdsymv pdsymv_
#endif

/* BLAS routines */
/**
 * \brief \f$\vec x \leftarrow \alpha \vec x\f$
 *
 * \c dscal() computes the product of a vector by a scalar \f$\vec x \leftarrow \alpha \vec x\f$. Uses
 * unrolled loops for increment equal to one.
 *
 * \f[\vec x \leftarrow \alpha \vec x\f]
 *
 * where \f$\alpha\f$ is a scalar and \f$\vec x\f$ is \f$n\f$ element vector.
 *
 * \param[in]     n     On entry, \c n specifies the order of the matrix \c a. \f$n \geq 0\f$.
 * \param[in]     alpha On entry, \c alpha specifies the scalar \f$alpha\f$.
 * \param[in,out] x     Array of dimension \f$ \geq 1 + (n -1)\cdot \|incx\|\f$.
 *                      - Before entry, the incremented array \c x must contain the \c n element vector
 *                        \f$\vec x\f$.
 *                      - On exit, \c x is overwritten by the updated vector \f$\vec x\f$ if \f$n \neq
 *                        0\f$. Otherwise it remains unchanged.
 * \param[in]     incx  On entry, \c incx specifies the increment for the elements of \c x. \f$incx \neq 0\f$.
 *
 * For a complete documentation the reader should refer to \cite BLAS_webpage.
 */
void dscal( int *n, double *alpha, double *x, int *incx );

/**
 * \brief \f$\vec y \leftarrow \vec x\f$
 *
 * \c dcopy() copies a vector,\f$\vec x\f$, to another vector, \f$\vec y\f$ \f$\vec y \leftarrow x\f$. Uses
 * unrolled loops for increments equal to one.
 *
 * \f[\vec y \leftarrow \vec x\f]
 *
 * where \f$\vec x\f$ and \f$\vec y\f$ are \f$n\f$ element vectors.
 *
 * \param[in]  n    On entry, \c n specifies the order of the matrix \c a. \f$n \geq 0\f$.
 * \param[in]  x    Array of dimension \f$ \geq  1 + ( n - 1 )\cdot\|incx\|\f$. Before entry, the incremented
    *               array \c x must contain the \c n element vector \f$\vec x\f$.
 * \param[in]  incx On entry, \c incx specifies the increment for the elements of \c x. \f$incx \neq 0\f$.
 * \param[out] y    Array of dimension \f$ \geq  1 + ( n - 1 )\cdot\|incy\|\f$. If \f$n > 0\f$ it contains a
 *                  copy of the vector \f$\vec x\f$. Otherwise, elements are unaltered.
 * \param[in]  incy On entry, \c incy specifies the increment for the elements of \c y. \f$incy \geq 0\f$.
 *
 * For a complete documentation the reader should refer to \cite BLAS_webpage.
 */
void dcopy( int *n, double *x, int *incx, double *y, int *incy );

/**
 * \brief \f$\vec y \leftarrow \alpha \vec x + \vec y\f$
 *
 * \c daxpy() computes a vector-scalar product, \f$\alpha \vec x\f$ and adds the result to another vector
 * \f$\vec y\f$. Uses unrolled loops for increments equal to one.
 *
 * \f[\vec y \leftarrow \alpha \vec x + \vec y \f]
 *
 * where \f$\alpha\f$ is a scalar and \f$\vec x\f$ and \f$\vec y\f$ are \f$n\f$ element vectors.
 *
 * \param[in]     n     On entry, \c n specifies the order of the matrix \c a. \f$n \geq 0\f$.
 * \param[in]     alpha On entry, \c alpha specifies the scalar \f$alpha\f$.
 * \param[in]     x     Array of dimension \f$ \geq  1 + ( n - 1 )\cdot\|incx\|\f$. Before entry, the
 *                      incremented array \c x must contain the \c n element vector \f$\vec x\f$.
 * \param[in]     incx  On entry, \c incx specifies the increment for the elements of \c x. \f$incx \neq 0\f$.
 * \param[in,out] y     Array of dimension \f$ \geq  1 + ( n - 1 )\cdot\|incy\|\f$.
 *                      - Before entry, the incremented array \c y must contain the \c n element vector
 *                        \f$\vec y\f$.
 *                      - On exit, \c y is overwritten by the updated vector \f$\vec y\f$.
 * \param[in]     incy  On entry, \c incy specifies the increment for the elements of \c y. \f$incy \geq 0\f$.
 *
 * For a complete documentation the reader should refer to \cite BLAS_webpage.
 */
void daxpy( int *n, double *alpha, double *x, int *incx, double *y, int *incy );

/**
 * \brief \f$\vec y \leftarrow \alpha \mathcal A\vec x + \beta \vec y\f$.
 *
 * \c dsymv() Performs the matrix-vector operation:
 *
 * \f[\vec y \leftarrow \alpha \mathcal A\vec x + \beta \vec y \f]
 *
 * where \f$\alpha\f$ and \f$beta\f$ are scalars, \f$\vec x\f$ and \f$\vec y\f$ are \f$n\f$ element vectors
 * and \f$\mathcal A\f$ is an \f$n\f$ by \f$n\f$ symmetric matrix.
 *
 * \param[in]     uplo  On entry uplo On entry, \c uplo specifies whether the upper or lower triangular part
 *                      of the array \c a is to be referenced as follows:
 *                      - \c uplo = 'U' or 'u'   Only the upper triangular part of \c a is to be referenced.
 *                      - \c uplo = 'L' or 'l'   Only the lower triangular part of \c a is to be referenced.
 * \param[in]     n     On entry, \c n specifies the order of the matrix \c a. \f$n \geq 0\f$.
 * \param[in]     alpha On entry, \c alpha specifies the scalar \f$alpha\f$.
 * \param[in]     a     Array of dimension (\c lda,\c n).
 *                      - Before entry with \c uplo = 'U' or 'u', the leading \c n by \c n upper triangular
 *                        part of the array \c a must contain the upper triangular part of the symmetric
 *                        matrix and the strictly lower triangular part of \c a is not referenced.
 *                      - Before entry with \c uplo = 'L' or 'l', the leading \c n by \c n lower triangular
 *                        part of the array \c a must contain the lower triangular part of  the symmetric
 *                        matrix and the strictly upper triangular part of \c a is not referenced.
 * \param[in]     lda     On entry, \c lda specifies the first dimension of \c a as declared in the calling
 *                      (sub) program. \f$lda \geq max(1,n)\f$.
 * \param[in]     x     Array of dimension \f$ \geq  1 + ( n - 1 )\cdot\|incx\|\f$. Before entry, the
 *                      incremented array \c x must contain the \c n element vector \f$\vec x\f$.
 * \param[in]     incx  On entry, \c incx specifies the increment for the elements of \c x. \f$incx \neq 0\f$.
 * \param[in]     beta  On entry, \c beta specifies the scalar \f$\beta\f$. When \c beta is supplied as zero
 *                      then \c y need not be set on input.
 * \param[in,out] y     Array of dimension \f$ \geq 1 + ( n - 1 )\cdot\|incy\|\f$.
 *                      - Before entry, the incremented array \c y must contain the \c n element vector
 *                        \f$\vec y\f$.
 *                      - On exit, \c y is overwritten by the updated vector \f$\vec y\f$.
 * \param[in]     incy  On entry, \c incy specifies the increment for the elements of \c y. \f$incy \geq 0\f$.
 *
 * For a complete documentation the reader should refer to \cite BLAS_webpage.
 */
void dsymv( char *uplo, int *n, double *alpha, double *a, int *lda, double *x, int *incx, double *beta,
	    double *y, int *incy );

/**
 * \brief \f$\vec y \leftarrow \alpha \mathcal A\vec x + \beta \vec y\f$ or \f$\vec y \leftarrow \alpha
 * \mathcal A^T\vec x + \beta \vec y\f$.
 *
 * \c dgemv()  performs one of the matrix-vector operations
 * - \f$\vec y \leftarrow \alpha \mathcal A\vec x + \beta \vec y\f$,   or
 * - \f$\vec y \leftarrow \alpha \mathcal A^T\vec x + \beta \vec y\f$.
 *
 * where \f$\alpha\f$ and \f$beta\f$ are scalars, \f$\vec x\f$ and \f$\vec y\f$ are vectors and \f$\mathcal
 * A\f$ is an \f$m\f$ by \f$n\f$ matrix.
 *
 * \param[in]     trans On entry, \c trans specifies the operation to be performed as follows:
 *                      - \c trans = 'N' or 'n' \f$\vec y \leftarrow \alpha \mathcal A\vec x + \beta \vec
 *                        y\f$.
 *                      - \c trans = 'T' or 't' \f$\vec y \leftarrow \alpha \mathcal A^T\vec x + \beta \vec
 *                        y\f$.
 *                      - \c trans = 'C' or 'c' \f$\vec y \leftarrow \alpha \mathcal A^T\vec x + \beta \vec
 *                        y\f$.
 * \param[in]     m     On entry, \c m specifies the number of rows of the matrix \c a. \f$M \geq 0\f$.
 * \param[in]     n     On entry, \c n specifies the number of columns of the matrix \c a. \f$n \geq 0\f$.
 * \param[in]     alpha On entry, \c alpha specifies the scalar \f$alpha\f$.
 * \param[in]     a     Array of dimension (\c lda,\c n ). Before entry, the leading \c m by \c n part of the
 *                      array \c a must contain the matrix of coefficients.
 * \param[in]     lda   On entry, \c lda specifies the first dimension of \c a as declared in the calling
 *                      (sub) program. \f$(lda \geq max( 1, m )\f$.
 * \param[in]     x     Array of dimension
 *                      - \f$ \geq 1 + ( n - 1 )\cdot\|incx\|\f$ when \c trans = 'N' or 'n' or,
 *                      - \f$ \geq 1 + ( m - 1 )\cdot\|incx\|\f$ otherwise.
 *                      Before entry, the incremented array \c x must contain the vector \f$\vec x\f$.
 * \param[in]     incx  On entry, \c incx specifies the increment for the elements of \c x. \f$incx \neq 0\f$.
 * \param[in]     beta  On entry, \c beta specifies the scalar \f$\beta\f$. When \c beta is supplied as zero
 *                      then \c y need not be set on input.
 * \param[in,out] y     Array of dimension
 *                      - \f$\geq 1 + ( m - 1 )\cdot\|incy\|\f$ when \c trans = 'N' or 'n' or
 *                      - \f$\geq 1 + ( n - 1 )\cdot\|incy\|\f$ otherwise.
 *                      - Before entry with \c beta non-zero, the incremented array \c y must contain the
 *                        vector \f$\vec y\f$.
 *                      - On exit, \c y is overwritten by the updated vector \f$\vec y\f$.
 * \param[in]     incy  On entry, \c incy specifies the increment for the elements of \c y. \f$incy \neq 0\f$.
 *
 * For a complete documentation the reader should refer to \cite BLAS_webpage.
 */
void dgemv_( char *trans, int *m, int *n, double *alpha, double *a, int *lda, double *x, int *incx,
	     double *beta, double *y, int *incy );

/* LAPACK routines */
/**
 * \brief Copies all or part of a two-dimensional matrix \c a to another.
 *
 * dlacpy_() copies all or part of a two-dimensional matrix \c a to another.
 *
 * \param[in]  uplo On entry, \c uplo specifies whether the upper or lower triangular part of the array \c a
 *                  is to be referenced as follows:
 *                  - \c uplo = 'U': Upper triangular part.
 *                  - \c uplo = 'L': Lower triangular part.
 *                  - Otherwise:  All of the matrix \c a
 * \param[in]  m    The number of rows of the matrix \c a. \f$m \geq 0\f$.
 * \param[in]  n    The number of columns of the matrix \c a.  \f$(n \geq 0\f$.
 * \param[in]  a    Array of dimension (\c lda,\c n). The \c m by \c n matrix \c a.  If \c uplo = 'U', only
 *                  the upper triangle or trapezoid is accessed; if \c uplo = 'L', only the lower triangle or
 *                  trapezoid is accessed.
 * \param[in]  lda  The leading dimension of the array \c a.  \f$lda \geq max(1,m)\f$.
 * \param[out] b    Array of dimension (\c ldb,\c n). On exit, \f$\mathcal B = \mathcal A\f$ in the locations
 *                  specified by \c uplo.
 * \param[in]  ldb  The leading dimension of the array B.  \f$ldb \geq max(1,M)\f$.
 *
 * For a complete documentation the reader should refer to \cite LAPACK_webpage \cite LUG.
 */
void dlacpy_( char *uplo, int *m, int *n, double *a, int *lda, double *b, int *ldb );

/**
 * \brief Multiplies the M by N real matrix \c a by a real scalar.
 *
 * \c dlascl_() multiplies the M by N real matrix \c a by the real scalar \f$cto/cfrom\f$. This is done
 * without over/underflow as long as the final result <tt>cto*a(i,j)/cfrom</tt> does not over/underflow. \c
 * type specifies that \c a may be full, upper triangular, lower triangular, upper Hessenberg, or banded.
 *
 * \param[in]     type  Indices the storage type of the input matrix.
 *                      - \c type = 'G': \c a is a full matrix.
 *                      - \c type = 'L': \c a is a lower triangular matrix.
 *                      - \c type = 'U': \c a is an upper triangular matrix.
 *                      - \c type = 'H': \c a is an upper Hessenberg matrix.
 *                      - \c type = 'B': \c a is a symmetric band matrix with lower bandwidth \c KL and upper
 *                        bandwidth \c KU and with the only the lower half stored.
 *                      - \c type = 'Q': \c a is a symmetric band matrix with lower bandwidth \c KL and upper
 *                        bandwidth \c KU and with the only the upper half stored.
 *                      - \c type = 'Z': \c a is a band matrix with lower bandwidth \c kl and upper bandwidth
 *                        \c KU.
 * \param[in]     kl    The lower bandwidth of \c a.  Referenced only if \c type = 'B', 'Q' or 'Z'.
 * \param[in]     ku    The upper bandwidth of \c a.  Referenced only if \c type = 'B', 'Q' or 'Z'.
 * \param[in]     cfrom The divisor. It must be nonzero.
 * \param[in]     cto   The matrix \c a is multiplied by \f$cto/cfrom\f$. \c a(i,j) is computed without
 *                      over/underflow if the final result \f$cto*a(i,j)/cfrom\f$ can be represented without
 *                      over/underflow. \c cfrom must be nonzero.
 * \param[in]     m     The number of rows of the matrix \c a.  \f$m \geq 0\f$.
 * \param[in]     n     The number of columns of the matrix \c a.  \f$n \geq 0\f$.
 * \param[in,out] a     Array of dimension (\c lda,\c n). The matrix to be multiplied by \f$cto/cfrom\f$. See
 *                      \c type for the storage type.
 * \param[in]     lda   The leading dimension of the array \c a.  \f$lda \geq max(1,m)\f$.
 * \param[out]    info  if:
 *                      - \f$info = 0\f$: successful exit.
 *                      - \f$info < 0\f$: if \f$info = -i\f$, the \em i-th argument had an illegal value.
 *
 * For a complete documentation the reader should refer to \cite LAPACK_webpage \cite LUG.
 */
void dlascl_( char *type, int *kl, int *ku, double *cfrom, double *cto, int *m, int *n, double *a, int *lda,
	      int *info );

/**
 * \brief Computes the Cholesky factorization of a real symmetric positive definite matrix \c a.
 *
 * \c dpotrf_() Computes the Cholesky factorization of a real symmetric positive definite matrix \c a.
 *
 * The factorization has the form
 * - \f$\mathcal A = \mathcal U^T\cdot \mathcal U\f$,  if \c uplo = 'U', or
 * - \f$\mathcal A = \mathcal L\cdot \mathcal L^T\f$,  if \c uplo = 'L',
 * where \f$\mathcal U\f$ is an upper triangular matrix and \f$\mathcal L\f$ is lower triangular.
 *
 * This is the block version of the algorithm, calling Level 3 BLAS.
 *
 * \param[in]     uplo Specifies the stored part of the matrix.
 *                     - = 'U':  Upper triangle of \c a is stored;
 *                     - = 'L':  Lower triangle of \c a is stored.
 * \param[in]     n    The order of the matrix \c a.  \f$n \geq 0\f$.
 * \param[in,out] a    Array of dimension (\c lda,\c n).
 *                     - On entry, the symmetric matrix \c a.  If \c uplo = 'U', the leading N-by-N upper
 *                       triangular part of \c a contains the upper triangular part of the matrix \c a, and
 *                       the strictly lower triangular part of \c a is not referenced.  If \c uplo = 'L', the
 *                       leading n-by-n lower triangular part of \c a contains the lower triangular part of
 *                       the matrix \c a, and the strictly upper triangular part of \c a is not referenced.
 *                     - On exit, if \f$info = 0\f$, the factor U or L from the Cholesky factorization
 *                       \f$\mathcal A = \mathcal U^T\cdot U\f$ or \f$\mathcal A = \mathcal L\cdot \mathcal
 *                       L^T\f$.
 * \param[in]     lda  The leading dimension of the array \c a. \f$lda \geq max(1,n)\f$.
 * \param[out]    info if:
 *                     - \f$info = 0\f$: successful exit.
 *                     - \f$info < 0\f$: if \f$info = -i\f$, the \em i-th argument had an illegal value.
 *                     - \f$info > 0\f$: if \f$info = i\f$·, the leading minor of order \em i is not positive
 *                       definite, and the factorization could not be completed.
 *
 * For a complete documentation the reader should refer to \cite LAPACK_webpage \cite LUG.
 */
void dpotrf_( char *uplo, int *n, double *a, int *lda, int *info );

/**
 * \brief Computes the Cholesky factorization of a real symmetric positive definite matrix \c a stored in
 * packed format.
 *
 * \c dpptrf_() computes the Cholesky factorization of a real symmetric positive definite matrix \f$\mathcal
 * A\f$ stored in packed format.
 *
 * The factorization has the form:
 * - \f$\mathcal A = \mathcal U^T\cdot \mathcal U\f$, if \c uplo = 'U', or
 * - \f$\mathcal A = \mathcal L\cdot L^T\f$,  if \c uplo = 'L',
 * where \f$\mathcal U\f$ is an upper triangular matrix and \f$\mathcal L\f$ is lower triangular.
 *
 * \param[in]     uplo Specifies the stored part of the matrix.
 *                     - = 'U':  Upper triangle of \c a is stored;
 *                     - = 'L':  Lower triangle of \c a is stored.
 * \param[in]     n    The order of the matrix \c A.  \f$n \geq 0\f$.
 * \param[in,out] ap   Array of dimension \f$(n*(n+1)/2)\f$.
 *                     - On entry, the upper or lower triangle of the symmetric matrix \c a, packed columnwise
 *                       in a linear array. The \em j-th column of \c a is stored in the array \c ap as
 *                       follows:
 *                       - if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
 *                       - if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
 *                     - On exit, if \c info = 0, the triangular factor U or L from the Cholesky factorization
 *                       A = U**T*U or A = L*L**T, in the same storage format as A.
 * \param[out]    info if:
 *                     - \f$info = 0\f$: successful exit.
 *                     - \f$info < 0\f$: if \f$info = -i\f$, the \em i-th argument had an illegal value.
 *                     - \f$info > 0\f$: if \f$info = i\f$·, the leading minor of order \em i is not positive
 *                       definite, and the factorization could not be completed.
 *
 * For a complete documentation the reader should refer to \cite LAPACK_webpage \cite LUG.
 */
void dpptrf_( char *uplo, int *n, double *ap, int *info );

/**
 * \brief Computes the inverse of a real symmetric positive definite matrix \c a using the Cholesky factorization.
 *
 * \c dpotri_() computes the inverse of a real symmetric positive definite matrix \c a using the Cholesky
 * factorization \f$\mathcal A = \mathcal U**T*\mathcal U\f$ or \f$\mathcal A = \mathcal L*\mathcal L**T\f$
 * computed by \c dpotrf_().
 *
 * \param[in]     uplo On entry, \c uplo specifies whether the upper or lower triangular part of the array \c
 *                     a is to be referenced as follows:
 *                     - \c uplo = 'U':  Upper triangle of \c a is stored;
 *                     - \c uplo = 'L':  Lower triangle of \c a is stored.
 * \param[in]     n    The order of the matrix \c a. \f$n \geq 0\f$.
 * \param[in,out] a    Array of dimension (\c lda,\c n).
 *                     - On entry, the triangular factor U or L from the Cholesky factorization \f$\mathcal A
 *                       = \mathcal U^T\cdot \mathcal U\f$ or \f$\mathcal A = \mathcal L\cdot \mathcal L^T\f$,
 *                       as computed by dpotrf_().
 *                     - On exit, the upper or lower triangle of the (symmetric) inverse of \c a, overwriting
 *                       the input factor U or L.
 * \param[in]     lda  The leading dimension of the array \c a. \f$lda \geq max(1,n)\f$.
 * \param[out]    info if
 *                     - \f$info = 0\f$: successful exit
 *                     - \f$info < 0\f$: if \f$info = -i\f$, the \em i-th argument had an illegal value.
 *                     - \f$info < 0\f$: if \f$info = i\f$, the \f$(i,i)\f$ element of the factor U or L is
 *                       zero, and the inverse could not be computed.
 *
 * For a complete documentation the reader should refer to \cite LAPACK_webpage \cite LUG.
 */
void dpotri_( char *uplo, int *n, double *a, int *lda, int *info );

/**
 * \brief Computes the inverse of a packed symmetric (Hermitian) positive-definite matrix.
 *
 * dpptri_() Computes the inverse of a packed symmetric (Hermitian) positive-definite matrix \c a using the
 * Cholesky factorization \f$\mathcal A = \mathcal U**T*\mathcal U\f$ or \f$\mathcal A = \mathcal L*\mathcal
 * L**T\f$ computed by dpptrf_().
 *
 * \param[in]     uplo On entry, \c uplo specifies whether the upper or lower triangular part of the array \c
 *                     a is to be referenced as follows:
 *                     - \c uplo = 'U': Upper triangular part.
 *                     - \c uplo = 'L': Lower triangular part.
 * \param[in]     n    The order of the matrix. \f$n \geq 0\f$.
 * \param[in,out] ap   Array of dimension at least \f$n(n+1)/2\f$.
 *                     - On entry, the triangular factor U or L from the Cholesky factorization \f$\mathcal A
 *                       = \mathcal U^T\cdot \mathcal U\f$ or \f$\mathcal A = \mathcal L\cdot \mathcal L^T\f$,
 *                       packed columnwise as a linear array. The \em j-th column of \f$\mathcal U\f$ or
 *                       \f$\mathcal L\f$ is stored in the array \c ap as follows:
 *                       - if \c uplo = 'U': \f$\mathcal A (i + (j-1)*j/2) = \mathcal U (i,j)\f$ for \f$1 \leq
 *                         i \leq j\f$
 *                       - if \c uplo = 'L': \f$\mathcal A (i + (j-1)*(2n-j)/2) = \mathcal L (i,j)\f$ for \f$j \leq
 *                         i \leq n\f$
 *                     - On exit, the upper or lower triangle of the (symmetric) inverse of \f$\mathcal A\f$ ,
 *                       overwriting the input factor U or L.
 * \param[out]    info if
 *                     - \f$info = 0\f$: successful exit
 *                     - \f$info < 0\f$: if \f$info = -i\f$, the \em i-th argument had an illegal value.
 *                     - \f$info < 0\f$: if \f$info = i\f$, the \em i-th diagonal element of the Cholesky
 *                       factor U or L is zero, and the inverse could not be computed.
 *
 * For a complete documentation the reader should refer to \cite LAPACK_webpage \cite LUG.
 */
void dpptri_( char *uplo, int *n, double *ap, int *info );


#if _MPI_

/* PBLAS routines */

/**
 * \brief Scales an \em n-element distributed vector \f$sub(\vec x)\f$ by a scalar \f$\alpha\f$ (PBLAS
 * routine).
 *
 * pdscal() scales an \em n-element distributed vector \f$sub(\vec x)\f$ by a scalar \f$\alpha\f$, where
 * \f$sub(\vec x)\f$ denotes:
 * - \f$x(ix,jx:jx+n-1)\f$ if \c incx = \f$m_x\f$.
 * - \f$x(ix:ix+n-1,jx)\f$ if \c incx = 1 and \c incx <> \f$m_x\f$.
 *
 * \param[in]     n     (globalx) The number of components of the distributed vector \f$sub(\vec x)\f$. \f$n
 *                      \geq 0\f$.
 * \param[in]     alpha (global) The scalar used to multiply each component of \f$sub(\vec x)\f$.
 * \param[in,out] x     (local) Array of dimension \f$(jx - 1)m_x + ix + (n-1)|incx|\f$. It contains the
 *                      entries of the distributed vector \f$sub(\vec x)\f$. On output its values are scaled
 *                      by \f$\alpha\f$.
 * \param[in]     ix    (global) The global row index of the submatrix of the distributed matrix \f$\mathcal
 *                      X\f$ to operate on.
 * \param[in]     jx    (global) The global column index of the submatrix of the distributed matrix \f$\mathcal
 *                      X\f$ to operate on.
 * \param[in]     descx The array descripor of the distributed matrix \f$\mathcal X\f$.
 * \param[in]     incx  (global) The global increment for the elements of \f$\mathcal X\f$. Only two values
 *                      are supported in this version, namely 1 and \f$m_x\f$.
 *
 * For a complete documentation the reader should refer to \cite PBLAS_webpage \cite SLUG.
 */
void pdscal( int *n, double *alpha, double *x, int *ix, int *jx, int *descx, int *incx );

/**
 * \brief Copies one distributed vector into another (PBLAS routine).
 *
 * pdcopy() copies one distributed vector into another \f$sub(\vec y) = sub(\vec x)\f$, where:
 * - \f$sub(\vec x)\f$:
 *   - \f$x(ix,jx:jx+n-1)\f$ if \c incx = \f$m_x\f$.
 *   - \f$x(ix:ix+n-1,jx)\f$ if \c incx = 1 and \c incx <> \f$m_x\f$.
 * - \f$sub(\vec y)\f$:
 *   - \f$y(iy,jy:jy+n-1)\f$ if \c incy = \f$m_y\f$.
 *   - \f$y(iy:iy+n-1,jy)\f$ if \c incy = 1 and \c incy <> \f$m_y\f$.
 *
 * \param[in]     n     (global) The length of the distributed vectors to be copied. \f$n \geq 0\f$.
 * \param[in]     x     (local) Array of dimension \f$(jx - 1)m_x + ix + (n-1)|incx|\f$. It contains the
 *                      entries of the distributed vector \f$sub(\vec x)\f$.
 * \param[in]     ix    (global) The global row index of the submatrix of the distributed matrix \f$\mathcal
 *                      X\f$ to operate on.
 * \param[in]     jx    (global) The global column index of the submatrix of the distributed matrix \f$\mathcal
 *                      X\f$ to operate on.
 * \param[in]     descx The array descripor of the distributed matrix \f$\mathcal X\f$.
 * \param[in]     incx  (global) The global increment for the elements of \f$\mathcal X\f$. Only two values
 *                      are supported in this version, namely 1 and \f$m_x\f$.
 * \param[in,out] y     (local) Array of dimension \f$(jy - 1)m_y + iy + (n-1)|incy|\f$. It contains the
 *                      entries of the distributed vector \f$sub(\vec y)\f$. On exit they are overwritten by
 *                      the entries in \f$sub(\vec x)\f$.
 * \param[in]     iy    (global) The global row index of the submatrix of the distributed matrix \f$\mathcal
 *                      Y\f$ to operate on.
 * \param[in]     jy    (global) The global column index of the submatrix of the distributed matrix \f$\mathcal
 *                      Y\f$ to operate on.
 * \param[in]     descy The array descripor of the distributed matrix \f$\mathcal Y\f$.
 * \param[in]     incy  (global) The global increment for the elements of \f$\mathcal Y\f$. Only two values
 *                      are supported in this version, namely 1 and \f$m_y\f$.
 * 
 * For a complete documentation the reader should refer to \cite PBLAS_webpage \cite SLUG.
 */
void pdcopy(int *n, double *x, int *ix, int *jx, int *descx, int *incx, double *y, int *iy, int *jy, int *descy, int *incy );

/**
 * \brief Adds one distributed vector to another (PBLAS routine).
 *
 * pdaxpy() adds one distributed vector to another \f$sub(\vec y) = sub(\vec y) + \alpha \cdot sub(\vec x)\f$, where
 * - \f$sub(\vec x)\f$:
 *   - \f$x(ix,jx:jx+n-1)\f$ if \c incx = \f$m_x\f$.
 *   - \f$x(ix:ix+n-1,jx)\f$ if \c incx = 1 and \c incx <> \f$m_x\f$.
 * - \f$sub(\vec y)\f$:
 *   - \f$y(iy,jy:jy+n-1)\f$ if \c incy = \f$m_y\f$.
 *   - \f$y(iy:iy+n-1,jy)\f$ if \c incy = 1 and \c incy <> \f$m_y\f$.
 *
 * \param[in]     n     (global) The length of the distributed vectors to be added. \f$n \geq 0\f$.
 * \param[in]     alpha (global) The scalar used to multiply each component of \f$sub(\vec x)\f$.
 * \param[in]     x     (local) Array of dimension \f$(jx - 1)m_x + ix + (n-1)|incx|\f$. It contains the
 *                      entries of the distributed vector \f$sub(\vec x)\f$.
 * \param[in]     ix    (global) The global row index of the submatrix of the distributed matrix \f$\mathcal
 *                      X\f$ to operate on.
 * \param[in]     jx    (global) The global column index of the submatrix of the distributed matrix \f$\mathcal
 *                      X\f$ to operate on.
 * \param[in]     descx The array descripor of the distributed matrix \f$\mathcal X\f$.
 * \param[in]     incx  (global) The global increment for the elements of \f$\mathcal X\f$. Only two values
 *                      are supported in this version, namely 1 and \f$m_x\f$.
 * \param[in,out] y     (local) Array of dimension \f$(jy - 1)m_y + iy + (n-1)|incy|\f$. It contains the
 *                      entries of the distributed vector \f$sub(\vec y)\f$. On exit they are overwritten by
 *                      the aforementioned operation.
 * \param[in]     iy    (global) The global row index of the submatrix of the distributed matrix \f$\mathcal
 *                      Y\f$ to operate on.
 * \param[in]     jy    (global) The global column index of the submatrix of the distributed matrix \f$\mathcal
 *                      Y\f$ to operate on.
 * \param[in]     descy The array descripor of the distributed matrix \f$\mathcal Y\f$.
 * \param[in]     incy  (global) The global increment for the elements of \f$\mathcal Y\f$. Only two values
 *                      are supported in this version, namely 1 and \f$m_y\f$. 
 *
 * For a complete documentation the reader should refer to \cite PBLAS_webpage \cite SLUG.
 */
void pdaxpy( int *n, double *alpha, double *x, int *ix, int *jx, int *descx, int *incx, double *y, int *iy,
	     int *jy, int *descy, int *incy);

/**
 * \brief Adds a matrix to another (PBLAS routine).
 *
 * pdgeadd_() adds a general matrix to another:
 *
 * \f[sub(\mathcal C) = \beta\cdot sub(\mathcal C) + \alpha \cdot op(sub(\mathcal A))\f]
 * 
 * where:
 * - \f$sub(\mathcal C)\f$ denotes \f$\mathcal C (ic:ic+m-1,jc:jc+n -1)\f$;
 * - \f$op(\mathcal X)\f$ is either \f$op(\mathcal X) = \mathcal X\f$ or \f$op(\mathcal X) = \mathcal X^T\f$.
 *
 * \param[in]     trans (global) It specifies the operation to be performed as follows:
 *                      - \c trans = 'N' or 'n': \f$op(sub(\mathcal A)) = sub(\mathcal A)\f$.
 *                      - \c trans = 'T' or 't': \f$op(sub(\mathcal A)) = sub(\mathcal A)^T\f$.
 *                      - \c trans = 'C' or 'c': \f$op(sub(\mathcal A)) = sub(\mathcal A)^T\f$.
 * \param[in]     m     (global) The number of rows of the distributed submatrix \f$sub(\mathcal C)\f$ and the
 *                      number of columns of the submatrix \f$sub(\mathcal A)\f$. \f$m \geq 0\f$. 
 * \param[in]     n     (global) The number of columns of the distributed submatrix \f$sub(\mathcal
 *                      X)\f$ and the number of rows of the submatrix \f$sub(\mathcal A)\f$. \f$n \geq 0\f$. 
 * \param[in]     alpha (global) The scalar \f$\alpha\f$. If it is zero, then the entries of the submatrix
 *                      \f$sub(\mathcal A)\f$ neet not be set on input.
 * \param[in]     a     (local) Array of dimension \f$(lld_a, LOCq(ja+m-1))\f$ and contains the local entries
 *                      of the distributed matrix \f$\mathcal A\f$ to operate on.
 * \param[in]     ia    (global) The global row index of the submatrix of the distributed matrix \f$\mathcal
 *                      A\f$ to operate on.
 * \param[in]     ja    (global) The global column index of the submatrix of the distributed matrix \f$\mathcal
 *                      A\f$ to operate on.
 * \param[in]     desca Array descriptor of the distributed matrix \f$\mathcal A\f$.
 * \param[in]     beta  (global) The scalar \f$\beta\f$. If it is zero, then the entries of the submatrix
 *                      \f$sub(\mathcal C)\f$ neet not be set on input.
 * \param[in,out] c     (local) Array of dimension \f$(lld_c, LOCq(jc+n-1))\f$ and contains the local entries
 *                      of the distributed matrix \f$\mathcal C\f$ to operate on. On output the entries
 *                      contain the updated values according to the performed operation.
 * \param[in]     ic    (global) The global row index of the submatrix of the distributed matrix \f$\mathcal
 *                      C\f$ to operate on.
 * \param[in]     jc    (global) The global column index of the submatrix of the distributed matrix \f$\mathcal
 *                      C\f$ to operate on.
 * \param[in]     descc The array descripor of the distributed matrix \f$\mathcal C\f$.

 *
 * For a complete documentation the reader should refer to \cite PBLAS_webpage \cite SLUG.
 */
void pdgeadd_(char *trans, int *m, int *n, double *alpha, double *a, int *ia, int *ja, int *desca, double *beta,
	      double *c, int *ic, int *jc, int *descc);

/**
 * \brief Adds a trapezoidal matrix to another (PBLAS routine).
 *
 * pdtradd_() adds a trapezoidal matrix to another:
 *
 * \f[sub(\mathcal C) = \beta\cdot sub(\mathcal C) + \alpha \cdot op(sub(\mathcal A))\f]
 * 
 * where:
 * - \f$sub(\mathcal C)\f$ denotes \f$\mathcal C (ic:ic+m-1,jc:jc+n -1)\f$;
 * - \f$op(\mathcal X)\f$ is either \f$op(\mathcal X) = \mathcal X\f$ or \f$op(\mathcal X) = \mathcal X^T\f$.
 *
 * Therefore \f$op(sub(\mathcal A))\f$ denotes:
 * - \f$\mathcal A (ia:ia+m-1,ja:ja+n-1)\f$ if \c trans = 'N' or 'n',
 * - \f$\mathcal A (ia:ia+n-1,ja:ja+m-1)\f$ if \c trans = 'T' or 't',
 * - \f$\mathcal A (ia:ia+n-1,ja:ja+m-1)\f$ if \c trans = 'C' or 'c',
 *
 * \f$sub(\mathcal C)\f$ and \f$op(sub(\mathcal A))\f$ are \em m by \em n  upper or lower trapezoidal
 * submatrices.
 *
 * \param[in]     uplo  (global) Specifies whether the upper or lower triangular part of the distributed
 *                      matrix \f$\mathcal C\f$ is to be referenced.
 *                      - \c uplo = 'U' or 'u': Only the upper triangular part of the triangular submatrix
 *                        \f$sub(\mathcal C)\f$ is to be referenced.
 *                      - \c uplo = 'L' or 'l': Only the lower triangular part of the triangular submatrix
 *                        \f$sub(\mathcal C)\f$ is to be referenced.
 * \param[in]     trans (global) It specifies the operation to be performed as follows:
 *                      - \c trans = 'N' or 'n': \f$op(sub(\mathcal A)) = sub(\mathcal C)\f$.
 *                      - \c trans = 'T' or 't': \f$op(sub(\mathcal A)) = sub(\mathcal C)^T\f$.
 *                      - \c trans = 'C' or 'c': \f$op(sub(\mathcal A)) = sub(\mathcal C)^T\f$.
 * \param[in]     m     (global) The number of rows of the distributed submatrix \f$sub(\mathcal C)\f$ and the
 *                      number of columns of the submatrix \f$sub(\mathcal A)\f$. \f$m \geq 0\f$. 
 * \param[in]     n     (global) The number of columns of the distributed submatrix \f$sub(\mathcal
 *                      X)\f$ and the number of rows of the submatrix \f$sub(\mathcal A)\f$. \f$n \geq 0\f$. 
 * \param[in]     alpha (global) The scalar \f$\alpha\f$. If it is zero, then the entries of the submatrix
 *                      \f$sub(\mathcal A)\f$ neet not be set on input.
 * \param[in]     a     (local) Array of dimension \f$(lld_a, LOCq(ja+m-1))\f$ and contains the local entries
 *                      of the distributed matrix \f$\mathcal A\f$ to operate on.
 * \param[in]     ia    (global) The global row index of the submatrix of the distributed matrix \f$\mathcal
 *                      A\f$ to operate on.
 * \param[in]     ja    (global) The global column index of the submatrix of the distributed matrix \f$\mathcal
 *                      A\f$ to operate on.
 * \param[in]     desca Array descriptor of the distributed matrix \f$\mathcal A\f$.
 * \param[in]     beta  (global) The scalar \f$\beta\f$. If it is zero, then the entries of the submatrix
 *                      \f$sub(\mathcal C)\f$ neet not be set on input.
 * \param[in,out] c     (local) Array of dimension \f$(lld_c, LOCq(jc+n-1))\f$ and contains the local entries
 *                      of the distributed matrix \f$\mathcal C\f$ to operate on. On output the entries
 *                      contain the updated values according to the performed operation.
 * \param[in]     ic    (global) The global row index of the submatrix of the distributed matrix \f$\mathcal
 *                      C\f$ to operate on.
 * \param[in]     jc    (global) The global column index of the submatrix of the distributed matrix \f$\mathcal
 *                      C\f$ to operate on.
 * \param[in]     descc The array descripor of the distributed matrix \f$\mathcal C\f$.
 *
 * For a complete documentation the reader should refer to \cite PBLAS_webpage \cite SLUG.
 */
void pdtradd_( char *uplo, char *trans, int *m, int *n, double *alpha, double *a, int *ia, int *ja, int *desca,
	       double *beta, double *c, int *ic, int *jc, int *descc );

/**
 * \brief Performs a distributed matrix-vector operation (PBLAS routine).
 *
 * pdgemv_() performs one of the distributed matrix-vector operations:
 * - \f$sub(\vec y) = \alpha \cdot sub(\mathcal A)*sub(\vec x) + \beta \cdot sub(\vec y)\f$, or
 * - \f$sub(\vec y) = \alpha \cdot sub(\mathcal A)^T*sub(\vec x) + \beta \cdot sub(\vec y)\f$
 *
 * where:
 * - \f$sub(\mathcal A)\f$ denotes \f$\mathcal A (ia:ia+m-1,ja:ja+n -1)\f$;
 * - \f$sub(\vec x)\f$ denotes if \c trans = 'N';
 *   - \f$x(ix,jx:jx+n-1)\f$ if \c incx = \f$m_x\f$;
 *   - \f$x(ix:ix+n-1,jx)\f$ if \c incx = 1 and \c incx <> \f$m_x\f$.
 *
 *   else
 *   - \f$x(ix,jx:jx+m-1)\f$ if \c incx = \f$m_x\f$;
 *   - \f$x(ix:ix+m-1,jx)\f$ if \c incx = 1 and \c incx <> \f$m_x\f$.
 * - \f$sub(\vec y)\f$ denotes if \c trans = 'N';
 *   - \f$y(iy,jy:jy+n-1)\f$ if \c incy = \f$m_x\f$;
 *   - \f$y(iy:iy+n-1,jy)\f$ if \c incy = 1 and \c incy <> \f$m_y\f$.
 *
 *   else
 *   - \f$y(iy,jy:jy+m-1)\f$ if \c incy = \f$m_y\f$.
 *   - \f$y(iy:iy+m-1,jy)\f$ if \c incy = 1 and \c incx <> \f$m_y\f$.
 *
 * \param[in]     trans (global) It specifies the operation to be performed as follows:
 *                      - \c trans = 'N' or 'n': \f$sub(\vec y) = \alpha \cdot sub(\mathcal A)*sub(\vec x) +
 *                        \beta \cdot sub(\vec y)\f$
 *                      - \c trans = 'T' or 't': \f$sub(\vec y) = \alpha \cdot sub(\mathcal A)^T*sub(\vec x) +
 *                        \beta \cdot sub(\vec y)\f$
 *                      - \c trans = 'C' or 'c': \f$sub(\vec y) = \alpha*conj(sub(\mathcal A))^T*sub(\vec x) +
 *                        \beta \cdot sub(\vec y)\f$
 * \param[in]     m     (global) The number of rows of the distributed submatrix \f$sub(\mathcal A)\f$. \f$m
 *                      \geq 0\f$. 
 * \param[in]     n     (global) The number of columns of the distributed submatrix \f$sub(\mathcal
 *                      A)\f$. \f$n \geq 0\f$. 
 * \param[in]     alpha (global) The scalar \f$\alpha\f$.
 * \param[in]     a     (local) Array of dimension \f$(lld, LOCq(ja+n-1))\f$ and contains the local entries
 *                      of the distributed matrix \f$\mathcal A\f$ to operate on.
 * \param[in]     ia    (global) The global row index of the submatrix of the distributed matrix \f$\mathcal
 *                      A\f$ to operate on.
 * \param[in]     ja    (global) The global column index of the submatrix of the distributed matrix \f$\mathcal
 *                      A\f$ to operate on.
 * \param[in]     desca Array descriptor of the distributed matrix \f$\mathcal A\f$.
 * \param[in]     x     (local) Array of dimension at least:
 *                      - \f$(jx - 1)m_x + ix + (n-1)|incx|\f$ if \c trans = 'N';
 *                      - \f$(jx - 1)m_x + ix + (n-1)|incx|\f$ otherwise.
 *
 *                      It contains the entries of the distributed vector \f$sub(\vec x)\f$.
 * \param[in]     ix    (global) The global row index of the submatrix of the distributed matrix \f$\mathcal
 *                      X\f$ to operate on.
 * \param[in]     jx    (global) The global column index of the submatrix of the distributed matrix \f$\mathcal
 *                      X\f$ to operate on.
 * \param[in]     descx The array descripor of the distributed matrix \f$\mathcal X\f$.
 * \param[in]     incx  (global) The global increment for the elements of \f$\mathcal X\f$. Only two values
 *                      are supported in this version, namely 1 and \f$m_x\f$.
 * \param[in]     beta  (global) The scalar \f$\beta\f$.
 * \param[in,out] y     (local) Array of dimension  at least:
 *                      - \f$(jy - 1)m_y + iy + (n-1)|incy|\f$ if \c trans = 'N';
 *                      - \f$(jy - 1)m_y + iy + (m-1)|incy|\f$ otherwise.
 *
 *                      It contains the entries of the distributed vector \f$sub(\vec y)\f$. On exit they are
 *                      overwritten by the aforementioned operation.
 * \param[in]     iy    (global) The global row index of the submatrix of the distributed matrix \f$\mathcal
 *                      Y\f$ to operate on.
 * \param[in]     jy    (global) The global column index of the submatrix of the distributed matrix \f$\mathcal
 *                      Y\f$ to operate on.
 * \param[in]     descy The array descripor of the distributed matrix \f$\mathcal Y\f$.
 * \param[in]     incy  (global) The global increment for the elements of \f$\mathcal Y\f$. Only two values
 *                      are supported in this version, namely 1 and \f$m_y\f$. 
 *
 * For a complete documentation the reader should refer to \cite PBLAS_webpage \cite SLUG.
 */
void pdgemv_( char *trans, int *m, int *n, double *alpha, double *a, int *ia, int *ja, int *desca, double *x,
	      int *ix, int *jx, int *descx, int *incx, double *beta, double *y, int *iy, int *jy, int *descy,
	      int *incy );

/**
 * \brief Performs a distributed matrix-vector operation (symmetric, PBLAS routine).
 *
 * pdsymv() performs the distributed matrix-vector operation \f$sub(\vec y) = \alpha \cdot sub(\mathcal A)*sub(\vec
 * x) + \beta \cdot sub(\vec y)\f$
 *
 * where:
 * - \f$sub(\mathcal A)\f$ denotes \f$\mathcal A (ia:ia+m-1,ja:ja+n -1)\f$;
 * - \f$sub(\vec x)\f$ denotes:
 *   - \f$x(ix,jx:jx+n-1)\f$ if \c incx = \f$m_x\f$;
 *   - \f$x(ix:ix+n-1,jx)\f$ if \c incx = 1 and \c incx <> \f$m_x\f$.
 * - \f$sub(\vec y)\f$ denotes:
 *   - \f$y(iy,jy:jy+n-1)\f$ if \c incy = \f$m_x\f$;
 *   - \f$y(iy:iy+n-1,jy)\f$ if \c incy = 1 and \c incy <> \f$m_y\f$.
 *
 * \param[in]     uplo  (global) Specifies whether the upper or lower triangular part of the distributed
 *                      matrix \f$\mathcal A\f$ is to be referenced.
 *                      - \c uplo = 'U' or 'u': Only the upper triangular part of \f$sub(\mathcal A)\f$ is to
 *                        be referenced.
 *                      - \c uplo = 'L' or 'l': Only the lower triangular part of \f$sub(\mathcal A)\f$ is to
 *                        be referenced.
 * \param[in]     n     (global) The order of the distributed submatrix \f$sub(\mathcal A)\f$. \f$n \geq 0\f$. 
 * \param[in]     alpha (global) The scalar \f$\alpha\f$.
 * \param[in]     a     (local) Array of dimension \f$(lld, LOCq(ja+n-1))\f$ and contains the local entries
 *                      of the distributed matrix \f$\mathcal A\f$ to operate on.
 * \param[in]     ia    (global) The global row index of the submatrix of the distributed matrix \f$\mathcal
 *                      A\f$ to operate on.
 * \param[in]     ja    (global) The global column index of the submatrix of the distributed matrix \f$\mathcal
 *                      A\f$ to operate on.
 * \param[in]     desca Array descriptor of the distributed matrix \f$\mathcal A\f$.
 * \param[in]     x     (local) Array of dimension at least \f$(jx - 1)m_x + ix + (n-1)|incx|\f$. It
 *                      contains the entries of the distributed vector \f$sub(\vec x)\f$.
 * \param[in]     ix    (global) The global row index of the submatrix of the distributed matrix \f$\mathcal
 *                      X\f$ to operate on.
 * \param[in]     jx    (global) The global column index of the submatrix of the distributed matrix \f$\mathcal
 *                      X\f$ to operate on.
 * \param[in]     descx The array descripor of the distributed matrix \f$\mathcal X\f$.
 * \param[in]     incx  (global) The global increment for the elements of \f$\mathcal X\f$. Only two values
 *                      are supported in this version, namely 1 and \f$m_x\f$.
 * \param[in]     beta  (global) The scalar \f$\beta\f$.
 * \param[in,out] y     (local) Array of dimension  at least \f$(jy - 1)m_y + iy + (n-1)|incy|\f$. It
 *                      contains the entries of the distributed vector \f$sub(\vec y)\f$. On exit they are
 *                      overwritten by the aforementioned operation.
 * \param[in]     iy    (global) The global row index of the submatrix of the distributed matrix \f$\mathcal
 *                      Y\f$ to operate on.
 * \param[in]     jy    (global) The global column index of the submatrix of the distributed matrix \f$\mathcal
 *                      Y\f$ to operate on.
 * \param[in]     descy The array descripor of the distributed matrix \f$\mathcal Y\f$.
 * \param[in]     incy  (global) The global increment for the elements of \f$\mathcal Y\f$. Only two values
 *                      are supported in this version, namely 1 and \f$m_y\f$. 
 *
 * For a complete documentation the reader should refer to \cite PBLAS_webpage \cite SLUG.
 */
void pdsymv( char *uplo, int *n, double *alpha, double *a, int *ia, int *ja, int *desca, double *x, int *ix,
	     int *jx, int *descx, int *incx, double *beta, double *y, int *iy, int *jy, int *descy, int *incy );

/* ScaLAPACK routines */
/**
 * \brief Copies all or part of a distributed matrix \f$\mathcal A\f$ to another distributed matrix
 * \f$\mathcal B\f$ without communication. (ScaLAPACK routine).
 *
 * pdlacpy_() copies all or part of a distributed matrix \f$\mathcal A\f$ to another distributed matrix
 * \f$\mathcal B\f$ without communication. The copy is performed locally \f$sub(\mathcal A) = sub(\mathcal
 * B)\f$, where:
 * - \f$sub(\mathcal A)\f$ denotes \f$\mathcal A (ia:ia+m-1,ja:ja+n -1)\f$;
 * - \f$sub(\mathcal B)\f$ denotes \f$\mathcal B (ib:ib+m-1,jb:jb+n -1)\f$;
 *
 * \param[in] uplo  (global) Specifies whether the upper or lower triangular part of the distributed matrix
 *                  \f$\mathcal A\f$ is to be referenced.
 *                   - \c uplo = 'U' or 'u': Only the upper triangular part of \f$sub(\mathcal A)\f$ is to be
 *                    referenced.
 *                   - \c uplo = 'L' or 'l': Only the lower triangular part of \f$sub(\mathcal A)\f$ is to be
 *                    referenced.
 * \param[in]  m     (global) The number of rows ofthe distributed submatrix \f$sub(\mathcal A)\f$. \f$m \geq
 *                  0\f$. 
 * \param[in]  n     (global) The number of columns ofthe distributed submatrix \f$sub(\mathcal A)\f$. \f$n
 *                  \geq 0\f$. 
 * \param[in]  a     (local) Array of dimension \f$(lld_a, LOCc(ja+n-1))\f$ and contains the local entries of
 *                   the distributed matrix \f$\mathcal A\f$.
 * \param[in]  ia    (global) The global row index of the submatrix of the distributed matrix \f$\mathcal A\f$
 *                   to operate on.
 * \param[in]  ja    (global) The global column index of the submatrix of the distributed matrix \f$\mathcal
 *                   A\f$ to operate on.
 * \param[in]  desca Array descriptor of the distributed matrix \f$\mathcal A\f$.
 * \param[out] b     (local) Array of dimension \f$(lld_b, LOCc(jb+n-1))\f$ and contains on exit the local
 *                   entries of the distributed matrix \f$\mathcal B\f$ as follows:
 *                   - if \c uplo = 'U' or 'u': \f$\mathcal B (ib+i-1,jb+j-1) = \mathcal A (ia+i-1,ja+j-1)\f$,
 *                     \f$1 \leq i \leq j\f$, \f$1 \leq j \leq N\f$;
 *                   - if \c uplo = 'L' or 'l', \f$\mathcal B (ib+i-1,jb+j-1) = \mathcal A (ia+i-1,ja+j-1)\f$,
 *                     \f$j\leq i \leq M\f$, \f$1 \leq j \leq N\f$;
 *                   - otherwise, \f$\mathcal B (ib+i-1,jb+j-1) = \mathcal A (ia+i-1,ja+j-1)\f$, \f$j\leq i
 *                     \leq M\f$, \f$1 \leq j \leq N\f$;
 * \param[in]  ib    (global) The global row index of the submatrix of the distributed matrix \f$\mathcal
 *                   B\f$.
 * \param[in]  jb    (global) The global column index of the submatrix of the distributed matrix \f$\mathcal
 *                   B\f$.
 * \param[in]  descb The array descripor of the distributed matrix \f$\mathcal X\f$.
 *
 * For a complete documentation the reader should refer to \cite ScaLAPACK_webpage \cite SLUG.
 */
void pdlacpy_( char *uplo, int *m, int *n, double *a, int *ia, int *ja, int *desca, double *b, int *ib, int *jb,
	       int *descb);

/**
 * \brief Multiplies a distributed matrix \f$sub(\mathcal A)\f$ with a scalar (ScaLAPACK routine).
 *
 * pdlascl_() multiplies a distributed matrix \f$sub(\mathcal A)\f$ with a scalar \c cto/cfrom without
 * over/underflow as long as the result \f$cto \cdot \mathcal A(i,j)/cfrom\f$ does not overflow. \c type specifies that
 * \f$sub(\mathcal A)\f$ may be full, upper triangular, lower triangular or upper Hessenberg.
 * 
 * \param[in]     type  (global) Indices the storage type of the input matrix.
 *                      - \c type = 'G': \f$sub(\mathcal A)\f$ is a full matrix.
 *                      - \c type = 'L': \f$sub(\mathcal A)\f$ is a lower triangular matrix.
 *                      - \c type = 'U': \f$sub(\mathcal A)\f$ is an upper triangular matrix.
 *                      - \c type = 'H': \f$sub(\mathcal A)\f$ is an upper Hessenberg matrix.
 * \param[in]     uplo  (global) Specifies whether the upper or lower triangular part of the distributed
 *                      matrix \f$\mathcal A\f$ is to be referenced.
 *                      - \c uplo = 'U' or 'u': Only the upper triangular part of \f$sub(\mathcal A)\f$ is to
 *                        be referenced.
 *                      - \c uplo = 'L' or 'l': Only the lower triangular part of \f$sub(\mathcal A)\f$ is to
 *                        be referenced.

 * \param[in]     cfrom (global) The divisor. It must be nonzero.
 * \param[in]     cto   (global) The distributed matrix \f$sub(\mathcal A)\f$a is multiplied by
 *                      \f$cto/cfrom\f$. \f$\mathcal A(i,j)\f$ is computed without over/underflow if the final
 *                      result \f$cto\cdot \mathcal A(i,j)/cfrom\f$ can be represented without
 *                      over/underflow.
 * \param[in]     m     (global) The number of rows ofthe distributed submatrix \f$sub(\mathcal A)\f$. \f$m
 *                      \geq 0\f$. 
 * \param[in]     n     (global) The number of columns ofthe distributed submatrix \f$sub(\mathcal A)\f$. \f$n
 *                      \geq 0\f$. 
 * \param[in,out] a     (local) Array of dimension \f$(lld_a, LOCc(ja+n-1))\f$ and contains the local entries
 *                      of the distributed matrix \f$\mathcal A\f$. On exit it contains the entries multiplied
 *                      by \c cto/cfrom.
 * \param[in]     ia    (global) The global row index of the submatrix of the distributed matrix \f$\mathcal
 *                      A\f$ to operate on.
 * \param[in]     ja    (global) The global column index of the submatrix of the distributed matrix \f$\mathcal
 *                      A\f$ to operate on.
 * \param[in]     desca Array descriptor of the distributed matrix \f$\mathcal A\f$.
 * \param[out]    info  (global) if:
 *                      - \f$info = 0\f$: successful exit.
 *                      - \f$info < 0\f$: if \f$info = -i\f$, the \em i-th argument had an illegal value.
 *
 * For a complete documentation the reader should refer to \cite ScaLAPACK_webpage \cite SLUG.
 */
void pdlascl_( char *type, double *cfrom, double *cto, int *m, int *n, double *a, int *ia, int *ja,
	       int *desca, int *info);

/**
 * \brief Computes the Cholesky factorisation of an n-by-n real symmetric positive distributed
 * matrix.(ScaLAPACK routine).
 *
 * pdpotrf_() computes the Cholesky factorisation of an n-by-n real symmetric positive distributed
 * matrix \f$sub(\mathcal A)\f$ denoted by \f$\mathcal A (ia:ia+m-1,ja:ja+n -1)\f$. The factorization has the
 * form:
 * - \f$sub(\mathcal A) = \mathcal U^T\cdot \mathcal U\f$ if \c uplo = 'U';
 * - \f$sub(\mathcal A) = \mathcal L\cdot \mathcal L^T\f$ if \c uplo = 'L';
 *
 * where \f$\mathcal U\f$ is an upper triangular matrix and \f$\mathcal L\f$ is lower triangular.
 *
 * \param[in]     uplo  (global) Specifies whether the upper or lower triangular part of the distributed
 *                      matrix \f$\mathcal A\f$ is to be referenced.
 *                      - \c uplo = 'U' or 'u': Only the upper triangular part of \f$sub(\mathcal A)\f$ is to
 *                        be referenced.
 *                      - \c uplo = 'L' or 'l': Only the lower triangular part of \f$sub(\mathcal A)\f$ is to
 *                        be referenced.
 * \param[in]     n     (global) The order of the distributed submatrix \f$sub(\mathcal A)\f$. \f$n \geq
 *                      0\f$. 
 * \param[in,out] a     (local) Array of dimension \f$(lld_a, LOCc(ja+n-1))\f$ and contains the local entries
 *                      of the distributed matrix \f$\mathcal A\f$ to be factored.
 *                     - On entry, the symmetric matrix \c a.  If \c uplo = 'U', the leading N-by-N upper
 *                       triangular part of \c a contains the upper triangular part of the matrix \c a, and
 *                       the strictly lower triangular part of \c a is not referenced.  If \c uplo = 'L', the
 *                       leading n-by-n lower triangular part of \c a contains the lower triangular part of
 *                       the matrix \c a, and the strictly upper triangular part of \c a is not referenced.
 *                     - On exit, if \f$info = 0\f$, the factor U or L from the Cholesky factorization
 *                       \f$sub(\mathcal A) = \mathcal U^T\cdot U\f$ or \f$sub(\mathcal A) = \mathcal L\cdot \mathcal
 *                       L^T\f$.
 * \param[in]     ia    (global) The global row index of the submatrix of the distributed matrix \f$\mathcal
 *                      A\f$ to operate on.
 * \param[in]     ja    (global) The global column index of the submatrix of the distributed matrix \f$\mathcal
 *                      A\f$ to operate on.
 * \param[in]     desca Array descriptor of the distributed matrix \f$\mathcal A\f$.
 * \param[out]    info  (global) if:
 *                     - \f$info = 0\f$: successful exit.
 *                     - \f$info < 0\f$: if \f$info = -i\f$, the \em i-th argument had an illegal value.
 *                     - \f$info > 0\f$: if \f$info = k\f$·, the leading minor of order \em k, \f$\mathcal A
 *                       (ia:ia+k-1,ja:ja+k-1)\f$ is not positive definite, and the factorization could not be
 *                       completed.
 *
 * For a complete documentation the reader should refer to \cite ScaLAPACK_webpage \cite SLUG.
 */
void pdpotrf_( char *uplo, int *n, double *a, int *ia, int *ja, int *desca, int *info );

/**
 * \brief Computes the inverse off a real symmetric positive definite distribured matrix (ScaLAPACK routine).
 *
 * pdpotri_() computes the inverse off a real symmetric positive definite distribured matrix \f$sub(\mathcal
 * A)\f$ denoted by \f$\mathcal A (ia:ia+m-1,ja:ja+n -1)\f$ using the Cholesky factorization \f$sub(\mathcal
 * A) = \mathcal U**T*\mathcal U\f$ or \f$sub(\mathcal A)= \mathcal L*\mathcal L**T\f$ computed by pdpotrf_().
 *
 * \param[in]     uplo  (global) Specifies whether the upper or lower triangular part of the distributed
 *                      matrix \f$\mathcal A\f$ is to be referenced.
 *                      - \c uplo = 'U' or 'u': Only the upper triangular part of \f$sub(\mathcal A)\f$ is to
 *                        be referenced.
 *                      - \c uplo = 'L' or 'l': Only the lower triangular part of \f$sub(\mathcal A)\f$ is to
 *                        be referenced.
 * \param[in]     n     (global) The order of the distributed submatrix \f$sub(\mathcal A)\f$. \f$n \geq
 *                      0\f$. 
 * \param[in,out] a     (local) Array of dimension \f$(lld_a, LOCc(ja+n-1))\f$ and contains the local entries
 *                      of the distributed matrix \f$\mathcal A\f$ to be factored. On entry, the local pieces
 *                      of the triangular factor \f$\mathcal U\f$ or \f$\mathcal L\f$ from the Cholesky
 *                      factorization o the distributed matrix \f$sub(\mathcal A) = \mathcal U**T*\mathcal U\f$
 *                      or \f$sub(\mathcal A)= \mathcal L*\mathcal L**T\f$ as computed by pdpotrf_(). On exit
 *                      they are overwritten by the inverse of \f$sub(\mathcal A)\f$.
 * \param[in]     ia    (global) The global row index of the submatrix of the distributed matrix \f$\mathcal
 *                      A\f$ to operate on.
 * \param[in]     ja    (global) The global column index of the submatrix of the distributed matrix \f$\mathcal
 *                      A\f$ to operate on.
 * \param[in]     desca Array descriptor of the distributed matrix \f$\mathcal A\f$.
 * \param[out]    info  (global) if:
 *                     - \f$info = 0\f$: successful exit.
 *                     - \f$info < 0\f$: if \f$info = -i\f$, the \em i-th argument had an illegal value.
 *                     - \f$info > 0\f$: if \f$info = i\f$, the \f$(i,i)\f$ element of the factor \f$\mathcal
 *                      U\f$ or \f$\mathcal L\f$ is zero, and the inverse could not be computed.
 *
 * For a complete documentation the reader should refer to \cite ScaLAPACK_webpage \cite SLUG.
 */
void pdpotri_( char *uplo, int *n, double *a, int *ia, int *ja, int *desca, int *info );

#endif /* _MPI_ */

#endif /* NETILB_H_ */
