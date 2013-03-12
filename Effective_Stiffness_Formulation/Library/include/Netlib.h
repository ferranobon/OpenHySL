/**
 * \file Netlib.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 10th of August 2011
 *
 * \brief Prototypes of the function for the BLAS and LAPACK routines.
 *
 * This file contains the prototypes of the functions for the BLAS and LAPACK routines. They are only used in the case of
 * the reference libraries or the GOTOBLAS implementation.
 */

#ifndef NETLIB_H_
#define NETLIB_H_

#if _MPI_
#include <mpi.h>
#endif

/* BLAS routines */
/**
 * \brief \f$\vec x \leftarrow \alpha \vec x\f$
 *
 * \c dscal computes the product of a vector by a scalar \f$\vec x \leftarrow \alpha \vec x\f$. Uses unrolled loops for increment equal to one.
 *
 * where \f$\alpha\f$ is a scalar and \f$\vec x\f$ is \f$n\f$ element vector.
 *
 * \param[in] n On entry, \c n specifies the order of the matrix \c a. \f$n \geq 0\f$.
 * \param[in] alpha On entry, \c alpha specifies the scalar \f$alpha\f$.
 * \param[in,out] x Array of dimension \f$ \geq 1 + (n -1)\cdot \|incx\|\f$.
 * - Before entry, the incremented array \c x must contain the \c n element vector \f$\vec x\f$.
 * - On exit, \c x is overwritten by the updated vector \f$\vec x\f$ if \f$n \neq 0\f$. Otherwise it remains unchanged.
 * \param[in] incx On entry, \c incx specifies the increment for the elements of \c x. \f$incx \neq 0\f$.
 */
void dscal( int *n, double *alpha, double *x, int *incx );

/**
 * \brief \f$\vec y \leftarrow x\f$
 *
 * \c dcopy copies a vector,\f$\vec x\f$, to another vector, \f$\vec y\f$ \f$\vec y \leftarrow x\f$. Uses unrolled loops for increments equal to one.
 *
 * where \f$\vec x\f$ and \f$\vec y\f$ are \f$n\f$ element vectors.
 *
 * \param[in] n On entry, \c n specifies the order of the matrix \c a. \f$n \geq 0\f$.
 * \param[in] x Array of dimension \f$ \geq  1 + ( n - 1 )\cdot\|incx\|\f$. Before entry, the incremented array \c x must contain the \c n element vector \f$\vec x\f$.
 * \param[in] incx On entry, \c incx specifies the increment for the elements of \c x. \f$incx \neq 0\f$.
 * \param[out] y Array of dimension \f$ \geq  1 + ( n - 1 )\cdot\|incy\|\f$.
 * - If \f$n > 0\f$ it contains a copy of the vector \f$\vec x\f$. Otherwise, elements are unaltered.
 * \param[in] incy On entry, \c incy specifies the increment for the elements of \c y. \f$incy \geq 0\f$.
 */
void dcopy( int *n, double *x, int *incx, double *y, int *incy );

/**
 * \brief \f$\vec y \leftarrow \alpha \vec x + \vec y\f$
 *
 * \c daxpy computes a vector-scalar product, \f$\alpha \vec x\f$ and adds the result to another vector \f$\vec y\f$. Uses unrolled loops for increments equal to one.
 *
 * \f[y \leftarrow \alpha \vec x + \vec y \f]
 *
 * where \f$\alpha\f$ is a scalar and \f$\vec x\f$ and \f$\vec y\f$ are \f$n\f$ element vectors.
 *
 * \param[in] n On entry, \c n specifies the order of the matrix \c a. \f$n \geq 0\f$.
 * \param[in] alpha On entry, \c alpha specifies the scalar \f$alpha\f$.
 * \param[in] x Array of dimension \f$ \geq  1 + ( n - 1 )\cdot\|incx\|\f$. Before entry, the incremented array \c x must contain the \c n element vector \f$\vec x\f$.
 * \param[in] incx On entry, \c incx specifies the increment for the elements of \c x. \f$incx \neq 0\f$.
 * \param[in,out] y Array of dimension \f$ \geq  1 + ( n - 1 )\cdot\|incy\|\f$.
 * - Before entry, the incremented array \c y must contain the \c n element vector \f$\vec y\f$.
 * - On exit, \c y is overwritten by the updated vector \f$\vec y\f$.
 * \param[in] incy On entry, \c incy specifies the increment for the elements of \c y. \f$incy \geq 0\f$.
 */
void daxpy( int *n, double *alpha, double *x, int *incx, double *y, int *incy );

/**
 * \brief \f$\vec y \leftarrow \alpha \mathcal A\vec x + \beta \vec y\f$.
 *
 * \c dsymv Performs the matrix-vector operation:
 *
 * \f[\vec y \leftarrow \alpha \mathcal A\vec x + \beta \vec y \f]
 *
 * where \f$\alpha\f$ and \f$beta\f$ are scalars, \f$\vec x\f$ and \f$\vec y\f$ are \f$n\f$ element vectors and \f$\mathcal A\f$ is an \f$n\f$ by \f$n\f$ symmetric matrix.
 *
 * \param[in] uplo On entry uplo On entry, \c uplo specifies whether the upper or lower triangular part of the
 * array \c a is to be referenced as follows:
 * - \c uplo = 'U' or 'u'   Only the upper triangular part of \c a is to be referenced.
 * - \c uplo = 'L' or 'l'   Only the lower triangular part of \c a is to be referenced.
 * \param[in] n On entry, \c n specifies the order of the matrix \c a. \f$n \geq 0\f$.
 * \param[in] alpha On entry, \c alpha specifies the scalar \f$alpha\f$.
 * \param[in] a Array of dimension (\c lda,\c n).
 * - Before entry with \c uplo = 'U' or 'u', the leading \c n by \c n upper triangular part of the array \c a must contain the upper triangular part of
 *  the symmetric matrix and the strictly lower triangular part of \c a is not referenced.
 * - Before entry with \c uplo = 'L' or 'l', the leading \c n by \c n lower triangular part of the array \c a must contain the lower triangular part of
 *  the symmetric matrix and the strictly upper triangular part of \c a is not referenced.
 * \param[in] lda On entry, \c lda specifies the first dimension of \c a as declared in the calling (sub) program. \f$lda \geq max(1,n)\f$.
 * \param[in] x Array of dimension \f$ \geq  1 + ( n - 1 )\cdot\|incx\|\f$. Before entry, the incremented array \c x must contain the \c n element vector \f$\vec x\f$.
 * \param[in] incx On entry, \c incx specifies the increment for the elements of \c x. \f$incx \neq 0\f$.
 * \param[in] beta On entry, \c beta specifies the scalar \f$\beta\f$. When \c beta is supplied as zero then \c y need not be set on input.
 * \param[in,out] y Array of dimension \f$ \geq 1 + ( n - 1 )\cdot\|incy\|\f$.
 * - Before entry, the incremented array \c y must contain the \c n element vector \f$\vec y\f$.
 * - On exit, \c y is overwritten by the updated vector \f$\vec y\f$.
 * \param[in] incy On entry, \c incy specifies the increment for the elements of \c y. \f$incy \geq 0\f$.
 */
void dsymv( char *uplo, int *n, double *alpha, double *a, int *lda, double *x, int *incx, double *beta, double *y, int *incy );

/**
 * \brief \f$\vec y \leftarrow \alpha \mathcal A\vec x + \beta \vec y\f$ or \f$\vec y \leftarrow \alpha \mathcal A^Tx + beta y\f$,
 *
 * \c dgemv  performs one of the matrix-vector operations
 * - \f$\vec y \leftarrow \alpha \mathcal A\vec x + \beta \vec y\f$,   or
 * - \f$\vec y \leftarrow \alpha \mathcal A^T\vec x + \beta \vec y\f$.
 *
 * where \f$\alpha\f$ and \f$beta\f$ are scalars, \f$\vec x\f$ and \f$\vec y\f$ are vectors and \f$\mathcal A\f$ is an \f$m\f$ by \f$n\f$ matrix.

 * \param[in] trans On entry, \c trans specifies the operation to be performed as follows:
 * - \c trans = 'N' or 'n' \f$\vec y \leftarrow \alpha \mathcal A\vec x + \beta \vec y\f$.
 * - \c trans = 'T' or 't' \f$\vec y \leftarrow \alpha \mathcal A^T\vec x + \beta \vec y\f$.
 * - \c trans = 'C' or 'c' \f$\vec y \leftarrow \alpha \mathcal A^T\vec x + \beta \vec y\f$.
 * \param[in] m On entry, \c m specifies the number of rows of the matrix \c a. \f$M \geq 0\f$.
 * \param[in] n On entry, \c n specifies the number of columns of the matrix \c a. \f$n \geq 0\f$.
 * \param[in] alpha On entry, \c alpha specifies the scalar \f$alpha\f$.
 * \param[in] a Array of dimension (\c lda,\c n ). Before entry, the leading \c m by \c n part of the array \c a must contain the matrix of coefficients.
 * \param[in] lda On entry, \c lda specifies the first dimension of \c a as declared in the calling (sub) program. \f$(lda \geq max( 1, m )\f$.
 * \param[in] x Array of dimension
 * - \f$ \geq 1 + ( n - 1 )\cdot\|incx\|\f$ when \c trans = 'N' or 'n' or,
 * - \f$ \geq 1 + ( m - 1 )\cdot\|incx\|\f$ otherwise.
 * Before entry, the incremented array \c x must contain the vector \f$\vec x\f$.
 * \param[in] incx On entry, \c incx specifies the increment for the elements of \c x. \f$incx \neq 0\f$.
 * \param[in] beta On entry, \c beta specifies the scalar \f$\beta\f$. When \c beta is supplied as zero then \c y need not be set on input.
 * \param[in,out] y Array of dimension
 * - \f$\geq 1 + ( m - 1 )\cdot\|incy\|\f$ when \c trans = 'N' or 'n' or
 * - \f$\geq 1 + ( n - 1 )\cdot\|incy\|\f$ otherwise.
 * - Before entry with \c beta non-zero, the incremented array \c y must contain the vector \f$\vec y\f$.
 * - On exit, \c y is overwritten by the updated vector \f$\vec y\f$.
 * \param[in] incy On entry, \c incy specifies the increment for the elements of \c y. \f$incy \neq 0\f$.
 */
void dgemv_( char *trans, int *m, int *n, double *alpha, double *a, int *lda, double *x, int *incx, double *beta, double *y, int *incy );

/* LAPACK routines */
/**
 * \brief Copies all or part of a two-dimensional matrix \c a to another.
 *
 * dlacpy_ copies all or part of a two-dimensional matrix \c a to another.
 *
 * \param[in] uplo On entry, \c uplo specifies whether the upper or lower triangular part of the array \c a is to be referenced as follows:
 * - \c uplo = 'U': Upper triangular part.
 * - \c uplo = 'L': Lower triangular part.
 * - Otherwise:  All of the matrix \c a
 * \param[in] m The number of rows of the matrix \c a. \f$m \geq 0\f$.
 * \param[in] n The number of columns of the matrix \c a.  \f$(n \geq 0\f$.
 * \param[in] a Array of dimension (\c lda,\c n). The \c m by \c n matrix \c a.  If \c uplo = 'U', only the upper triangle or trapezoid is accessed;
 *  if \c uplo = 'L', only the lower triangle or trapezoid is accessed.
 * \param[in] lda The leading dimension of the array \c a.  \f$lda \geq max(1,m)\f$.
 * \param[out] b Array of dimension (\c ldb,\c n). On exit, \f$\mathcal B = \mathcal A\f$ in the locations specified by \c uplo.
 * \param[in] ldb The leading dimension of the array B.  \f$ldb \geq max(1,M)\f$. *
 */
void dlacpy_( char *uplo, int *m, int *n, double *a, int *lda, double *b, int *ldb );

/**
 * \brief Multiplies the M by N real matrix \c a by a real scalar.
 *
 * \c dlascl multiplies the M by N real matrix \c a by the real scalar \f$cto/cfrom\f$.  This is done without over/underflow as long as the final
 * result <tt>cto*a(i,j)/cfrom</tt> does not over/underflow. \c type specifies that \c a may be full, upper triangular, lower triangular, upper
 * Hessenberg, or banded.
 *
 * \param[in] type Indices the storage type of the input matrix.
 * - = 'G':  \c a is a full matrix.
 * - = 'L':  \c a is a lower triangular matrix.
 * - = 'U':  \c a is an upper triangular matrix.
 * - = 'H':  \c a is an upper Hessenberg matrix.
 * - = 'B':  \c a is a symmetric band matrix with lower bandwidth \c KL and upper bandwidth \c KU and with the only the lower half stored.
 * - = 'Q':  \c a is a symmetric band matrix with lower bandwidth \c KL and upper bandwidth \c KU and with the only the upper half stored.
 * - = 'Z':  \c a is a band matrix with lower bandwidth \c kl and upper bandwidth \c KU.
 * \param[in] kl The lower bandwidth of \c a.  Referenced only if \c type = 'B', 'Q' or 'Z'.
 * \param[in] ku The upper bandwidth of \c a.  Referenced only if \c type = 'B', 'Q' or 'Z'.
 * \param[in] cfrom The divisor.
 * \param[in] cto The matrix \c a is multiplied by \f$cto/cfrom\f$. \c a(i,j) is computed without over/underflow if the final result
 * \f$cto*a(i,j)/cfrom\f$ can be represented without over/underflow. \c cfrom must be nonzero.
 * \param[in] m The number of rows of the matrix \c a.  \f$m \geq 0\f$.
 * \param[in] n The number of columns of the matrix \c a.  \f$n \geq 0\f$.
 * \param[in,out] a Array of dimension (\c lda,\c n). The matrix to be multiplied by \f$cto/cfrom\f$. See \c type for the storage type.
 * \param[in] lda The leading dimension of the array \c a.  \f$lda \geq max(1,m)\f$.
 * \param[out] info if:
 * - \f$info = 0\f$: successful exit.
 * - \f$info < 0\f$: if \f$info = -i\f$, the i-th argument had an illegal value.
 */
void dlascl_( char *type, int *kl, int *ku, double *cfrom, double *cto, int *m, int *n, double *a, int *lda, int *info );

/**
 * \brief Computes the Cholesky factorization of a real symmetric positive definite matrix \c a.
 *
 * \c dpotrf Computes the Cholesky factorization of a real symmetric positive definite matrix \c a.
 *
 * The factorization has the form
 * - \f$\mathcal A = \mathcal U^T\cdot \mathcal U\f$,  if \c uplo = 'U', or
 * - \f$\mathcal A = \mathcal L\cdot \mathcal L^T\f$,  if \c uplo = 'L',
 * where \f$\mathcal U\f$ is an upper triangular matrix and \f$\mathcal L\f$ is lower triangular.
 *
 * This is the block version of the algorithm, calling Level 3 BLAS.
 *
 * \param[in] uplo Specifies the stored part of the matrix.
 * - = 'U':  Upper triangle of \c a is stored;
 * - = 'L':  Lower triangle of \c a is stored.
 * \param[in] n The order of the matrix \c a.  \f$n \geq 0\f$.
 * \param[in,out] a Array of dimension (\c lda,\c n).
 * - On entry, the symmetric matrix \c a.  If \c uplo = 'U', the leading N-by-N upper triangular part of \c a contains the upper
 * triangular part of the matrix \c a, and the strictly lower triangular part of \c a is not referenced.  If \c uplo = 'L', the leading n-by-n lower triangular part
 *  of \c a contains the lower triangular part of the matrix \c a, and the strictly upper triangular part of \c a is not referenced.
 * - On exit, if \f$info = 0\f$, the factor U or L from the Cholesky factorization \f$\mathcal A = \mathcal U^T\cdot
 * U\f$ or \f$\mathcal A = \mathcal L\cdot \mathcal L^T\f$.
 * \param[in] lda The leading dimension of the array \c a. \f$lda \geq max(1,n)\f$.
 * \param[out] info if:
 * - \f$info = 0\f$: successful exit.
 * - \f$info < 0\f$: if \f$info = -i\f$, the i-th argument had an illegal value.
 * - \f$info > 0\f$: if \f$info = i\f$·, the leading minor of order i is not positive definite, and the factorization could not be completed.
 */
void dpotrf_( char *uplo, int *n, double *a, int *lda, int *info );

/**
 * \brief Computes the Cholesky factorization of a real symmetric positive definite matrix
 * \c a stored in packed format.
 *
 * \c dpptrf computes the Cholesky factorization of a real symmetric positive definite
 * matrix A stored in packed format.
 *
 * The factorization has the form:
 * - \f$\mathcal A = \mathcal U^T\cdot \mathcal U\f$, if \c uplo = 'U', or
 * - \f$\mathcal A = \mathcal L\cdot L^T\f$,  if \c uplo = 'L',
 * where \f$\mathcal U\f$ is an upper triangular matrix and \f$\mathcal L\f$ is lower triangular.
 *
 * \param[in] uplo Specifies the stored part of the matrix.
 * - = 'U':  Upper triangle of \c a is stored;
 * - = 'L':  Lower triangle of \c a is stored.
 * \param[in] n The order of the matrix \c A.  \f$n \geq 0\f$.
 * \param[in,out] ap Array of dimension \f$(n*(n+1)/2)\f$.
 * - On entry, the upper or lower triangle of the symmetric matrix \c a, packed columnwise
 * in a linear array. The \em j-th column of \c a is stored in the array \c ap as follows:
 *     -# if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
 *     -# if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
 * - On exit, if \c info = 0, the triangular factor U or L from the
 *          Cholesky factorization A = U**T*U or A = L*L**T, in the same
 *          storage format as A.
 * \param[out] info if:
 * - \f$info = 0\f$: successful exit.
 * - \f$info < 0\f$: if \f$info = -i\f$, the i-th argument had an illegal value.
 * - \f$info > 0\f$: if \f$info = i\f$·, the leading minor of order i is not positive
 * definite, and the factorization could not be completed.
 */
void dpptrf_( char *uplo, int *n, double *ap, int *info );

/**
 * \brief Computes the inverse of a real symmetric positive definite matrix \c a using the Cholesky factorization.
 *
 * \c dpotri computes the inverse of a real symmetric positive definite matrix \c a using the Cholesky
 * factorization \f$\mathcal A = \mathcal U**T*\mathcal U\f$ or \f$\mathcal A = \mathcal L*\mathcal L**T\f$ computed by \c dpotrf_().
 *
 * \param[in] uplo On entry, \c uplo specifies whether the upper or lower triangular part of the array \c a is to be referenced as follows:
 * - \c uplo = 'U':  Upper triangle of \c a is stored;
 * - \c uplo = 'L':  Lower triangle of \c a is stored.
 * \param[in] n The order of the matrix \c a. \f$n \geq 0\f$.
 * \param[in,out] a Array of dimension (\c lda,\c n).
 * - On entry, the triangular factor U or L from the Cholesky factorization \f$\mathcal A = \mathcal U^T\cdot \mathcal
 * U\f$ or \f$\mathcal A = \mathcal L\cdot \mathcal L^T\f$, as computed by dpotrf_().
 * - On exit, the upper or lower triangle of the (symmetric) inverse of \c a, overwriting the input factor U or L.
 * \param[in] lda The leading dimension of the array \c a. \f$lda \geq max(1,n)\f$.
 * \param[out] info if
 * - \f$info = 0\f$: successful exit
 * - \f$info < 0\f$: if \f$info = -i\f$, the i-th argument had an illegal value.
 * - \f$info < 0\f$: if \f$info = i\f$, the \f$(i,i)\f$ element of the factor U or L is zero, and the inverse could not be computed.
 */
void dpotri_( char *uplo, int *n, double *a, int *lda, int *info );

void dpptri_( char *uplo, int *n, double *a, int *info );


#if _MPI_

// BLACS routines
void Cblacs_get( int icontxt, int what, int *val );
void Cblacs_gridinit(int *icontxt, char *order, int nprow, int npcol);
void Cblacs_gridinfo( int icontxt, int *nprow, int *npcol, int *myrow, int *mycol );
void Cblacs_gridmap( int *icontxt, int *usermap, int ldumap, int nprow, int npcol );
int Cblacs_pnum( int icontxt, int prow, int pcol );

int Csys2blacs_handle(MPI_Comm SysCtxt);
int Cfree_blacs_system_handle(int BlacsHandle);
void Cblacs_barrier(int icontxt, char *scope);
int Cblacs_gridexit( int icontxt );

// PBLAS routines
void pdscal( int *n, double *alpha, double *x, int *ix, int *jx, int *descx, int *incx );
void pdcopy(int *n, double *x, int *ix, int *jx, int *descx, int *incx, double *y, int *iy, int *jy, int *descy, int *incy );
void pdaxpy( int *n, double *alpha, double *x, int *ix, int *jx, int *descx, int *incx, double *y, int *iy, int *jy, int *descy, int *incy);

void pdgeadd_(char *trans, int *m, int *n, double *alpha, double *a, int *ia, int *ja, int *desca, double *beta, double *c, int *ic, int *jc, int *descc);
void pdtradd_( char *uplo, char *trans, int *m, int *n, double *alpha, double *a, int *ia, int *ja, int *desca, double *beta, double *c, int *ic, int *jc, int *descc );
void pdgemv_( char *trans, int *m, int *n, double *alpha, double *a, int *ia, int *ja, int *desca, double *x, int *ix, int *jx, int *descx, int *incx, double *beta, double *y, int *iy, int *jy, int *descy, int *incy );
void pdsymv( char *uplo, int *n, double *alpha, double *a, int *ia, int *ja, int *desca, double *x, int *ix, int *jx, int *descx, int *incx, double *beta, double *y, int *iy, int *jy, int *descy, int *incy );

// ScaLAPACK routines
void descinit_(int *desc, int *m, int *n, int *mb, int *nb, int *irsrc, int *icsrc, int *ictxt, int *lld, int *info);
int numroc_(int *n, int *nb, int *iproc, int *isrcproc, int *nprocs);
void infog2l_( int *grindx, int *gcindx, int *desc, int *nprow, int *npcol, int *myrow, int *mycol, int *lrindx, int *lcindx, int *rsrc, int *csrc );

void pdlacpy_( char *uplo, int *m, int *n, double *a, int *ia, int *ja, int *desca, double *b, int *ib, int *jb, int *descb);
void pdlascl_( char *uplo, double *cfrom, double *cto, int *m, int *n, double *a, int *ia, int *ja, int *desca, int *info);

void pdpotrf_( char *uplo, int *n, double *a, int *ia, int *ja, int *desca, int *info );
void pdpotri_( char *uplo, int *n, double *a, int *ia, int *ja, int *desca, int *info );

void pdgetri_( int *n, double *a, int *ia, int *ja, int *desca, int *ipiv, double *work, int *lwork, int *iwork, int *liwork, int *info );
void pdgetrf_( int *m, int *n, double *a, int *ia, int *ja, int *desca, int *ipiv, int *info );
#endif

#endif /* NETILB_H_ */
