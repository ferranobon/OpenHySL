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
 * \brief \f$ x \leftarrow \alpha x\f$
 *
 * \c dscal computes the product of a vector by a scalar \f$ x \leftarrow \alpha x\f$. Uses unrolled loops for increment equal to one.
 *
 * where \f$\alpha\f$ is a scalar and \f$x\f$ is \f$n\f$ element vector.
 *
 * \param[in] N On entry, \c N specifies the order of the matrix \c A. \f$N \geq 0\f$.
 * \param[in] ALPHA On entry, \c ALPHA specifies the scalar \f$alpha\f$.
 * \param[in,out] X Array of dimension \f$ \geq 1 + (N -1)\cdot \|INCX\|\f$.
 * - Before entry, the incremented array \c X must contain the \c N element vector \e x.
 * - On exit, \c X is overwritten by the updated vector \e x if \f$N \neq 0\f$. Otherwise it remains unchanged.
 * \param[in] INCX On entry, \c INCX specifies the increment for the elements of \c X. \f$INCX \neq 0\f$.
 */
void dscal_( int *N, double *ALPHA, double *X, int *INCX );

/**
 * \brief \f$ y \leftarrow x\f$
 *
 * \c dcopy copies a vector,\e x, to another vector, \e y \f$ y \leftarrow x\f$. Uses unrolled loops for increments equal to one.
 *
 * where \f$x\f$ and \f$y\f$ are \f$n\f$ element vectors.
 *
 * \param[in] N On entry, \c N specifies the order of the matrix \c A. \f$N \geq 0\f$.
 * \param[in] X Array of dimension \f$ \geq  1 + ( N - 1 )\cdot\|INCX\|\f$. Before entry, the incremented array \c X must contain the \c N element vector \e x.
 * \param[in] INCX On entry, \c INCX specifies the increment for the elements of \c X. \f$INCX \neq 0\f$.
 * \param[out] Y Array of dimension \f$ \geq  1 + ( N - 1 )\cdot\|INCY\|\f$.
 * - If \f$N > 0\f$ it contains a copy of the vector \e x. Otherwise, elements are unaltered.
 * \param[in] INCY On entry, \c INCY specifies the increment for the elements of \c Y. \f$INCY \geq 0\f$.
 */
void dcopy_( int *N, double *X, int *INCX, double *Y, int *INCY );

/**
 * \brief \f$y \leftarrow \alpha x + y\f$
 *
 * \c daxpy computes a vector-scalar product, \f$\alpha x\f$ and adds the result to another vector \f$y\f$. Uses unrolled loops for increments equal to one.
 *
 * \f[
 * y \leftarrow \alpha x + y
 * \f]
 *
 * where \f$\alpha\f$ is a scalar and \f$x\f$ and \f$y\f$ are \f$n\f$ element vectors.
 *
 * \param[in] N On entry, \c N specifies the order of the matrix \c A. \f$N \geq 0\f$.
 * \param[in] ALPHA On entry, \c ALPHA specifies the scalar \f$alpha\f$.
 * \param[in] X Array of dimension \f$ \geq  1 + ( N - 1 )\cdot\|INCX\|\f$. Before entry, the incremented array \c X must contain the \c N element vector \e x.
 * \param[in] INCX On entry, \c INCX specifies the increment for the elements of \c X. \f$INCX \neq 0\f$.
 * \param[in,out] Y Array of dimension \f$ \geq  1 + ( N - 1 )\cdot\|INCY\|\f$.
 * - Before entry, the incremented array \c Y must contain the \c N element vector \e y.
 * - On exit, \c Y is overwritten by the updated vector \e y.
 * \param[in] INCY On entry, \c INCY specifies the increment for the elements of \c Y. \f$INCY \geq 0\f$.
 */
void daxpy_( int *N, double *ALPHA, double *X, int *INCX, double *Y, int *INCY );

/**
 * \brief \f$y \leftarrow \alpha Ax + \beta y\f$.
 *
 * \c dsymv Performs the matrix-vector operation
 * \f[
 * y \leftarrow \alpha Ax + \beta y
 * \f]
 *
 * where \f$\alpha\f$ and \f$beta\f$ are scalars, \f$x\f$ and \f$y\f$ are \f$n\f$ element vectors and \f$A\f$ is an \f$n\f$ by \f$n\f$ symmetric matrix.
 *
 * \param[in] UPLO On entry UPLO On entry, \c UPLO specifies whether the upper or lower triangular part of the array A is to be referenced as follows:
 * - \c UPLO = 'U' or 'u'   Only the upper triangular part of \c A is to be referenced.
 * - \c UPLO = 'L' or 'l'   Only the lower triangular part of \c A is to be referenced.
 * \param[in] N On entry, \c N specifies the order of the matrix \c A. \f$N \geq 0\f$.
 * \param[in] ALPHA On entry, \c ALPHA specifies the scalar \f$alpha\f$.
 * \param[in] A Array of dimension (\c LDA,\c N).
 * - Before entry with \c UPLO = 'U' or 'u', the leading \c N by \c N upper triangular part of the array \c A must contain the upper triangular part of
 *  the symmetric matrix and the strictly lower triangular part of A is not referenced.
 * - Before entry with \c UPLO = 'L' or 'l', the leading \c N by \c N lower triangular part of the array \c A must contain the lower triangular part of
 *  the symmetric matrix and the strictly upper triangular part of A is not referenced.
 * \param[in] LDA On entry, \c LDA specifies the first dimension of \c A as declared in the calling (sub) program. \f$LDA \geq max(1,N)\f$.
 * \param[in] X Array of dimension \f$ \geq  1 + ( N - 1 )\cdot\|INCX\|\f$. Before entry, the incremented array \c X must contain the \c N element vector \e x.
 * \param[in] INCX On entry, \c INCX specifies the increment for the elements of \c X. \f$INCX \neq 0\f$.
 * \param[in] BETA On entry, \c BETA specifies the scalar \f$\beta\f$. When \c BETA is supplied as zero then \c Y need not be set on input.
 * \param[in,out] Y Array of dimension \f$ \geq 1 + ( N - 1 )\cdot\|INCY\|\f$.
 * - Before entry, the incremented array \c Y must contain the \c N element vector \e y.
 * - On exit, \c Y is overwritten by the updated vector \e y.
 * \param[in] INCY On entry, \c INCY specifies the increment for the elements of \c Y. \f$INCY \geq 0\f$.
 */
void dsymv_( char *UPLO, int *N, double *ALPHA, double *A, int *LDA, double *X, int *INCX, double *BETA, double *Y, int *INCY );

/**
 * \brief \f$ y \leftarrow \alpha Ax + \beta y\f$ or \f$y \leftarrow \alpha A^Tx + beta y\f$,
 *
 * \c dgemv  performs one of the matrix-vector operations
 * - \f$y \leftarrow \alpha Ax + \beta y\f$,   or
 * - \f$y \leftarrow \alpha A^Tx + \beta y\f$.
 *
 * where \f$\alpha\f$ and \f$beta\f$ are scalars, \f$x\f$ and \f$y\f$ are vectors and \f$A\f$ is an \f$m\f$ by \f$n\f$ matrix.

 * \param[in] TRANS On entry, TRANS specifies the operation to be performed as follows:
 * - \c TRANS = 'N' or 'n' \f$y \leftarrow \alpha Ax + \beta y\f$.
 * - \c TRANS = 'T' or 't' \f$y \leftarrow \alpha A^Tx + \beta y\f$.
 * - \c TRANS = 'C' or 'c' \f$y \leftarrow \alpha A^Tx + \beta y\f$.
 * \param[in] M On entry, \c M specifies the number of rows of the matrix \c A. \f$M \geq 0\f$.
 * \param[in] N On entry, \c N specifies the number of columns of the matrix \c A. \f$N \geq 0\f$.
 * \param[in] ALPHA On entry, \c ALPHA specifies the scalar \f$alpha\f$.
 * \param[in] A Array of dimension (\c LDA,\c n ). Before entry, the leading \c M by \c N part of the array \c A must contain the matrix of coefficients.
 * \param[in] LDA On entry, \c LDA specifies the first dimension of \c A as declared in the calling (sub) program. \f$(LDA \geq max( 1, M )\f$.
 * \param[in] X Array of dimension
 * - \f$ \geq 1 + ( N - 1 )\cdot\|INCX\|\f$ when \c TRANS = 'N' or 'n' or,
 * - \f$ \geq 1 + ( M - 1 )\cdot\|INCX\|\f$ otherwise.
 * Before entry, the incremented array \c X must contain the vector \e x.
 * \param[in] INCX On entry, \c INCX specifies the increment for the elements of \c X. \f$INCX \neq 0\f$.
 * \param[in] BETA On entry, \c BETA specifies the scalar \f$\beta\f$. When \c BETA is supplied as zero then \c Y need not be set on input.
 * \param[in,out] Y Array of dimension
 * - \f$\geq 1 + ( M - 1 )\cdot\|INCY\|\f$ when \c TRANS = 'N' or 'n' or
 * - \f$\geq 1 + ( N - 1 )\cdot\|INCY\|\f$ otherwise.
 * - Before entry with \c BETA non-zero, the incremented array \c Y must contain the vector \e y.
 * - On exit, \c Y is overwritten by the updated vector \e y.
 * \param[in] INCY On entry, \c INCY specifies the increment for the elements of \c Y. \f$INCY \neq 0\f$.
 */
void dgemv_( char *TRANS, int *M, int *N, double *ALPHA, double *A, int *LDA, double *X, int *INCX, double *BETA, double *Y, int *INCY );

/* LAPACK routines */
#if !_SPARSE_
/**
 * \brief Copies all or part of a two-dimensional matrix A to another.
 *
 * dlacpy_ copies all or part of a two-dimensional matrix A to another.
 *
 * \param[in] UPLO On entry, \c UPLO specifies whether the upper or lower triangular part of the array A is to be referenced as follows:
 * - \c UPLO = 'U': Upper triangular part.
 * - \c UPLO = 'L': Lower triangular part.
 * - Otherwise:  All of the matrix A
 * \param[in] M The number of rows of the matrix A. \f$M \geq 0\f$.
 * \param[in] N The number of columns of the matrix A.  \f$(N \geq 0\f$.
 * \param[in] A Array of dimension (\c LDA,\c N). The m by n matrix A.  If \c UPLO = 'U', only the upper triangle or trapezoid is accessed;
 *  if \c UPLO = 'L', only the lower triangle or trapezoid is accessed.
 * \param[in] LDA The leading dimension of the array A.  \f$LDA \geq max(1,M)\f$.
 * \param[out] B Array of dimension (\c LDB,\c N). On exit, \f$B = A\f$ in the locations specified by \c UPLO.
 * \param[in] LDB The leading dimension of the array B.  \f$LDB \geq max(1,M)\f$. *
 */
void dlacpy_( char *UPLO, int *M, int *N, double *A, int *LDA, double *B, int *LDB );

/**
 * \brief Multiplies the M by N real matrix A by a real scalar.
 *
 * \c dlascl multiplies the M by N real matrix A by the real scalar  CTO/CFROM.  This is done without over/underflow as long as the final
 * result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that A may be full, upper triangular, lower triangular, upper
 * Hessenberg, or banded.
 *
 * \param[in] TYPE Indices the storage type of the input matrix.
 * - = 'G':  A is a full matrix.
 * - = 'L':  A is a lower triangular matrix.
 * - = 'U':  A is an upper triangular matrix.
 * - = 'H':  A is an upper Hessenberg matrix.
 * - = 'B':  A is a symmetric band matrix with lower bandwidth \c KL and upper bandwidth \c KU and with the only the lower half stored.
 * - = 'Q':  A is a symmetric band matrix with lower bandwidth \c KL and upper bandwidth \c KU and with the only the upper half stored.
 * - = 'Z':  A is a band matrix with lower bandwidth \c KL and upper bandwidth \c KU.
 * \param[in] KL The lower bandwidth of A.  Referenced only if \c TYPE = 'B', 'Q' or 'Z'.
 * \param[in] KU The upper bandwidth of A.  Referenced only if \c TYPE = 'B', 'Q' or 'Z'.
 * \param[in] CFROM The divisor.
 * \param[in] CTO The matrix A is multiplied by \f$CTO/CFROM\f$. A(I,J) is computed without over/underflow if the final result
 * \f$CTO*A(I,J)/CFROM\f$ can be represented without over/underflow. \c CFROM must be nonzero.
 * \param[in] M The number of rows of the matrix A.  \f$M \geq 0\f$.
 * \param[in] N The number of columns of the matrix A.  \f$N \geq 0\f$.
 * \param[in,out] A Array of dimension (\c LDA,\c N). The matrix to be multiplied by \f$CTO/CFROM\f$. See \c TYPE for the storage type.
 * \param[in] LDA The leading dimension of the array A.  \f$LDA \geq max(1,M)\f$.
 * \param[out] INFO if:
 * - \f$INFO = 0\f$: successful exit.
 * - \f$INFO < 0\f$: if \f$INFO = -i\f$, the i-th argument had an illegal value.
 */
void dlascl_( char *TYPE, int *KL, int *KU, double *CFROM, double *CTO, int *M, int *N, double *A, int *LDA, int *INFO );

/**
 * \brief Computes the Cholesky factorization of a real symmetric positive definite matrix A.
 *
 * \c dpotrf Computes the Cholesky factorization of a real symmetric positive definite matrix A.
 *
 * The factorization has the form
 * - \f$A = U^T\cdot U\f$,  if \c UPLO = 'U', or
 * - \f$A = L\cdot L^T\f$,  if \c UPLO = 'L',
 * where U is an upper triangular matrix and L is lower triangular.
 *
 * This is the block version of the algorithm, calling Level 3 BLAS.
 *
 * \param[in] UPLO Specifies the stored part of the matrix.
 * - = 'U':  Upper triangle of A is stored;
 * - = 'L':  Lower triangle of A is stored.
 * \param[in] N The order of the matrix A.  \f$N \geq 0\f$.
 * \param[in,out] A Array of dimension (\c LDA,\c N).
 * - On entry, the symmetric matrix \c A.  If \c UPLO = 'U', the leading N-by-N upper triangular part of \c A contains the upper
 * triangular part of the matrix \c A, and the strictly lower triangular part of A is not referenced.  If \c UPLO = 'L', the leading N-by-N lower triangular part
 *  of \c A contains the lower triangular part of the matrix \c A, and the strictly upper triangular part of \c A is not referenced.
 * - On exit, if \f$INFO = 0\f$, the factor U or L from the Cholesky factorization \f$A = U^T\cdot U\f$ or \f$A = L\cdot L^T\f$.
 * \param[in] LDA The leading dimension of the array A. \f$LDA \geq max(1,N)\f$.
 * \param[out] INFO if:
 * - \f$INFO = 0\f$: successful exit.
 * - \f$INFO < 0\f$: if \f$INFO = -i\f$, the i-th argument had an illegal value.
 * - \f$INFO > 0\f$: if \f$INFO = i\f$·, the leading minor of order i is not positive definite, and the factorization could not be completed.
 */
void dpotrf_( char *UPLO, int *N, double *A, int *LDA, int *INFO );

/**
 * \brief Computes the inverse of a real symmetric positive definite matrix A using the Cholesky factorization.
 *
 * \c dpotri computes the inverse of a real symmetric positive definite matrix A using the Cholesky factorization A = U**T*U or A = L*L**T computed by \c dpotrf_().
 *
 * \param[in] UPLO On entry, \c UPLO specifies whether the upper or lower triangular part of the array A is to be referenced as follows:
 * - \c UPLO = 'U':  Upper triangle of A is stored;
 * - \c UPLO = 'L':  Lower triangle of A is stored.
 * \param[in] N The order of the matrix A. \f$N \geq 0\f$.
 * \param[in,out] A Array of dimension (\c LDA,\c N).
 * - On entry, the triangular factor U or L from the Cholesky factorization \f$A = U^T\cdot U\f$ or \f$A = L\cdot L^T\f$, as computed by dpotrf_().
 * - On exit, the upper or lower triangle of the (symmetric) inverse of \c A, overwriting the input factor U or L.
 * \param[in] LDA The leading dimension of the array \c A. \f$LDA \geq max(1,N)\f$.
 * \param[out] INFO if
 * - \f$INFO = 0\f$: successful exit
 * - \f$INFO < 0\f$: if \f$INFO = -i\f$, the i-th argument had an illegal value.
 * - \f$INFO < 0\f$: if \f$INFO = i\f$, the \f$(i,i)\f$ element of the factor U or L is zero, and the inverse could not be computed.
 */
void dpotri_( char *UPLO, int *N, double *A, int *LDA, int *INFO );
#endif

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
void pdscal_( int *n, double *alpha, double *x, int *ix, int *jx, int *descx, int *incx );
void pdcopy_(int *n, double *x, int *ix, int *jx, int *descx, int *incx, double *y, int *iy, int *jy, int *descy, int *incy );
void pdaxpy_( int *n, double *alpha, double *x, int *ix, int *jx, int *descx, int *incx, double *y, int *iy, int *jy, int *descy, int *incy);

void pdgeadd_(char *trans, int *m, int *n, double *alpha, double *a, int *ia, int *ja, int *desca, double *beta, double *c, int *ic, int *jc, int *descc);
void pdtradd_( char *uplo, char *trans, int *m, int *n, double *alpha, double *a, int *ia, int *ja, int *desca, double *beta, double *c, int *ic, int *jc, int *descc );
void pdgemv_( char *trans, int *m, int *n, double *alpha, double *a, int *ia, int *ja, int *desca, double *x, int *ix, int *jx, int *descx, int *incx, double *beta, double *y, int *iy, int *jy, int *descy, int *incy );
void pdsymv_( char *uplo, int *n, double *alpha, double *a, int *ia, int *ja, int *desca, double *x, int *ix, int *jx, int *descx, int *incx, double *beta, double *y, int *iy, int *jy, int *descy, int *incy );

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
