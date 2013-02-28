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
void dscal( int *n, double *alpha, double *x, int *incx );
void dcopy( int *n, double *x, int *incx, double *y, int *incy );
void daxpy( int *n, double *alpha, double *x, int *incx, double *y, int *incy );
void dsymv( char *uplo, int *n, double *alpha, double *a, int *lda, double *x, int *incx, double *beta, double *y, int *incy );
void dgemv_( char *trans, int *m, int *n, double *alpha, double *a, int *lda, double *x, int *incx, double *beta, double *y, int *incy );
/* lapack routines */
void dlacpy_( char *uplo, int *m, int *n, double *a, int *lda, double *b, int *ldb );
void dlascl_( char *type, int *kl, int *ku, double *cfrom, double *cto, int *m, int *n, double *a, int *lda, int *info );
void dpotrf_( char *uplo, int *n, double *a, int *lda, int *info );
void dpotri_( char *uplo, int *n, double *a, int *lda, int *info );


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
