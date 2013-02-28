/**
 * \file Init_Rayleigh.c
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 9th of February 2013
 * 
 * \todo Add support for packed storage to reduce memory use.
 *
 * \brief Rayleigh damping routines.
 *
 * Routines for calculating the proportional viscous damping matrix using Rayleigh Damping. The routines
 * make use of the BLAS library to perform the linear algebra operations and they support both single and
 * double precision. Sparse BLAS operations are supported through the Intel MKL library.
 */
#include <stdlib.h>        /* For exit() */

#include "Initiation.h"     /* Header files for the initiation phase */
#include "MatrixVector.h"   /* MatrixVector definition */
#include "Print_Messages.h" /* For Print_Message() */

#include "Auxiliary_Math.h" /* For Max() */

#if _MKL_
#include <mkl_blas.h>
#include <mkl_lapack.h>
#else
#include "Netlib.h"
#endif

void Rayleigh_Damping( const MatrixVector_t *const Mass, const MatrixVector_t *const Stiff, MatrixVector_t *const Damp,
		       const Rayleigh_t *const Rayleigh )
{
     char uplo;
     int ione;
     int incx, incy;
     int info, lda, ldb;
     double done;

     int Rows, Cols, Length;
     double alpha, beta;

     int i;

     ione = 1; done = 1.0;
     incx = 1; incy = 1;
     uplo = 'L';  /* The lower part of the matrix will be used and the upper part will strictly
		   * not be referenced */
     Rows = Mass->Rows;
     Cols = Mass->Cols;
     alpha = Rayleigh->Alpha;
     beta = Rayleigh->Beta;

     lda = Max( 1, Rows );
     ldb = Max( 1, Damp->Rows );

     /* LAPACK: C = M */
     dlacpy_( &uplo, &Rows, &Cols, Mass->Array, &lda, Damp->Array, &ldb );

     /* LAPACK: C = Rayleigh.alpha*M = Rayleigh.alpha*C */
     dlascl_( &uplo, &ione, &ione, &done, &alpha, &Damp->Rows, &Damp->Cols, Damp->Array, &lda, &info );

     if ( info < 0 ){
	  Print_Message( ERROR, 3, STRING, "dlascl: The ", INT, -info, STRING, "-th argument had an illegal value." );
	  exit( EXIT_FAILURE );
     }

     /* BLAS: C = alpha*M + beta*K = C + beta*K. Only half of the matrix is calculated */
     for ( i = 0; i < Damp->Rows; i++ ){
	  Length = Damp->Rows - i;
	  daxpy( &Length, &beta, &Stiff->Array[i*Stiff->Rows + i], &incx, &Damp->Array[i*Damp->Rows +i], &incy);
     }

     Print_Message( SUCCESS, 1, STRING, "Damping matrix successfully calculated." );
}
