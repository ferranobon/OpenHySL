#include <stdio.h>          /* For printf(), fprintf() */
#include <stdlib.h>         /* For exit() */

#include "Auxiliary_Math.h" /* For Max() */
#include "MatrixVector.h"   /* MatrixVector definition */
#include "Print_Messages.h" /* For Print_Header() */
#include "Rayleigh.h"       /* Rayleigh damping routines */

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
	  Print_Header( ERROR );
	  fprintf( stderr, "dlascl: The %d-th argument had an illegal value.\n", -info );
	  exit( EXIT_FAILURE );
     }

     /* BLAS: C = alpha*M + beta*K = C + beta*K. Only half of the matrix is calculated */
     for ( i = 0; i < Damp->Rows; i++ ){
	  Length = Damp->Rows - i;
	  daxpy( &Length, &beta, &Stiff->Array[i*Stiff->Rows + i], &incx, &Damp->Array[i*Damp->Rows +i], &incy);
     }

     Print_Header( SUCCESS );
     printf("Damping matrix successfully calculated.\n" );
}

void Rayleigh_Damping_PS( const MatrixVector_t *const Mass, const MatrixVector_t *const Stiff, MatrixVector_t *const Damp,
			  const Rayleigh_t *const Rayleigh )
{
     int incx, incy;
     int Length;
     double alpha, beta;

     incx = 1; incy = 1;

     alpha = Rayleigh->Alpha;
     beta = Rayleigh->Beta;

     Length = (Damp->Rows*Damp->Cols + Damp->Rows)/2;

     /* BLAS: C = M */
     dcopy( &Length, Mass->Array, &incx, Damp->Array, &incy );

     /* BLAS: C = Rayleigh.alpha*M = Rayleigh.alpha*C */
     dscal( &Length, &alpha, Damp->Array, &incx );

     /* BLAS: C = alpha*M + beta*K = C + beta*K. Only half of the matrix is calculated */
     daxpy( &Length, &beta, Stiff->Array, &incx, Damp->Array, &incy);
     
     Print_Header( SUCCESS );
     printf("Damping matrix successfully calculated.\n" );
}
