/**
 * \file Init_GainMatrix.c
 * \author Ferran Obón Santacana
 * \version 1.0 PARDISO solver for matrix inversion is deprecated.
 * \date 9th of February 2013
 * 
 * \todo Add support for packed storage to reduce memory use.
 *
 * \brief Routines to compute the gain matrix.
 *
 * Routines for calculating the gain matrix. The routines make use of the BLAS and LAPACK libraries to perform the linear
 * algebra operations, including the matrix inversion in single and double precision. Sparse BLAS operations are supported
 * through the Intel MKL library, but since matrix inversion is a dense operations they still rely on LAPACK for this
 * operation. The PARDISO solver is no longer supported since it requires more memory and time than the LAPACK equivalent
 * routines.
 */
#include <stdio.h>          /* For printf(), fprintf() */
#include <stdlib.h>         /* For exit() */

#include "Print_Messages.h" /* For Print_Header() */
#include "Initiation.h"
#include "MatrixVector.h"
#include "MatrixVector_PS.h"

#include "Auxiliary_Math.h"

#if _MKL_
#include <mkl_blas.h>
#include <mkl_lapack.h>
#else 
#include "Netlib.h"
#endif

void IGainMatrix( MatrixVector_t *const IGain, const MatrixVector_t *const Mass, const MatrixVector_t *const Damp,
		  const MatrixVector_t *const Stiff, const Scalars_t Const )
{

     char uplo;    
     int lda, info; /* Leading dimension of the matrix, LAPACK error handling variable */

     int ione = 1;     /* Integer of value one */
     double one = 1.0; /* Double precision one for dlascl() parameter cfrom */
     double Scalar;    /* Double precision scalar */

     uplo = 'L';  /* The lower part of the matrix will be used and the upper part will strictly
		   * not be referenced */
     lda = Max( 1, IGain->Rows );

     /* IGain = Const.Alpha*M + Const.Beta*C + Const.Gamma*K */
     MatrixVector_Add3Mat( Mass, Damp, Stiff, Const, IGain );

     /* LAPACK: Compute the Cholesky factorization */
     dpotrf_( &uplo, &IGain->Rows, IGain->Array, &lda, &info );

     if ( info == 0 ){
	  Print_Header( SUCCESS );
	  printf( "Cholesky factorization successfully completed.\n" );
     }
     else if (info < 0){
	  Print_Header( ERROR );
	  fprintf( stderr, "Cholesky factorization: the %d-th argument has an illegal value.\n", -info );
	  exit( EXIT_FAILURE );
     } else if (info > 0){
	  Print_Header( ERROR );
	  fprintf( stderr, "Cholesky factorization: the leading minor of order %d is not positive definite,", info );
	  fprintf( stderr, " and the factorization could not be completed.\n" );
	  exit( EXIT_FAILURE );
     }

     /* LAPACK: Compute the inverse of the IGain matrix using the Cholesky factorization computed
      * by dpotrf_( ) */
     dpotri_( &uplo, &IGain->Rows, IGain->Array, &lda, &info );

     if ( info == 0 ){
	  Print_Header( SUCCESS );
	  printf( "Matrix Inversion successfully completed.\n" );
     } else if (info < 0){
	  Print_Header( ERROR );
	  fprintf( stderr, "Matrix inversion: the %d-th argument has an illegal value.\n", -info );
	  exit( EXIT_FAILURE );
     } else if (info > 0){
	  Print_Header( ERROR );
	  fprintf( stderr, "Matrix Inversion: the (%d,%d) element of the factor U or L is zero,", info, info );
	  fprintf( stderr, "and the inverse could not be computed.\n" );
	  exit( EXIT_FAILURE );
     }

     Scalar = Const.Lambda;
     dlascl_( &uplo, &ione, &ione, &one, &Scalar, &IGain->Rows, &IGain->Cols, IGain->Array, &lda, &info );

     if( info < 0 ){
	  Print_Header( ERROR );
	  fprintf( stderr, "Gain matrix: the %d-th argument has an illegal value and cannot be scaled.\n", -info );
	  exit( EXIT_FAILURE );
     } else {
	  Print_Header( SUCCESS );
	  printf( "Gain matrix successfully calculated.\n" );
     }

}

void IGainMatrix_PS( MatrixVector_t *const IGain, const MatrixVector_t *const Mass, const MatrixVector_t *const Damp,
		  const MatrixVector_t *const Stiff, const Scalars_t Const )
{
     char uplo;
     int info;         /* LAPACK error handling variable */
     double Scalar;    /* Double precision scalar */
     int incx, Length;

     incx = 1;
     uplo = 'L';  /* The lower part of the matrix will be used and the upper part will
		   * strictly not be referenced */

     /* IGain = Const.Alpha*M + Const.Beta*C + Const.Gamma*K */
     MatrixVector_Add3Mat_PS( Mass, Damp, Stiff, Const, IGain );

     /* LAPACK: Compute the Cholesky factorization */
     dpptrf_( &uplo, &IGain->Rows, IGain->Array, &info );

     if ( info == 0 ){
	  Print_Header( SUCCESS );
	  printf( "Cholesky factorization successfully completed.\n" );
     }
     else if (info < 0){
	  Print_Header( ERROR );
	  fprintf( stderr, "Cholesky factorization: the %d-th argument has an illegal value.\n", -info );
	  exit( EXIT_FAILURE );
     } else if (info > 0){
	  Print_Header( ERROR );
	  fprintf( stderr, "Cholesky factorization: the leading minor of order %d is not positive definite,", info );
	  fprintf( stderr, " and the factorization could not be completed.\n" );
	  exit( EXIT_FAILURE );
     }

     /* LAPACK: Compute the inverse of the IGain matrix using the Cholesky factorization
      * computed by dpptrf_( ) */
     dpptri_( &uplo, &IGain->Rows, IGain->Array, &info );

     if ( info == 0 ){
	  Print_Header( SUCCESS );
	  printf( "Matrix Inversion successfully completed.\n" );
     } else if (info < 0){
	  Print_Header( ERROR );
	  fprintf( stderr, "Matrix inversion: the %d-th argument has an illegal value.\n", -info );
	  exit( EXIT_FAILURE );
     } else if (info > 0){
	  Print_Header( ERROR );
	  fprintf( stderr, "Matrix Inversion: the (%d,%d) element of the factor U or L is zero,", info, info );
	  fprintf( stderr, "and the inverse could not be computed.\n" );
	  exit( EXIT_FAILURE );
     }

     Length = (IGain->Rows*IGain->Cols + IGain->Rows)/2;
     Scalar = Const.Lambda;
     dscal( &Length, &Scalar, IGain->Array, &incx );
     
     Print_Header( SUCCESS );
     printf( "Gain matrix successfully calculated.\n" );
}
