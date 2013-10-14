#include <stdio.h>          /* For printf(), fprintf() */
#include <stdlib.h>         /* For exit() */

#include "GainMatrix.h"
#include "MatrixVector.h"
#include "MatrixVector_PS.h"
#include "Print_Messages.h"  /* For Print_Header() */

#include "Auxiliary_Math.h"

#if _MKL_
#include <mkl_blas.h>
#include <mkl_lapack.h>
#else 
#include "Netlib.h"
#endif

void IGainMatrix( MatrixVector_t *const IGain, const MatrixVector_t *const Mass, 
		  const MatrixVector_t *const Damp, const MatrixVector_t *const Stiff,
		  const Scalars_t Const )
{

     char uplo;    
     int lda, info;    /* Leading dimension of the matrix, LAPACK error handling
			* variable */

     int ione = 1;     /* Integer of value one */
     double one = 1.0; /* Double precision one for dlascl() parameter cfrom */
     double Scalar;    /* Double precision scalar */

     uplo = 'L';       /* The lower part of the matrix will be used and the upper part will strictly not be
			* referenced */
     lda = Max( 1, IGain->Rows );

     /* IGain = Const.Alpha*M + Const.Beta*C + Const.Gamma*K */
     MatrixVector_Add3Mat( Mass, Damp, Stiff, Const, IGain );

     /* LAPACK: Compute the Cholesky factorization */
     dpotrf_( &uplo, &IGain->Rows, IGain->Array, &lda, &info );

     if ( info == 0 ){
	  Print_Header( SUCCESS );
	  printf( "Cholesky factorization successfully completed.\n" );
     } else if (info < 0){
	  Print_Header( ERROR );
	  fprintf( stderr, "Cholesky factorization: the %d-th argument has an illegal value.\n", -info );
	  exit( EXIT_FAILURE );
     } else if (info > 0){
	  Print_Header( ERROR );
	  fprintf( stderr, "Cholesky factorization: the leading minor of order %d is not positive definite,", info );
	  fprintf( stderr, " and the factorization could not be completed.\n" );
	  exit( EXIT_FAILURE );
     }

     /* LAPACK: Compute the inverse of the IGain matrix using the Cholesky factorization computed by
      * dpotrf() */
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
     dlascl_( &uplo, &ione, &ione, &one, &Scalar, &IGain->Rows, &IGain->Cols,
	      IGain->Array, &lda, &info );

     if( info < 0 ){
	  Print_Header( ERROR );
	  fprintf( stderr, "Gain matrix: the %d-th argument has an illegal value and cannot be scaled.\n", -info );
	  exit( EXIT_FAILURE );
     } else {
	  Print_Header( SUCCESS );
	  printf( "Gain matrix successfully calculated.\n" );
     }

}

void IGainMatrix_PS( MatrixVector_t *const IGain, const MatrixVector_t *const Mass,
		     const MatrixVector_t *const Damp, const MatrixVector_t *const Stiff,
		     const Scalars_t Const )
{
     char uplo;
     int info;         /* LAPACK error handling variable */
     double Scalar;    /* Double precision scalar */
     int incx, Length;

     incx = 1;
     uplo = 'L';       /* The lower part of the matrix will be used and the upper part will strictly not be
			* referenced */

     /* IGain = Const.Alpha*M + Const.Beta*C + Const.Gamma*K */
     MatrixVector_Add3Mat_PS( Mass, Damp, Stiff, Const, IGain );

     /* LAPACK: Compute the Cholesky factorization */
     dpptrf_( &uplo, &IGain->Rows, IGain->Array, &info );

     if ( info == 0 ){
	  Print_Header( SUCCESS );
	  printf( "Cholesky factorization successfully completed.\n" );
     } else if (info < 0){
	  Print_Header( ERROR );
	  fprintf( stderr, "Cholesky factorization: the %d-th argument has an illegal value.\n", -info );
	  exit( EXIT_FAILURE );
     } else if (info > 0){
	  Print_Header( ERROR );
	  fprintf( stderr, "Cholesky factorization: the leading minor of order %d is not positive definite,", info );
	  fprintf( stderr, " and the factorization could not be completed.\n" );
	  exit( EXIT_FAILURE );
     }

     /* LAPACK: Compute the inverse of the IGain matrix using the Cholesky factorization computed by
      *	dpptrf() */
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

void IGainMatrix_Float2Double( MatrixVector_t *const IGain, const MatrixVector_t *const Mass, 
		  const MatrixVector_t *const Damp, const MatrixVector_t *const Stiff,
		  const Scalars_t Const )
{
     char uplo;    
     int lda, info;    /* Leading dimension of the matrix, LAPACK error handling
			* variable */

     int ione = 1;     /* Integer of value one */
     double one = 1.0; /* Double precision one for dlascl() parameter cfrom */
     double Scalar;    /* Double precision scalar */

     double *TempMatDouble = NULL;  /* Temporal double matrix for matrix inversion purposes. */
     size_t i;
     
     uplo = 'L';       /* The lower part of the matrix will be used and the upper part will strictly not be
			* referenced */
     lda = Max( 1, IGain->Rows );
     
     /* IGain = Const.Alpha*M + Const.Beta*C + Const.Gamma*K */
     MatrixVector_Add3Mat( Mass, Damp, Stiff, Const, IGain );
     
     TempMatDouble = (double *) calloc ( (size_t) IGain->Rows*(size_t) IGain->Cols, sizeof(double) );

     if( TempMatDouble == NULL ){
	  Print_Header( ERROR );
	  fprintf( stderr, "IGainMatrix_Float2Double(): Out of memory.\n" );
	  exit( EXIT_FAILURE );
     }

#pragma omp parallel for
     for( i = 0; i < (size_t) IGain->Rows*(size_t) IGain->Cols; i++ ){
	  TempMatDouble[i] = (double) IGain->Array[i];
     }

     /* LAPACK: Compute the Cholesky factorization */
     dpotrf_( &uplo, &IGain->Rows, IGain->Array, &lda, &info );

     if ( info == 0 ){
	  Print_Header( SUCCESS );
	  printf( "Cholesky factorization successfully completed.\n" );
     } else if (info < 0){
	  Print_Header( ERROR );
	  fprintf( stderr, "Cholesky factorization: the %d-th argument has an illegal value.\n", -info );
	  exit( EXIT_FAILURE );
     } else if (info > 0){
	  Print_Header( ERROR );
	  fprintf( stderr, "Cholesky factorization: the leading minor of order %d is not positive definite,", info );
	  fprintf( stderr, " and the factorization could not be completed.\n" );
	  exit( EXIT_FAILURE );
     }

     /* LAPACK: Compute the inverse of the IGain matrix using the Cholesky factorization computed by
      * dpotrf() */
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
     dlascl_( &uplo, &ione, &ione, &one, &Scalar, &IGain->Rows, &IGain->Cols,
	      IGain->Array, &lda, &info );

#pragma omp parallel for
     for( i = 0; i < (size_t) IGain->Rows*(size_t) IGain->Cols; i++ ){
	  IGain->Array[i] = (float) TempMatDouble[i];
     }

     /* Free the temporal double matrix */
     free( TempMatDouble );

     if( info < 0 ){
	  Print_Header( ERROR );
	  fprintf( stderr, "Gain matrix: the %d-th argument has an illegal value and cannot be scaled.\n", -info );
	  exit( EXIT_FAILURE );
     } else {
	  Print_Header( SUCCESS );
	  printf( "Gain matrix successfully calculated.\n" );
     }

}


void IGainMatrix_Float2Double_PS( MatrixVector_t *const IGain, const MatrixVector_t *const Mass,
		     const MatrixVector_t *const Damp, const MatrixVector_t *const Stiff,
		     const Scalars_t Const )
{

     char uplo;
     int info;         /* LAPACK error handling variable */
     double Scalar;    /* Double precision scalar */
     int incx;
     double *TempMatDouble = NULL;  /* Temporal double matrix for matrix inversion purposes. */
     size_t Length, i;

     incx = 1;
     uplo = 'L';       /* The lower part of the matrix will be used and the upper part will strictly not be
			* referenced */

     /* IGain = Const.Alpha*M + Const.Beta*C + Const.Gamma*K */
     MatrixVector_Add3Mat_PS( Mass, Damp, Stiff, Const, IGain );

     Length = ((size_t) IGain->Rows*(size_t) IGain->Cols + (size_t)IGain->Rows)/2;
     TempMatDouble = (double *) calloc ( Length, sizeof(double) );
     if( TempMatDouble == NULL ){
	  Print_Header( ERROR );
	  fprintf( stderr, "IGainMatrix_Float2Double_PS(): Out of memory.\n" );
	  exit( EXIT_FAILURE );
     }

#pragma omp parallel for
     for( i = 0; i < Length; i++ ){
	  TempMatDouble[i] = (double) IGain->Array[i];
     }

     /* LAPACK: Compute the Cholesky factorization */
     dpptrf_( &uplo, &IGain->Rows, IGain->Array, &info );

     if ( info == 0 ){
	  Print_Header( SUCCESS );
	  printf( "Cholesky factorization successfully completed.\n" );
     } else if (info < 0){
	  Print_Header( ERROR );
	  fprintf( stderr, "Cholesky factorization: the %d-th argument has an illegal value.\n", -info );
	  exit( EXIT_FAILURE );
     } else if (info > 0){
	  Print_Header( ERROR );
	  fprintf( stderr, "Cholesky factorization: the leading minor of order %d is not positive definite,", info );
	  fprintf( stderr, " and the factorization could not be completed.\n" );
	  exit( EXIT_FAILURE );
     }

     /* LAPACK: Compute the inverse of the IGain matrix using the Cholesky factorization computed by
      *	dpptrf() */
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

     Scalar = Const.Lambda;
     dscal( (int *) &Length, &Scalar, IGain->Array, &incx );

#pragma omp parallel for
     for( i = 0; i < Length; i++ ){
	  IGain->Array[i] = (float) TempMatDouble[i];
     }

     /* Free the temporal double matrix */
     free( TempMatDouble );
     
     Print_Header( SUCCESS );
     printf( "Gain matrix successfully calculated.\n" );
}
