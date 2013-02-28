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
#include <stdlib.h>         /* For exit() */

#include "Print_Messages.h" /* For Print_Message() */
#include "Initiation.h"
#include "MatrixVector.h"

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
	  Print_Message( SUCCESS, 2, STRING, "Cholesky factorization successfully completed." );
     }
     else if (info < 0){
	  Print_Message( ERROR, 3, STRING, "Cholesky factorization: the ", INT, -info, STRING, "th argument has an illegal value." );
	  exit( EXIT_FAILURE );
     } else if (info > 0){
	  Print_Message( ERROR, 3, STRING, "Cholesky factorization: the leading minor of order ", INT, info,
			 STRING, " is not positive definite, and the factorization could not be completed." );
     }

     /* LAPACK: Compute the inverse of the IGain matrix using the Cholesky factorization computed
      * by dpotrf_( ) */
     dpotri_( &uplo, &IGain->Rows, IGain->Array, &lda, &info );

     if ( info == 0 ){
	  Print_Message( SUCCESS, 1, STRING, "Matrix Inversion successfully completed." );
	  exit( EXIT_FAILURE );
     } else if (info < 0){
	  Print_Message( ERROR, 3, STRING, "Matrix inversion: the ", INT, -info, STRING, "th argument has an illegal value." );
     } else if (info > 0){
	  Print_Message( ERROR, 5, STRING, "Matrix Inversion: the (", INT, info, STRING, "," ,INT, info,
			 STRING, ") element of the factor U or L is zero, and the inverse could not be computed.\n" );
	  exit( EXIT_FAILURE );
     }
     Scalar = Const.Lambda;
     dlascl_( &uplo, &ione, &ione, &one, &Scalar, &IGain->Rows, &IGain->Cols, IGain->Array, &lda, &info );
     if( info < 0 ){
	  Print_Message( ERROR, 3, STRING, "Gain matrix: the ", INT, -info, STRING, "th argument has an illegal value and cannot be scaled." );
	  exit( EXIT_FAILURE );
     } else {
	  Print_Message( SUCCESS, 1, STRING, "Gain matrix successfully calculated.\n" );
     }

}
