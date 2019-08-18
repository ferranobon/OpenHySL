#include <stdio.h>          /* For printf(), fprintf() */
#include <stdlib.h>         /* For exit() */

#include "GainMatrix.h"
#include "MatrixVector_MPI.h"
#include "Print_Messages.h"  /* For Print_Header() */

#include "Auxiliary_Math.h"
#include "Definitions.h"

#if _MKL_
#include <mkl_scalapack.h>
#else 
#include "Netlib.h"
#endif

void IGainMatrix_MPI( PMatrixVector_t *const IGain, PMatrixVector_t *const Mass,
		      PMatrixVector_t *const Damp, PMatrixVector_t *const Stiff,
		      const Scalars_t Const )
{
     
     char uplo;
     int ione, info;
     hysl_float_t Scalar;
     hysl_float_t one = 1.0;

     ione = 1;
     uplo = 'L';  /* The lower part of the matrix will be used; the upper part will strictly not be
		   * referenced */

     /* IGain = Const.Alpha*M + Const.Beta*C + Const.Gamma*K */
     PMatrixVector_Add3Mat( Mass, Damp, Stiff, Const, IGain );
	
     /* ScaLAPACK: Compute the Cholesky factorization of the symmetric positive definite matrix IGain */
     hysl_ppotrf( &uplo, &IGain->GlobalSize.Row, IGain->Array, &ione, &ione, IGain->Desc, &info );
	
     if ( info == 0 ){
	  Print_Header( SUCCESS );
	  printf( "Cholesky factorization successfully completed.\n" );
     } else if (info < 0){
	  Print_Header( ERROR );
	  fprintf( stderr, "Cholesky factorization: the %d-th argument has an illegal value.\n", -info );
	  exit( EXIT_FAILURE );
     } else if (info > 0){
	  Print_Header( ERROR );
	  fprintf( stderr, "Cholesky factorization: the leading minor of order %d is not positive definite, and the factorization could not be completed.\n", info );
	  exit( EXIT_FAILURE );
     }
     /* SCALAPACK: Compute the inverse of Me using the Cholesky factorization computed by pdpotrf_() */
     hysl_ppotri( &uplo, &IGain->GlobalSize.Row, IGain->Array, &ione, &ione, IGain->Desc, &info );

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
     hysl_plascl( &uplo, &one, &Scalar, &IGain->GlobalSize.Row, &IGain->GlobalSize.Col, IGain->Array, &ione, &ione, IGain->Desc, &info );
     
     if( info < 0 ){
	  Print_Header( ERROR );
	  fprintf( stderr, "Gain matrix: the %d-th argument has an illegal value and cannot be scaled.\n", -info );
	  exit( EXIT_FAILURE );
     } else {
	  Print_Header( SUCCESS );
	  printf( "Gain matrix successfully calculated.\n" );
     }


}
