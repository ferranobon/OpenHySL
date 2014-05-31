#include <stdio.h>          /* For printf(), fprintf() */
#include <stdlib.h>         /* For exit() */
#include <math.h>

#include "Auxiliary_Math.h" /* For Max() */
#include "MatrixVector.h"   /* MatrixVector definition */
#include "Print_Messages.h" /* For Print_Header() */
#include "Modal_Damping.h"       /* Rayleigh damping routines */
#include "Definitions.h"

#if _MKL_
#include <mkl_blas.h>
#include <mkl_lapack.h>
#else
#include "Netlib.h"
#endif

void Modal_Damping( const MatrixVector_t *const Mass, const MatrixVector_t *const Stiff, MatrixVector_t *const Damp,
		    double DampFactor )
{
     char uplo, transa, transb;
     int ione;
     int incx, incy;
     int info, lda, ldb;
     HYSL_FLOAT done;

     int Rows, Cols, Length;
     HYSL_FLOAT alpha, beta;

     MatrixVector_t EVectors;
     MatrixVector_t EValues;
     MatrixVector_t temp, temp1;

     int i;
     Print_Header( WARNING );
     fprintf( stderr, "Modal_Damping(): Untested routine." );

     ione = 1; done = 1.0;
     incx = 1; incy = 1;
     uplo = 'L';  /* The lower part of the matrix will be used and the upper part will strictly
		   * not be referenced */
     Rows = Mass->Rows;
     Cols = Mass->Cols;

     MatrixVector_Create( Rows, Cols, &EVectors );
     MatrixVector_Create( Rows, 1, &EValues );
     MatrixVector_Create( Rows, Cols, &temp );
     MatrixVector_Create( Rows, Cols, &temp1 );

     for( i = 0; i < Rows; i++ ){
	  if (Mass->Array[i*Rows +i] == 0.0){
	       temp1.Array[i*Rows + i] = 1E-12;
	  } else {
	       temp1.Array[i*Rows + i] = Mass->Array[i*Rows + i];
	  }
	  printf( "%lE %lE\n", Mass->Array[i*Rows + i], temp1.Array[i*Rows + i] );
     }

     Compute_Eigenvalues_Eigenvectors( Stiff, &temp1, &EValues, &EVectors );
     for( i = 0; i < Rows; i++ ){
	  Damp->Array[i*Rows + i] = 2.0*sqrt( EValues.Array[i] )*DampFactor;
     }

     alpha = 1.0;
     beta = 0.0;
     transa = 'T'; transb = 'N';
     dgemm_( &transa, &transb, &Rows, &Cols, &Rows, &alpha, EVectors.Array, &Rows, Damp->Array, &Rows, &beta, temp.Array, &Rows );
     transa = 'N'; transb = 'T';
     dgemm_( &transa, &transb, &Rows, &Cols, &Rows, &alpha, temp.Array, &Rows, EVectors.Array, &Rows, &beta, Damp->Array, &Rows );

     MatrixVector_ToFile( Damp, "Damp.txt" );

     lda = Max( 1, Rows );
     ldb = Max( 1, Damp->Rows );

     MatrixVector_Destroy( &EVectors );
     MatrixVector_Destroy( &EValues );
     MatrixVector_Destroy( &temp );
     MatrixVector_Destroy( &temp1 );

     Print_Header( SUCCESS );
     printf("Damping matrix successfully calculated.\n" );
}
