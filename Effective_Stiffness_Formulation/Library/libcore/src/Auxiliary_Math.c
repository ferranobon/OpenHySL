#include <stdio.h>              /* For printf(), fprintf() */
#include <stdlib.h>             /* For exit() */

#include "Auxiliary_Math.h"
#include "Print_Messages.h"

#if _MKL_
#include <mkl_blas.h>
#include <mkl_lapack.h>
#else
#include "Netlib.h"
#endif

int Max ( const int a, const int b )
{
     if ( a >= b ){
	  return a;
     } else return b;
}

int Min ( const int a, const int b )
{
     if ( a <= b ){
	  return a;
     } else return b;
}

MatrixVector_t Generate_IdentityMatrix( int Rows, int Cols )
{
     MatrixVector_t Identity;
     unsigned short int i;

     if( Rows != Cols ){
	  Print_Header( ERROR );
	  fprintf( stderr, "Generate_IdentityMatrix: The number of rows and columns must be the same.\n" );
	  exit( EXIT_FAILURE );
     }

     MatrixVector_Create( Rows, Cols, &Identity );
     
     for( i = 0; i < Rows; i++ ){
	  Identity.Array[i + Rows*i] = 1.0;
     }

     return Identity;
}

unsigned int MatrixVector_ReturnIndex_UPS( const unsigned int RowIndex, const unsigned int ColIndex, const int n )
{
     unsigned int Index;

     if( RowIndex >= ColIndex ){
	  Index = RowIndex + (2*(unsigned int) n - ColIndex)*(ColIndex - 1)/2 - 1;
     } else {
	  Index = ColIndex + (2*(unsigned int) n - RowIndex)*(RowIndex - 1)/2 - 1;
     }
     return Index;
}

unsigned int MatrixVector_ReturnIndex_LPS( const unsigned int RowIndex, const unsigned int ColIndex )
{
     unsigned int Index;

     if( ColIndex >= RowIndex ){
	  Index = RowIndex + ColIndex*(ColIndex - 1)/2 - 1;
     } else {
	  Index = ColIndex + RowIndex*(RowIndex - 1)/2 - 1;
     }
     return Index;
}


void Compute_Eigenvalues_Eigenvectors ( MatrixVector_t *const MatrixA, MatrixVector_t *const MatrixB, MatrixVector_t *const EigenValues, MatrixVector_t *const EigenVectors )
{
     int i, j;   /* Counters */
     int one = 1, Length;
     int lda, ldb, info;
     int lwork; /* Dimension of the array work */
     double *work, *TempMat, temp;

     if( MatrixA->Rows != MatrixB->Rows || MatrixA->Cols != MatrixB->Cols ){
	  Print_Header( ERROR );
	  fprintf( stderr, "Compute_Eigenvalues_Eigenvectors: The matrices must be identical.\n" );
	  exit( EXIT_FAILURE );
     }

     Length = MatrixA->Rows;
     lwork = 3*Length - 1;
     lda = Max (1, Length);
     ldb = lda;

     Length = MatrixA->Rows*MatrixA->Cols;
     TempMat = (double*) calloc( (size_t) Length, sizeof (double) );
     work = (double*) calloc( (size_t) lwork, sizeof (double) );

     /* DSYGV_:On Entry EigenVectors must contain the Matrix A */
     dcopy( &Length, MatrixA->Array, &one, EigenVectors->Array, &one );
     dcopy( &Length, MatrixB->Array, &one, TempMat, &one );

     Length = MatrixA->Rows;
     dsygv_( &one, "V", "L", &Length, EigenVectors->Array, &lda, TempMat, &ldb, EigenValues->Array, work, &lwork, &info );

     if( info == 0 ){
	  Print_Header( SUCCESS );
	  printf( "Successfully calculated the eigenvalues and eigenvectors\n." );
     } else if ( info < 0 ){
	  Print_Header( ERROR );
	  fprintf( stderr, "Compute_Eigenvalues_Eigenvectors: the %d-th argument of the function dsygv_() had an illegal value", info );
	  exit( EXIT_FAILURE );
     } else if ( info > 0 ){
	  if ( info <= EigenVectors->Rows ){
	       Print_Header( ERROR );
	       fprintf( stderr, "Compute_Eigenvalues_Eigenvectors: %d off-diagonal elements of an intermediate tridiagonal form did not converge to zero.\n", info );
	       exit( EXIT_FAILURE );
	  } else {
	       Print_Header( ERROR );
	       fprintf( stderr, "Compute_Eigenvalues_Eigenvectors: the leading minor of order %d of MatrixB (TempMat) is not positive definite. The factorization of MatrixB could not be completed and no eigenvalues or eigenvectors were computed.\n", info - MatrixB->Rows );
	       exit( EXIT_FAILURE );
	  }
     }

     for ( i = 0; i < Length - 1; i++){
	  
	  /* Order the Eigenvalues and eigenvectors in ascendent order */
	  if ( EigenValues->Array[i] > EigenValues->Array[i+1] ){
	       /* Swap Eigenvalues */
	       temp = EigenValues->Array[i];
	       EigenValues->Array[i] = EigenValues->Array[i+1];
	       EigenValues->Array[i+1] = temp;
	       /* Now Swap Eigenvectors */
	       for( j = 0; j < Length; j++ ){
		    temp = EigenVectors->Array[Length*(i+1) + j];
		    EigenVectors->Array[Length*i + j] = EigenVectors->Array[Length*(i+1) + j];
		    EigenVectors->Array[Length*(i+1) + j] = temp;
	       }
	  }
     }

     /* Free the dynamically allocated memory */
     free( TempMat );
     free( work );
}
