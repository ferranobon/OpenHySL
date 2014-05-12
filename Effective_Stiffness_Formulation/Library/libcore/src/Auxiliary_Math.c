#include <stdio.h>              /* For printf(), fprintf() */
#include <stdlib.h>             /* For exit() */
#include <stdbool.h>            /* For bool, true and false */
#include <math.h>               /* For sqrt(), log() */

#include "Auxiliary_Math.h"
#include "Print_Messages.h"

#include "Definitions.h"

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
     HYSL_FLOAT *work, *TempMat, temp;

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
     TempMat = (HYSL_FLOAT*) calloc( (size_t) Length, sizeof (HYSL_FLOAT) );
     work = (HYSL_FLOAT*) calloc( (size_t) lwork, sizeof (HYSL_FLOAT) );

     /* DSYGV_:On Entry EigenVectors must contain the Matrix A */
     hysl_copy( &Length, MatrixA->Array, &one, EigenVectors->Array, &one );
     hysl_copy( &Length, MatrixB->Array, &one, TempMat, &one );

     Length = MatrixA->Rows;
     hysl_sygv( &one, "V", "L", &Length, EigenVectors->Array, &lda, TempMat, &ldb, EigenValues->Array, work, &lwork, &info );

     if( info == 0 ){
	  Print_Header( SUCCESS );
	  printf( "Successfully calculated the eigenvalues and eigenvectors.\n" );
     } else if ( info < 0 ){
	  Print_Header( ERROR );
	  fprintf( stderr, "Compute_Eigenvalues_Eigenvectors: the %d-th argument of the function hysl_sygv() had an illegal value", info );
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

/* Routine based on ran1 from Numerical Receipes in C */
HYSL_FLOAT RandomNumber( long int *const idum )
{
     int j;
     long int k;
     static long int iy = 0;
     static long iv[NTAB];
     HYSL_FLOAT temp;

     if( *idum <= 0 || !iy ){
	  if( -(*idum) < 1 ){
	       *idum = 1;
	  } else {
	       *idum = -(*idum);
	  }

	  for( j = NTAB + 7; j>=0; j-- ){
	       k = (*idum)/IQ;
	       *idum = IA*(*idum-k*IQ) - IR*k;
	       if( *idum < 0 ){
		    *idum += IM;
	       }
	       if( j < NTAB ){
		    iv[j] = *idum;
	       }
	  }
	  iy = iv[0];
     }
  
     k = (*idum)/(long int)IQ;
     *idum = IA*(*idum - k*IQ) - IR*k;
     if( *idum < 0 ){
	  *idum += IM;
     }
  
     j = iy/(long int)NDIV;
     iy = iv[j];
     iv[j] = *idum;
     if( (temp = (HYSL_FLOAT) AM*(HYSL_FLOAT) iy) > (HYSL_FLOAT) RNMX ){
	  return (HYSL_FLOAT) RNMX;
     } else {
	  return temp;
     }

}

/* Routine based on Gasdev() from Numerical Receipes in C */
HYSL_FLOAT Gaussian_Deviate( const HYSL_FLOAT *const mu, const HYSL_FLOAT *const sigma, long int *const idum )
{
     static bool iset;
     static HYSL_FLOAT GD_value;
     HYSL_FLOAT fac, rsq, v1, v2;

     if( *idum < 0 ){   /* Reinitialise */
	  iset = false;
     }
  
     if( iset == false ){
	  do{
#if _FLOAT_
	       v1 = 2.0*RandomNumber( idum ) - 1.0;
	       v2 = 2.0*RandomNumber( idum ) - 1.0;
#else
	       v1 = 2.0*RandomNumber( idum ) - 1.0;
	       v2 = 2.0*RandomNumber( idum ) - 1.0;
#endif
	       rsq = v1*v1 + v2*v2;
	  } while ( rsq >= 1.0 || rsq == 0.0 );

#if _FLOAT_
	  fac = sqrtf(-2.0*logf(rsq)/rsq);
#else
	  fac = sqrt(-2.0*log(rsq)/rsq);
#endif

	  GD_value = (*mu) + (*sigma)*v1*fac;
	  iset = true;
	  return (*mu) + (*sigma)*v2*fac;
     } else {
	  iset = false;
	  return GD_value;
     }
}

