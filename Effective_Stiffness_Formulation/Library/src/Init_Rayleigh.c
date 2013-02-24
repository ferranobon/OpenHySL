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

#if _SPARSE_
#include <mkl.h>
#else
#include "Netlib.h"
#endif

void Rayleigh_Damping( const MatrixVector *const Mass, const MatrixVector *const Stif, MatrixVector *const Damp, const RayleighConst *const Rayleigh )
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
	  daxpy_( &Length, &beta, Stif->Array[i*Stif->Rows + i], &incx, Damp->Array[i*Damp->Rows +i], &incy);
     }

     Print_Message( SUCCESS, 1, STRING "Damping matrix successfully calculated." );
}

#if _SPARSE_
void Rayleigh_Damping_Sp( const MatrixVector_Sp *const Mass, const MatrixVector_Sp *const Stif, MatrixVector_Sp *const Damp, const RayleighConst *const Rayleigh )
{
     MatrixVector_Sp Temp;  /* Temporal matrix */
     int i;                 /* A counter */
     int Length;
     int incx, incy;        /* Stride in the operations */
     double alpha, beta;    /* Constants */
     char trans;
     int job, sort, info;   /* MKL variables */

     alpha = Rayleigh->Alpha;
     beta = Rayleigh->Beta;

     Init_MatrixVector_Sp( &Temp, Mass->Rows, Mass->Cols, Mass->Num_Nonzero );

     incx = 1; incy = 1;
     Length = Temp.Num_Nonzero;
     dcopy_( &Length, Mass->Values, &incx, Temp.Values, &incy );

     /* Copy the column array */
#pragma omp parallel for
     for (i = 0; i < Length; i++ ){
	  Temp.Columns[i] = Mass->Columns[i];
     }

     /* Copy the RowIndex array */
     Length = Temp.Rows + 1;
#pragma omp parallel for
     for (i = 0; i < Length; i++ ){
	  Temp.RowIndex[i] = Mass->RowIndex[i];
     }

     /* Scal the Values array */
     dscal_( &Length, &alpha, Temp.Values, &incx );

     trans = 'N';  /* The transpose of the matrix is not used */
     job = 0;      /* The routine computes the addition */
     sort = 0;     /* The routine does not perform any reordering */
     mkl_dcsradd( &trans, &job, &sort, &Temp.Rows, &Temp.Cols, Temp.Values, Temp.Columns, Temp.RowIndex,
		  &beta, Stif->Values, Stif->Columns, Stif->RowIndex,
		  Damp->Values, Damp->Columns, Damp->RowIndex, &Damp->Num_Nonzero, &info );

     /* Delete the previously allocated sparse matrix */
     Destroy_MatrixVector_Sparse( &Temp );

     if ( info > 0){
	  Print_Message( ERROR, 1, STRING, "Number of elements exceeded while calculating the Damping matrix." );
     } else if ( info < 0 ){
	  Print_Message( ERROR, 1, STRING, "I do not understand." );
     } else {
	  Print_Message( SUCCESS, 1, STRING, "Damping matrix successfully calculated." );
     }
}
#endif
