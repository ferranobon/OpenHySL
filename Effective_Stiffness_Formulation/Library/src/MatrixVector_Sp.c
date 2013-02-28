#include <stdio.h>
#include <stdlib.h>
#include <assert.h>  /* For assert() */

#include "MatrixVector_Sp.h"
#include "Print_Messages.h"
#include "Auxiliary_Math.h" /* For Max() */

#include <mkl_blas.h>
#include <mkl_spblas.h>

void MatrixVector_SetRowsCols_Sp( const int Rows, const int Cols, MatrixVector_Sp_t *const MatVec_Sp )
{
     
     if ( Rows > 0 && Cols > 0 ){ 
	  MatVec_Sp->Rows = Rows;
	  MatVec_Sp->Cols = Cols;
	  /* Set the number of non-zero elements to 0 to indicate
	   * that the memory is not yet allocated. */
	  MatVec_Sp->Num_Nonzero = 0;   
     } else {
	  Print_Message( ERROR, 1, STRING, "MatrixVector_SetRowsCols_Sp: The number of rows or columns must be greater than zero." );
	  exit( EXIT_FAILURE );
     }

}

void MatrixVector_AllocateSpace_Sp( const int nnz, MatrixVector_Sp_t *const MatVec_Sp )
{

     /* Only perform this operation if the number of non-zero elements is equal
      * to zero. A different number would mean that the matrix is already
      * initialised */
     if( MatVec_Sp->Num_Nonzero != 0 ){
	  Print_Message( ERROR, 1, STRING, "MatrixVector_AllocateSpace_Sp: Error when initialising a sparse matrix since it is already initialised." );
     } else {
	  MatVec_Sp->Num_Nonzero = nnz;

	  /* Allocate the memory for the Values and Columns matrices */
	  MatVec_Sp->Values = (double *) calloc( (size_t) MatVec_Sp->Num_Nonzero, sizeof(double) );
	  MatVec_Sp->Columns = (int *) calloc( (size_t) MatVec_Sp->Num_Nonzero, sizeof(int) );
	  /* Allocate the RowIndex matrix. Length = Rows + 1 */
	  MatVec_Sp->RowIndex = (int *) calloc( (size_t) MatVec_Sp->Rows + 1, sizeof(int) );
     }

}

void MatrixVector_Create_Sp( const int Rows, const int Cols, const int nnz, MatrixVector_Sp_t *const MatVec_Sp )
{
     MatrixVector_SetRowsCols_Sp( Rows, Cols, MatVec_Sp );
     MatrixVector_AllocateSpace_Sp( nnz, MatVec_Sp );
}

int MatrixVector_CountNNZ_GE( const double *const Matrix, const int Rows, const int Cols )
{
     int i, j;
     int Count;

     Count = 0;
     for ( i = 0; i < Rows; i++ ){
	  for ( j = 0; j < Cols; j++ ){
	       if ( Matrix[i*Cols + j] != 0.0f ){
		    Count = Count + 1;
	       }
	  }
     }

     return Count;
}

int MatrixVector_CountNNZ_SY( const double *const Sym_Matrix, const int Rows )
{
     int i, j;
     int Count;

     Count = 0;
     for ( i = 0; i < Rows; i++ ){
	  for ( j = 0; j < Rows; j++ ){
	       if ( Sym_Matrix[i*Rows + j] != 0.0 ){
		    Count = Count + 1;
	       }
	  }
     }

     return Count;
}

#if _MKL_
void MatrixVector_Dense2CSR( const MatrixVector_t *const MatVec, const int Operation, MatrixVector_Sp_t *const MatVec_Sp )
{
     int job[6], lda;  /* Variables required by mkl_sdnscsr_( ). lda = leading dimension of
			 * the matrix lda = max(1,Num_Rows). */
     int info;

     /* Copy the number of Rows and Columns */
     MatrixVector_SetRowsCols_Sp( MatVec->Rows, MatVec->Cols, MatVec_Sp );

     /* Count the number of non-zero elements */
     if ( Operation == 0 ){      /* Count the non-zero elements of a symmetric matrix */
	  MatVec_Sp->Num_Nonzero = MatrixVector_CountNNZ_SY( MatVec->Array, MatVec->Rows );
     } else if ( Operation == 1 ){ /* Count the non-zero elements of a general matrix */       
	  MatVec_Sp->Num_Nonzero = MatrixVector_CountNNZ_GE( MatVec->Array, MatVec->Rows, MatVec->Cols );
     } else assert( Operation < 0 || Operation > 1 );

     /* Allocate the necessary space for the Value and Columns arrays */
     MatrixVector_AllocateSpace_Sp( MatVec_Sp->Num_Nonzero, MatVec_Sp );

     /* MKL: Transform the dense matrix into a CSR-three array variation matrix */
     job[0] = 0; /* The matrix is converted to CSR format. */
     job[1] = 0; /* Zero-based indexing is used for the dense matrix. */
     job[2] = 1; /* One-based indexing for the sparse matrix is used. */

     if ( Operation == 0 ){ /* Symmetric matrix */
	  job[3] = 1; /* Values will contain the upper triangular part of the dense matrix. */
     } else if ( Operation == 1 ){ /* General matrix */
	  job[3] = 2; /* All the elements of the dense matrix will be considered */
     }

     job[4] = MatVec_Sp->Num_Nonzero; /* Maximum number of non-zero elements allowed. */
     job[5] = 1; /* Values, Columns and RowIndex arrays are generated. */
     lda = Max( 1, MatVec_Sp->Rows );
     mkl_ddnscsr( job, &MatVec_Sp->Rows, &MatVec_Sp->Cols, MatVec->Array, &lda, MatVec_Sp->Values, MatVec_Sp->Columns, MatVec_Sp->RowIndex, &info );

     if (info != 0 ){
	  Print_Message( ERROR, 1, STRING, "MatrixVector_Dense2CSR: An error occurred during the mkl_ddnscsr() operation." );
	  exit( EXIT_FAILURE );
     }
}

void MatrixVector_CSR2Dense( const MatrixVector_Sp_t *const MatVec_Sp,  const int MatVec_Type, MatrixVector_t *const MatVec )
{
     int job[6], lda;  /* Variables required by mkl_sdnscsr_( ). lda = leading dimension of
			 * the matrix lda = max(1,Num_Rows). */
     int info;

     if( MatVec_Sp->Rows != MatVec->Rows || MatVec_Sp->Cols != MatVec->Cols ){
	  Print_Message( ERROR, 1, STRING, "MatrixVector_CSR2Dense: The dimensions of the matrices must match." );
	  exit( EXIT_FAILURE );
     }

     /* MKL: Transform the dense matrix into a CSR-three array variation matrix */
     job[0] = 1; /* The matrix is converted to dense format. */
     job[1] = 0; /* Zero-based indexing is used for the dense matrix. */
     job[2] = 1; /* One-based indexing for the sparse matrix is used. */

     if ( MatVec_Type == 0 ){ /* Symmetric matrix */
	  job[3] = 1; /* Values will contain the upper triangular part of the sparse matrix. */
     } else if ( MatVec_Type == 1 ){ /* General matrix */
	  job[3] = 2; /* All the elements of the sparse matrix will be considered */
     } else assert( MatVec_Type == 0 || MatVec_Type == 1 );

     job[4] = MatVec_Sp->Num_Nonzero; /* Maximum number of non-zero elements allowed. */
     job[5] = 1; /* Values, Columns and RowIndex arrays are generated. */
     lda = Max( 1, MatVec->Rows );
     mkl_ddnscsr( job, &MatVec->Rows, &MatVec->Cols, MatVec->Array, &lda, MatVec_Sp->Values, MatVec_Sp->Columns, MatVec_Sp->RowIndex, &info );

     if (info != 0 ){
	  Print_Message( ERROR, 1, STRING, "MatrixVector_CSR2Dense: An error occurred during the mkl_ddnscsr() operation." );
	  exit( EXIT_FAILURE );
     }
}

void MatrixVector_Add3Mat_Sp( const MatrixVector_Sp_t *const MatA, const MatrixVector_Sp_t *const MatB, const MatrixVector_Sp_t *const MatC,
			      const Scalars_t Const, MatrixVector_Sp_t *const MatY )
{

     MatrixVector_Sp_t Temp;
     int Length, i;
     int incx, incy;
     double alpha, beta, gamma;
     char trans;
     int job, sort, info;

     alpha = Const.Alpha;
     beta = Const.Beta;
     gamma = Const.Gamma;

     MatrixVector_Create_Sp( MatA->Rows, MatA->Cols, MatA->Num_Nonzero, &Temp );

     incx = 1; incy = 1;
     Length = Temp.Num_Nonzero;
     dcopy( &Length, MatA->Values, &incx, Temp.Values, &incy );

#pragma omp parallel for
     for (i = 0; i < Length; i++ ){
	  Temp.Columns[i] = MatA->Columns[i];
     }

     /* Scal the Values array */
     dscal( &Length, &alpha, Temp.Values, &incx );

     /* Copy the RowIndex array */
     Length = Temp.Rows + 1;
#pragma omp parallel for
     for (i = 0; i < Length; i++ ){
	  Temp.RowIndex[i] = MatA->RowIndex[i];
     }

     trans = 'N';  /* The operation C = Temp + beta*B is performed */
     job = 0;      /* The routine computes the addition */
     sort = 0;     /* The routine does not perform any reordering */
     mkl_dcsradd( &trans, &job, &sort, &Temp.Rows, &Temp.Cols, Temp.Values, Temp.Columns, Temp.RowIndex,
		  &beta, MatB->Values, MatB->Columns, MatB->RowIndex,
		  MatY->Values, MatY->Columns, MatY->RowIndex, &MatY->Num_Nonzero, &info );

     if (info != 0 ){
	  Print_Message( ERROR, 1, STRING, "MatrixVector_Add3Mat_Sp: An error occurred during the mkl_dcsradd() operation." );
	  exit( EXIT_FAILURE );
     }

     /* Delete the previously allocated sparse matrix */
     MatrixVector_Destroy_Sp( &Temp );

     mkl_dcsradd( &trans, &job, &sort, &MatY->Rows, &MatY->Cols, MatY->Values, MatY->Columns, MatY->RowIndex,
		  &gamma, MatC->Values, MatC->Columns, MatC->RowIndex,
		  MatY->Values, MatY->Columns, MatY->RowIndex, &MatY->Num_Nonzero, &info );

     if (info != 0 ){
	  Print_Message( ERROR, 1, STRING, "MatrixVector_Add3Mat_Sp: An error occurred during the mkl_dcsradd() operation." );
	  exit( EXIT_FAILURE );
     }
     
}
#endif

void MatrixVector_Destroy_Sp( MatrixVector_Sp_t *const MatVec_Sp )
{

     /* Set the number of rows, columns and non-zero elements to 0 */
     MatVec_Sp->Rows = 0;
     MatVec_Sp->Cols = 0;
     MatVec_Sp->Num_Nonzero = 0;

     /* Free the memory */
     free( MatVec_Sp->Values );
     free( MatVec_Sp->Columns );
     free( MatVec_Sp->RowIndex );
}
