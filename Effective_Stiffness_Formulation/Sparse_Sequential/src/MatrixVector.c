/**
 * \file MatrixVector.c
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 10th of August 2011
 * \todo Add support for packaged storages to decrease the memory use
 *
 * \brief Source code of MatrixVector creation and manipulation functions.
 *
 * This file contains the source code of those functions involved in creating/destroying matrices and vectors. Essential
 * matrix/vector manipulations are also contemplated.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ErrorHandling.h"
#include "MatrixVector.h"
#include "Netlib.h"

#include <mkl_spblas.h>

void Init_Dense_MatrixVector( Dense_MatrixVector *const Mat, const int Rows, const int Cols )
{
     if ( Rows >= 0 ){ 
	  (*Mat).Rows = Rows;
	  if ( Rows > 0 ){
	       (*Mat).Cols = Cols;
	  } else {
	       (*Mat).Cols = 0;
	  }
     } else {
	  PrintErrorAndExit( "The number of rows must be equal or greater than zero" );
     }
	     
	(*Mat).Array = calloc( (*Mat).Rows*(*Mat).Cols, sizeof(float));
}

void Dense_to_CSR_SY( const Dense_MatrixVector *const Mat, Sp_MatrixVector *const Sp_Mat )
{
     int job[6], lda;  /* Variables required by mkl_sdnscsr_( ). lda = leading dimension of
			 * the matrix lda = max(1,Num_Rows). */

     /* Copy the number of Rows and Columns */
     Sp_Mat->Rows = Mat->Rows;
     Sp_Mat->Cols = Mat->Cols;

     /* Count the number of non-zero elements */
     Sp_Mat->Num_Nonzero = Count_Nonzero_Elements_SY( Mat->Array, Mat->Rows, Mat->Cols );
     
     /* Allocate the necessary space for the Value and Columns arrays */
     Sp_Mat->Values = (float *) calloc( Sp_Mat->Num_Nonzero, sizeof(float) );
     Sp_Mat->Columns = (int *) calloc( Sp_Mat->Num_Nonzero, sizeof(int) );

     /* Allocate memory for the RowIndex array */
     Sp_Mat->RowIndex = (int *) calloc( Mat->Rows, sizeof( int ) );

     /* MKL: Transform the dense matrix into a CSR-three array variation matrix */
     job[0] = 0; /* The matrix is converted to CSR format. */
     job[1] = 1; /* One-based indexing is used for the dense matrix. */
     job[2] = 1; /* One-based indexing for the sparse matrix is used. */
     job[3] = 1; /* Values will contain the upper triangular part of the dense matrix. */
     job[4] = Sp_Mat->Num_Nonzero; /* Maximum number of non-zero elements allowed. */
     job[5] = 1; /* Values, Columns and RowIndex arrays are generated. */
     lda = Max( 1, Sp_Mat->Rows );
     mkl_sdnscsr_( job, &Sp_Mat->Rows, &Sp_Mat->Cols, Mat->Array, &lda, Sp_Mat->Values, Sp_Mat->Columns, Sp_Mat->RowIndex );
}

int Count_Nonzero_Elements_SY( const float *const Sym_Matrix, const int Rows, const int Cols )
{
     int i, j;
     int Count;

     Count = 0;
     for ( i = 0; i < Rows; i++ ){
	  for ( j = i; j < Cols; j++ ){
	       if ( Sym_Matrix[i*Cols + j] != 0.0 ){
		    Count = Count + 1;
	       }
	  }
     }

     return Count;
}

void Dense_MatrixVector_From_File( Dense_MatrixVector *const Mat, const char *Filename )
{

	FILE *InFile;

	InFile = fopen( Filename, "r" );

	if ( InFile != NULL ){

		int i;
		for ( i = 0; i < (*Mat).Rows*(*Mat).Cols; i++ ){
			fscanf(InFile,"%f", &(*Mat).Array[i]);
		}
		fclose( InFile );
	} else ErrorFileAndExit( "It is not possible to read data because it was not possible to open: ", Filename );
}

void Set2Value( Dense_MatrixVector *const Mat, const float Value )
{

  static int incx = 0;
  static int incy = 1;
  static int Length;
  static float Val;

  Length = (*Mat).Rows*(*Mat).Cols;
  Val = Value;

  scopy_( &Length, &Val, &incx, (*Mat).Array, &incy );
}

void Modify_Element_DM( Dense_MatrixVector *const Mat, const int RowIndex, const int ColIndex, const float Value, const char *Operation )
{

  const char *OpSet = "Set";
  const char *OpAdd = "Add";
  const char *OpMult = "Multiply";
  const char *OpDiv = "Divide";

  if ( strcmp( Operation, OpSet ) == 0 ){
      (*Mat).Array[(RowIndex - 1)*(*Mat).Cols + (ColIndex - 1)] = Value;
  }else if ( strcmp( Operation, OpAdd ) == 0 ){
      (*Mat).Array[(RowIndex - 1)*(*Mat).Cols + (ColIndex - 1)] = (*Mat).Array[(RowIndex - 1)*(*Mat).Cols + (ColIndex - 1)] + Value;
  } else if ( strcmp( Operation, OpMult ) == 0 ){
      (*Mat).Array[(RowIndex - 1)*(*Mat).Cols + (ColIndex - 1)] = (*Mat).Array[(RowIndex - 1)*(*Mat).Cols + (ColIndex - 1)]*Value;
  } else if ( strcmp( Operation, OpDiv ) == 0 ){
	  (*Mat).Array[(RowIndex - 1)*(*Mat).Cols + (ColIndex - 1)] = (*Mat).Array[(RowIndex - 1)*(*Mat).Cols + (ColIndex - 1)]/Value;
  } else {
      fprintf( stderr, "%s: Operation not identified. Valid operations are:\n", Operation );
      fprintf( stderr, "\t1) %s\n", OpSet );
      fprintf( stderr, "\t1) %s\n", OpAdd );
      fprintf( stderr, "\t1) %s\n", OpMult );
      fprintf( stderr, "\t1) %s\n", OpDiv );
      exit( EXIT_FAILURE );
  }
}

void Dense_Add3Mat( Dense_MatrixVector *const MatY, const Dense_MatrixVector *const MatA, const Dense_MatrixVector *const MatB, const Dense_MatrixVector *const MatC, const Scalars Const )
{

	int i;           /* A counter */
	int Length;      /* Variable to store the size of the vector to multiply. */
	int incx, incy;  /* Stride in the vectors for BLAS library */
	int lda, ldb;
	int ione;
	float done, Alpha;     /* Constant to use in the BLAS library */
	char uplo;       /* BLAS & LAPACK: Character to specify which part of the matrix has been referenced. */
	int info;        /* LAPACK: Variable to inform if the operations of Cholesky factorization and inverse were successful or not */

	ione = 1;
	incx = 1; incy = 1;
	uplo = 'L';      /* Character defining that the lower part of the symmetric matrix is referenced (see man dpotrf) */
	done = 1.0;

	lda = Max( 1, (*MatA).Rows );
	ldb = Max( 1, (*MatY).Rows );

	/* LAPACK: Y = A */
	slacpy_( &uplo, &(*MatY).Rows, &(*MatY).Cols, (*MatA).Array, &lda, (*MatY).Array, &ldb );

	/* LAPACK: Calculates Y = A  */
	Alpha = Const.Alpha;
	slascl_( &uplo, &ione, &ione, &done, &Alpha, &(*MatY).Rows, &(*MatY).Cols, (*MatY).Array, &ldb, &info);

	if ( info < 0 ){
		LAPACKPErrorAndExit( "dlascl: The ", -info, "-th argument had an illegal value" );
	}

	/* BLAS: Calculates Y = Beta*B + Y. Only computes half of the matrix */
	Alpha = Const.Beta;
	for (i = 0; i < (*MatY).Rows; i++){
		Length = (*MatY).Rows - i;
		saxpy_( &Length, &Alpha, &(*MatB).Array[i*(*MatB).Rows + i], &incx, &(*MatY).Array[i*(*MatY).Rows +i], &incy);
	}

	/* BLAS: Calculates Y = Gamma*C + Y. Only computes half of the matrix */
	Alpha = Const.Gamma;
	for (i = 0; i < (*MatY).Rows; i++){
		Length = (*MatY).Rows - i;
		saxpy_( &Length, &Alpha, &(*MatC).Array[i*(*MatC).Rows + i], &incx, &(*MatY).Array[i*(*MatY).Rows +i], &incy);
	}
}

void Dense_MatrixVector_To_File( const Dense_MatrixVector *const Mat, const char *Filename )
{
	FILE *OutFile;

	OutFile = fopen( Filename, "w" );

	if ( OutFile != NULL ){

		int i;
		int j;
		for ( i = 0; i < (*Mat).Rows; i++){
			for( j = 0; j < (*Mat).Cols; j++ ){
				fprintf(OutFile,"%e\t", (*Mat).Array[i + j*(*Mat).Rows]);
			}
			fprintf( OutFile, "\n" );
		}
		fclose( OutFile );
	} else ErrorFileAndExit( "It is not possible to read data because it was not possible to open: ", Filename );
}

void Sp_MatrixVector_To_File_SY( const Sp_MatrixVector *const Sp_Mat, const char *Filename )
{
     int i, j;          /* Counters */
     int lda, job[6];   /* Variables required by mkl_sdnscsr_( ). lda = leading dimension of
			 * the matrix lda = max(1,Num_Rows). */
     float *Temp_Array; /* Temporal full array */
     FILE *OutFile;

     OutFile = fopen( Filename, "w" );

     if ( OutFile != NULL ){

	  /* Allocate memory for the temporal array */
	  Temp_Array = (float *) calloc( Sp_Mat->Rows*Sp_Mat->Cols, sizeof( float ) );

	  /* MKL: Transform the dense matrix into a CSR-three array variation matrix */
	  job[0] = 1; /* The matrix is restored from CSR format. */
	  job[1] = 1; /* One-based indexing is used for the dense matrix. */
	  job[2] = 1; /* One-based indexing for the sparse matrix is used. */
	  job[3] = 1; /* Values will contain the upper triangular part of the dense matrix. */
	  job[4] = Sp_Mat->Num_Nonzero; /* Maximum number of non-zero elements allowed. */
	  job[5] = 1; /* Values, Columns and RowIndex arrays are generated. */
	  lda = Max( 1, Sp_Mat->Rows );
	  mkl_sdnscsr_( job, &Sp_Mat->Rows, &Sp_Mat->Cols, Temp_Array, &lda, Sp_Mat->Values, Sp_Mat->Columns, Sp_Mat->RowIndex );

	  for ( i = 0; i < Sp_Mat->Rows; i++ ){
	       for ( j = 0; j < Sp_Mat->Cols; j++ ){
		    fprintf( OutFile, "%e\t", Temp_Array[i*Sp_Mat->Rows + j] );
	       }
	       fprintf( OutFile, "\n" );
	  }

	  /* Deallocate memory */
	  free( Temp_Array );

	  /* Close file */
	  fclose( OutFile );
	       
     } else ErrorFileAndExit( "It is not possible to read data because it was not possible to open: ", Filename );

}

void Destroy_Dense_MatrixVector( Dense_MatrixVector *const Mat)
{
	/* Set the number of rows and columns to 0 */
	(*Mat).Rows = 0;
	(*Mat).Cols = 0;
	free( (*Mat).Array );
}

void Destroy_Sparse_MatrixVector( Sp_MatrixVector *const Sp_Mat )
{

     /* Set the number of rows, columns and non-zero elements to 0 */
     Sp_Mat->Rows = 0;
     Sp_Mat->Cols = 0;
     Sp_Mat->Num_Nonzero = 0;

     /* Free the memory */
     free( Sp_Mat->Values );
     free( Sp_Mat->Columns );
     free( Sp_Mat->RowIndex );
}

int Max ( const int a, const int b )
{
	if ( a >= b ){
		return a;
	} else return b;
}
