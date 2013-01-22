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
#include <assert.h>

#include "Netlib.h"
#include "ErrorHandling.h"
#include "MatrixVector.h"
#include "Colors.h"

/* MatrixMarket format */
#include "mmio.h"

#if _SPARSE_
#include <mkl.h>
#include <mkl_spblas.h>
#endif

void Init_MatrixVector( MatrixVector *const Mat, const int Rows, const int Cols )
{
     if ( Rows >= 0 ){ 
	  Mat->Rows = Rows;
	  if ( Rows > 0 ){
	       Mat->Cols = Cols;
	  } else {
	       Mat->Cols = 0;
	  }
     } else {
	  PrintErrorAndExit( "The number of rows must be equal or greater than zero" );
     }
     Mat->Array = NULL;
     Mat->Array = (double *) calloc( (size_t) Mat->Rows*(size_t) Mat->Cols, sizeof(double));
     if ( Mat->Array == NULL ){
	  PrintErrorAndExit( "Out of memory\n");
     }
}

void Init_MatrixVector_Sp( Sp_MatrixVector *const Mat, const int Rows, const int Cols, const int nnz )
{
     if ( Rows >= 0 ){ 
	  Mat->Rows = Rows;
	  if ( Rows > 0 ){
	       Mat->Cols = Cols;
	  } else {
	       Mat->Cols = 0;
	  }
     } else {
	  PrintErrorAndExit( "The number of rows must be equal or greater than zero" );
     }
     
     if ( nnz == 0 ){
	  PrintErrorAndExit( "The number of non-zero elements must be greater that zero" );
     } else {
	  Mat->Num_Nonzero = nnz;
     }
     
     /* Allocate the memory for the Values and Columns matrices */
     Mat->Values = (double *) calloc( (size_t) Mat->Num_Nonzero, sizeof(double) );
     Mat->Columns = (int *) calloc( (size_t) Mat->Num_Nonzero, sizeof(int) );

     /* Allocate the RowIndex matrix. Length = Rows + 1 */
     Mat->RowIndex = (int *) calloc( (size_t) Mat->Rows + 1, sizeof(int) );
}

#if _SPARSE_
void Dense_to_CSR( const MatrixVector *const Mat, Sp_MatrixVector *const Sp_Mat, const int Operation )
{
     int job[6], lda;  /* Variables required by mkl_sdnscsr_( ). lda = leading dimension of
			 * the matrix lda = max(1,Num_Rows). */
     int info;

     /* Copy the number of Rows and Columns */
     Sp_Mat->Rows = Mat->Rows;
     Sp_Mat->Cols = Mat->Cols;

     /* Count the number of non-zero elements */
     if ( Operation == 0 ){      /* Count the non-zero elements of a symmetric matrix */
	  Sp_Mat->Num_Nonzero = Count_Nonzero_Elements_SY( Mat->Array, Mat->Rows );
     } else if ( Operation == 1 ){ /* Count the non-zero elements of a general matrix */       
	  Sp_Mat->Num_Nonzero = Count_Nonzero_Elements_GE( Mat->Array, Mat->Rows, Mat->Cols );
     } else assert( Operation < 0 || Operation > 1 );

     /* Allocate the necessary space for the Value and Columns arrays */
     Sp_Mat->Values = (double *) calloc( (size_t) Sp_Mat->Num_Nonzero, sizeof(double) );
     Sp_Mat->Columns = (int *) calloc( (size_t) Sp_Mat->Num_Nonzero, sizeof(int) );

     /* Allocate memory for the RowIndex array */
     Sp_Mat->RowIndex = (int *) calloc( (size_t) Mat->Rows + 1, sizeof( int ) );

     /* MKL: Transform the dense matrix into a CSR-three array variation matrix */
     job[0] = 0; /* The matrix is converted to CSR format. */
     job[1] = 0; /* Zero-based indexing is used for the dense matrix. */
     job[2] = 1; /* One-based indexing for the sparse matrix is used. */

     if ( Operation == 0 ){ /* Symmetric matrix */
	  job[3] = 1; /* Values will contain the upper triangular part of the dense matrix. */
     } else if ( Operation == 1 ){ /* General matrix */
	  job[3] = 2; /* All the elements of the dense matrix will be considered */
     }

     job[4] = Sp_Mat->Num_Nonzero; /* Maximum number of non-zero elements allowed. */
     job[5] = 1; /* Values, Columns and RowIndex arrays are generated. */
     lda = Max( 1, Sp_Mat->Rows );
     mkl_ddnscsr( job, &Sp_Mat->Rows, &Sp_Mat->Cols, Mat->Array, &lda, Sp_Mat->Values, Sp_Mat->Columns, Sp_Mat->RowIndex, &info );

     if (info != 0 ){
	  printf("ERROROORRRR\n");
     }
}
#endif

int Count_Nonzero_Elements_SY( const double *const Sym_Matrix, const int Rows )
{
     int i, j;
     int Count;

     Count = 0;
     for ( i = 0; i < Rows; i++ ){
	  for ( j = i; j < Rows; j++ ){
	       if ( Sym_Matrix[i*Rows + j] != 0.0f ){
		    Count = Count + 1;
	       }
	  }
     }

     return Count;
}

int Count_Nonzero_Elements_GE( const double *const Matrix, const int Rows, const int Cols )
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

void MatrixVector_From_File( MatrixVector *const Mat, const char *Filename )
{

     FILE *InFile;

     InFile = fopen( Filename, "r" );

     if ( InFile != NULL ){

	  int i;
	  for ( i = 0; i < Mat->Rows*Mat->Cols; i++ ){
	       fscanf(InFile,"%lf", &Mat->Array[i]);
	  }
	  fclose( InFile );
     } else ErrorFileAndExit( "It is not possible to read data because it was not possible to open: ", Filename );


}

void MatrixVector_From_File_Sp2Dense( MatrixVector *const Mat, const char *Filename )
{
     FILE *InFile;          /* Input file */
     MM_typecode matcode;   /* MatrixMarket: type of the matrix (symmetric, dense, complex, ...)  */
     int return_code;       /* MatrixMarket: return code for the functions */
     int i, j;              /* Indexes of the position within the matrix of the readen value */
     double Value;           /* Value to be saved in the position (i,j) of the matrix */
     int Rows, Cols;        /* Number of Rows and Columns */
     int nnz;               /* Number of non-zero elements */
     int innz;              /* Counter for the number of non-zero elements */

     /* Open the file */
     InFile = fopen( Filename, "r" );
     if ( InFile == NULL) {
	  ErrorFileAndExit( "Could not read the Load Vector Form. Failed to open: ",
			    Filename );
     }

     /* Read the banner and identify which type of matrix is in the file */
     if( mm_read_banner( InFile, &matcode ) != 0 ){
	  ErrorFileAndExit( "Could not process Market Matrix banner in ", Filename );
     }
     
     /* Only sparse matrices are accepted */
     if ( !mm_is_sparse(matcode) ){
	  fprintf( stderr, "[ " RED "ERROR" RESET " ] The Load vector form should be of");
	  fprintf( stderr, "  type sparse or dense for this application to work\n" );
	  fprintf( stderr, "[ " RED "ERROR" RESET " ] Market Market type: [%s]\n",
		   mm_typecode_to_str(matcode));
	  exit( EXIT_FAILURE );
     }
     
     /* Get the sizes */
     if ( (return_code = mm_read_mtx_crd_size( InFile, &Rows, &Cols, &nnz)) !=0){
	  exit( EXIT_FAILURE );
     }

     /* Check if the dimensions of the matrices are the same */
     if ( Rows != Mat->Rows || Cols != Mat->Cols ){
	  fprintf( stderr, "[ " RED "ERROR" RESET " ] The sizes of the load vector (%d, %d)", Rows, Cols ); 
	  fprintf( stderr, "do not match with the specified ones in the configuration file (%d, %d).n",
		   Mat->Rows, Mat->Cols );
	  exit( EXIT_FAILURE );
     }

     /* Read the values. The MatrixMarket format imposes that the file should contain only the
      * lower part of the matrix in 1-based index. Since C and FORTRAN use row-major and column-major
      * ordering respectively, the matrices will be stored as upper part in the C so that when
      * calling the FORTRAN routines from BLAS they access the lower part of the matrix without
      * requiring transposing it.
      */
     for( innz = 0; innz < nnz; innz++ ){
	  fscanf( InFile, "%d %d %lf", &i, &j, &Value );
	  Mat->Array[(j-1)*Mat->Cols + (i-1)] = Value;
     }
}

#if _SPARSE_
void MatrixVector_From_File_Sp( Sp_MatrixVector *const Mat, const char *Filename )
{
     FILE *InFile;  /* Input file */
     int i, j;      /* Indexes of the position within the matrix of the readen value */
     char d;        /* Dump character between two values */
     double Value;   /* Value to be saved in the position (i,j) of the matrix */
     int Position;  /* Counter for the Values and columns array */
     int Pos_RI;    /* Counter for the RowIndex array */
     int innz;      /* Counter for the number of non-zeros */
     int Rows, Cols, nnz;       /* Number of rows, columns and non-zero elements */

     InFile = fopen( Filename, "r" );

     if( InFile != NULL ){
	  /* Read the number of rows, columns and non-zero elements */
	  fscanf( InFile, "%d %d %d", &Rows, &Cols, &nnz );

	  Init_MatrixVector_Sp( Mat, Rows, Cols, nnz );
	  Position = 0;
	  Pos_RI = 0;
	  innz = 1;
	  Mat->RowIndex[Pos_RI] = innz;

	  while( innz <= Mat->Num_Nonzero ) { 
	       fscanf( InFile, "%i%c%i%c%le", &i, &d, &j, &d, &Value );
	       if ( j >= i ){    /* Consider only the upper part */
		    innz = innz + 1;

		    Mat->Values[Position] = Value;
		    Mat->Columns[Position] = j + 1;  /* One based index */
		    Position = Position + 1;
		    if ( i > Pos_RI ){
			 while ( Pos_RI < i ){
			      Pos_RI = Pos_RI + 1;
			      Mat->RowIndex[Pos_RI] = innz - 1;
			 }
		    }

	       }        
	  }

	  /* Add the number of non-zero elements at the final position of the
	   * RowIndex array */
	  Pos_RI = Pos_RI + 1;
	  while( Pos_RI <= Rows ){
	       Mat->RowIndex[Pos_RI] = innz;
	       Pos_RI = Pos_RI + 1;
	  }

	  /* The program has reached the end of the file */
	  fclose( InFile );

     } else ErrorFileAndExit( "It is not possible to read data because it was not possible to open: ", Filename );
}
#endif

void Set2Value( MatrixVector *const Mat, const double Value )
{

     static int incx = 0;
     static int incy = 1;
     static int Length;
     static double Val;

     Length = Mat->Rows*Mat->Cols;
     Val = Value;

     dcopy_( &Length, &Val, &incx, Mat->Array, &incy );
}

void Modify_Element( MatrixVector *const Mat, const int RowIndex, const int ColIndex, const double Value, const char *Operation )
{

     const char *OpSet = "Set";
     const char *OpAdd = "Add";
     const char *OpMult = "Multiply";
     const char *OpDiv = "Divide";

     if ( strcmp( Operation, OpSet ) == 0 ){
	  Mat->Array[(RowIndex - 1)*Mat->Cols + (ColIndex - 1)] = Value;
     }else if ( strcmp( Operation, OpAdd ) == 0 ){
	  Mat->Array[(RowIndex - 1)*Mat->Cols + (ColIndex - 1)] = Mat->Array[(RowIndex - 1)*Mat->Cols + (ColIndex - 1)] + Value;
     } else if ( strcmp( Operation, OpMult ) == 0 ){
	  Mat->Array[(RowIndex - 1)*Mat->Cols + (ColIndex - 1)] = Mat->Array[(RowIndex - 1)*Mat->Cols + (ColIndex - 1)]*Value;
     } else if ( strcmp( Operation, OpDiv ) == 0 ){
	  Mat->Array[(RowIndex - 1)*Mat->Cols + (ColIndex - 1)] = Mat->Array[(RowIndex - 1)*Mat->Cols + (ColIndex - 1)]/Value;
     } else {
	  fprintf( stderr, "%s: Operation not identified. Valid operations are:\n", Operation );
	  fprintf( stderr, "\t1) %s\n", OpSet );
	  fprintf( stderr, "\t1) %s\n", OpAdd );
	  fprintf( stderr, "\t1) %s\n", OpMult );
	  fprintf( stderr, "\t1) %s\n", OpDiv );
	  exit( EXIT_FAILURE );
     }
}

void Add3Mat( MatrixVector *const MatY, const MatrixVector *const MatA, const MatrixVector *const MatB, const MatrixVector *const MatC, const Scalars Const )
{

     int i;           /* A counter */
     int Length;      /* Variable to store the size of the vector to multiply. */
     int incx, incy;  /* Stride in the vectors for BLAS library */
     int lda, ldb;
     int ione;
     double done, Alpha;     /* Constant to use in the BLAS library */
     char uplo;       /* BLAS & LAPACK: Character to specify which part of the matrix has been referenced. */
     int info;        /* LAPACK: Variable to inform if the operations of Cholesky factorization and inverse were successful or not */

     ione = 1;
     incx = 1; incy = 1;
     uplo = 'L';      /* Character defining that the lower part of the symmetric matrix is referenced (see man dpotrf) */
     done = 1.0;

     lda = Max( 1, (*MatA).Rows );
     ldb = Max( 1, (*MatY).Rows );

     /* LAPACK: Y = A */
     dlacpy_( &uplo, &(*MatY).Rows, &(*MatY).Cols, (*MatA).Array, &lda, (*MatY).Array, &ldb );

     /* LAPACK: Calculates Y = A  */
     Alpha = Const.Alpha;
     dlascl_( &uplo, &ione, &ione, &done, &Alpha, &(*MatY).Rows, &(*MatY).Cols, (*MatY).Array, &ldb, &info);

     if ( info < 0 ){
	  LAPACKPErrorAndExit( "dlascl: The ", -info, "-th argument had an illegal value" );
     }

     /* BLAS: Calculates Y = Beta*B + Y. Only computes half of the matrix */
     Alpha = Const.Beta;
     for (i = 0; i < (*MatY).Rows; i++){
	  Length = (*MatY).Rows - i;
	  daxpy_( &Length, &Alpha, &(*MatB).Array[i*(*MatB).Rows + i], &incx, &(*MatY).Array[i*(*MatY).Rows +i], &incy);
     }

     /* BLAS: Calculates Y = Gamma*C + Y. Only computes half of the matrix */
     Alpha = Const.Gamma;
     for (i = 0; i < (*MatY).Rows; i++){
	  Length = (*MatY).Rows - i;
	  daxpy_( &Length, &Alpha, &(*MatC).Array[i*(*MatC).Rows + i], &incx, &(*MatY).Array[i*(*MatY).Rows +i], &incy);
     }
}

void Add3Mat_Sparse( Sp_MatrixVector *const MatY, const Sp_MatrixVector *const MatA, const Sp_MatrixVector *const MatB, const Sp_MatrixVector *const MatC, const Scalars Const )
{

     Sp_MatrixVector Temp;
     int Length, i;
     int incx, incy;
     double alpha, beta, gamma;
     char trans;
     int job, sort, info;

     alpha = Const.Alpha;
     beta = Const.Beta;
     gamma = Const.Gamma;

     Init_MatrixVector_Sp( &Temp, MatA->Rows, MatA->Cols, MatA->Num_Nonzero );

     incx = 1; incy = 1;
     Length = Temp.Num_Nonzero;
     dcopy_( &Length, MatA->Values, &incx, Temp.Values, &incy );

#pragma omp parallel for
     for (i = 0; i < Length; i++ ){
	  Temp.Columns[i] = MatA->Columns[i];
     }

     /* Scal the Values array */
     dscal_( &Length, &alpha, Temp.Values, &incx );

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

     /* Delete the previously allocated sparse matrix */
     Destroy_MatrixVector_Sparse( &Temp );

     mkl_dcsradd( &trans, &job, &sort, &MatY->Rows, &MatY->Cols, MatY->Values, MatY->Columns, MatY->RowIndex,
		  &gamma, MatC->Values, MatC->Columns, MatC->RowIndex,
		  MatY->Values, MatY->Columns, MatY->RowIndex, &MatY->Num_Nonzero, &info );

     
}

void MatrixVector_To_File( const MatrixVector *const Mat, const char *Filename )
{
     FILE *OutFile;

     OutFile = fopen( Filename, "w" );

     if ( OutFile != NULL ){

	  int i;
	  int j;
	  for ( i = 0; i < Mat->Rows; i++){
	       for( j = 0; j < Mat->Cols; j++ ){
		    fprintf(OutFile,"%le\t", Mat->Array[i + j*Mat->Rows]);
	       }
	       fprintf( OutFile, "\n" );
	  }
	  fclose( OutFile );
     } else ErrorFileAndExit( "It is not possible to read data because it was not possible to open: ", Filename );
}

void MatrixVector_To_File_Sparse( const Sp_MatrixVector *const Sp_Mat, const char *Filename )
{
     int i;          /* Counter */

     FILE *OutFile;

     OutFile = fopen( Filename, "w" );

     if ( OutFile != NULL ){
	  
	  fprintf( OutFile, "Values: Nonzero elements.\n" );
	  for ( i = 0; i < Sp_Mat->Num_Nonzero; i++ ){
	       fprintf( OutFile, "%lf\t", Sp_Mat->Values[i] );
	  }
	  fprintf( OutFile, "\n" );

	  fprintf( OutFile, "Columns array.\n" );
	  for ( i = 0; i < Sp_Mat->Num_Nonzero; i++ ){
	       fprintf( OutFile, "%i\t", Sp_Mat->Columns[i] );
	  }
	  fprintf( OutFile, "\n" );

	  fprintf( OutFile, "RowIndex array.\n" );
	  for ( i = 0; i < Sp_Mat->Rows + 1; i++ ){
	       fprintf( OutFile, "%i\t", Sp_Mat->RowIndex[i] );
	  }
	  fprintf( OutFile, "\n" );	  

	  /* Close file */
	  fclose( OutFile );
	       
     } else ErrorFileAndExit( "It is not possible to read data because it was not possible to open: ", Filename );

}

void Destroy_MatrixVector( MatrixVector *const Mat)
{
     /* Set the number of rows and columns to 0 */
     Mat->Rows = 0;
     Mat->Cols = 0;
     free( Mat->Array );
}

void Destroy_MatrixVector_Sparse( Sp_MatrixVector *const Sp_Mat )
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
