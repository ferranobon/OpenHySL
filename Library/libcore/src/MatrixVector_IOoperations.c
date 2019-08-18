#include <stdio.h>
#include <stdlib.h>

#include "Auxiliary_Math.h" /* For Max() */
#include "MatrixVector.h"
#include "Print_Messages.h" /* For Print_Header() */
#include "Definitions.h"

#if _MATRIXMARKET_
#include "mmio.h"
#endif

void MatrixVector_FromFile( const char *Filename, MatrixVector_t *const MatVec )
{

     FILE *InFile;
     int i;         /* A counter */

     InFile = fopen( Filename, "r" );

     if ( InFile == NULL ){
	  Print_Header( ERROR );
	  fprintf( stderr, "MatrixVector_FromFile: It is not possible to open %s.\n", Filename );
	  exit( EXIT_FAILURE );

     }
     
     for ( i = 0; i < MatVec->Rows*MatVec->Cols; i++ ){
#if _FLOAT_
	  fscanf( InFile,"%f", &MatVec->Array[i] );
#else
	  fscanf( InFile,"%lf", &MatVec->Array[i] );
#endif
     }
     fclose( InFile );
     
     Print_Header( SUCCESS );
     printf( "MatrixVector_FromFile: Contents of %s successfully readen.\n", Filename );
}

#if _MATRIXMARKET_
void MatrixVector_FromFile_MM( const char *Filename, MatrixVector_t *const MatVec )
{
     FILE *InFile;          /* Input file */
     MM_typecode matcode;   /* MatrixMarket: type of the matrix (symmetric, dense, complex, ...)  */
     int return_code;       /* MatrixMarket: return code for the functions */
     int i, j;              /* Indexes of the position within the matrix of the readen value */
     hysl_float_t Value;    /* Value to be saved in the position (i,j) of the matrix */
     int Rows, Cols;        /* Number of Rows and Columns */
     int nnz;               /* Number of non-zero elements */
     int innz;              /* Counter for the number of non-zero elements */

     /* Open the file */
     InFile = fopen( Filename, "r" );
     if ( InFile == NULL) {
	  Print_Header( ERROR );
	  fprintf( stderr, "MatrixVector_FromFile_MM: It is not possible to open %s.\n", Filename );
	  exit( EXIT_FAILURE );
     }

     /* Read the banner and identify which type of matrix is in the file */
     if( mm_read_banner( InFile, &matcode ) != 0 ){
	  Print_Header( ERROR );
	  fprintf( stderr, "MatrixVector_FromFile_MM: Could not process Market Matrix banner in %s\n.", Filename );
	  exit( EXIT_FAILURE );
     }
     
     /* Only sparse matrices are accepted */
     if ( !mm_is_sparse(matcode) ){
	  Print_Header( ERROR );
	  fprintf( stderr, "MatrixVector_FromFile_MM: the matrix or vector should be of type sparse for this application to work.\n" );
	  Print_Header( ERROR );
	  fprintf( stderr, "Specified Matrix Market type: %s.\n", mm_typecode_to_str(matcode) );
	  exit( EXIT_FAILURE );
     }
     
     /* Get the sizes */
     if ( (return_code = mm_read_mtx_crd_size( InFile, &Rows, &Cols, &nnz)) !=0){
	  exit( EXIT_FAILURE );
     }

     /* Check if the dimensions of the matrices are the same */
     if ( Rows != MatVec->Rows || Cols != MatVec->Cols ){
	  Print_Header( ERROR );
	  fprintf( stderr, "MatrixVector_From_File_MM: The sizes of the matrix or vector (%d,%d) ", Rows, Cols );
	  fprintf( stderr, "do not match with the specified ones in the configuration file (%d,%d)\n", MatVec->Rows, MatVec->Cols );
	  exit( EXIT_FAILURE );
     }

     /* Read the values. The MatrixMarket format imposes that the file should contain only the
      * lower part of the matrix in 1-based index. Since C and FORTRAN use row-major and column-major
      * ordering respectively, the matrices will be stored as upper part in the C so that when
      * calling the FORTRAN routines from BLAS they access the lower part of the matrix without
      * requiring transposing it.
      */
     for( innz = 0; innz < nnz; innz++ ){
#if _FLOAT_
	  fscanf( InFile, "%d %d %E", &i, &j, &Value );
#else
	  fscanf( InFile, "%d %d %lE", &i, &j, &Value );
#endif
	  MatVec->Array[(j-1)*MatVec->Cols + (i-1)] = Value;
     }

     Print_Header( SUCCESS );
     printf( "MatrixVector_FromFile_MM: Contents of %s successfully readen.\n", Filename );
}

void MatrixVector_ToFile_MM( const MatrixVector_t *const MatVec, const char *Filename ){

     int i, j, Num_Nonzero;      /* Counters */
     FILE *OutFile;

     OutFile = fopen( Filename, "w" );

     if ( OutFile == NULL ){
	  Print_Header( ERROR );
	  fprintf( stderr, "MatrixVector_ToFile_MM: It is not possible to open %s.\n", Filename );
	  exit( EXIT_FAILURE );
     }

     fprintf( OutFile, "%%%%MatrixMarket matrix coordinate real symmetric\n" );

     Num_Nonzero = 0;
     /* Count the non-zero elements */
     for( i = 0; i < MatVec->Rows; i++ ){
	  for( j = i; j < MatVec->Cols; j++ ){
	       if( MatVec->Array[i*MatVec->Cols + j] != 0.0 ){
		    Num_Nonzero = Num_Nonzero + 1;
	       }
	  }
     }
     fprintf( OutFile, "%d %d %d\n", MatVec->Rows, MatVec->Cols, Num_Nonzero );

     /* Print the elements (lower part = Transpose) */
     for( i = 0; i < MatVec->Rows; i++ ){
	  for( j = i; j < MatVec->Cols; j++ ){
	       if( MatVec->Array[i*MatVec->Cols + j] != 0.0 ){
		    fprintf( OutFile, "%d %d %.8lE\n", j+1, i+1, MatVec->Array[i*MatVec->Cols + j] );
	       }
	  }
     }

     fclose( OutFile );
     
     Print_Header( SUCCESS );
     printf( "MatrixVector_ToFile_MM: Matrix successfully saved to %s.\n", Filename );
}
#endif /* _MATRIXMARKET_ */

void MatrixVector_ToFile( const MatrixVector_t *const MatVec, const char *Filename )
{
     int i, j;      /* Counters */
     FILE *OutFile;

     OutFile = fopen( Filename, "w" );

     if ( OutFile == NULL ){
	  Print_Header( ERROR );
	  fprintf( stderr, "MatrixVector_ToFile: It is not possible to open %s.\n", Filename );
	  exit( EXIT_FAILURE );
     }

     for ( i = 0; i < MatVec->Rows; i++){
	  for( j = 0; j < MatVec->Cols; j++ ){
#if _FLOAT_
	       fprintf( OutFile,"%E\t", MatVec->Array[i + j*MatVec->Rows] );
#else
	       fprintf( OutFile,"%lE\t", MatVec->Array[i + j*MatVec->Rows] );
#endif
	  }
	  fprintf( OutFile, "\n" );
     }
     fclose( OutFile );
     
     Print_Header( SUCCESS );
     printf( "MatrixVector_ToFile: Matrix successfully saved to %s.\n", Filename );
}
