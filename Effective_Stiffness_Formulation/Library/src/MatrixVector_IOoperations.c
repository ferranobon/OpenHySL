#include <stdio.h>
#include <stdlib.h>

#include "MatrixVector.h"
#include "Print_Messages.h" /* For Print_Message() */

void MatrixVector_FromFile( const char *Filename, MatrixVector_t *const MatVec )
{

     FILE *InFile;

     InFile = fopen( Filename, "r" );

     if ( InFile != NULL ){

	  int i;
	  for ( i = 0; i < MatVec->Rows*MatVec->Cols; i++ ){
	       fscanf(InFile,"%lf", &MatVec->Array[i]);
	  }
	  fclose( InFile );
     } else{
	  PrintMessage( ERROR, 3, STRING, "MatrixVector_FromFile: It is not possible to open: ", STRING, Filename, STRING, "." );
	  exit( EXIT_FAILURE );
     }

}

#if _MATRIXMARKET_
void MatrixVector_FromFile_MM( const char *Filename, MatrixVector_t *const MatVec )
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
	  Print_Message( ERROR, 3, STRING, "MatrixVector_FromFIle_Sp2Dense: It is not possible to open: ", STRING, Filename,
			 STRING, "." );
	  exit( EXIT_FAILURE );
     }

     /* Read the banner and identify which type of matrix is in the file */
     if( mm_read_banner( InFile, &matcode ) != 0 ){
	  Print_Message( ERROR, 3, STRING, "MatrixVector_FromFile_Sp2Dense: Could not process Market Matrix banner in ", STRING, Filename,
			 STRING, "." );
	  exit( EXIT_FAILURE );
     }
     
     /* Only sparse matrices are accepted */
     if ( !mm_is_sparse(matcode) ){
	  Print_Message( ERROR, 1, STRING, "MatrixVector_FromFile_Sp2Dense: the matrix or vector should be of type sparse or dense for this application to work." );
	  Print_Message( ERROR, 3, STRING, "Specified Matrix Market type: [", STRING, mm_typecode_to_str(matcode), STRING, "]." );
	  exit( EXIT_FAILURE );
     }
     
     /* Get the sizes */
     if ( (return_code = mm_read_mtx_crd_size( InFile, &Rows, &Cols, &nnz)) !=0){
	  exit( EXIT_FAILURE );
     }

     /* Check if the dimensions of the matrices are the same */
     if ( Rows != MatVec->Rows || Cols != MatVec->Cols ){
	  Print_Message( ERROR, 10, STRING, "The sizes of the load vector (", INT, Rows, STRING, ",", INT, Cols, STRING, ") do not match with the"
			 STRING, " specified ones in the configuration file (", INT, MatVec->Rows, STRING, ",", INT, MatVec->Cols, STRING, ")." );
	  exit( EXIT_FAILURE );
     }

     /* Read the values. The MatrixMarket format imposes that the file should contain only the
      * lower part of the matrix in 1-based index. Since C and FORTRAN use row-major and column-major
      * ordering respectively, the matrices will be stored as upper part in the C so that when
      * calling the FORTRAN routines from BLAS they access the lower part of the matrix without
      * requiring transposing it.
      */
     for( innz = 0; innz < nnz; innz++ ){
	  fscanf( InFile, "%d %d %lE", &i, &j, &Value );
	  MatVec->Array[(j-1)*MatVec->Cols + (i-1)] = Value;
     }
}
#endif /* _MATRIXMARKET_ */

void MatrixVector_ToFile( const char *Filename, const MatrixVector_t *const MatVec )
{
     int i, j;      /* Counters */
     FILE *OutFile;

     OutFile = fopen( Filename, "w" );

     if ( OutFile != NULL ){
	  for ( i = 0; i < MatVec->Rows; i++){
	       for( j = 0; j < MatVec->Cols; j++ ){
		    fprintf(OutFile,"%le\t", MatVec->Array[i + j*MatVec->Rows]);
	       }
	       fprintf( OutFile, "\n" );
	  }
	  fclose( OutFile );
     } else{
	  Print_Message( ERROR, 3, STRING, "MatrixVector_ToFile: It is not possible to open: ", STRING,
			 Filename, STRING, "." );
     }
}
