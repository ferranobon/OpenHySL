#include <stdio.h>
#include <stdlib.h>

#include "MatrixVector_Sp.h"
#include "Print_Messages.h"
#include "Definitions.h"

#if _MATRIXMARKET_
#include "mmio.h"
#endif

#if _MATRIXMARKET_
void MatrixVector_FromFile_MM_Sp( const char *Filename, MatrixVector_Sp_t *const MatVec_Sp )
{
     FILE *InFile;          /* Input file */
     MM_typecode matcode;   /* MatrixMarket: type of the matrix (symmetric, dense, complex, ...)  */
     int return_code;       /* MatrixMarket: return code for the functions */
     int i, j;              /* Indexes of the position within the matrix of the readen value */
     HYSL_FLOAT Value;      /* Value to be saved in the position (i,j) of the matrix */
     int Rows, Cols;        /* Number of Rows and Columns */
     int nnz;               /* Number of non-zero elements */
     int innz;              /* Counter for the number of non-zero elements */
     int Pos_RI;            /* Counter for the RowIndex array */

     /* Open the file */
     InFile = fopen( Filename, "r" );
     if ( InFile == NULL) {
	  Print_Header( ERROR );
	  fprintf( stderr, "MatrixVector_FromFile_MM_Sp: It is not possible to open: %s.\n", Filename );
	  exit( EXIT_FAILURE );
     }

     /* Read the banner and identify which type of matrix is in the file */
     if( mm_read_banner( InFile, &matcode ) != 0 ){
	  Print_Header( ERROR );
	  fprintf( stderr, "MatrixVector_FromFile_MM_Sp: Could not process Market Matrix banner in %s\n.", Filename );
	  exit( EXIT_FAILURE );
     }

     /* Only sparse matrices are accepted */     
     if ( !mm_is_sparse(matcode) ){
	  Print_Header( ERROR );
	  fprintf( stderr, "MatrixVector_FromFile_MM_Sp: the matrix or vector should be of type sparse for this application to work.\n" );
	  Print_Header( ERROR );
	  fprintf( stderr, "Specified Matrix Market type: %s.\n", mm_typecode_to_str(matcode) );
	  exit( EXIT_FAILURE );
     }
     
     /* Get the sizes */
     if ( (return_code = mm_read_mtx_crd_size( InFile, &Rows, &Cols, &nnz)) !=0){
	  exit( EXIT_FAILURE );
     }

     /* Check if the dimensions of the matrices are the same */
     if ( Rows != MatVec_Sp->Rows || Cols != MatVec_Sp->Cols ){
	  Print_Header( ERROR );
	  fprintf( stderr, "MatrixVector_From_File_MM_Sp: The sizes of the matrix or vector (%d,%d) ", Rows, Cols );
	  fprintf( stderr, "do not match with the specified ones in the configuration file (%d,%d)\n", MatVec_Sp->Rows, MatVec_Sp->Cols );
	  exit( EXIT_FAILURE );
     }

     MatrixVector_AllocateSpace_Sp( nnz, MatVec_Sp );

     /* Read the values. The MatrixMarket format imposes that the file should contain only the
      * lower part of the matrix in 1-based index. Since C and FORTRAN use row-major and column-major
      * ordering respectively, the matrices will be stored as upper part in the C so that when
      * calling the FORTRAN routines from BLAS they access the lower part of the matrix without
      * requiring transposing it.
      */
     Pos_RI = 0;
     for( innz = 0; innz < nnz; innz++ ){
#if _FLOAT_
	  fscanf( InFile, "%d %d %E", &i, &j, &Value );
#else
	  fscanf( InFile, "%d %d %lE", &i, &j, &Value );
#endif
	  MatVec_Sp->Values[innz] = Value;
	  MatVec_Sp->Columns[innz] = i;            /* The read matrix is the lower triangular part but we store the upper triangular part */
	  if ( j > Pos_RI ){
	       while ( Pos_RI < j ){
		    MatVec_Sp->RowIndex[Pos_RI] = innz + 1;
		    Pos_RI = Pos_RI + 1;
	       }
	  }
     }

     while( Pos_RI <= MatVec_Sp->Rows ){
	  MatVec_Sp->RowIndex[Pos_RI] = innz + 1;
	  Pos_RI = Pos_RI + 1;
     }

     Print_Header( SUCCESS );
     printf( "MatrixVector_FromFile_MM_Sp: Contents of %s successfully readen.\n", Filename );
}
#endif /* _MATRIXMARKET_ */

void MatrixVector_ToFile_Sp( const MatrixVector_Sp_t *const MatVec_Sp, const char *Filename )
{
     int i;          /* Counter */

     FILE *OutFile;

     OutFile = fopen( Filename, "w" );

     if ( OutFile == NULL ){
	  Print_Header( ERROR );
	  fprintf( stderr, "MatrixVector_ToFile_Sp: It is not possible to open %s.\n", Filename );
     }

	  
     fprintf( OutFile, "Values: Nonzero elements.\n" );
     for ( i = 0; i < MatVec_Sp->Num_Nonzero; i++ ){
#if _FLOAT_
	  fprintf( OutFile, "%E\t", MatVec_Sp->Values[i] );
#else
	  fprintf( OutFile, "%lE\t", MatVec_Sp->Values[i] );
#endif
     }
     fprintf( OutFile, "\n" );

     fprintf( OutFile, "Columns array.\n" );
     for ( i = 0; i < MatVec_Sp->Num_Nonzero; i++ ){
	  fprintf( OutFile, "%i\t", MatVec_Sp->Columns[i] );
     }
     fprintf( OutFile, "\n" );

     fprintf( OutFile, "RowIndex array.\n" );
     for ( i = 0; i < MatVec_Sp->Rows + 1; i++ ){
	  fprintf( OutFile, "%i\t", MatVec_Sp->RowIndex[i] );
     }
     fprintf( OutFile, "\n" );	  

     /* Close file */
     fclose( OutFile );

     Print_Header( SUCCESS );
     printf( "MatrixVector_ToFile_Sp: Sparse matrix successfully saved to %s in CSR three array variation format.\n", Filename );

}
