#include <stdio.h>
#include <stdlib.h>

#include "Auxiliary_Math.h" /* For Max() */
#include "MatrixVector.h"
#include "Print_Messages.h" /* For Print_Header() */
#include "Definitions.h"

#if _MATRIXMARKET_
#include "mmio.h"
#endif

void MatrixVector_FromFile (const char *Filename, MatrixVector_t *const MatVec) {
    FILE *InFile = fopen(Filename, "r");

    if (InFile == NULL) {
        Print_Header( ERROR);
        fprintf( stderr, "MatrixVector_FromFile: It is not possible to open %s.\n", Filename);
        exit( EXIT_FAILURE);

    }

    for (int32_t idx = 0; idx < MatVec->Rows * MatVec->Cols; idx++) {
#if _FLOAT_
	  fscanf( InFile,"%f", &MatVec->Array[idx] );
#else
        fscanf(InFile, "%lf", &MatVec->Array[idx]);
#endif
    }
    fclose(InFile);

    Print_Header( SUCCESS);
    printf("MatrixVector_FromFile: Contents of %s successfully readen.\n", Filename);
}

#if _MATRIXMARKET_
void MatrixVector_FromFile_MM(const char *Filename, MatrixVector_t *const MatVec) {
     /* Open the file */
     FILE *InFile = fopen(Filename, "r");
     if (InFile == NULL) {
         Print_Header(ERROR);
         fprintf(stderr, "MatrixVector_FromFile_MM: It is not possible to open %s.\n", Filename);
         exit(EXIT_FAILURE);
     }

     MM_typecode matcode = "";   /* MatrixMarket: type of the matrix (symmetric, dense, complex, ...)  */
     /* Read the banner and identify which type of matrix is in the file */
     if (mm_read_banner(InFile, &matcode) != 0){
         Print_Header(ERROR);
         fprintf(stderr, "MatrixVector_FromFile_MM: Could not process Market Matrix banner in %s\n.", Filename);
         exit(EXIT_FAILURE);
     }
     
     /* Only sparse matrices are accepted */
     if (!mm_is_sparse(matcode)){
         Print_Header(ERROR);
         fprintf( stderr, "MatrixVector_FromFile_MM: the matrix or vector should be of type sparse for this application to work.\n" );
         Print_Header(ERROR);
         fprintf(stderr, "Specified Matrix Market type: %s.\n", mm_typecode_to_str(matcode));
         exit(EXIT_FAILURE);
     }
     
     /* Get the sizes */
     int32_t return_code = 0;   /* MatrixMarket: return code for the functions */
     int32_t Rows = 0;          /* Number of Rows */
     int32_t Cols = 0;          /* Number of Columns */
     int32_t nnz = 0;           /* Number of non-zero elements */
     if ((return_code = mm_read_mtx_crd_size( InFile, &Rows, &Cols, &nnz)) !=0) {
         exit(EXIT_FAILURE);
     }

     /* Check if the dimensions of the matrices are the same */
     if ((Rows != MatVec->Rows) || (Cols != MatVec->Cols)){
         Print_Header(ERROR);
         fprintf(stderr, "MatrixVector_From_File_MM: The sizes of the matrix or vector (%d,%d) ", Rows, Cols);
         fprintf(stderr, "do not match with the specified ones in the configuration file (%d,%d)\n", MatVec->Rows, MatVec->Cols);
         exit(EXIT_FAILURE);
     }

     /* Read the values. The MatrixMarket format imposes that the file should contain only the
      * lower part of the matrix in 1-based index. Since C and FORTRAN use row-major and column-major
      * ordering respectively, the matrices will be stored as upper part in the C so that when
      * calling the FORTRAN routines from BLAS they access the lower part of the matrix without
      * requiring transposing it.
      */
     for(int32_t innz = 0; innz < nnz; innz++ ){
         int32_t rowIdx = 0;
         int32_t colIdx = 0;
         hysl_float_t Value = 0.0;
#if _FLOAT_
         fscanf( InFile, "%d %d %E", &rowIdx, &colIdx, &Value );
#else
         fscanf( InFile, "%d %d %lE", &rowIdx, &colIdx, &Value );
#endif
         MatVec->Array[((colIdx - 1)*MatVec->Cols) + rowIdx - 1] = Value;
     }

     Print_Header(SUCCESS);
     printf("MatrixVector_FromFile_MM: Contents of %s successfully readen.\n", Filename);
}

void MatrixVector_ToFile_MM(const MatrixVector_t *const MatVec, const char *Filename){
    FILE *OutFile = fopen( Filename, "w" );

     if (OutFile == NULL){
         Print_Header(ERROR);
         fprintf(stderr, "MatrixVector_ToFile_MM: It is not possible to open %s.\n", Filename);
         exit(EXIT_FAILURE);
     }

     fprintf( OutFile, "%%%%MatrixMarket matrix coordinate real symmetric\n" );

     int32_t num_non_zero = 0;
     /* Count the non-zero elements */
     for (int32_t idx = 0; idx < MatVec->Rows; idx++){
         for (int32_t jdx = idx; jdx < MatVec->Cols; jdx++){
             if (MatVec->Array[(idx * MatVec->Cols) + jdx] != 0.0){
                 num_non_zero = num_non_zero + 1;
             }
         }
     }
     fprintf(OutFile, "%d %d %d\n", MatVec->Rows, MatVec->Cols, num_non_zero);

     /* Print the elements (lower part = Transpose) */
     for (int32_t idx = 0; idx < MatVec->Rows; idx++){
         for (int32_t jdx = idx; jdx < MatVec->Cols; jdx++){
             if( MatVec->Array[(idx * MatVec->Cols) + jdx] != 0.0){
                 fprintf( OutFile, "%d %d %.8lE\n", jdx + 1, idx + 1, MatVec->Array[(idx * MatVec->Cols) + jdx]);
             }
         }
     }

     fclose(OutFile);

     Print_Header( SUCCESS );
     printf( "MatrixVector_ToFile_MM: Matrix successfully saved to %s.\n", Filename );
}
#endif /* _MATRIXMARKET_ */

void MatrixVector_ToFile (const MatrixVector_t *const MatVec, const char *Filename) {
    FILE *OutFile = fopen(Filename, "w");

    if (OutFile == NULL) {
        Print_Header(ERROR);
        fprintf(stderr, "MatrixVector_ToFile: It is not possible to open %s.\n", Filename);
        exit(EXIT_FAILURE);
    }

    for (int32_t idx = 0; idx < MatVec->Rows; idx++) {
        for (int32_t jdx = 0; jdx < MatVec->Cols; jdx++) {
#if _FLOAT_
            fprintf(OutFile,"%E\t", MatVec->Array[idx + (jdx * MatVec->Rows)]);
#else
            fprintf(OutFile, "%lE\t", MatVec->Array[idx + (jdx * MatVec->Rows)]);
#endif
        }
        fprintf(OutFile, "\n");
    }
    fclose(OutFile);

    Print_Header( SUCCESS);
    printf("MatrixVector_ToFile: Matrix successfully saved to %s.\n", Filename);
}
