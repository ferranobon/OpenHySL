#include <stdio.h>
#include <stdlib.h>

#include "Auxiliary_Math.h"
#include "MatrixVector_MPI.h"
#include "Print_Messages.h"

#if _MATRIXMARKET_
#include "mmio.h"
#endif

#if _MKL_
#include "mkl_pblas.h"
#include "Cblacs.h"
#include "Scalapack_Aux.h"
#else
#include "Netlib.h"
#endif

void PMatrixVector_FromFile( const char *Filename, PMatrixVector_t *const Mat )
{

     int i, j;		/* Counters */
     int myrow, mycol;	/* Grid variables */
     int nprow, npcol;
     
     char trans = 'N';
     int izero = 0, ione = 1;
     double done = 1.0, dzero = 0.0;
     int lld, info;
     
     double *LocalMatrix;
     int descLocal[9];
     
     FILE *InFile;
     
     
     Cblacs_gridinfo( Mat->Desc[1], &nprow, &npcol, &myrow, &mycol );
     
     lld = Max(1, Mat->GlobalSize.Row);
     
     descinit_( descLocal, &Mat->GlobalSize.Row, &Mat->GlobalSize.Col, &Mat->GlobalSize.Row,
		&Mat->GlobalSize.Col, &izero, &izero, &Mat->Desc[1], &lld, &info );
     
     if ( info < 0 ){
	  Print_Header( ERROR );
	  fprintf( stderr, "descinit: The %d-th argument had an illegal value.\n", -info );
	  exit( EXIT_FAILURE );
     }
     
     if ( myrow == 0 && mycol == 0 ){
	  LocalMatrix = calloc( (size_t) Mat->GlobalSize.Row * (size_t)Mat->GlobalSize.Col, sizeof(double) );
	  if( LocalMatrix == NULL ){
	       Print_Header( ERROR );
	       fprintf( stderr, "PMatrixVector_FromFile: Out of memory.\n");
	       exit( EXIT_FAILURE );
	  }

	  InFile = fopen( Filename, "r" );
	  if( InFile == NULL ){
	       Print_Header( ERROR );
	       fprintf( stderr, "PMatrixVector_FromFile: Could not open %s.\n", Filename );
	       exit( EXIT_FAILURE );
	  }

	  for ( i = 0; i < Mat->GlobalSize.Row; i++ ){
	       for ( j = 0; j < Mat->GlobalSize.Col; j++ ){
		    fscanf( InFile, "%lE", &LocalMatrix[i + j*Mat->GlobalSize.Row] );
	       }
	  }
	  fclose( InFile );
     } else {
	  LocalMatrix = NULL;
     }
     
     pdgeadd_( &trans, &Mat->GlobalSize.Row, &Mat->GlobalSize.Col, &done, LocalMatrix, &ione, &ione,
	       descLocal, &dzero, Mat->Array, &ione, &ione, Mat->Desc );
     
     if ( myrow == 0 && mycol == 0 ){
	  Print_Header( SUCCESS );
	  printf( "PMatrixVector_FromFile: Contents of %s successfully readen.\n", Filename );     
	  free( LocalMatrix );
     }     
}

#if _MATRIXMARKET_
void PMatrixVector_FromFile_MM( const char *Filename, PMatrixVector_t *const Mat )
{

     int i, j, innz;				/* Counters */
     double Value;
     int myrow, mycol;		/* Grid variables */
     int nprow, npcol;
     
     MM_typecode matcode;   /* MatrixMarket: type of the matrix (symmetric, dense, complex, ...)  */
     int return_code;       /* MatrixMarket: return code for the functions */

     char trans = 'N';
     int izero = 0, ione = 1;
     double done = 1.0, dzero = 0.0;
     int lld, info;

     int Rows, Cols, nnz;
     
     double *LocalMatrix;
     int descLocal[9];
     
     FILE *InFile;
     
     
     Cblacs_gridinfo( Mat->Desc[1], &nprow, &npcol, &myrow, &mycol );
     
     lld = Max(1, Mat->GlobalSize.Row );
     
     descinit_( descLocal, &Mat->GlobalSize.Row, &Mat->GlobalSize.Col, &Mat->GlobalSize.Row,
		&Mat->GlobalSize.Col, &izero, &izero, &Mat->Desc[1], &lld, &info );
     
     if ( info < 0 ){
	  Print_Header( ERROR );
	  fprintf( stderr, "descinit: The %d-th argument had an illegal value.\n", -info );
	  exit( EXIT_FAILURE );
     }
     
     if ( myrow == 0 && mycol == 0 ){
	  LocalMatrix = calloc( (size_t) Mat->GlobalSize.Row * (size_t)Mat->GlobalSize.Col, sizeof(double) );
	  if( LocalMatrix == NULL ){
	       Print_Header( ERROR );
	       fprintf( stderr, "PMatrixVector_FromFile_MM: Out of memory.\n");
	       exit( EXIT_FAILURE );
	  }
	  
	  InFile = fopen( Filename, "r" );
	  if( InFile == NULL ){
	       Print_Header( ERROR );
	       fprintf( stderr, "PMatrixVector_FromFile_MM: Could not open %s.\n", Filename );
	  }
	  
	  
	  /* Read the banner and identify which type of matrix is in the file */
	  if( mm_read_banner( InFile, &matcode ) != 0 ){
	       Print_Header( ERROR );
	       fprintf( stderr, "PMatrixVector_FromFile_MM: Could not process Market Matrix banner in %s\n.",
			Filename );
	       exit( EXIT_FAILURE );
	  }
     
	  /* Only sparse matrices are accepted */
	  if ( !mm_is_sparse(matcode) ){
	       Print_Header( ERROR );
	       fprintf( stderr, "PMatrixVector_FromFile_MM: the matrix or vector should be of type sparse for this application to work.\n" );
	       Print_Header( ERROR );
	       fprintf( stderr, "Specified Matrix Market type: %s.\n", mm_typecode_to_str(matcode) );
	       exit( EXIT_FAILURE );
	  }
	  
	  /* Get the sizes */
	  if ( (return_code = mm_read_mtx_crd_size( InFile, &Rows, &Cols, &nnz)) !=0){
	       exit( EXIT_FAILURE );
	  }

	  /* Check if the dimensions of the matrices are the same */
	  if ( Rows != Mat->GlobalSize.Row || Cols != Mat->GlobalSize.Col ){
	       Print_Header( ERROR );
	       fprintf( stderr, "MatrixVector_From_File_MM: The sizes of the matrix or vector (%d,%d) ",
			Rows, Cols );
	       fprintf( stderr, "do not match with the specified ones in the configuration file (%d,%d)\n",
			Mat->GlobalSize.Row, Mat->GlobalSize.Col );
	       exit( EXIT_FAILURE );
	  }

	  /* Read the values. The MatrixMarket format imposes that the file should contain only the lower part
	   * of the matrix in 1-based index. Since C and FORTRAN use row-major and column-major ordering
	   * respectively, the matrices will be stored as upper part in the C so that when calling the FORTRAN
	   * routines from BLAS they access the lower part of the matrix without requiring transposing it.
	   */
	  for( innz = 0; innz < nnz; innz++ ){
	       fscanf( InFile, "%d %d %lE", &i, &j, &Value );
	       LocalMatrix[(j-1)*Mat->GlobalSize.Col + (i-1)] = Value;
	  }
     } else {
	  LocalMatrix = NULL;
     }
     
     pdgeadd_( &trans, &Mat->GlobalSize.Row, &Mat->GlobalSize.Col, &done, LocalMatrix, &ione, &ione, descLocal,
	       &dzero, Mat->Array, &ione, &ione, Mat->Desc );
     
     if ( myrow == 0 && mycol == 0 ){
	  Print_Header( SUCCESS );
	  printf( "PMatrixVector_FromFile_MM: Contents of %s successfully readen.\n", Filename );     
	  free( LocalMatrix );
     }
}
#endif /* _MATRIXMARKET_ */

void PMatrixVector_ToFile( PMatrixVector_t *const Mat, const char *Filename )
{

     int i, j;			/* Counters */
     int myrow, mycol;		/* Grid variables */
     int nprow, npcol;

     char trans = 'N';
     int izero = 0, ione = 1;
     double done = 1.0, dzero = 0.0;
     int lld, info;

     double *LocalMatrix;
     int descLocal[9];

     FILE *OutFile;

     Cblacs_gridinfo( Mat->Desc[1], &nprow, &npcol, &myrow, &mycol );

     if ( Mat->GlobalSize.Col > 1 ){
	  lld = Max(1, numroc_( &Mat->GlobalSize.Row, &Mat->GlobalSize.Col, &myrow, &izero, &nprow ) );
     } else {
	  lld = Mat->GlobalSize.Row;
     }
     descinit_( descLocal, &Mat->GlobalSize.Row, &Mat->GlobalSize.Col, &Mat->GlobalSize.Row, &Mat->GlobalSize.Col,
		&izero, &izero, &Mat->Desc[1], &lld, &info );
     if ( info < 0 ){
	  Print_Header( ERROR );
	  fprintf( stderr, "descinit: The %d-th argument had an illegal value.\n", -info );
	  exit( EXIT_FAILURE );
     }

     if ( myrow == 0 && mycol == 0 ){
	  LocalMatrix = calloc( (size_t) (Mat->GlobalSize.Row*Mat->GlobalSize.Col), sizeof(double) );
	  if( LocalMatrix == NULL ){
	       Print_Header( ERROR );
	       fprintf( stderr, "PMatrixVector_ToFile: Out of memory.\n");
	       exit( EXIT_FAILURE );
	  }
     } else {
	  LocalMatrix = NULL;
     }

     pdgeadd_( &trans, &Mat->GlobalSize.Row, &Mat->GlobalSize.Col, &done, Mat->Array, &ione, &ione, Mat->Desc,
	       &dzero, LocalMatrix, &ione, &ione, descLocal );

     if ( myrow == 0 && mycol == 0 ){

	  OutFile = fopen( Filename, "w" );
	  if ( OutFile == NULL ){
	       Print_Header( ERROR );
	       fprintf( stderr, "PMatrixVector_ToFile: Could not open %s.\n", Filename );
	       exit( EXIT_FAILURE );
	  }

	  for ( i = 0; i < Mat->GlobalSize.Row; i++ ){
	       for ( j = 0; j < Mat->GlobalSize.Col; j++ ){
		    fprintf( OutFile, "%lE\t", LocalMatrix[i + Mat->GlobalSize.Row*j] );  /* Print in row major order */
	       }
	       fprintf( OutFile, "\n" );
	  }
	  fclose( OutFile );		
     }

     if ( myrow == 0 && mycol == 0 ){
	  free( LocalMatrix );
	  Print_Header( SUCCESS );
	  printf( "PMatrixVector_ToFile: Matrix successfully saved to %s.\n", Filename );
     }

}
