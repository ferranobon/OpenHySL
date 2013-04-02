#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "MatrixVector.h"
#include "MatrixVector_MPI.h"
#include "Print_Messages.h"
#include "Auxiliary_Math.h"

#if _MKL_
#include "mkl_blas.h"
#include "mkl_pblas.h"
#include "mkl_blacs.h"
#include "mkl_scalapack.h"
#else
#include "Netlib.h"
#endif

void PMatrixVector_Create( int icntxt, const int NumRows, const int NumCols, const int BlRows, int const BlCols,
			   PMatrixVector_t *const Mat )
{

	int myrow, mycol;  /* Variables to store row and column in the process grid */
	int nprow, npcol;

	int izero = 0;     /* Zero of type integer. Used by numroc_( ) */
	int lld, info;

	Mat->GlobalSize.Row = NumRows;   /* Number of rows in the global array */
	Mat->GlobalSize.Col = NumCols;   /* Number of columns in the global array */
	
	Mat->BlockSize.Row = BlRows;  /* Block size in the vertical dimension */
	Mat->BlockSize.Col = BlCols;  /* Block size in the horizontal dimension */

	Cblacs_gridinfo( icntxt, &nprow, &npcol, &myrow, &mycol );  /* Get information about the grid */

	/* Compute the size of the local matrices */
	Mat->LocalSize.Row = numroc_( &Mat->GlobalSize.Row, &Mat->BlockSize.Row, &myrow, &izero, &nprow );
	Mat->LocalSize.Col = numroc_( &Mat->GlobalSize.Col, &Mat->BlockSize.Col, &mycol, &izero, &npcol );
	lld = Max( 1, Mat->LocalSize.Row );
	
	descinit_( Mat->Desc, &Mat->GlobalSize.Row, &Mat->GlobalSize.Col, &Mat->BlockSize.Row,
		   &Mat->BlockSize.Col, &izero, &izero, &icntxt, &lld, &info );
	if ( info < 0 ){
	     Print_Header( ERROR );
	     fprintf( stderr, "descinit: The %d-th argument had an illegal value.\n", -info );
	     exit( EXIT_FAILURE );
	}
	
	/* Allocate memory for the local array */
	Mat->Array = (double *) calloc( (size_t) Mat->LocalSize.Row*(size_t) Mat->LocalSize.Col,
					sizeof(double) );
	if ( Mat->Array == NULL ){
	     Print_Header( ERROR );
	     fprintf( stderr, "PMatrixVector_Create: Out of memory.\n" );
	     exit( EXIT_FAILURE );
	}
}

void PMatrixVector_Set2Value( const double Value, PMatrixVector_t *const Mat )
{

     int incx, incy;
     int Length;
     double Val;

     incx = 0; incy = 1;
     
     Length = Mat->LocalSize.Row*Mat->LocalSize.Col;
     Val = Value;
     
     dcopy( &Length, &Val, &incx, Mat->Array, &incy );
}

void PMatrixVector_ModifyElement( const int GRowIndex, const int GColIndex, const double Value,
				  const char *Operation, PMatrixVector_t *const Mat )
{

	int myrow, mycol, nprow, npcol;

	const char *OpSet = "Set";
	const char *OpAdd = "Add";
	const char *OpMult = "Multiply";
	const char *OpDiv = "Divide";

	int LRowIndex, LColIndex;
	int RowProcess, ColProcess;

	/* Get grid info */
	Cblacs_gridinfo( Mat->Desc[1], &nprow, &npcol, &myrow, &mycol );

	/* Given the global index of an element (GRowIndex, GColIndex) returns the local index of the element
	 * (LRowIndex, LColIndex) and the coordinates of the process (Row Process, ColProcess) */
	infog2l_( &GRowIndex, &GColIndex, Mat->Desc, &nprow, &npcol, &myrow, &mycol, &LRowIndex, &LColIndex,
		  &RowProcess, &ColProcess );

	/* Modify the value in the local matrix. -1 is substracted because C starts the array indexes at 0,
	 * while FORTRAN starts them at 1 (infog2l is a FORTRAN routine) */
	if ( myrow == RowProcess && mycol == ColProcess ){

		if ( strcmp( Operation, OpSet ) == 0 ){
			Mat->Array[(LRowIndex - 1)*Mat->LocalSize.Col + (LColIndex -1)] = Value;
		}else if ( strcmp( Operation, OpAdd ) == 0 ){
			Mat->Array[(LRowIndex - 1)*Mat->LocalSize.Col + (LColIndex -1)] =
			     Mat->Array[(LRowIndex - 1)*Mat->LocalSize.Col + (LColIndex -1)] + Value;
		} else if ( strcmp( Operation, OpMult ) == 0 ){
			Mat->Array[(LRowIndex - 1)*Mat->LocalSize.Col + (LColIndex -1)] =
			     Mat->Array[(LRowIndex - 1)*Mat->LocalSize.Col + (LColIndex -1)]*Value;
		} else if ( strcmp( Operation, OpDiv ) == 0 ){
			Mat->Array[(LRowIndex - 1)*Mat->LocalSize.Col + (LColIndex -1)] =
			     Mat->Array[(LRowIndex - 1)*Mat->LocalSize.Col + (LColIndex -1)]/Value;
		} else {
		     Print_Header( ERROR );
		     fprintf( stderr, "PMatrixVector_ModifyElement: Operation '%s' not identified. Valid operations are:\n" , Operation );
		     fprintf( stderr, "[......] 1) %s.\n", OpSet );
		     fprintf( stderr, "[......] 2) %s.\n", OpAdd );
		     fprintf( stderr, "[......] 3) %s.\n", OpMult );
		     fprintf( stderr, "[......] 4) %s.\n", OpDiv );
		     exit( EXIT_FAILURE );
		}
	}

}

void PMatrixVector_Add3Mat( PMatrixVector_t *const MatA, PMatrixVector_t *const MatB,
			    PMatrixVector_t *const MatC, const Scalars_t Const, PMatrixVector_t *const MatY )
{

	char uplo, trans;
	int ione;
	double ScalarA, ScalarB;

	ione = 1;
	trans = 'N'; /* The operation will not use the transpose matrix */
	uplo = 'L';  /* The lower part of the matrix will be used; the upper part will strictly not be
		      * referenced */
	/* ScaLAPACK: Perform Y = A (locally. There is no communication) */
	pdlacpy( &uplo, &MatY->GlobalSize.Row, &MatY->GlobalSize.Col, MatA->Array, &ione, &ione,
		 MatA->Desc, MatY->Array, &ione, &ione, MatY->Desc );

	/* ScaLAPACK: Perform Y = beta*B + alpha*A = beta*B + alpha*Y */
	ScalarA = Const.Alpha; ScalarB = Const.Beta;
	pdtradd_( &uplo, &trans, &MatY->GlobalSize.Row, &MatY->GlobalSize.Col, &ScalarB, MatB->Array,
		  &ione, &ione, MatB->Desc, &ScalarA, MatY->Array, &ione, &ione, MatY->Desc );

	/* ScaLAPACK: Perform Y = gamma*C + beta*B = gamma*C + 1.0*Y */
	ScalarA = 1.0; ScalarB = Const.Gamma;
	pdtradd_( &uplo, &trans, &MatY->GlobalSize.Row, &MatY->GlobalSize.Col, &ScalarB, MatC->Array,
		  &ione, &ione, MatC->Desc, &ScalarA, MatY->Array, &ione, &ione, MatY->Desc );
}

void PMatrixVector_Destroy( PMatrixVector_t *const Mat )
{

	/* Set Global and local sizes to 0 */
	Mat->GlobalSize.Row = 0;
	Mat->GlobalSize.Col = 0;

	Mat->LocalSize.Row = 0;
	Mat->LocalSize.Col = 0;

	/* Deallocate memory */
	free( Mat->Array );
}
