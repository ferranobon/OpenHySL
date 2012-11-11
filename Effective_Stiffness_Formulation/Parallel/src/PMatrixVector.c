/*
 * PMatrixVector.c
 *
 *  Created on: 22/07/2011
 *      Author: ferran
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* For strcmp( ) */

#include "ErrorHandling.h"
#include "Netlib.h"
#include "PMatrixVector.h"

void CreateDistMatrix( int icntxt, PMatrixVector *const Mat, const int NumRows, const int NumCols, const int BlRows, int const BlCols )
{

	int myrow, mycol;  /* Variables to store row and column in the process grid */
	int nprow, npcol;

	int izero = 0;     /* Zero of type integer. Used by numroc_( ) */
	int lld, info;
	static int counter = 1;

	Mat->GlobalSize.Row = NumRows;   /* Number of rows in the global array */
	Mat->GlobalSize.Col = NumCols;   /* Number of columns in the global array */
	
	Mat->BlockSize.Row = BlRows;  /* Block size in the vertical dimension */
	Mat->BlockSize.Col = BlCols;  /* Block size in the horizontal dimension */

	Cblacs_gridinfo( icntxt, &nprow, &npcol, &myrow, &mycol );  // Get information about the grid

	/* Compute the size of the local matrices */
	Mat->LocalSize.Row = numroc_( &Mat->GlobalSize.Row, &Mat->BlockSize.Row, &myrow, &izero, &nprow );
	Mat->LocalSize.Col = numroc_( &Mat->GlobalSize.Col, &Mat->BlockSize.Col, &mycol, &izero, &npcol );
	lld = Max( 1, Mat->LocalSize.Row );
	
	descinit_( Mat->Desc, &Mat->GlobalSize.Row, &Mat->GlobalSize.Col, &Mat->BlockSize.Row, &Mat->BlockSize.Col, &izero, &izero, &icntxt, &lld, &info );
	if ( info < 0 ){
		LAPACKPErrorAndExit( "descinit: The ", -info, "-th argument had an illegal value" );
		
	}

	/* Allocate memory for the local array */
	Mat->Array = (float *) calloc( (size_t) Mat->LocalSize.Row*Mat->LocalSize.Col, sizeof(float) );
	counter = counter + 1;
}

void Set2Value( PMatrixVector *const Mat, const float Value )
{

	static int incx = 0;
	static int incy = 1;
	static int Length;
	static float Val;

	Length = Mat->LocalSize.Row*Mat->LocalSize.Col;
	Val = Value;

	scopy_( &Length, &Val, &incx, Mat->Array, &incy );
}
void ModifyElement( PMatrixVector *const Mat, int GRowIndex, int GColIndex, const float Value, const char *Operation )
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

	/* Given the global index of an element (GRowIndex, GColIndex) returns the local index of the element (LRowIndex, LColIndex) and the coordinates of the process (Row Process, ColProcess) */
	infog2l_( &GRowIndex, &GColIndex, Mat->Desc, &nprow, &npcol, &myrow, &mycol, &LRowIndex, &LColIndex, &RowProcess, &ColProcess );

	/* Modify the value in the local matrix. -1 is substracted because C starts the array indexes at 0, while FORTRAN starts them at 1 (infog2l is a FORTRAN routine) */
	if ( myrow == RowProcess && mycol == ColProcess ){

		if ( strcmp( Operation, OpSet ) == 0 ){
			Mat->Array[(LRowIndex - 1)*Mat->LocalSize.Col + (LColIndex -1)] = Value;
		}else if ( strcmp( Operation, OpAdd ) == 0 ){
			Mat->Array[(LRowIndex - 1)*Mat->LocalSize.Col + (LColIndex -1)] = Mat->Array[(LRowIndex - 1)*Mat->LocalSize.Col + (LColIndex -1)] + Value;
		} else if ( strcmp( Operation, OpMult ) == 0 ){
			Mat->Array[(LRowIndex - 1)*Mat->LocalSize.Col + (LColIndex -1)] = Mat->Array[(LRowIndex - 1)*Mat->LocalSize.Col + (LColIndex -1)]*Value;
		} else if ( strcmp( Operation, OpDiv ) == 0 ){
			Mat->Array[(LRowIndex - 1)*Mat->LocalSize.Col + (LColIndex -1)] = Mat->Array[(LRowIndex - 1)*Mat->LocalSize.Col + (LColIndex -1)]/Value;
		} else {
			fprintf( stderr, "%s: Operation not identified. Valid operations are:\n", Operation );
			fprintf( stderr, "\t1) %s\n", OpSet );
			fprintf( stderr, "\t1) %s\n", OpAdd );
			fprintf( stderr, "\t1) %s\n", OpMult );
			fprintf( stderr, "\t1) %s\n", OpDiv );
			exit( EXIT_FAILURE );
		}
	}

}

void PAdd3Mat( PMatrixVector *const MatY, PMatrixVector *const MatA, PMatrixVector *const MatB, PMatrixVector *const MatC, Scalars Const )
{

	char uplo, trans;
	int ione;
	float done;

	ione = 1;
	trans = 'N'; /* The operation will not use the transpose matrix */
	uplo = 'L';  /* The lower part of the matrix will be used; the upper part will strictly not be referenced */
	done = 1.0;

	/* ScaLAPACK: Perform Y = A (locally. There is no communication) */
	pslacpy_( &uplo, &(*MatY).GlobalSize.Row, &(*MatY).GlobalSize.Col, (*MatA).Array, &ione, &ione, (*MatA).Desc, (*MatY).Array, &ione, &ione, (*MatY).Desc );

	/* ScaLAPACK: Perform Y = beta*B + alpha*A = beta*B + alpha*Y */
	pstradd_( &uplo, &trans, &(*MatY).GlobalSize.Row, &(*MatY).GlobalSize.Col, &Const.Beta, (*MatB).Array, &ione, &ione, (*MatB).Desc, &Const.Alpha, (*MatY).Array, &ione, &ione, (*MatY).Desc );

	/* ScaLAPACK: Perform Y = gamma*C + beta*B = gamma*C + 1.0*Y */
	pstradd_( &uplo, &trans, &(*MatY).GlobalSize.Row, &(*MatY).GlobalSize.Col, &Const.Gamma, (*MatC).Array, &ione, &ione, (*MatC).Desc, &done, (*MatY).Array, &ione, &ione, (*MatY).Desc );
}


void DistMatrixFromFile( PMatrixVector *const Mat, const char *Filename )
{

     int i, j;				/* Counters */
     int myrow, mycol;		/* Grid variables */
     int nprow, npcol;
     
     char trans = 'N';
     int izero = 0, ione = 1;
     float done = 1.0, dzero = 0.0;
     int lld, info;
     
     float *LocalMatrix;
     int descLocal[9];
     
     FILE *InFile;
     
     
     Cblacs_gridinfo( Mat->Desc[1], &nprow, &npcol, &myrow, &mycol );
     
     lld = Mat->GlobalSize.Row; //Max(1, numroc_( &Mat->GlobalSize.Row, &Mat->GlobalSize.Col, &myrow, &izero, &nprow ) );
     
     descinit_( descLocal, &Mat->GlobalSize.Row, &Mat->GlobalSize.Col, &Mat->GlobalSize.Row, &Mat->GlobalSize.Col, &izero, &izero, &Mat->Desc[1], &lld, &info );
     
     if ( info < 0 ){
	  LAPACKPErrorAndExit( "descinit: The ", -info, "-th argument had an illegal value" );
     }
     
     if ( myrow == 0 && mycol == 0 ){
	  LocalMatrix = calloc( Mat->GlobalSize.Row*Mat->GlobalSize.Col, sizeof(float) );
	  InFile = fopen( Filename, "r" );
	  
	  if ( InFile != NULL ){
	       for ( i = 0; i < Mat->GlobalSize.Row; i++ ){
		    for ( j = 0; j < Mat->GlobalSize.Col; j++ ){
			 fscanf( InFile, "%e", &LocalMatrix[i + j*Mat->GlobalSize.Row] );
		    }
	       }
	       fclose( InFile );
	  } else {
	       ErrorFileAndExit( "It is not possible to read data because it was not possible to open: ", Filename );
	  }
     } else {
	  LocalMatrix = NULL;
     }
     
     psgeadd_( &trans, &Mat->GlobalSize.Row, &Mat->GlobalSize.Col, &done, LocalMatrix, &ione, &ione, descLocal, &dzero, Mat->Array, &ione, &ione, Mat->Desc );
     
     if ( myrow == 0 && mycol == 0 ){
	  free( LocalMatrix );
     }
     
}
/*
void DistVectorFromFile( PMatrixVector *const Mat, const char *Filename )
{
     int i, j;
*/


/**************************************************************
 *****                                                    *****
 *****                  WriteMatrixFile                   *****
 *****   Writes a distributed matrix into a binary file   *****
 *****                                                    *****
 *************************************************************/

void DistMatrixToFile( PMatrixVector *const Mat, const char *Filename )
{

	int i, j;				// Counters
	int myrow, mycol;		// Grid variables
	int nprow, npcol;

	char trans = 'N';
	int izero = 0, ione = 1;
	float done = 1.0, dzero = 0.0;
	int lld, info;

	float *LocalMatrix;
	int descLocal[9];

	FILE *OutFile;

	Cblacs_gridinfo( Mat->Desc[1], &nprow, &npcol, &myrow, &mycol );

//	if ( Mat->GlobalSize.Col > 1 ){
//		lld = Max(1, numroc_( &Mat->GlobalSize.Row, &Mat->GlobalSize.Col, &myrow, &izero, &nprow ) );
//	} else {
		lld = Mat->GlobalSize.Row;
//	}
	descinit_( descLocal, &Mat->GlobalSize.Row, &Mat->GlobalSize.Col, &Mat->GlobalSize.Row, &Mat->GlobalSize.Col, &izero, &izero, &Mat->Desc[1], &lld, &info );
	if ( info < 0 ){
		LAPACKPErrorAndExit( "descinit: The ", -info, "-th argument had an illegal value" );
	}

	if ( myrow == 0 && mycol == 0 ){
		LocalMatrix = calloc( Mat->GlobalSize.Row*Mat->GlobalSize.Col, sizeof(float) );
	} else {
		LocalMatrix = NULL;
	}

	psgeadd_( &trans, &Mat->GlobalSize.Row, &Mat->GlobalSize.Col, &done, Mat->Array, &ione, &ione, Mat->Desc, &dzero, LocalMatrix, &ione, &ione, descLocal );

	if ( myrow == 0 && mycol == 0 ){

		OutFile = fopen( Filename, "w" );
		if ( OutFile != NULL ){
			for ( i = 0; i < Mat->GlobalSize.Row; i++ ){
				for ( j = 0; j < Mat->GlobalSize.Col; j++ ){
					fprintf( OutFile, "%e\t", LocalMatrix[i + Mat->GlobalSize.Row*j] );  //*Print in row major order */
				}
				fprintf( OutFile, "\n" );
			}
			fclose( OutFile );
		} else {
			ErrorFileAndExit( "It is not possible to read data because it was not possible to open: ", Filename );
		}
	}

	if ( myrow == 0 && mycol == 0 ){
		free( LocalMatrix );
	}

}

void DestroyDistMatrix( PMatrixVector *const Mat )
{

	/* Set Global and local sizes to 0 */
	Mat->GlobalSize.Row = 0;
	Mat->GlobalSize.Col = 0;

	Mat->LocalSize.Row = 0;
	Mat->LocalSize.Col = 0;

	/* Deallocate memory */
	free( Mat->Array );
}

int Max ( const int a, const int b )
{
	if ( a > b ){
		return a;
	} else return b;
}
