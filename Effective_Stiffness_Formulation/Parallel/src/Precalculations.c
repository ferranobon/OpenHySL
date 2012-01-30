/*
 * Precalculations.c
 *
 *  Created on: 22/07/2011
 *      Author: ferran
 */


#include <stdio.h>
#include <stdlib.h>
#include <assert.h>		/* For assert( ) */
#include <mpi.h>

#include "ErrorHandling.h"
#include "Netlib.h"
#include "PMatrixVector.h"
#include "Precalculations.h"

void ReadDataEarthquake( float *displacement, float *velocity, float *acceleration, const int NumSteps, const char *Filename )
{

	int i;					/* A counter */
	float unnecessary;		/* Variable to store unnecessary data */
	float temp1, temp2, temp3;
	FILE *InFile;

	InFile = fopen( Filename, "r" );

	if ( InFile != NULL ){
		for ( i = 0; i < NumSteps; i++ ){
			/* The first column contains the number of steps (not required) */
			fscanf( InFile, "%f %f %f %f", &unnecessary, &temp1, &temp2, &temp3 );
			acceleration[i] = temp1/10.0;
			velocity[i] = temp2/10.0;
			displacement[i] = temp3/10.0;
		}
		/* Close File */
		fclose( InFile );
	} else {
		ErrorFileAndExit( "It is not possible to read data because it was not possible to open: ", Filename );
	}
}

void CopyDiagonalValues( MPI_Comm Comm, PMatrixVector *const Mat, PMatrixVector *const Vec )
{

	int nprow, npcol, myrow, mycol;       /* Variables required by BLACS routines Cblacs_grindinfo( ) and infog2l( ) */
	int nprowV, npcolV, myrowV, mycolV;

	int GIndex;                           /* Global Index */

	int LRowIndex, LColIndex;             /* Local Indexes */
	int LVecIndex, LVecIndexC;

	int RowProcess, ColProcess;
	int VRowProcess, VColProcess;

	MPI_Status status;

	int ione = 1;

	Cblacs_gridinfo( (*Mat).Desc[1], &nprow, &npcol, &myrow, &mycol );
	Cblacs_gridinfo( (*Vec).Desc[1], &nprowV, &npcolV, &myrowV, &mycolV );

	for ( GIndex = 1; GIndex <= (*Mat).GlobalSize.Row; GIndex++ ){

		infog2l_( &GIndex, &GIndex, (*Mat).Desc, &nprow, &npcol, &myrow, &mycol, &LRowIndex, &LColIndex, &RowProcess, &ColProcess );
		infog2l_( &GIndex, &ione, (*Vec).Desc, &nprowV, &npcolV, &myrowV, &mycolV, &LVecIndex, &LVecIndexC, &VRowProcess, &VColProcess );

		if ( Cblacs_pnum( (*Vec).Desc[1], VRowProcess, VColProcess) == Cblacs_pnum( (*Mat).Desc[1], RowProcess, ColProcess) && myrow == RowProcess && mycol==ColProcess && myrowV==VRowProcess && mycolV == VColProcess ){
			(*Vec).Array[(LVecIndex - 1) + (LVecIndexC - 1)*(*Vec).LocalSize.Row] = (*Mat).Array[(LRowIndex - 1) + (LColIndex - 1)*(*Mat).LocalSize.Row];
		} else {
			if ( myrow == RowProcess && mycol == ColProcess ){
				MPI_Send( &(*Mat).Array[(LRowIndex - 1) + (LColIndex - 1)*(*Mat).LocalSize.Row], 1, MPI_FLOAT, Cblacs_pnum( (*Vec).Desc[1], VRowProcess, VColProcess ), 0, Comm );
			}

			if ( myrowV == VRowProcess && mycolV == VColProcess ){
				MPI_Recv( &(*Vec).Array[(LVecIndex - 1) + (LVecIndexC - 1)*(*Vec).LocalSize.Row], 1, MPI_FLOAT, Cblacs_pnum( (*Mat).Desc[1], RowProcess, ColProcess ), 0, Comm, &status );
			}
		}
	}

}

void Calc_Input_Load( PMatrixVector *const InLoad, PMatrixVector *const Stif, PMatrixVector *const Damp, PMatrixVector *const Mass, PMatrixVector *const DiagM, PMatrixVector *const D, PMatrixVector *const V, PMatrixVector *const A )
{

  static int incx, incy;       /* Stride in the vectors for PBLAS library */
  static int ione;             /* Integer variable of value 1 for PBLAS library */
  static float Alpha, Beta;    /* Constants to use in the PBLAS library */
  static char uplo;            /* Character to use in the PBLAS library */

  incx = 1; incy = 1;
  Alpha = 1.0; Beta = 0.0;
  uplo = 'L';     /* Character defining that the lower part of the symmetric matrix is referenced (see man ssymv) */

  pssymv_( &uplo, &(*InLoad).GlobalSize.Row, &Alpha, (*Stif).Array, &ione, &ione, (*Stif).Desc,
	   (*D).Array, &ione, &ione, (*D).Desc, &incx, 
	   &Beta, (*InLoad).Array, &ione, &ione, (*InLoad).Desc, &incy );
  Beta = 1.0;
  pssymv_( &uplo, &(*InLoad).GlobalSize.Row, &Alpha, (*Damp).Array, &ione, &ione, (*Damp).Desc,
	   (*V).Array, &ione, &ione, (*V).Desc, &incx, 
	   &Beta, (*InLoad).Array, &ione, &ione, (*InLoad).Desc, &incy );
  pssymv_( &uplo, &(*InLoad).GlobalSize.Row, &Alpha, (*Mass).Array, &ione, &ione, (*Mass).Desc,
	   (*A).Array, &ione, &ione, (*A).Desc, &incx, 
	   &Beta, (*InLoad).Array, &ione, &ione, (*InLoad).Desc, &incy );

  Alpha = -(*A).Array[0];
  psaxpy_( &(*InLoad).GlobalSize.Row, &Alpha, (*DiagM).Array, &ione, &ione, (*DiagM).Desc, &incx, (*InLoad).Array, &ione, &ione, (*InLoad).Desc, &incy );

}
