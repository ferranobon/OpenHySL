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

void ReadDataEarthquake_AbsValues( float *Velocity, float *Displacement, const unsigned int NumSteps, const char *Filename )
{

     unsigned int i;					/* A counter */
     float unnecessary;		/* Variable to store unnecessary data */
     float temp1, temp2, temp3;
     FILE *InFile;

     InFile = fopen( Filename, "r" );

     if ( InFile != NULL ){
	  for ( i = 0; i < NumSteps; i++ ){
	       fscanf( InFile, "%E %E %E %E", &unnecessary, &temp1, &temp2, &temp3 );
	       Velocity[i] = temp2/1000.0f;
	       Displacement[i] = temp3/1000.0f;
	  }

	  /* Close File */
	  fclose( InFile );
     } else {
	  ErrorFileAndExit( "The earthquake data cannot be read because it was not possible to open ", Filename );
     }
}

void ReadDataEarthquake_RelValues( float *Acceleration, const unsigned int NumSteps, const char *Filename )
{

     unsigned int i;					/* A counter */
     float unnecessary;		/* Variable to store unnecessary data */
     float temp1, temp2, temp3;
     FILE *InFile;

     InFile = fopen( Filename, "r" );


     if ( InFile != NULL ){
	  for ( i = 0; i < NumSteps; i++ ){
	       fscanf( InFile, "%E %E %E %E", &unnecessary, &temp1, &temp2, &temp3 );
	       Acceleration[i] = temp1/1000.0f;
	  }

	  /* Close File */
	  fclose( InFile );
     } else {
	  ErrorFileAndExit( "The earthquake data cannot be read because it was not possible to open ", Filename );
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
void Calc_Input_Load_AbsValues( PMatrixVector *const InLoad, const PMatrixVector *const Stif, const PMatrixVector *const Damp, const PMatrixVector *const D, const PMatrixVector *const V )
{

     static int incx, incy;     /* Stride in the vectors for PBLAS library */
     static int ione;           /* Integer variable of value 1 for PBLAS library */
     static float Alpha, Beta;  /* Constants to use in the PBLAS library */
     static char uplo;          /* Character to use in the PBLAS library */

     incx = 1; incy = 1;
     Alpha = 1.0f; Beta = 0.0f;
     uplo = 'L';
     ione = 1;

     /* {r} is the load form vector */
     /* li = K*{r}*ug */
     pssymv_( &uplo, &InLoad->GlobalSize.Row, &Alpha, Stif->Array, &ione, &ione,
	      Stif->Desc, D->Array, &ione, &ione, D->Desc, &incx, 
	      &Beta, InLoad->Array, &ione, &ione, InLoad->Desc, &incy );

     /* li = K*{r}*ug + C*{r}*vg = li + C*{r}*vg */
     Beta = 1.0f;
     pssymv_( &uplo, &InLoad->GlobalSize.Row, &Alpha, Damp->Array, &ione, &ione,
	      Damp->Desc, V->Array, &ione, &ione, V->Desc, &incx, 
	      &Beta, InLoad->Array, &ione, &ione, InLoad->Desc, &incy );
}

void Calc_Input_Load_RelValues( PMatrixVector *const InLoad, const PMatrixVector *const Mass, const PMatrixVector *const A )
{
     static int incx, incy;     /* Stride in the vectors for PBLAS library */
     static int ione;           /* Integer variable of value 1 for PBLAS library */
     static float Alpha, Beta;  /* Constants to use in the PBLAS library */
     static char uplo;          /* Character to use in the PBLAS library */

     incx = 1; incy = 1;
     Alpha = 1.0f; Beta = 0.0f;
     uplo = 'L';

     pssymv_( &uplo, &InLoad->GlobalSize.Row, &Alpha, Mass->Array, &ione, &ione,
	      Mass->Desc, A->Array, &ione, &ione, A->Desc, &incx, 
	      &Beta, InLoad->Array, &ione, &ione, InLoad->Desc, &incy );
}

void Apply_LoadVectorForm( PMatrixVector *const Vector, const PMatrixVector *const LoadForm, const float Value )
{
     static int incx = 1;
     static int incy = 1;
     static int ione = 1;
     static float Scalar;

     Scalar = Value;


     pscopy_( &Vector->GlobalSize.Row, LoadForm->Array, &ione, &ione, LoadForm->Desc,
	      &ione, Vector->Array, &ione, &ione, Vector->Desc, &ione );
     psscal_( &Vector->GlobalSize.Row, &Scalar, Vector->Array, &ione, &ione,
	      Vector->Desc, &incx );
}
