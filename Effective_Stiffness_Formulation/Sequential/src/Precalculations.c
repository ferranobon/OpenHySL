/**
 * \file Precalculations.c
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 10th of August 2011
 *
 * \brief Source code of the functions used in the precalculations phase.
 *
 * This file contains the source code of the functions that are called during the precalculations phase. This includes,
 * reading the earthquake data from a file, and perform some operations like Copying the diagonal values of the Mass matrix
 * and calculate the input load.
 */

#include <stdio.h>
#include <stdlib.h>

#include "ErrorHandling.h"
#include "MatrixVector.h"
#include "Netlib.h"  /* BLAS and LAPACK prototypes. */
#include "Precalculations.h"

void ReadDataEarthquake( float *Acceleration, float *Velocity, float *Displacement, const int NumSteps, const char *Filename )
{

  int i;					/* A counter */
  float unnecessary;		/* Variable to store unnecessary data */
  float temp1, temp2, temp3;
  FILE *InFile;

  InFile = fopen( Filename, "r" );


  if ( InFile != NULL ){
      for ( i = 0; i < NumSteps; i++ ){
    	  /* EFAST The first column contains the number of steps (not required) */
	   //fscanf( InFile, "%f%f%f%f", &unnecessary, &temp1, &temp2, &temp3 );
	   //fscanf( InFile, "%f%f%f%f", &unnecessary, &temp1, &temp2, &temp3 );
	   //fscanf( InFile, "%f%f%f%f", &unnecessary, &temp1, &temp2, &temp3 );

	   fscanf( InFile, "%e %e %e %e", &unnecessary, &temp1, &temp2, &temp3 );
	   Acceleration[i] = temp1/1000.0;
	   Velocity[i] = temp2/1000.0;
	   Displacement[i] = temp3/1000.0;
      }

      /* Close File */
      fclose( InFile );
  } else {
       ErrorFileAndExit( "The earthquake data cannot be read because it was not possible to open ", Filename );
  }
}

/* Stores in the variable Vec, the diagonal values of the Matrix Mat. */
void CopyDiagonalValues( const MatrixVector *const Mat, MatrixVector *const Vec )
{

  int incx, incy;  /* Stride in the vectors for the BLAS library */

  incx = (*Mat).Rows + 1;
  incy = 1;

  /* BLAS routine */
  scopy_( &(*Vec).Rows, (*Mat).Array, &incx, (*Vec).Array, &incy );

}

void Calc_Input_Load( MatrixVector *const InLoad, const MatrixVector *const Stif, const MatrixVector *const Damp, const MatrixVector *const Mass, const MatrixVector *const DiagM, const MatrixVector *const D, const MatrixVector *const V, const MatrixVector *const A )
{

  static int incx, incy;       /* Stride in the vectors for BLAS library */
  static float Alpha, Beta;    /* Constants to use in the BLAS library */
  static char uplo;            /* Character to use in the BLAS library */

  incx = 1; incy = 1;
  Alpha = 1.0; Beta = 0.0;
  uplo = 'L';     /* Character defining that the lower part of the symmetric matrix is referenced (see man dsymv) */

  ssymv_( &uplo, &(*InLoad).Rows, &Alpha, (*Stif).Array, &(*InLoad).Rows, (*D).Array, &incx, &Beta, (*InLoad).Array, &incy );
  Beta = 1.0;
  ssymv_( &uplo, &(*InLoad).Rows, &Alpha, (*Damp).Array, &(*InLoad).Rows, (*V).Array, &incx, &Beta, (*InLoad).Array, &incy );
  ssymv_( &uplo, &(*InLoad).Rows, &Alpha, (*Mass).Array, &(*InLoad).Rows, (*A).Array, &incx, &Beta, (*InLoad).Array, &incy );

  Alpha = -(*A).Array[0];
  saxpy_( &(*InLoad).Rows, &Alpha, (*DiagM).Array, &incx, (*InLoad).Array, &incy );

}
