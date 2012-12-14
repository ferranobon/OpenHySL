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

#if _SPARSE_
#include <mkl_spblas.h>
#endif

#include "ErrorHandling.h"
#include "MatrixVector.h"
#include "Netlib.h"  /* BLAS and LAPACK prototypes. */
#include "Precalculations.h"

void ReadDataEarthquake_AbsValues( float *Acceleration, float *Velocity, float *Displacement, const unsigned int NumSteps, const char *Filename )
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

/* Stores in the variable Vec, the diagonal values of the Matrix Mat. */
void CopyDiagonalValues( const MatrixVector *const Mat, MatrixVector *const Vec )
{

     int incx, incy;  /* Stride in the vectors for the BLAS library */

     incx = (*Mat).Rows + 1;
     incy = 1;

     /* BLAS routine */
     scopy_( &(*Vec).Rows, (*Mat).Array, &incx, (*Vec).Array, &incy );

}

void Calc_Input_Load_AbsValues( MatrixVector *const InLoad, const MatrixVector *const Stif, const MatrixVector *const Damp, const MatrixVector *const D, const MatrixVector *const V )
{

     static int incx, incy;       /* Stride in the vectors for BLAS library */
     static float Alpha, Beta;    /* Constants to use in the BLAS library */
     static char uplo;            /* Character to use in the BLAS library */

     incx = 1; incy = 1;
     Alpha = 1.0f; Beta = 0.0f;
     uplo = 'L';     /* Character defining that the lower part of the symmetric matrix is referenced (see man dsymv) */

     /* {r} is the load form vector */
     /* li = K*{r}*ug */
     ssymv_( &uplo, &(*InLoad).Rows, &Alpha, (*Stif).Array, &(*InLoad).Rows, (*D).Array, &incx, &Beta, (*InLoad).Array, &incy );
     MatrixVector_To_File( InLoad, "InLoad" );
     /* li = K*{r}*ug + C*{r}*vg = li + C*{r}*vg */
     Beta = 1.0;
     ssymv_( &uplo, &(*InLoad).Rows, &Alpha, (*Damp).Array, &(*InLoad).Rows, (*V).Array, &incx, &Beta, (*InLoad).Array, &incy );
}

void Calc_Input_Load_RelValues( MatrixVector *const InLoad, const MatrixVector *const Mass, const MatrixVector *const A )
{

     static int incx = 1, incy = 1;         /* Stride in the vectors for BLAS library */
     static float Alpha = -1.0, Beta = 0.0;  /* Constants to use in the BLAS library */
     static char uplo = 'L';                /* Character defining that the lower part of the symmetric matrix is referenced (see man ssymv) */

     ssymv_( &uplo, &(*InLoad).Rows, &Alpha, (*Mass).Array, &(*InLoad).Rows, (*A).Array, &incx, &Beta, (*InLoad).Array, &incy );
}

void Apply_LoadVectorForm ( MatrixVector *const Vector, const MatrixVector *const LoadForm, const float Value )
{
     static int incx = 1;
     static int incy = 1;
     static float Scalar;

     Scalar = Value;
     scopy_( &Vector->Rows, LoadForm->Array, &incx, Vector->Array, &incy );
     sscal_( &Vector->Rows, &Scalar, Vector->Array, &incx );
}

#if _SPARSE_
void Calc_Input_Load_AbsValues_Sparse( MatrixVector *const InLoad, const Sp_MatrixVector *const Stif, const Sp_MatrixVector *const Damp, const MatrixVector *const D, const MatrixVector *const V )
{

     static float Alpha, Beta;    /* Constants to use in the BLAS library */
     static char trans = 'N';
     static char matdescra[6] = {'S', 'U', 'N', 'F'};

     Alpha = 1.0f; Beta = 0.0f;

     mkl_scsrmv( &trans, &InLoad->Rows, &InLoad->Rows, &Alpha, matdescra, Stif->Values, Stif->Columns, Stif->RowIndex, &(Stif->RowIndex[1]), D->Array, &Beta, InLoad->Array );
     MatrixVector_To_File( InLoad, "InLoad" );
     Beta = 1.0;
     mkl_scsrmv( &trans, &InLoad->Rows, &InLoad->Rows, &Alpha, matdescra, Damp->Values, Damp->Columns, Damp->RowIndex, &(Damp->RowIndex[1]), V->Array, &Beta, InLoad->Array );
}

void Calc_Input_Load_RelValues_Sparse( MatrixVector *const InLoad, const Sp_MatrixVector *const Mass, const MatrixVector *const A )
{

     static float Alpha, Beta;    /* Constants to use in the BLAS library */
     static char trans = 'N';
     static char matdescra[6] = {'S', 'U', 'N', 'F'};

     Alpha = 1.0; Beta = 0.0;

     mkl_scsrmv( &trans, &InLoad->Rows, &InLoad->Rows, &Alpha, matdescra, Mass->Values, Mass->Columns, Mass->RowIndex, &Mass->RowIndex[1], A->Array, &Beta, InLoad->Array );
}

#endif
