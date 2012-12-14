#include <stdio.h>
#include <stdlib.h>

#include <petscmat.h>

#include "ErrorHandling.h"

void ReadDataEarthquake_AbsValues( PetscScalar *Velocity, PetscScalar *Displacement, PetscInt NumSteps, const char *Filename )
{

     PetscInt i;					/* A counter */
     PetscScalar unnecessary;		/* Variable to store unnecessary data */
     PetscScalar temp1, temp2, temp3;
     FILE *InFile;

     InFile = fopen( Filename, "r" );

     if ( InFile != NULL ){
	  for ( i = 0; i < NumSteps; i++ ){
	       fscanf( InFile, "%E %E %E %E", (float *) &unnecessary, (float *) &temp1, (float *) &temp2, (float *) &temp3 );
	       Velocity[i] = temp2/1000.0f;
	       Displacement[i] = temp3/1000.0f;
	  }

	  /* Close File */
	  fclose( InFile );
     } else {
	  ErrorFileAndExit( "The earthquake data cannot be read because it was not possible to open ", Filename );
     }
}

void ReadDataEarthquake_RelValues( PetscScalar *Acceleration, PetscInt NumSteps, const char *Filename )
{

     PetscInt i;					/* A counter */
     PetscScalar unnecessary;		/* Variable to store unnecessary data */
     PetscScalar temp1, temp2, temp3;
     FILE *InFile;

     InFile = fopen( Filename, "r" );


     if ( InFile != NULL ){
	  for ( i = 0; i < NumSteps; i++ ){
	       fscanf( InFile, "%E %E %E %E", (float *) &unnecessary, (float *) &temp1, (float *) &temp2, (float *) &temp3 );
	       Acceleration[i] = temp1/1000.0f;
	  }

	  /* Close File */
	  fclose( InFile );
     } else {
	  ErrorFileAndExit( "The earthquake data cannot be read because it was not possible to open ", Filename );
     }
}

void Calc_Input_Load_AbsValues( Vec InLoad, Mat Stif, const Mat Damp, Vec D, Vec V )
{

     /* {r} is the load form vector */
     /* li = K*{r}*ug */
     MatMult( Stif, D, InLoad );
     /* li = K*{r}*ug + C*{r}*vg = li + C*{r}*vg */
     MatMultAdd( Damp, V, InLoad, InLoad );
}

void Calc_Input_Load_RelValues( Vec InLoad, Mat Mass, Vec A )
{
     
     MatMult( Mass, A, InLoad );
     VecScale( InLoad, -1.0 );
}

void Apply_LoadVectorForm ( Vec Vector, Vec LoadForm, PetscScalar Value )
{
     VecCopy( LoadForm, Vector );
     VecScale( Vector, Value );
}
