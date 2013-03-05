#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Print_Messages.h"
#include "Substructure_SimMeasured.h"

void Substructure_SimMeasured_Init( const char *FileName, const unsigned int NSteps, const unsigned int NSubsteps, const char *Description, MeasuredSim_t *const Sub )
{
     unsigned int i;   /* A counter */
     FILE *InFile;

     Sub->Description = strdup( Description );

     Sub->Length = NSteps*NSubsteps;

     Sub->Values = (double *) calloc( (size_t) Sub->Length, sizeof(double) );
     if ( Sub->Values == NULL ){
	  Print_Header( ERROR );
	  fprintf(stderr, "Substructure_SimMeasured_Init: Out of memory.\n" );
	  exit( EXIT_FAILURE );
     }

     InFile = fopen( FileName, "r" );
     if ( InFile == NULL ){
	  Print_Header( ERROR );
	  fprintf(stderr, "Substructure_SimMeasured_Init: Could not open the file with measured forces '%s'.\n",
		  FileName );
	  exit( EXIT_FAILURE );
     }

     i = 0;
     while( i < Sub->Length && (fscanf( InFile, "%lE", &Sub->Values[i] ) != EOF) ){  
	  i = i + 1;
     }

     if( Sub->Length != i){
	  Print_Header( ERROR );
	  fprintf( stderr, "Substructure_SimMeasured_Init: End of file reached before the required number of input data in %s.", FileName );
	  fprintf( stderr, " The amount of data present in the file should be NSteps*NSubsteps: %d.\n", Sub->Length );
	  exit( EXIT_FAILURE );
     }
}


void Substructure_SimMeasured( const MeasuredSim_t *const Sub, double *const fc )
{
     static unsigned int i = 0;

     (*fc) = Sub->Values[i];

     i = i + 1;
}

void Substructure_SimMeasured_Destroy( MeasuredSim_t *const Sub )
{
     free( Sub->Description );
}
