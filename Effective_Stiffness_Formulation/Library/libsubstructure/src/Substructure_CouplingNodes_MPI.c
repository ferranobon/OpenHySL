#include <stdlib.h>

#include "mpi.h"

#include "Print_Messages.h"
#include "Substructure.h"
#include "Substructure_CouplingNodes.h"

#include "Cblacs.h"


void Substructure_BroadCastCouplingNodes( CouplingNode_t *const CNodes )
{

     int rank;

     MPI_Comm_rank( MPI_COMM_WORLD, &rank );

     MPI_Bcast( &CNodes->Order, 1, MPI_INT, 0, MPI_COMM_WORLD );

     if( rank != 0 ){
	  CNodes->Array = (int *) calloc( (size_t) CNodes->Order, sizeof(int) );
     }

     MPI_Bcast( CNodes->Array, CNodes->Order, MPI_INT, 0, MPI_COMM_WORLD );
}

void Substructure_InfoLocation_Init( const int Length, int *Position, int *const Descriptor,
				     InfoLocation_t *const InfoLoc )
{

     int ione = 1, i;

     int nprow, npcol, myrow, mycol;

     InfoLoc->LRowIndex = (int *) calloc( (size_t) Length, sizeof(int) );
     InfoLoc->LColIndex = (int *) calloc( (size_t) Length, sizeof(int) );
     InfoLoc->RowProcess = (int *) calloc( (size_t) Length, sizeof(int) );
     InfoLoc->ColProcess = (int *) calloc( (size_t) Length, sizeof(int) );

     /* Get grid information */
     Cblacs_gridinfo( Descriptor[1], &nprow, &npcol, &myrow, &mycol );

     /* Get the local row, column indeces and the coordinates in the process grid */
     for( i = 0; i < Length; i++ ){
	  infog2l_( &Position[i], &ione, Descriptor, &nprow, &npcol, &myrow, &mycol, &InfoLoc->LRowIndex[i],
		    &InfoLoc->LColIndex[i], &InfoLoc->RowProcess[i], &InfoLoc->ColProcess[i] );
     }
}

void Substructure_InfoLocation_Destroy( InfoLocation_t *const InfoLoc )
{
     free( InfoLoc->LRowIndex );
     free( InfoLoc->LColIndex );
     free( InfoLoc->RowProcess );
     free( InfoLoc->ColProcess );
}
