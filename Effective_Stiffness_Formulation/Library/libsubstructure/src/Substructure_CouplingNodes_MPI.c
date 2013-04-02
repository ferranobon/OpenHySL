#include <stdlib.h>

#include "mpi.h"

#include "Print_Messages.h"
#include "Substructure.h"
#include "Substructure_CouplingNodes.h"


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
