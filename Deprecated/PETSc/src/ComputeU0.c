#include <stdio.h>
#include <stdlib.h>

#include <petscmat.h>

#include "ComputeU0.h"
#include "Initiation.h"

void EffK_Calc_Effective_Force( Mat Mass, Mat Damp, Vec Disp, Vec Vel, Vec Acc, Vec Tempvec, PetscScalar a0, PetscScalar a1, PetscScalar a2, PetscScalar a3, PetscScalar a4, PetscScalar a5, Vec Eff_Force )
{
     /* tempvec = Disp */
     VecCopy( Disp, Tempvec );
     /* tempvec = a0*Disp + a2*Vel = a0*tempvec + a2*Vel */
     VecAXPBY( Tempvec, a2, a0, Vel );
     /* tempvec = a0*Disp + a2*Vel + a3*Acc = tempvec + a3*Acc */
     VecAXPY( Tempvec, a3, Acc );
     /* Eff_Force = Mass*(a0*Disp + a2*Vel + a3*Acc) = Mass*tempvec */
     MatMult( Mass, Tempvec, Eff_Force );

     /* tempvec = Disp */
     VecCopy( Disp, Tempvec );
     /* tempvec = a1*Disp + a4*Vel = a0*tempvec + a4*Vel */
     VecAXPBY( Tempvec, a4, a1, Vel );
     /* tempvec = a1*Disp + a4*Vel + a5*Acc = tempvec + a5*Acc */
     VecAXPY( Tempvec, a5, Acc );
     /* Eff_Force = Mass*(a0*Disp + a2*Vel + a3*Acc) + Damp*(a1*Disp + a4*Vel + a5*Acc) = Eff_Force + Damp*tempvec */
     MatMultAdd( Damp, Tempvec, Eff_Force, Eff_Force );
}


void EffK_ComputeU0( Vec Eff_Force, Vec In_Load, Vec Err_Force, PetscScalar PID_P, Mat Keinv,
		     Vec Tempvec, Vec Disp0 )
{
     /* Tempvec = Eff_Forcce + LoadTdT */
     VecWAXPY( Tempvec, 1.0, Eff_Force, In_Load );
     /* Tempvec =  Eff_Forcce + LoadTdT + Err_Force = tempvec + Err_Force. The sign of Err_Force was already applied when calculating it. */
     VecAXPY( Tempvec, PID_P, Err_Force );
     /* Disp0 = Keinv*(Eff_Force + LoadTdT + Err_Force) = Keinv*Tempvec */
     MatMult( Keinv, Tempvec, Disp0 );  
}

void CreateVectorXm( MPI_Comm Comm, Vec VecX, Vec VecXm, Coupling_Node *const CNodes )
{

     PetscInt i, irow, *ix;
     PetscInt Start_Position, Position;
     PetscInt Iam, rstart, rend, Rows, RowsXm;
     PetscScalar *Values;
     PetscMPIInt rank;

     MPI_Comm_rank( Comm, &rank );
     VecGetSize( VecX, &Rows );
     VecGetSize( VecXm, &RowsXm );
     PetscMalloc( RowsXm*sizeof(PetscInt), &ix );

     Start_Position = 0;
     irow = 0;
     for( i = 0; i < CNodes->Order; i++ ){
	  Position = CNodes->Array[i] - 1;
	  while ( Start_Position < Position ){
	       ix[irow] = Start_Position;
	       Start_Position = Start_Position + 1;
	       irow = irow + 1;
	  }
	  Start_Position = Start_Position + 1;
     }

     while( Start_Position < Rows ){
	  ix[irow] = Start_Position;
	  irow = irow + 1;
	  Start_Position = Start_Position + 1;
     }
     Values = GetVec_Values( Comm, VecX, RowsXm, ix );
     if( rank == 0 ){
	  for( i = 0; i < 503; i++ ){
	       //      printf("%f \t", Values[i] );
	  }
     }

     /* Transform to Xm coordinates */
     for( i = 0; i < RowsXm; i++ ){
	  ix[i] = i;
     }

     VecSetValues( VecXm, RowsXm, ix, Values, INSERT_VALUES );

     VecAssemblyBegin( VecXm );
     VecAssemblyEnd( VecXm );

     PetscFree( Values );
     PetscFree( ix );

}

void CreateVectorXc( MPI_Comm Comm, Vec VecX, PetscScalar *VecXc, Coupling_Node *const CNodes )
{

     PetscInt i;
     PetscInt Position;
     PetscInt Iam, rstart, rend, GRow;
     PetscScalar *array, Value;

     MPI_Status status;
     PetscMPIInt rank;

     MPI_Comm_rank( Comm, &rank );
     VecGetOwnershipRange( VecX, &rstart, &rend );
     PetscMalloc( 200*sizeof(PetscScalar), &array );
     VecGetArray( VecX, &array );

     for ( i = 0; i < CNodes->Order; i++ ){
	  Position = CNodes->Array[i] - 1;
	  Iam = GetOwner_Position_Vector( Comm, VecX, Position );

	  printf( "Iam %d rank %d, rstart %d, rend %d\n", Iam, rank, rstart, rend );
	  if( rstart <= Position && Position < rend ){
	       Value = array[Position - rstart];
	       for( i = 0; i < 31; i++ ){
		    printf("%e\t", array[i] );
	       }
	       printf("\nValue %e\n", Value );
	  }
	  if ( Iam == rank && rank == 0 ){
	       VecXc[i] = Value;
	  } else {
	       if ( rank == Iam ){
		    MPI_Send( &Value, 1, MPIU_SCALAR, 0, 1, Comm );
	       }

	       if ( rank == 0 ){
		    MPI_Recv( &VecXc[i], 1, MPIU_SCALAR, Iam, 1, Comm, &status );
	       }
	  }
     }
     VecRestoreArray( VecX, &array );

}


PetscInt GetOwner_Position_Vector( MPI_Comm Comm, Vec Vector, PetscInt Row )
{
     PetscErrorCode ierr;
     PetscMPIInt size, rank;
     PetscInt rstart, rend;
     PetscInt i, Found, *Iams, Iam;

     Iams = NULL;
     MPI_Comm_size( Comm, &size );
     MPI_Comm_rank( Comm, &rank );

     PetscMalloc( size*sizeof(PetscInt), &Iams );
     Iam = -1;

     VecGetOwnershipRange( Vector, &rstart, &rend );

     if ( rstart <= Row && Row < rend ){
	  Iam = rank;
     }
	       
     MPI_Allgather( &Iam, 1, MPIU_INT, Iams, 1, MPIU_INT, Comm );

     i = 0;
     Found = 0;
     while( (i < size) && !Found ){
	  if( Iams[i] != -1 ){
	       Found = 1;
	       Iam = Iams[i];
	  } else {
	       i = i + 1;
	  }
     }
     PetscFree( Iams );
     return Iam;
}

PetscScalar* GetVec_Values( MPI_Comm Comm, Vec Vector, PetscInt NumRows, PetscInt *Rows )
{
     PetscMPIInt size, rank;
     PetscInt i,j, Iam;
     PetscInt rstart, rend;
     PetscScalar *Value;
     
     PetscMalloc( NumRows*sizeof(PetscScalar), &Value );
     VecGetOwnershipRange( Vector, &rstart, &rend );

     for( i = 0; i < NumRows; i++ ){
	  Iam = GetOwner_Position_Vector( Comm, Vector, Rows[i] );
	  VecAssemblyBegin( Vector );
	  VecAssemblyEnd( Vector );
	  if( (rstart <= Rows[i]) && (Rows[i] < rend) ){
	       VecGetValues( Vector, 1, &Rows[i], &Value[i] );
	  }
	  MPI_Bcast( &Value[i], 1, MPIU_SCALAR, Iam, Comm );
     }

     return Value;
     
}
