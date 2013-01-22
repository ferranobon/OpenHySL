#include <stdio.h>
#include <stdlib.h>

#include <petscmat.h>

#include "EndingStep.h"
#include "Initiation.h"
#include "ComputeU0.h"

void JoinNonCouplingPart( MPI_Comm Comm, Vec VecXm, Mat MatrixXm, Vec fcprevsub, Vec VecX, Coupling_Node *const CNodes )
{

     PetscInt i, icoup, *ix = NULL;
     PetscInt PosX, PosXm, Length;
     PetscInt Rows;
     PetscScalar *Values = NULL;


     MatMultAdd( MatrixXm, fcprevsub, VecXm, VecXm );

     PosX = 0; PosXm = 0;
     for( icoup = 0; icoup < CNodes->Order; icoup++ ){
	  Length = CNodes->Array[icoup] - PosX - 1;
	  if( Length > 0 ){
	       PetscMalloc( Length*sizeof(PetscInt), &ix );
	       for ( i = 0; i < Length; i++ ){
		    ix[i] = i + PosXm;
	       }
	       Values = GetVec_Values( Comm, VecXm, Length, ix );
	       
	       for( i = 0; i < Length; i++ ){
		    ix[i] = i + PosX;
	       }
	       
	       VecSetValues( VecX, Length, ix, Values, INSERT_VALUES );
	       VecAssemblyBegin( VecX );
	       VecAssemblyEnd( VecX );
	       
	       PosXm = PosXm + Length;
	       PetscFree( Values ); Values = NULL;
	       PetscFree( ix ); ix = NULL;
	  }
	  PosX = CNodes->Array[icoup];
     }

     VecGetSize( VecX, &Rows );
     Length = Rows - CNodes->Array[CNodes->Order - 1];  
     PetscMalloc( Length*sizeof(PetscInt), &ix );
     for ( i = 0; i < Length; i++ ){
	  ix[i] = i + PosXm;
     }

     Values = GetVec_Values( Comm, VecXm, Length, ix );
	       
     for( i = 0; i < Length; i++ ){
	  ix[i] = i + PosX;
     }
     
     VecSetValues( VecX, Length, ix, Values, INSERT_VALUES );
     VecAssemblyBegin( VecX );
     VecAssemblyEnd( VecX );
     
     PetscFree( Values );
     PetscFree( ix );
}

void Compute_Acceleration( Vec DispTdT, Vec DispT, Vec VelT, Vec AccT,
			   PetscScalar a0, PetscScalar a2, PetscScalar a3,
			   Vec AccTdT )
{
    
     /* AccTdT = DispTdT - DispT */
     VecWAXPY( AccTdT, -1.0, DispT, DispTdT );
     /* AccTdT = a0*(DispTdT - DispT) - a2*VelT = a0*AccTdT - a2*VelT */
     VecAXPBY( AccTdT, -1.0*a2, a0, VelT );
     /* AccTdT = a0*(DispTdT - DispT) - a2*VelT - a3*AccT = AccTdT - a3*AccT */
     VecAXPY( AccTdT, -1.0*a3, AccT );

}

void Compute_Velocity( Vec VelT, Vec AccT, Vec AccTdT,
		       PetscScalar a6, PetscScalar a7, 
		       Vec VelTdT )
{
     /* VelTdT = VelT + a6*AccT = */
     VecWAXPY( VelTdT, a6, AccT, VelT );
     /* VelTdT = VelT + a6*AccT + a7*AccTdT = VelTdT + a7*AccTdT */
     VecAXPY( VelTdT, a7, AccTdT );
}

void Compute_Force_Error( Mat Mass, Mat Damp, Mat Stiff, Vec AccTdT, Vec VelTdT, Vec DispTdT,
			  Vec fc, Vec LoadTdT, Vec fu )
{
     
     /* fu = Mass*AccTdT */
     MatMult( Mass, AccTdT, fu );
     /* fu = Mass*AccTdT + Damp*VelTdT = fu + Damp*VelTdT */
     MatMultAdd( Damp, VelTdT, fu, fu );
     /* fu = Mass*AccTdT + Damp*VelTdT + Stiff*DispTdT = fu + Stiff*DispTdT */
     MatMultAdd( Stiff, DispTdT, fu, fu );

     /* fu = fc -(Mass*AccTdT + Damp*VelTdT + Stiff*DispTdT) = fc -fu */
     VecAXPBY( fu, 1.0, -1.0, fc );
     /* fu = LoatTdT + fc -(Mass*AccTdT + Damp*VelTdT + Stiff*DispTdT) = LoadTdT + fu */
     VecAXPY( fu, 1.0, LoadTdT );
}
 
