#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <petscmat.h>

#include "MatrixVector.h"
#include "Initiation.h"
#include "Precalculations.h"
#include "ComputeU0.h"
#include "Send_Receive_Data.h"
#include "EndingStep.h"


static char help[] = "Reads in a Symmetric matrix in MatrixMarket format. Writes\n\
it using the PETSc sparse format. It also adds a Vector set to random values to the\n\
output file. Input parameters are:\n\
  -fin <filename> : input file\n\
  -fout <filename> : output file\n\n";


int main( int argc, char **argv )
{

     AlgConst InitCnt;
     Scalars Constants;

     PetscInt istep;

     PetscScalar *AccAll, *VelAll, *DispAll;
     PetscErrorCode ierr;
     PetscMPIInt size, rank;

     Mat Mass, Stiff, Damp;
     Mat Keinv;
     Mat Keinv_m; PetscScalar *Keinv_c;


     Vec EffT;
     Vec DispT, DispTdT0, DispTdT;
     Vec DispTdT0_m;
     Vec tempvec;

     Vec VelT, VelTdT;
     Vec AccT, AccTdT;

     Vec LoadVectorForm, LoadTdT, LoadTdT1;

     Vec fc, fcprevsub;
     Vec fu;

     Vec ErrCompT, ErrCompTdT;

     Vec Disp, Vel, Acc;

     Coupling_Node CNodes;
     int Socket;
     PetscScalar *DispTdT0_c;

     PetscInitialize( &argc, &argv, (char *) 0, help );

     ierr = MPI_Comm_size( PETSC_COMM_WORLD, &size ); CHKERRQ( ierr );
     ierr = MPI_Comm_rank( PETSC_COMM_WORLD, &rank ); CHKERRQ( ierr );


     if( rank == 0 ){
	  /* Constants definitions. */
	  InitConstants( &InitCnt, "ConfFile.conf" );
	  
	  /* Read the coupling nodes from a file */
	  Read_Coupling_Nodes( &CNodes, InitCnt.FileCNodes );
     }

     BroadcastConfFile( &InitCnt );
     BroadCast_Coupling_Nodes( &CNodes );

     /* Load the Mass matrix */
     PETSc_LoadMatrix_FromFile( PETSC_COMM_WORLD, InitCnt.FileM, &Mass,
		       InitCnt.Order, InitCnt.Order );

     /* Load the Stiffness matrix */
     PETSc_LoadMatrix_FromFile( PETSC_COMM_WORLD, InitCnt.FileK, &Stiff,
		       InitCnt.Order, InitCnt.Order );

     /* Allocate memory and compute damping matrix */
     MatDuplicate( Stiff, MAT_DO_NOT_COPY_VALUES, &Damp );
     CalculateMatrixC( &Mass, &Stiff, &Damp, &InitCnt.Rayleigh );

     MatCreate( PETSC_COMM_WORLD, &Keinv );
     MatSetSizes( Keinv, PETSC_DECIDE, PETSC_DECIDE, InitCnt.Order, InitCnt.Order );
     MatSetType( Keinv, MATMPIDENSE );
     MatSetUp( Keinv );

     Constants.Alpha = 1.0f;
     Constants.Beta = InitCnt.a0;
     Constants.Gamma = InitCnt.a1;

     //EffK_Calculate_Keinv( Mass, Stiff, Damp, Keinv, Constants );
     PETSc_LoadMatrix_FromFile( PETSC_COMM_WORLD, "Keinv.bin", &Keinv,
			   InitCnt.Order, InitCnt.Order );

     if( rank == 0 ){
	  PetscMalloc( CNodes.Order*CNodes.Order*sizeof(PetscScalar), &Keinv_c );
	  PetscMalloc( CNodes.Order*sizeof(PetscScalar), &DispTdT0_c );
     }

     BuildMatrixXc( MPI_COMM_WORLD, Keinv, Keinv_c, &CNodes );
     if (rank == 0 ){
	  for ( istep = 0; istep < CNodes.Order*CNodes.Order; istep++ ){
	       printf("%lf\t", Keinv_c[istep] );
	  }
     }

     MatCreate( PETSC_COMM_WORLD, &Keinv_m );
     MatSetSizes( Keinv_m, PETSC_DECIDE, PETSC_DECIDE, InitCnt.Order - CNodes.Order, CNodes.Order );
     MatSetType( Keinv_m, MATMPIDENSE );
     MatSetUp( Keinv_m );

     BuildMatrixXcm( MPI_COMM_WORLD, Keinv, Keinv_m, &CNodes );

    /* Send the coupling part of the effective matrix */
     if( rank == 0 ){
	  Send_Effective_Matrix( (float *) Keinv_c, (unsigned int) CNodes.Order, &Socket, InitCnt.Remote );
     }

     /* Allocate the memory to store the earthquake recorded values */
     if( InitCnt.Use_Absolute_Values ){
	  PetscMalloc( InitCnt.Nstep*sizeof(PetscScalar), &VelAll );
	  PetscMalloc( InitCnt.Nstep*sizeof(PetscScalar), &DispAll );
	  AccAll = NULL;
	  if ( rank == 0 ){
	       ReadDataEarthquake_AbsValues( VelAll, DispAll, InitCnt.Nstep, InitCnt.FileData );
	  }
     } else {
	  PetscMalloc( InitCnt.Nstep*sizeof(PetscScalar), &AccAll );
	  VelAll = NULL;
	  DispAll = NULL;
	  if( rank == 0 ){
	       ReadDataEarthquake_RelValues( AccAll, InitCnt.Nstep, InitCnt.FileData );
	  }
     }

     if ( InitCnt.Use_Absolute_Values ){
	  MPI_Bcast( DispAll, InitCnt.Nstep, MPIU_SCALAR, 0, MPI_COMM_WORLD );
	  MPI_Bcast( VelAll, InitCnt.Nstep, MPIU_SCALAR, 0, MPI_COMM_WORLD );
     } else {
	  MPI_Bcast( AccAll, InitCnt.Nstep, MPIU_SCALAR, 0, MPI_COMM_WORLD );
     }


     istep = 1;

     PETSc_LoadVector_FromFile( MPI_COMM_WORLD, InitCnt.FileLVector, &LoadVectorForm, InitCnt.Order );

     VecDuplicate( LoadVectorForm, &Disp );
     VecDuplicate( LoadVectorForm, &Vel );
     VecDuplicate( LoadVectorForm, &Acc );

     VecDuplicate( LoadVectorForm, &LoadTdT ); VecSet( LoadTdT, 0.0 );

     if( InitCnt.Use_Absolute_Values ){
	  Apply_LoadVectorForm( Disp, LoadVectorForm, DispAll[istep - 1] );
	  Apply_LoadVectorForm( Vel, LoadVectorForm, VelAll[istep - 1] );
	  Calc_Input_Load_AbsValues( LoadTdT, Stiff, Damp, Disp, Vel );
     } else {
	  Apply_LoadVectorForm( Acc, LoadVectorForm, AccAll[istep - 1] );
	  Calc_Input_Load_RelValues( LoadTdT, Mass, Acc );
     }


     VecDuplicate( LoadVectorForm, &EffT ); VecSet( EffT, 0.0 );

     VecDuplicate( LoadVectorForm, &DispT ); VecSet( DispT, 0.0 );
     VecDuplicate( LoadVectorForm, &DispTdT0 ); VecSet( DispTdT0, 0.0 );
     VecSetUp(DispTdT0 );
     VecView( DispTdT0, PETSC_VIEWER_STDOUT_WORLD );

     VecDuplicate( LoadVectorForm, &DispTdT ); VecSet( DispTdT, 0.0 );

     VecDuplicate( LoadVectorForm, &tempvec ); VecSet( tempvec, 0.0 );

     VecDuplicate( LoadVectorForm, &VelT ); VecSet( VelT, 0.0 );
     VecDuplicate( LoadVectorForm, &VelTdT ); VecSet( VelTdT, 0.0 );

     VecDuplicate( LoadVectorForm, &AccT ); VecSet( AccT, 0.0 );
     VecDuplicate( LoadVectorForm, &AccTdT ); VecSet( AccTdT, 0.0 );

     VecDuplicate( LoadVectorForm, &LoadTdT1 ); VecSet( LoadTdT1, 0.0 );

     VecDuplicate( LoadVectorForm, &fc ); VecSet( fc, 0.0 );
     VecDuplicate( LoadVectorForm, &fu ); VecSet( fu, 0.0 );

     VecDuplicate( LoadVectorForm, &ErrCompT ); VecSet( ErrCompT, 0.0 );
     VecDuplicate( LoadVectorForm, &ErrCompTdT ); VecSet( ErrCompTdT, 0.0 );

     VecCreate( MPI_COMM_WORLD, &DispTdT0_m );
     VecSetSizes( DispTdT0_m, PETSC_DECIDE, InitCnt.Order - CNodes.Order );
     VecSetType( DispTdT0_m, VECMPI );
     VecSetUp( DispTdT0_m ); VecSet( DispTdT0_m, 0.0 );

     VecCreate( MPI_COMM_WORLD, &fcprevsub );
     VecSetSizes( fcprevsub, PETSC_DECIDE, CNodes.Order );
     VecSetType( fcprevsub, VECMPI );
     VecSetUp( fcprevsub ); VecSet( fcprevsub, 0.0 );

     if( rank == 0 ){
	  printf( "Starting stepping process\n" );
     }

     while( istep <= 500 ){
	  printf("istep: %d\n", istep );
	  /* Calculate the effective force vector
	     Fe = M*(a0*u + a2*v + a3*a) + C*(a1*u + a4*v + a5*a) */
	  EffK_Calc_Effective_Force( Mass, Damp, DispT, VelT, AccT, tempvec,
				     InitCnt.a0, InitCnt.a1, InitCnt.a2, InitCnt.a3, InitCnt.a4, InitCnt.a5,
				     EffT );

	  /* Compute the new Displacement u0 */
	  EffK_ComputeU0( EffT, LoadTdT, fu, InitCnt.PID.P, Keinv, tempvec, DispTdT0 );
	  
	  /* Split DispTdT into coupling and non-coupling part */
	  //CreateVectorXm( MPI_COMM_WORLD, DispTdT0, DispTdT0_m, &CNodes );
	  //CreateVectorXc( MPI_COMM_WORLD, DispTdT0, DispTdT0_c, &CNodes );
	  
	  if(istep == 500 ){
      
	       if(rank == 0){
		    printf("DispTdT0_c %e\n",DispTdT0_c[0] );
 	       }
	  }

	  /* Perform substepping */
	  Do_Substepping( MPI_COMM_WORLD, (float *) DispTdT0_c, DispTdT, fcprevsub, fc, InitCnt.Remote.Type, InitCnt.Delta_t*(float) istep, Socket, (unsigned int) CNodes.Order, (unsigned int *) CNodes.Array  );

	  if ( istep < InitCnt.Nstep ){
	       /* Calculate the input load for the next step during the
		  sub-stepping process. */
	       if( InitCnt.Use_Absolute_Values ){
		    Apply_LoadVectorForm( Disp, LoadVectorForm, DispAll[istep] );
		    Apply_LoadVectorForm( Vel, LoadVectorForm, VelAll[istep] );
		    Calc_Input_Load_AbsValues( LoadTdT1, Stiff, Damp, Disp, Vel );
		    VecView( LoadTdT1, PETSC_VIEWER_STDOUT_WORLD );
	       } else {
		    Apply_LoadVectorForm( Acc, LoadVectorForm, AccAll[istep] );
		    Calc_Input_Load_RelValues( LoadTdT1, Mass, Acc );
	       }
	  }

	  /* Join the non-coupling part. DispTdT_m = Keinv_m*fc + DispTdT0_m. Although DispTdT0_m is what has been received from the other computer,
	     it has the name of DispTdT_m to avoid further operations if using the NETLIB libraries. */
	  //JoinNonCouplingPart( MPI_COMM_WORLD, DispTdT0_m, Keinv_m, fcprevsub, DispTdT, &CNodes );

	  /* Compute acceleration ai1 = a0*(ui1 -ui) - a2*vi -a3*ai */
	  Compute_Acceleration( DispTdT, DispT, VelT, AccT, InitCnt.a0, InitCnt.a2, InitCnt.a3, AccTdT );

	  /* Compute Velocity. vi = vi + a6*ai +a7*ai */
	  Compute_Velocity( VelT, AccT, AccTdT, InitCnt.a6, InitCnt.a7, VelTdT );

	  /* Error Compensation. fu = LoatTdT + fc -(Mass*AccTdT + Damp*VelTdT + Stiff*DispTdT) */
	  Compute_Force_Error( Mass, Damp, Stiff, AccTdT, VelTdT, DispTdT, fc, LoadTdT, fu );

	  /* Backup vectors */
	  VecCopy(  LoadTdT1, LoadTdT );
	  VecCopy( DispTdT, DispT );
	  VecCopy( VelTdT, VelT );
	  VecCopy( AccTdT, AccT );

	  istep = istep + 1;
     }

     /* Close the Connection */
     Close_Connection( &Socket, InitCnt.Remote.Type, (unsigned int) CNodes.Order, InitCnt.Nstep, 4 );
 

     if( rank == 0 ){
	  PetscFree( Keinv_c );
	  free( DispTdT0_c );
     }

     free( CNodes.Array );

     MatDestroy( &Mass );
     MatDestroy( &Stiff );
     MatDestroy( &Damp );

     MatDestroy( &Keinv );
     MatDestroy( &Keinv_m );

     VecDestroy( &EffT );

     VecDestroy( &DispT );
     VecDestroy( &DispTdT );
     VecDestroy( &DispTdT0 );
     VecDestroy( &DispTdT0_m );

     VecDestroy( &tempvec );

     VecDestroy( &VelT );
     VecDestroy( &VelTdT );

     VecDestroy( &AccT );
     VecDestroy( &AccTdT );

     VecDestroy( &LoadVectorForm );

     VecDestroy( &LoadTdT );
     VecDestroy( &LoadTdT1 );

     VecDestroy( &fc );
     VecDestroy( &fcprevsub );

     VecDestroy( &fu );

     VecDestroy( &ErrCompT );
     VecDestroy( &ErrCompTdT );

     VecDestroy( &Disp );
     VecDestroy( &Vel );
     VecDestroy( &Acc );

     PetscFinalize( );
     return 0;
}

