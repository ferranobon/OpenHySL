#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <petscksp.h>

#include "ErrorHandling.h"
#include "Initiation.h"
#include "Send_Receive_Data.h"
#include "Conf_Parser.h"

void InitConstants( AlgConst *const InitConst, const char *FileName )
{

     /* Use Relative or absolute values */
     InitConst->Use_Absolute_Values = 1;
     if ( InitConst->Use_Absolute_Values != 0 && InitConst->Use_Absolute_Values != 1 ){
	  PrintErrorAndExit( "Invalid option for Use_Absolute_Values" );
     }

     /* Order of the matrices */
     InitConst->Order = 504;
     if ( InitConst->Order <= 0 ){
	  PrintErrorAndExit( "Invalid option for the order of the matrices" );
     }

     /* Number of steps and Time step */
     InitConst->Nstep = 4096;
     if ( InitConst->Nstep <= 0 ){
	  PrintErrorAndExit( "Invalid number of steps" );
     }

     InitConst->Delta_t = 0.01f;
     if ( InitConst->Delta_t <= 0.0f ){
	  PrintErrorAndExit( "Invalid time step" );
     }

     /* Rayleigh values */
     InitConst->Rayleigh.Alpha = 1.4f;
     InitConst->Rayleigh.Beta = 0.0004f;

     /* Newmark integration constants */
     InitConst->Newmark.Gamma = 0.5f;
     InitConst->Newmark.Beta = 0.25f;

     /* PID Constants */
     InitConst->PID.P = 0.95f;
     InitConst->PID.I = 0.0f;
     InitConst->PID.D = 0.0f;

     /* Several constants to multiply the vectors */
     InitConst->Const1 = InitConst->Newmark.Beta*InitConst->Delta_t*InitConst->Delta_t;
     InitConst->Const2 = (0.5f - 2.0f*InitConst->Newmark.Beta + InitConst->Newmark.Gamma)*InitConst->Delta_t*InitConst->Delta_t;
     InitConst->Const3 = (0.5f + InitConst->Newmark.Beta - InitConst->Newmark.Gamma)*InitConst->Delta_t*InitConst->Delta_t;

     /* Constants for Ending Step */
     InitConst->a0 = 1.0f/(InitConst->Newmark.Beta*InitConst->Delta_t*InitConst->Delta_t);
     InitConst->a1 = InitConst->Newmark.Gamma/(InitConst->Newmark.Beta*InitConst->Delta_t);
     InitConst->a2 = 1.0f/(InitConst->Newmark.Beta*InitConst->Delta_t);
     InitConst->a3 = 1.0f/(2.0f*InitConst->Newmark.Beta) - 1.0f;
     InitConst->a4 = InitConst->Newmark.Gamma/InitConst->Newmark.Beta - 1.0f;
     InitConst->a5 = (InitConst->Delta_t/2.0f)*(InitConst->Newmark.Gamma/InitConst->Newmark.Beta - 2.0f);
     InitConst->a6 = (1.0f - InitConst->Newmark.Gamma)*InitConst->Delta_t;
     InitConst->a7 = InitConst->Newmark.Gamma*InitConst->Delta_t;

     /* File Names */
     InitConst->FileM =  "504M.bin";
     InitConst->FileK = "504K.bin";
     InitConst->FileC = "Damp.txt";
     InitConst->FileLVector = "504LV.bin";
     InitConst->FileCNodes = "Couple_Nodes.txt";
     InitConst->FileData = "GroundMovement.txt";
     InitConst->FileOutput = "Out.txt";

     InitConst->Remote.Type = Get_Type_Protocol( "TCPCustom" );
     InitConst->Remote.IP = "192.168.1.53";
     InitConst->Remote.Port = "3333";
     InitConst->Remote.Account_Name = "NotValid.txt";
     InitConst->Remote.Account_Password = "NotValid.txt";
}

void BroadcastConfFile( AlgConst *const InitConst )
{

	/* MPI Variables */
	int rank;

	int	LengthArrays;
	int i;     /* A counter */

	/* Setup three blocks */
	int          blockcounts[5] = {3, 19, 0, 0, 1};
	MPI_Datatype types[5];
	MPI_Aint     displs[5];
	MPI_Datatype InfoFile;

	MPI_Comm_rank( MPI_COMM_WORLD, &rank );

	if ( rank == 0 ){
	     LengthArrays = strlen( (*InitConst).FileM ) + 1;
	     LengthArrays = LengthArrays + strlen( (*InitConst).FileK ) + 1;
	     LengthArrays = LengthArrays + strlen( (*InitConst).FileC ) + 1;
	     LengthArrays = LengthArrays + strlen( (*InitConst).FileLVector ) + 1;
	     LengthArrays = LengthArrays + strlen( (*InitConst).FileCNodes ) + 1;
	     LengthArrays = LengthArrays + strlen( (*InitConst).FileData ) + 1;
	     LengthArrays = LengthArrays + strlen( (*InitConst).FileOutput ) + 1;
	}
	
	MPI_Bcast( &LengthArrays, 1, MPI_INT, 0, MPI_COMM_WORLD );

	blockcounts[2] = LengthArrays;

	if( rank == 0 ){
	     LengthArrays = strlen( (*InitConst).Remote.IP ) + 1;
	     LengthArrays = LengthArrays + strlen( (*InitConst).Remote.Port ) + 1;
	     LengthArrays = LengthArrays + strlen( (*InitConst).Remote.Account_Name ) + 1;
	     LengthArrays = LengthArrays + strlen( (*InitConst).Remote.Account_Password ) + 1;
	}

	MPI_Bcast( &LengthArrays, 1, MPI_INT, 0, MPI_COMM_WORLD );

	blockcounts[3] = LengthArrays;


	/* Initialize types and displs with addresses anof items */
	MPI_Address( &(*InitConst).Order, &displs[0] );
	MPI_Address( &(*InitConst).Delta_t,   &displs[1] );
	MPI_Address( &(*InitConst).FileM, &displs[2] );
	MPI_Address( &(*InitConst).Remote, &displs[3] );
	MPI_Address( &(*InitConst).Remote.Type, &displs[4] );


	types[0] = MPIU_INT;
	types[1] = MPIU_SCALAR;
	types[2] = MPI_CHAR;
	types[3] = MPI_CHAR;
	types[4] = MPI_INT;

	/* Adjust the displacement array so that the displacements are offsets from the beginning of the structure */
	for (i = 4; i >=0; i--){
		displs[i] -= displs[0];
	}

	MPI_Type_create_struct( 5, blockcounts, displs, types, &InfoFile );
	MPI_Type_commit( &InfoFile );

	MPI_Bcast( &(*InitConst), 1, InfoFile, 0, MPI_COMM_WORLD );
}


void Delete_InitConstants( AlgConst *const InitConst )
{

     free( InitConst->FileM );
     free( InitConst->FileK );
     if( InitConst->FileC != NULL ){
	  free( InitConst->FileC );
     }
     free( InitConst->FileLVector );
     free( InitConst->FileCNodes );
     free( InitConst->FileData );
     free( InitConst->FileOutput );

     Delete_ServerInformation( &InitConst->Remote );
}

void Read_Coupling_Nodes( Coupling_Node *const CNodes, const char *Filename )
{
     FILE *InFile;
     PetscInt i;
     
     InFile = fopen( Filename, "r" );

     if ( InFile != NULL ){
	  /* The first value should be the number of Coupling nodes */
	  fscanf( InFile, "%i", &CNodes->Order );
	  
	  /* Allocate the necessary memory */
	  PetscMalloc( CNodes->Order*sizeof(PetscInt), &CNodes->Array );
	  
	  /* Read the contents of the file */
	  for( i = 0; i < CNodes->Order; i++ ){
	       fscanf( InFile, "%i", &CNodes->Array[i] );
	  }

	  /* Close the file */
	  fclose( InFile );
     } else ErrorFileAndExit( "It is not possible to read data because it was not possible to open: ", Filename );

}

void BroadCast_Coupling_Nodes( Coupling_Node *const CNodes )
{

     PetscMPIInt rank;

     MPI_Comm_rank( MPI_COMM_WORLD, &rank );

     MPI_Bcast( &CNodes->Order, 1, MPIU_INT, 0, MPI_COMM_WORLD );

     if( rank != 0 ){
	  PetscMalloc( CNodes->Order*sizeof(PetscInt), &CNodes->Array );
     }

     MPI_Bcast( CNodes->Array, CNodes->Order, MPIU_INT, 0, MPI_COMM_WORLD );
}


void CalculateMatrixC( const Mat *const Mass, const Mat *const Stiff, Mat *const Damp, const RayleighConst *const Rayleigh )
{

     MatAXPY( (*Damp), Rayleigh->Alpha, (*Mass), SUBSET_NONZERO_PATTERN );
     MatAXPY( (*Damp), Rayleigh->Beta, (*Stiff), SAME_NONZERO_PATTERN );
}

void EffK_Calculate_Keinv( Mat Mass, Mat Stiff, Mat Damp, Mat Keinv, const Scalars Const )
{

     KSP ksp;
     PC pc;
     PetscScalar Value = 1.0f, *Array;
     PetscInt rstart, rend, lsize, j;
     Mat Temp;
     Vec x, b;
     PetscInt IndexR, IndexC, Index;

     MatDuplicate( Stiff, MAT_DO_NOT_COPY_VALUES, &Temp );
     MatAXPY( Temp, Const.Alpha, Stiff, SAME_NONZERO_PATTERN );
     MatAXPY( Temp, Const.Beta, Mass, SUBSET_NONZERO_PATTERN );
     MatAXPY( Temp, Const.Gamma, Damp, SAME_NONZERO_PATTERN );
     MatGetSize( Temp, &IndexR, &IndexC );

     VecCreate( PETSC_COMM_WORLD, &x );
     VecSetSizes( x, PETSC_DECIDE, IndexR );
     VecSetType( x, VECMPI );

     VecDuplicate( x, &b );
     VecSet( b, 0.0 );

     for( Index = 0; Index < IndexR; Index++ ){
	  KSPCreate( PETSC_COMM_WORLD, &ksp);
	  
	  printf("Index: %d\n", Index );
	  /* 
	     Set operators. Here the matrix that defines the linear system
	     also serves as the preconditioning matrix.
	  */
	  KSPSetOperators(ksp,Temp,Temp,DIFFERENT_NONZERO_PATTERN);
	  
	  KSPGetPC(ksp,&pc);
	  PCSetType(pc,PCJACOBI);
	  KSPSetTolerances(ksp,1e-5f,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
	  KSPSetFromOptions(ksp);
	  
	  KSPSetUp( ksp );

	  if( Index > 0 ){
	       VecSetValue( b, Index -1, 0.0, INSERT_VALUES );
	  }
	  VecSetValue( b, Index, Value, INSERT_VALUES );

	  VecAssemblyBegin( b );
	  VecAssemblyEnd( b );


	  KSPSolve(ksp,b,x);
	  VecGetOwnershipRange( x, &rstart, &rend );
	  VecGetLocalSize( x, &lsize );
	  VecGetArray( x, &Array );
	  for( j = 0; j < lsize; j++ ){
	       MatSetValue( Keinv, rstart + j, Index, Array[j], INSERT_VALUES );
	  }
	  MatAssemblyBegin( Keinv, MAT_FINAL_ASSEMBLY );
	  MatAssemblyEnd( Keinv, MAT_FINAL_ASSEMBLY );
	  
	  VecRestoreArray( x, &Array );
	  KSPDestroy( &ksp );
   }

     MatView( Keinv, PETSC_VIEWER_STDOUT_WORLD );

     MatDestroy( &Temp );
     VecDestroy( &b );
     VecDestroy( &x );
}

void BuildMatrixXc( MPI_Comm Comm, Mat Matrix, PetscScalar *MatCouple, const Coupling_Node *const CNodes )
{

     PetscInt icoup;    /* Counter for the coupling nodes */
     PetscInt jcoup;    /* Another counter */
     PetscInt GRow, GCol, rstart, rend;
     PetscScalar Value;
     MPI_Status status;
     PetscInt Iam;
     PetscMPIInt size, rank;

     MPI_Comm_size( Comm, &size );
     MPI_Comm_rank( Comm, &rank );
     
     for( icoup = 0; icoup < CNodes->Order; icoup++ ){
	  
	  for (jcoup = icoup; jcoup < CNodes->Order; jcoup++){
	       GRow = CNodes->Array[icoup] -1;
	       GCol = CNodes->Array[jcoup] -1;
	       
	       MatGetOwnershipRange( Matrix, &rstart, &rend );
	       Iam = GetOwner_Position( Comm, Matrix, GRow ); 

	       if( rstart <= GRow && GRow < rend ){
		    MatGetValues( Matrix, 1, &GRow, 1, &GCol, &Value );
	       }

	       if ( Iam == rank && rank == 0 ){
		    MatCouple[icoup*CNodes->Order + jcoup] = Value;
		    if( icoup != jcoup ){
			 MatCouple[jcoup*CNodes->Order + icoup] = Value;
		    }
	       } else {
		    if ( rank == Iam ){
			 MPI_Send( &Value, 1, MPIU_SCALAR, 0, 1, Comm );
		    }
		    if ( rank == 0 ){
			 MPI_Recv( &MatCouple[icoup*CNodes->Order + jcoup], 1, MPIU_SCALAR, Iam, 1, Comm, &status );
			 if( icoup != jcoup ){
			      MatCouple[jcoup*CNodes->Order + icoup] = MatCouple[icoup*CNodes->Order + jcoup];
			 }
		    }
	       }
	  }
     }
}

void BuildMatrixXcm( MPI_Comm Comm, Mat Matrix, Mat MatXcm, const Coupling_Node *const CNodes )
{

     PetscInt Start_Position, Position;
     PetscInt Rows, Cols;
     PetscInt RowsXcm, ColsXcm;
     PetscInt i, *indxm, *indxn, irow;
     PetscScalar *Values;
     PetscInt Size;

     MatGetSize( Matrix, &Rows, &Cols );
     MatGetSize( MatXcm, &RowsXcm, &ColsXcm );

     Size = RowsXcm*ColsXcm;

     PetscMalloc( RowsXcm*sizeof(PetscInt), &indxm );
     PetscMalloc( ColsXcm*sizeof(PetscInt), &indxn );

     /* Initialise the index of columns */
     for( i = 0; i < CNodes->Order; i++ ){
	  indxn[i] = CNodes->Array[i] - 1;
     }

     Start_Position = 0;
     irow = 0;
     for( i = 0; i < CNodes->Order; i++ ){
	  Position = CNodes->Array[i] - 1;
	  while ( Start_Position < Position ){
	       indxm[irow] = Start_Position;
	       Start_Position = Start_Position + 1;
	       irow = irow + 1;
	  }
	  Start_Position = Start_Position + 1;
     }

     while( Start_Position < Rows ){
	  indxm[irow] = Start_Position;
	  irow = irow + 1;
	  Start_Position = Start_Position + 1;
     }

     Values = GetMat_Value( Comm, Matrix, RowsXcm, indxm, ColsXcm, indxn );

     /* Transform to Xcm coordinates */
     for( i = 0; i < RowsXcm; i++ ){
	  indxm[i] = i;
     }
     for( i = 0; i < ColsXcm; i++ ){
	  indxn[i] = i;
     }

     MatSetValues( MatXcm, RowsXcm, indxm, ColsXcm, indxn, Values, INSERT_VALUES );
     MatAssemblyBegin( MatXcm, MAT_FINAL_ASSEMBLY );
     MatAssemblyEnd( MatXcm, MAT_FINAL_ASSEMBLY );

     PetscFree( Values );
     PetscFree( indxm );
     PetscFree( indxn );
	  
}

PetscScalar* GetMat_Value( MPI_Comm Comm, Mat Matrix, PetscInt NumRows, PetscInt *Rows, PetscInt NumCols, PetscInt *Cols )
{
     PetscMPIInt size, rank;
     PetscInt i,j, Iam;
     PetscInt rstart, rend;
     PetscScalar *Value;

     PetscMalloc( NumRows*NumCols*sizeof(PetscScalar), &Value );
     MatGetOwnershipRange( Matrix, &rstart, &rend );
     
     for( i = 0; i < NumRows; i++ ){
	  for( j = 0; j < NumCols; j++ ){
	       Iam = GetOwner_Position( Comm, Matrix, Rows[i] );
	       if( rstart <= Rows[i] && Rows[i] < rend ){
		    MatGetValues( Matrix, 1, &Rows[i], 1, &Cols[j], &Value[i*NumCols + j] );
	       }
	       MPI_Bcast( &Value[i*NumCols + j], 1, MPIU_SCALAR, Iam, Comm );
	  }
     }
     return Value;
     
}

PetscInt GetOwner_Position( MPI_Comm Comm, Mat Matrix, PetscInt Row )
{
     PetscMPIInt size, rank;
     PetscInt rstart, rend;
     PetscInt i, Found, *Iams, Iam;

     Iams = NULL;
     MPI_Comm_size( Comm, &size );
     MPI_Comm_rank( Comm, &rank );

     PetscMalloc( size*sizeof(PetscInt), &Iams );
     Iam = -1;

     MatGetOwnershipRange( Matrix, &rstart, &rend );

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
