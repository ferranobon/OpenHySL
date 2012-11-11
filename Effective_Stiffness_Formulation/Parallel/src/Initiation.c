/*
 * Initiation.c
 *
 *  Created on: 22/07/2011
 *      Author: ferran
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>  /* For strlen( ) */
#include <math.h>    /* For sqrt( ) */
#include <assert.h>  /* For assert( ) */

#include <mpi.h>

#include "ErrorHandling.h"
#include "Initiation.h"
#include "Netlib.h"
#include "PMatrixVector.h"
#include "Send_Receive_Data.h"


void InitConstants( AlgConst *const InitConst, const int size )
{

     /* Grid information */
     (*InitConst).ProcessGrid.Rows = (int) ( sqrtf( (float) size ) );
     (*InitConst).ProcessGrid.Cols = size/(*InitConst).ProcessGrid.Rows;
     
     /* Block size */
     (*InitConst).BlockSize.Rows = 2;
     (*InitConst).BlockSize.Cols = 2;

     /* Use Relative or absolute values */
     InitConst->Use_Absolute_Values = 1;
     if ( InitConst->Use_Absolute_Values != 0 && InitConst->Use_Absolute_Values != 1 ){
	  PrintErrorAndExit( "Invalid option for Use_Absolute_Values" );
     }

     /* Order of the matrices */
     InitConst->Order = 33;
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
     InitConst->FileM =  "33M.txt";
     InitConst->FileK = "33K.txt";
     InitConst->FileC = "33C.txt";
     InitConst->FileLVector = "33LV.txt";
     InitConst->FileCNodes = "Couple_Nodes.txt";
     InitConst->FileData = "GroundMovement.txt";
     InitConst->FileOutput = "Out.txt";

     InitConst->Remote.Type = Get_Type_Protocol( "TCPCustom" );
     InitConst->Remote.IP = "192.168.1.53";
     InitConst->Remote.Port = "44000";
     InitConst->Remote.Account_Name = "NotValid.txt";
     InitConst->Remote.Account_Password = "NotValid.txt";
}

void BroadcastConfFile( AlgConst *const InitConst )
{

     /* MPI Variables */
     int rank;
     
     int LengthArrays;
     int i;     /* A counter */
     
     /* Setup three blocks */
     int          blockcounts[5] = {7, 19, 0, 0, 1};
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
     MPI_Address( &InitConst->ProcessGrid, &displs[0] );
     MPI_Address( &InitConst->Delta_t,   &displs[1] );
     MPI_Address( &InitConst->FileM, &displs[2] );
     MPI_Address( &(*InitConst).Remote, &displs[3] );
     MPI_Address( &(*InitConst).Remote.Type, &displs[4] );


     types[0] = MPI_INT;
     types[1] = MPI_FLOAT;
     types[2] = MPI_CHAR;
     types[3] = MPI_CHAR;
     types[4] = MPI_INT;

     /* Adjust the displacement array so that the displacements are offsets from the beginning of the structure */
     for (i = 4; i >=0; i--){
	  displs[i] -= displs[0];
     }

     MPI_Type_create_struct(5, blockcounts, displs, types, &InfoFile );
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
     int i;
     
     InFile = fopen( Filename, "r" );

     if ( InFile != NULL ){
	  /* The first value should be the number of Coupling nodes */
	  fscanf( InFile, "%i", &CNodes->Order );
	  
	  /* Allocate the necessary memory */
	  CNodes->Array = (int *) calloc( (size_t) CNodes->Order, sizeof(int) );
	  
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

     int rank;

     MPI_Comm_rank( MPI_COMM_WORLD, &rank );

     MPI_Bcast( &CNodes->Order, 1, MPI_INT, 0, MPI_COMM_WORLD );

     if( rank != 0 ){
	  CNodes->Array = (int *) calloc( (size_t) CNodes->Order, sizeof(int) );
     }

     MPI_Bcast( CNodes->Array, CNodes->Order, MPI_INT, 0, MPI_COMM_WORLD );
}

void CalculateMatrixC( PMatrixVector *const Mass, PMatrixVector *const Stif, PMatrixVector *const Damp, RayleighConst Rayleigh )
{

	char trans, uplo;
	int ione;

	ione = 1;
	trans = 'N'; /* The operation will not use the transpose matrix */
	uplo = 'L';  /* The lower part of the matrix will be used; the upper part will strictly not be referenced */

	/* ScaLAPACK: Perform C = M (locally. There is no communication) */
	pslacpy_( &uplo, &(*Damp).GlobalSize.Row, &(*Damp).GlobalSize.Col, (*Mass).Array, &ione, &ione, (*Mass).Desc, (*Damp).Array, &ione, &ione, (*Damp).Desc );

	/* ScaLAPACK: Perform C = alpha*M + beta*K = alpha*C + beta*K */
	pstradd_( &uplo, &trans, &(*Damp).GlobalSize.Row, &(*Damp).GlobalSize.Col, &Rayleigh.Beta, (*Stif).Array, &ione, &ione, (*Stif).Desc, &Rayleigh.Alpha, (*Damp).Array, &ione, &ione, (*Damp).Desc );

}

void CalculateMatrixKeinv( PMatrixVector *const Keinv, PMatrixVector *const Mass, PMatrixVector *const Damp, PMatrixVector *const Stif, Scalars Const )
{

	char uplo;
	int ione, info;

	ione = 1;
	uplo = 'L';  /* The lower part of the matrix will be used; the upper part will strictly not be referenced */

	/* Perform Keinv = [M + gamma*Delta_t*C + beta*Delta_t^2*K] */
	PAdd3Mat( &(*Keinv), &(*Stif), &(*Mass), &(*Damp), Const );

	/* ScaLAPACK: Compute the Cholesky factorization of the symmetric positive definite matrix Keinv */
	pspotrf_( &uplo, &(*Keinv).GlobalSize.Row, (*Keinv).Array, &ione, &ione, (*Keinv).Desc, &info );

	if ( info == 0 ){
		printf( "Cholesky factorization successfully completed.\n" );
	}
	else if (info < 0){
		LAPACKPErrorAndExit( "Cholesky factorization: the ", -info, "th argument has an illegal value." );
	} else if (info > 0){
		LAPACKPErrorAndExit( "Cholesky factorization: the leading minor of order ", info, " is not positive definite, and the factorization could not be completed.\n" );
	} else assert( 0 );

	// SCALAPACK: Compute the inverse of Me using the Cholesky factorization computed by pdpotrf_( )
	pspotri_( &uplo, &(*Keinv).GlobalSize.Row, (*Keinv).Array, &ione, &ione, (*Keinv).Desc, &info );

	if (info == 0){
		printf( "Matrix Inversion successfully completed.\n" );
	} else if (info < 0){
		LAPACKPErrorAndExit( "Matrix Inversion: the ", -info, "th argument has an illegal value." );
	} else if (info > 0){
		fprintf( stderr, "Matrix Inversion: the (%d,%d) element of the factor U or L is zero, and the inverse could not be computed.\n", info, info );
		fprintf( stderr, "Exiting program" );

		Cblacs_gridexit( (*Keinv).Desc[1] );
		MPI_Finalize( );
		exit( EXIT_FAILURE );
	} else assert( 0 );


}

void CalculateMatrixG( PMatrixVector *const Gain, PMatrixVector *const Keinv, float Const )
{

	int ione, info;
	char uplo;
	float cfrom;

	ione = 1; cfrom = 1.0;
	uplo = 'L';  /* The lower part of the matrix will be used; the upper part will strictly not be referenced */

	/* ScaLAPACK: Perform Gain = EffMassInf (locally. There is no communication) */
	pslacpy_( &uplo, &(*Gain).GlobalSize.Row, &(*Gain).GlobalSize.Col, (*Keinv).Array, &ione, &ione, (*Keinv).Desc, (*Gain).Array, &ione, &ione, (*Gain).Desc );

	/* ScaLAPACK: Perform Gain = beta*Delta_t^2*Gain = Const*Gain */
	pslascl_( &uplo, &cfrom, &Const, &(*Gain).GlobalSize.Row, &(*Gain).GlobalSize.Col, (*Gain).Array, &ione, &ione, (*Gain).Desc, &info );

	/* TODO: Implement error messages using info */
}

void BuildMatrixXc( MPI_Comm Comm, PMatrixVector *const Mat, float *MatCouple, const Coupling_Node *const CNodes )
{

     int icoup, jcoup;
     int nprow, npcol, myrow, mycol;
     int rank;
     int GRowIndex, GColIndex;
     int LRowIndex, LColIndex;
     int RowProcess, ColProcess;

     MPI_Status status;

     MPI_Comm_rank( Comm, &rank );

     /* Get grid info */
     Cblacs_gridinfo( (*Mat).Desc[1], &nprow, &npcol, &myrow, &mycol );
     
     for( icoup = 0; icoup < CNodes->Order; icoup++ ){
	  GColIndex = CNodes->Array[icoup];
	  for (jcoup = icoup; jcoup < CNodes->Order; jcoup++){
	       GRowIndex = CNodes->Array[jcoup];

	       /* Given the global index of an element (GRowIndex, GColIndex) returns the local index of the element
		  (LRowIndex, LColIndex) and the coordinates of the process (Row Process, ColProcess) */
	       infog2l_( &GRowIndex, &GColIndex, (*Mat).Desc, &nprow, &npcol, &myrow, &mycol, &LRowIndex, &LColIndex,
			 &RowProcess, &ColProcess );

	       if ( Cblacs_pnum( (*Mat).Desc[1], RowProcess, ColProcess ) == 0 && myrow == RowProcess && mycol == ColProcess ){ 
		    MatCouple[icoup*CNodes->Order + jcoup] = Mat->Array[(LRowIndex - 1) + Mat->LocalSize.Row*(LColIndex - 1)];

		    /* Now add the elements belonging to the same row as the current coupling
		     * node but also belonging to the same column as the rest of the coupling
		     * nodes */
		    MatCouple[jcoup*CNodes->Order + icoup] = MatCouple[icoup*CNodes->Order + jcoup];
	       } else {
		    if ( myrow == RowProcess && mycol == ColProcess ){			 
			 MPI_Send( &(*Mat).Array[(LRowIndex - 1)+(*Mat).LocalSize.Row*(LColIndex - 1)], 1, MPI_FLOAT, 0, 1, Comm );
		    }

		    if ( rank == 0 ){
			 /* Store the diagonal elements */
			 MPI_Recv( &MatCouple[icoup*CNodes->Order + jcoup], 1, MPI_FLOAT, Cblacs_pnum( (*Mat).Desc[1], RowProcess, ColProcess ), 1, Comm, &status );
			 /* Now add the elements belonging to the same row as the current coupling
			  * node but also belonging to the same column as the rest of the coupling
			  * nodes */
			 MatCouple[jcoup*CNodes->Order + icoup] = MatCouple[icoup*CNodes->Order + jcoup];
		    }
	       }
	  }

     }     
}

void BuildMatrixXcm( MPI_Comm Comm, PMatrixVector *const Mat, PMatrixVector *const VecXcm, const Coupling_Node *const CNodes )
{


     int Length, Accumulated_Length;
     int icoup;      /* Counter for the coupling nodes */
     int jcoup;
     int incx, incy;
     int PosXcm_Rows, PosXcm_Cols;
     int PosX_Rows, PosX_Cols;
     int rank;

     MPI_Comm_rank( Comm, &rank );

     incx = Mat->GlobalSize.Row;
     incy = 1;          /* The values in the Xm matrix are stored in columns. Therefore the stride
			 * has to be 1 because of the way FORTRAN handles arrays (major column ordering) */

     /* Since the matrix is symmetric and only a part of it is stored, this routine has to be splitted into two
      * parts. The first will copy the elements above the coupling node, while the second will focus on the
      * other part */

     /* Copy until the first coupling node */
     Length = CNodes->Array[0] - 1;
     Accumulated_Length = 0;
     for( icoup = 0; icoup < CNodes->Order; icoup++ ){
	  PosX_Rows = CNodes->Array[icoup];     /* 1 based index */
	  PosX_Cols = 1; /* 1 based index */
	  PosXcm_Rows = 1;
	  PosXcm_Cols = icoup + 1;
	  pscopy_( &Length, Mat->Array, &PosX_Rows, &PosX_Cols, Mat->Desc, &incx,
		   VecXcm->Array, &PosXcm_Rows, &PosXcm_Cols, VecXcm->Desc, &incy );
     }
     Accumulated_Length = Accumulated_Length + Length;
     for( icoup = 1; icoup < CNodes->Order; icoup++ ){
	  Length = CNodes->Array[icoup] - CNodes->Array[icoup - 1] - 1;
	  if ( Length > 0 ){
	       for( jcoup = icoup; jcoup < CNodes->Order; jcoup++ ){
		    
		    PosX_Rows = CNodes->Array[jcoup];     /* 1 based index */
		    PosX_Cols = CNodes->Array[icoup-1] + 1; /* 1 based index */
		    PosXcm_Rows = Accumulated_Length + 1;
		    PosXcm_Cols = jcoup + 1;		    
		    pscopy_( &Length, Mat->Array, &PosX_Rows, &PosX_Cols, Mat->Desc, &incx,
			     VecXcm->Array, &PosXcm_Rows, &PosXcm_Cols, VecXcm->Desc, &incy );
	       }
	       Accumulated_Length = Accumulated_Length + Length;
	  }
	  
     }

     incx = 1;
     Length = Mat->GlobalSize.Row - CNodes->Array[CNodes->Order - 1];
     if( Length > 0 ){
	  for ( icoup = CNodes->Order - 1; icoup >= 0; icoup = icoup -1 ){
	       PosX_Rows = CNodes->Array[CNodes->Order - 1] + 1;
	       PosX_Cols = CNodes->Array[icoup];
	       PosXcm_Rows = VecXcm->GlobalSize.Row - Length + 1;
	       PosXcm_Cols = icoup + 1;
	       pscopy_( &Length, Mat->Array, &PosX_Rows, &PosX_Cols, Mat->Desc, &incx,
			VecXcm->Array, &PosXcm_Rows, &PosXcm_Cols, VecXcm->Desc, &incy );
	  }
     }

     Accumulated_Length = Length;
     for( icoup = CNodes->Order -2; icoup >= 0; icoup = icoup -1 ){
	  Length = CNodes->Array[icoup + 1] - CNodes->Array[icoup] - 1;
	  if( Length > 0 ){
	       for( jcoup = icoup; jcoup >= 0; jcoup = jcoup - 1){
		    PosX_Rows = CNodes->Array[icoup] + 1;
		    PosX_Cols = CNodes->Array[jcoup];
		    PosXcm_Rows = VecXcm->GlobalSize.Row - Accumulated_Length - Length + 1;
		    PosXcm_Cols = jcoup + 1;	           
		    pscopy_( &Length, Mat->Array, &PosX_Rows, &PosX_Cols, Mat->Desc, &incx,
			     VecXcm->Array, &PosXcm_Rows, &PosXcm_Cols, VecXcm->Desc, &incy );
	       }
	       Accumulated_Length = Accumulated_Length + Length;
	  }

     }

}
