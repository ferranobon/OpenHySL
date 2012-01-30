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
	(*InitConst).ProcessGrid.Rows = (int) ( sqrt( (float) size ) );
	(*InitConst).ProcessGrid.Cols = size/(*InitConst).ProcessGrid.Rows;

	/* Block size */
	(*InitConst).BlockSize.Rows = 2;
	(*InitConst).BlockSize.Cols = 2;

	/* Order of the matrices and number of steps */
	(*InitConst).Order = 33;

	/* Order of the coupling nodes and their starting position */
	(*InitConst).OrderC = 1;
	(*InitConst).PosCouple = 31;

	/* Number of steps and Time step */
	(*InitConst).Nstep = 4096;
	(*InitConst).Delta_t = 0.01;

	/* Rayleigh values */
	(*InitConst).Rayleigh.Alpha = 0.1;
	(*InitConst).Rayleigh.Beta = 0.004;

	/* Newmark integration constants */
	(*InitConst).Newmark.Gamma = 0.5;
	(*InitConst).Newmark.Beta = 0.25;

	/* PID Constants */
	(*InitConst).PID.P = 0.75;
	(*InitConst).PID.I = 0.0;
	(*InitConst).PID.D = 0.0;

	/* Several constants to multiply the vectors */
	(*InitConst).Const1 = (*InitConst).Newmark.Beta*(*InitConst).Delta_t*(*InitConst).Delta_t;
	(*InitConst).Const2 = (0.5 - 2.0*(*InitConst).Newmark.Beta + (*InitConst).Newmark.Gamma)*(*InitConst).Delta_t*(*InitConst).Delta_t;
	(*InitConst).Const3 = (0.5 + (*InitConst).Newmark.Beta - (*InitConst).Newmark.Gamma)*(*InitConst).Delta_t*(*InitConst).Delta_t;

	/* Constants for Ending Step */
	(*InitConst).a0 = 1.0/((*InitConst).Newmark.Beta*(*InitConst).Delta_t*(*InitConst).Delta_t);
	(*InitConst).a1 = (*InitConst).Newmark.Gamma/((*InitConst).Newmark.Beta*(*InitConst).Delta_t);
	(*InitConst).a2 = 1.0/((*InitConst).Newmark.Beta*(*InitConst).Delta_t);
	(*InitConst).a3 = 1.0/(2.0*(*InitConst).Newmark.Beta) - 1.0;
	(*InitConst).a4 = (*InitConst).Newmark.Gamma/(*InitConst).Newmark.Beta - 1.0;
	(*InitConst).a5 = ((*InitConst).Delta_t/2.0)*((*InitConst).Newmark.Gamma/(*InitConst).Newmark.Beta - 2.0);
	(*InitConst).a6 = (1.0 - (*InitConst).Newmark.Gamma)*(*InitConst).Delta_t;
	(*InitConst).a7 = (*InitConst).Newmark.Gamma*(*InitConst).Delta_t;

	(*InitConst).FileM = "33M.txt";
	(*InitConst).FileK = "33K.txt";
	(*InitConst).FileC = "33C.txt";
	(*InitConst).FileData = "GroundMovement.txt";
}

void BroadcastConfFile( AlgConst *const InitConst )
{

	/* MPI Variables */
	int rank;

	int	LengthArrays;
	int i;     /* A counter */

	/* Setup three blocks */
	int          blockcounts[3] = {9, 18, 0};
	MPI_Datatype types[3];
	MPI_Aint     displs[3];
	MPI_Datatype InfoFile;

	MPI_Comm_rank( MPI_COMM_WORLD, &rank );

	if ( rank == 0 ){
		LengthArrays = strlen( (*InitConst).FileM ) + 1;
		LengthArrays = LengthArrays + strlen( (*InitConst).FileK ) + 1;
		LengthArrays = LengthArrays + strlen( (*InitConst).FileC ) + 1;
		LengthArrays = LengthArrays + strlen( (*InitConst).FileData ) + 1;
	}

	MPI_Bcast( &LengthArrays, 1, MPI_INT, 0, MPI_COMM_WORLD );

	blockcounts[2] = LengthArrays;

	/* Initialize types and displs with addresses anof items */
	MPI_Address( &(*InitConst).ProcessGrid, &displs[0] );
	MPI_Address( &(*InitConst).Delta_t,   &displs[1] );
	MPI_Address( &(*InitConst).FileM, &displs[2] );

	types[0] = MPI_INT;
	types[1] = MPI_FLOAT;
	types[2] = MPI_CHAR;

	/* Adjust the displacement array so that the displacements are offsets from the beginning of the structure */
	for (i = 2; i >=0; i--){
		displs[i] -= displs[0];
	}

	MPI_Type_create_struct(3, blockcounts, displs, types, &InfoFile );
	MPI_Type_commit( &InfoFile );

	MPI_Bcast( &(*InitConst), 1, InfoFile, 0, MPI_COMM_WORLD );
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

void BuildMatrixXc( MPI_Comm Comm, PMatrixVector *const Mat, float *MatCouple, const int PosCpl, const int orderc )
{

	int i, j, m;              /* Counters for the matrix Xc */
	int GRowIndex, GColIndex; /* Auxiliary variables to access the rows and columns of the Global matrix in the infog2l_( ) routine */

	/* Variables required by BLACS routines Cblacs_grindinfo( ) and infog2l( ) */
	int nprow, npcol, myrow, mycol;
	int LRowIndex, LColIndex;
	int RowProcess, ColProcess;

	/* Variables required by MPI routines */
	int rank;
	MPI_Status status;

	MPI_Comm_rank( Comm, &rank );

	/* Get grid info */
	Cblacs_gridinfo( (*Mat).Desc[1], &nprow, &npcol, &myrow, &mycol );

	j = 0;
	for (i = 0; i < orderc; i++){

		GRowIndex = PosCpl + i;
		GColIndex = PosCpl + j;

		/* Given the global index of an element (GRowIndex, GColIndex) returns the local index of the element (LRowIndex, LColIndex) and the coordinates of the process (Row Process, ColProcess) */
		infog2l_( &GRowIndex, &GColIndex, (*Mat).Desc, &nprow, &npcol, &myrow, &mycol, &LRowIndex, &LColIndex, &RowProcess, &ColProcess );

		/* Send the value to the process labeled as 0 only if the process are not the same. -1 is substracted because C starts the array indexes at 0, while FORTRAN starts them at 1 (infog2l is a FORTRAN routine) */
		if ( Cblacs_pnum( (*Mat).Desc[1], RowProcess, ColProcess ) == 0 && myrow == RowProcess && mycol == ColProcess ){
			MatCouple[i*orderc + j] = (*Mat).Array[(LRowIndex - 1) + (*Mat).LocalSize.Row*(LColIndex - 1)];
		} else {
			if ( myrow == RowProcess && mycol == ColProcess ){
				MPI_Send( &(*Mat).Array[(LRowIndex - 1) + (*Mat).LocalSize.Col*(LColIndex - 1)], 1, MPI_FLOAT, 0, 1, Comm );
			}

			if ( rank == 0 ){
				/* Store the diagonal elements */
				MPI_Recv( &MatCouple[i*orderc + j], 1, MPI_FLOAT, Cblacs_pnum( (*Mat).Desc[1], RowProcess, ColProcess ), 1, Comm, &status );
			}
		}

		for (m = 1; m < orderc -i; m++){

			GColIndex = PosCpl + j + m;
			infog2l_( &GRowIndex, &GColIndex, (*Mat).Desc, &nprow, &npcol, &myrow, &mycol, &LRowIndex, &LColIndex, &RowProcess, &ColProcess );

			if ( Cblacs_pnum( (*Mat).Desc[1], RowProcess, ColProcess ) == 0 && myrow == RowProcess && mycol == ColProcess ){
				MatCouple[i*orderc + j + m] = (*Mat).Array[(LRowIndex - 1) + (*Mat).LocalSize.Col*(LColIndex - 1)];
				MatCouple[i*orderc + j + m*orderc] = MatCouple[i*orderc + j + m];
			} else {
				if ( myrow == RowProcess && mycol == ColProcess ){
					MPI_Send( &(*Mat).Array[(LRowIndex - 1)*(*Mat).LocalSize.Col + (LColIndex - 1)], 1, MPI_FLOAT, 0, 1, Comm );
				}

				if ( rank == 0 ){
					/* Storing the triangular part of Xc */
					MPI_Recv( &MatCouple[i*orderc + j + m], 1, MPI_FLOAT, Cblacs_pnum( (*Mat).Desc[1], RowProcess, ColProcess ) , 1, Comm, &status );

					/* Storing the symmetric of Xc */
					MatCouple[i*orderc + j + m*orderc] = MatCouple[i*orderc + j + m];
				}
			}
		}
		j = j + 1;
	}
}

void BuildMatrixXcm( MPI_Comm Comm, PMatrixVector *const Mat, PMatrixVector *const VecXcm, const int PosCoupl, const int orderc )
{

	int i, j;  /* Counters */

	int GRowIndex, GColIndex; /* Auxiliary variables to access the rows and columns of the Global matrix in the infog2l_( ) routine */
	int GRowIndexXcm, GColIndexXcm;

	/* Variables required by BLACS routines Cblacs_grindinfo( ) and infog2l( ) */
	int nprow, npcol, myrow, mycol;
	int LRowIndex, LColIndex;
	int RowProcess, ColProcess;

	int nprowXcm, npcolXcm, myrowXcm, mycolXcm;
	int LRowIndexXcm, LColIndexXcm;
	int RowProcessXcm, ColProcessXcm;
	int aux;

	int rank;
	MPI_Status status;


	MPI_Comm_rank( Comm, &rank );

	/* Get grid info */
	Cblacs_gridinfo( (*Mat).Desc[1], &nprow, &npcol, &myrow, &mycol );
	Cblacs_gridinfo( (*VecXcm).Desc[1], &nprowXcm, &npcolXcm, &myrowXcm, &mycolXcm );


	for ( i = 0; i < orderc; i++ ){
		for( j = 0; j < PosCoupl -1; j++ ){
			GRowIndex = PosCoupl + i;
			GColIndex = j + 1;


			GRowIndexXcm = j + 1;
			GColIndexXcm = i + 1;

			/*
			 * Given the global index of an element (GRowIndex, GColIndex) returns the local index of the element (LRowIndex, LColIndex)
			 * and the coordinates of the process (Row Process, ColProcess)
			 */
			infog2l_( &GRowIndex, &GColIndex, (*Mat).Desc, &nprow, &npcol, &myrow, &mycol, &LRowIndex, &LColIndex, &RowProcess, &ColProcess );
			infog2l_( &GRowIndexXcm, &GColIndexXcm, (*VecXcm).Desc, &nprowXcm, &npcolXcm, &myrowXcm, &mycolXcm, &LRowIndexXcm, &LColIndexXcm, &RowProcessXcm, &ColProcessXcm );

			if ( Cblacs_pnum( (*VecXcm).Desc[1], RowProcessXcm, ColProcessXcm) == Cblacs_pnum( (*Mat).Desc[1], RowProcess, ColProcess) && myrow == RowProcess && mycol==ColProcess && myrowXcm==RowProcessXcm && mycolXcm == ColProcessXcm ){
				(*VecXcm).Array[(LColIndexXcm -1)*(*VecXcm).LocalSize.Row + (LRowIndexXcm -1)] = (*Mat).Array[(LRowIndex - 1) + (*Mat).LocalSize.Row*(LColIndex - 1)];
			} else {
				if ( myrow == RowProcess && mycol == ColProcess ){
					MPI_Send( &(*Mat).Array[(LRowIndex - 1) + (*Mat).LocalSize.Row*(LColIndex - 1)], 1, MPI_FLOAT, Cblacs_pnum( (*VecXcm).Desc[1], RowProcessXcm, ColProcessXcm ) , 1, Comm );
				}

				if ( myrowXcm == RowProcessXcm && mycolXcm == ColProcessXcm ){
					MPI_Recv( &(*VecXcm).Array[(LColIndexXcm -1)*(*VecXcm).LocalSize.Row + (LRowIndexXcm -1)], 1, MPI_FLOAT, Cblacs_pnum( (*Mat).Desc[1], RowProcess, ColProcess ) , 1, Comm, &status );
				}
			}
		}
	}

	aux = GRowIndexXcm;

	for ( j = 0; j < orderc; j++ ){

		GRowIndexXcm = aux;

		for ( i = 0; i < (*Mat).GlobalSize.Row - (PosCoupl + orderc - 1); i++ ){

			GRowIndex = PosCoupl + orderc + i;
			GColIndex = PosCoupl + j;

			GRowIndexXcm = GRowIndexXcm + 1;

			infog2l_( &GRowIndex, &GColIndex, (*Mat).Desc, &nprow, &npcol, &myrow, &mycol, &LRowIndex, &LColIndex, &RowProcess, &ColProcess );
			infog2l_( &GRowIndexXcm, &GColIndexXcm, (*VecXcm).Desc, &nprowXcm, &npcolXcm, &myrowXcm, &mycolXcm, &LRowIndexXcm, &LColIndexXcm, &RowProcessXcm, &ColProcessXcm );

			if ( Cblacs_pnum( (*VecXcm).Desc[1], RowProcessXcm, ColProcessXcm) == Cblacs_pnum( (*Mat).Desc[1], RowProcess, ColProcess) && myrow == RowProcess && mycol==ColProcess && myrowXcm==RowProcessXcm && mycolXcm == ColProcessXcm ){
				(*VecXcm).Array[(LColIndexXcm -1)*(*VecXcm).LocalSize.Row + (LRowIndexXcm -1)] = (*Mat).Array[(LRowIndex - 1) + (*Mat).LocalSize.Row*(LColIndex - 1)];
			} else {

				if ( myrow == RowProcess && mycol == ColProcess ){

					MPI_Send( &(*Mat).Array[(LRowIndex - 1) + (*Mat).LocalSize.Row*(LColIndex - 1)], 1, MPI_FLOAT, Cblacs_pnum( (*VecXcm).Desc[1], RowProcessXcm, ColProcessXcm ) , 1, Comm );
				}

				if ( myrowXcm == RowProcessXcm && mycolXcm == ColProcessXcm ){
					MPI_Recv( &(*VecXcm).Array[(LColIndexXcm -1)*(*VecXcm).LocalSize.Row + (LRowIndexXcm -1)], 1, MPI_FLOAT, Cblacs_pnum( (*Mat).Desc[1], RowProcess, ColProcess ) , 1, Comm, &status );
				}
			}
		}
		GColIndexXcm = GColIndexXcm + 1;
	}
}
