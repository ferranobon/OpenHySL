/**
 * \file Initiation.c
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 10th of August 2011
 * \todo Add support for packaged storages to decrease the memory use.
 * \todo Make the routines BuildMatrixXc() and BuildMatrixXcm() to be able to handle non consecutive coupling nodes.
 *
 * \brief Source code of the functions used in the Initiation phase.
 *
 * This file contains the source code of the routines used during the initiation phase of the substructure algorithm: construction of the Proportional Viscous
 * Damping Matrix, the Effective Mass Matrix and the Gain Matrix and their non-coupling and coupling part.
 */

#include <stdio.h>   /* For fprintf( ) and stderr */
#include <stdlib.h>  /* For exit( ) */
#include <string.h>  /* For strcmp( ) */
#include "ErrorHandling.h"
#include "Initiation.h"
#include "MatrixVector.h"
#include "Netlib.h"
#include "Send_Receive_Data.h"

void InitConstants( AlgConst *const InitConst )
{

	/* Order of the matrices */
	(*InitConst).Order = 1;

	/* Order of the coupling nodes and their starting position */
	(*InitConst).OrderC = 1;
	(*InitConst).PosCouple = 1;

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
	(*InitConst).PID.P = 0.95;
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

	/* File Names */
/*EFAST*/
	(*InitConst).FileM = "1M.txt";
	(*InitConst).FileK = "1K.txt";
	(*InitConst).FileC = "33C.txt";
	(*InitConst).FileData = "GroundMovement.txt";

	/* Get the communication protocol to be used */
	(*InitConst).Type_Protocol = Get_Type_Protocol( );
	if( (*InitConst).Type_Protocol == -1 ){
	     PrintErrorAndExit( "Invalid Protocol type." );
	}
}

int Get_Type_Protocol( )
{
     Remote_Machine_Info Remote;

     GetServerInformation( &Remote );

     if ( !strcmp( Remote.Type, "None" ) ){
	  return 0;
     } else if ( !strcmp( Remote.Type, "Custom" ) ){
	  return 1;
     } else if ( !strcmp( Remote.Type, "PNSE" ) ){
	  return 2; 
     } else if ( !strcmp( Remote.Type, "OpenFresco" ) ){
	  return 3; 
     } else {
	  return -1;
     }
}

void CalculateMatrixC( const Dense_MatrixVector *const Mass, const Dense_MatrixVector *const Stif, Dense_MatrixVector *const Damp, const RayleighConst *const Rayleigh )
{
	char uplo;
	int ione;
	int incx, incy;
	int info, lda, ldb;
	float done;

	int Rows, Cols, Length;
	float alpha, beta;

	int i;

	ione = 1; done = 1.0;
	incx = 1; incy = 1;
	uplo = 'L';  /* The lower part of the matrix will be used; the upper part will strictly not be referenced */
	Rows = (*Mass).Rows;
	Cols = (*Mass).Cols;
	alpha = (*Rayleigh).Alpha;
	beta = (*Rayleigh).Beta;

	lda = Max( 1, Rows );
	ldb = Max( 1, (*Damp).Rows );

	/* LAPACK: C = M */
	slacpy_( &uplo, &Rows, &Cols, (*Mass).Array, &lda, (*Damp).Array, &ldb );

	/* LAPACK: C = Rayleigh.alpha*M = Rayleigh.alpha*C */
	slascl_( &uplo, &ione, &ione, &done, &alpha, &(*Damp).Rows, &(*Damp).Cols, (*Damp).Array, &lda, &info );

	if ( info < 0 ){
		LAPACKPErrorAndExit( "dlascl: The ", -info, "-th argument had an illegal value" );
	}

	/* BLAS: C = alpha*M + beta*K = C + beta*K. Only half of the matrix is calculated */
	for ( i = 0; i < (*Damp).Rows; i++ ){
		Length = (*Damp).Rows - i;
		saxpy_( &Length, &beta, &(*Stif).Array[i*(*Stif).Rows + i], &incx, &(*Damp).Array[i*(*Damp).Rows +i], &incy);
	}

}

void CalculateMatrixKeinv( Dense_MatrixVector *const Keinv, const Dense_MatrixVector *const Mass, const Dense_MatrixVector *const Damp, const Dense_MatrixVector *const Stif, const Scalars Const )
{

	char uplo;
	int lda, info;

	uplo = 'L';  /* The lower part of the matrix will be used; the upper part will strictly not be referenced */
	lda = Max( 1, (*Keinv).Rows );

	/* Perform Meinv = [M + gamma*Delta_t*C + beta*Delta_t^2*K] */
	Dense_Add3Mat( &(*Keinv), &(*Stif), &(*Mass), &(*Damp), Const );

	/* LAPACK: Compute the Cholesky factorization of the symmetric positive definite matrix Meinv */
	spotrf_( &uplo, &(*Keinv).Rows, (*Keinv).Array, &lda, &info );

	if ( info == 0 ){
		printf( "Cholesky factorization successfully completed.\n" );
	}
	else if (info < 0){
		LAPACKPErrorAndExit( "Cholesky factorization: the ", -info, "th argument has an illegal value." );
	} else if (info > 0){
		LAPACKPErrorAndExit("Cholesky factorization: the leading minor of order ", info, " is not positive definite, and the factorization could not be completed." );
	}

	/* LAPACK: Compute the inverse of Me using the Cholesky factorization computed by pdpotrf_( ) */
	spotri_( &uplo, &(*Keinv).Rows, (*Keinv).Array, &lda, &info );

	if ( info == 0 ){
		printf( "Matrix Inversion successfully completed.\n" );
	} else if (info < 0){
		LAPACKPErrorAndExit( "Matrix Inversion: the ", -info, "th argument has an illegal value." );
	} else if (info > 0){
		fprintf( stderr, "Matrix Inversion: the (%d,%d) element of the factor U or L is zero, and the inverse could not be computed.\n", info, info );
		fprintf( stderr, "Exiting program.\n" );
		exit( EXIT_FAILURE );
	}

}

void BuildMatrixXc( const Dense_MatrixVector *const Mat, float *MatCouple, const int PosCpl, const int OrderC )
{

	int i, j, m;
	int GRowIndex, GColIndex;

	j = 0;
	for (i = 0; i < OrderC; i++){

		GRowIndex = PosCpl + i;
		GColIndex = PosCpl + j;

		MatCouple[i*OrderC + j] = (*Mat).Array[(GRowIndex - 1) + (*Mat).Rows*(GColIndex - 1)];

		for (m = 1; m < OrderC -i; m++){

			GColIndex = PosCpl + j + m;

			MatCouple[i*OrderC + j + m] = (*Mat).Array[(GRowIndex - 1) + (*Mat).Cols*(GColIndex - 1)];
			MatCouple[i*OrderC + j + m*OrderC] = MatCouple[i*OrderC + j + m];
		}
		j = j + 1;
	}
}

void BuildMatrixXcm( const Dense_MatrixVector *const Mat, Dense_MatrixVector *const VecXcm, const int PosCpl, const int OrderC )
{

	int i, j;  /* Counters */

	int GRowIndex, GColIndex;
	int GRowIndexXcm, GColIndexXcm;
	int aux;

	for ( i = 0; i < OrderC; i++ ){
		for( j = 0; j < PosCpl -1; j++ ){
			GRowIndex = PosCpl + i;
			GColIndex = j + 1;

			GRowIndexXcm = j + 1;
			GColIndexXcm = i + 1;

			(*VecXcm).Array[(GColIndexXcm -1)*(*VecXcm).Rows + (GRowIndexXcm -1)] = (*Mat).Array[(GRowIndex - 1) + (*Mat).Rows*(GColIndex - 1)];
		}
	}

	aux = GRowIndexXcm;

	for ( j = 0; j < OrderC; j++ ){

		GRowIndexXcm = aux;

		for ( i = 0; i < (*Mat).Rows - (PosCpl + OrderC - 1); i++ ){

			GRowIndex = PosCpl + OrderC + i;
			GColIndex = PosCpl + j;

			GRowIndexXcm = GRowIndexXcm + 1;

			(*VecXcm).Array[(GColIndexXcm -1)*(*VecXcm).Rows + (GRowIndexXcm -1)] = (*Mat).Array[(GRowIndex - 1) + (*Mat).Rows*(GColIndex - 1)];
		}
		GColIndexXcm = GColIndexXcm + 1;
	}
}
