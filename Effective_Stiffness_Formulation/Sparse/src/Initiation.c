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
#include "Conf_Parser.h"

void InitConstants( AlgConst *const InitConst, const char *FileName )
{

     ConfFile *Config;
     
     Config = ConfFile_Create( 35 );

     ConfFile_ReadFile( Config, FileName );

     /* Use Relative or absolute values */
     (*InitConst).Use_Absolute_Values = ConfFile_GetInt( Config, "General:Use_Absolute_Values" );
     if ( InitConst->Use_Absolute_Values != 0 && InitConst->Use_Absolute_Values != 1 ){
	  PrintErrorAndExit( "Invalid option for Use_Absolute_Values" );
     }

     /* Order of the matrices */
     (*InitConst).Order = ConfFile_GetInt( Config, "General:Order" );
     if ( InitConst->Order <= 0 ){
	  PrintErrorAndExit( "Invalid option for the order of the matrices" );
     }

     /* Number of steps and Time step */
     (*InitConst).Nstep = (unsigned int) ConfFile_GetInt( Config, "General:Num_Steps" );
     if ( InitConst->Nstep <= 0 ){
	  PrintErrorAndExit( "Invalid number of steps" );
     }

     (*InitConst).Delta_t = ConfFile_GetFloat( Config, "General:Delta" );
     if ( InitConst->Delta_t <= 0.0f ){
	  PrintErrorAndExit( "Invalid time step" );
     }

     /* Rayleigh values */
     (*InitConst).Rayleigh.Alpha = ConfFile_GetFloat( Config, "Rayleigh:Alpha" );
     (*InitConst).Rayleigh.Beta = ConfFile_GetFloat( Config, "Rayleigh:Beta" );

     /* Newmark integration constants */
     (*InitConst).Newmark.Gamma = ConfFile_GetFloat( Config, "Newmark:Gamma" );
     (*InitConst).Newmark.Beta = ConfFile_GetFloat( Config, "Newmark:Beta" );

     /* PID Constants */
     (*InitConst).PID.P = ConfFile_GetFloat( Config, "PID:P" );
     (*InitConst).PID.I = ConfFile_GetFloat( Config, "PID:I" );
     (*InitConst).PID.D = ConfFile_GetFloat( Config, "PID:D" );

     /* Several constants to multiply the vectors */
     (*InitConst).Const1 = (*InitConst).Newmark.Beta*(*InitConst).Delta_t*(*InitConst).Delta_t;
     (*InitConst).Const2 = (0.5f - 2.0f*(*InitConst).Newmark.Beta + (*InitConst).Newmark.Gamma)*(*InitConst).Delta_t*(*InitConst).Delta_t;
     (*InitConst).Const3 = (0.5f + (*InitConst).Newmark.Beta - (*InitConst).Newmark.Gamma)*(*InitConst).Delta_t*(*InitConst).Delta_t;

     /* Constants for Ending Step */
     (*InitConst).a0 = 1.0f/((*InitConst).Newmark.Beta*(*InitConst).Delta_t*(*InitConst).Delta_t);
     (*InitConst).a1 = (*InitConst).Newmark.Gamma/((*InitConst).Newmark.Beta*(*InitConst).Delta_t);
     (*InitConst).a2 = 1.0f/((*InitConst).Newmark.Beta*(*InitConst).Delta_t);
     (*InitConst).a3 = 1.0f/(2.0f*(*InitConst).Newmark.Beta) - 1.0f;
     (*InitConst).a4 = (*InitConst).Newmark.Gamma/(*InitConst).Newmark.Beta - 1.0f;
     (*InitConst).a5 = ((*InitConst).Delta_t/2.0f)*((*InitConst).Newmark.Gamma/(*InitConst).Newmark.Beta - 2.0f);
     (*InitConst).a6 = (1.0f - (*InitConst).Newmark.Gamma)*(*InitConst).Delta_t;
     (*InitConst).a7 = (*InitConst).Newmark.Gamma*(*InitConst).Delta_t;

     /* File Names */
/*EFAST*/
     (*InitConst).FileM = strdup( ConfFile_GetString( Config, "FileNames:Mass_Matrix" ) );
     (*InitConst).FileK = strdup( ConfFile_GetString( Config, "FileNames:Stiffness_Matrix" ) );
     (*InitConst).FileC = strdup( ConfFile_GetString( Config, "FileNames:Damping_Matrix" ) );
     (*InitConst).FileLVector = strdup( ConfFile_GetString( Config, "FileNames:FileLVector" ) );
     (*InitConst).FileCNodes = strdup( ConfFile_GetString( Config, "FileNames:Coupling_Nodes" ) );
     (*InitConst).FileData = strdup( ConfFile_GetString( Config, "FileNames:Ground_Motion" ) );
     (*InitConst).FileOutput = strdup( ConfFile_GetString( Config, "FileNames:OutputFile" ) );

     /* Get the communication protocol to be used */
     GetServerInformation( &InitConst->Remote, Config );
     ConfFile_Free( Config );
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

void CalculateMatrixC( const MatrixVector *const Mass, const MatrixVector *const Stif, MatrixVector *const Damp, const RayleighConst *const Rayleigh )
{
     char uplo;
     int ione;
     int incx, incy;
     int info, lda, ldb;
     float done;

     int Rows, Cols, Length;
     float alpha, beta;

     int i;

     ione = 1; done = 1.0f;
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

void CalculateMatrixKeinv( MatrixVector *const Keinv, const MatrixVector *const Mass, const MatrixVector *const Damp, const MatrixVector *const Stif, const Scalars Const )
{

     char uplo;
     int lda, info;

     uplo = 'L';  /* The lower part of the matrix will be used; the upper part will strictly not be referenced */
     lda = Max( 1, (*Keinv).Rows );

     /* Perform Meinv = [M + gamma*Delta_t*C + beta*Delta_t^2*K] */
     Add3Mat( &(*Keinv), &(*Stif), &(*Mass), &(*Damp), Const );

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

void CalculateMatrixKeinv_Pardiso( MatrixVector *const Keinv, const Sp_MatrixVector *const Mass, const Sp_MatrixVector *const Damp, const Sp_MatrixVector *const Stiff, const Scalars Const )
{

     int iparm[64];

     Add3Mat_Sparse( &(*Keinv), &(*Stif), &(*Mass), &(*Damp), Const );

     /* Setup the Pardiso control parameters */
     for (i = 0; i < 64; i++){
	  iparm[i] = 0;
     }
     iparm[0] = 1;	/* No solver default */
     iparm[1] = 2;	/* Fill-in reordering from METIS */
     iparm[3] = 0;	/* No iterative-direct algorithm */
     iparm[4] = 0;	/* No user fill-in reducing permutation */
     iparm[5] = 0;	/* Write solution into x */
     iparm[7] = 2;	/* Max numbers of iterative refinement steps */
     iparm[9] = 13;	/* Perturb the pivot elements with 1E-13 */
     iparm[10] = 1;	/* Use nonsymmetric permutation and scaling MPS */
     iparm[12] = 0;	/* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
     iparm[13] = 0;	/* Output: Number of perturbed pivots */
     iparm[17] = -1;	/* Output: Number of nonzeros in the factor LU */
     iparm[18] = -1;	/* Output: Mflops for LU factorization */
     iparm[19] = 0;	/* Output: Numbers of CG Iterations */
     iparm[28] = 1;     /* Input/output matrices are single precision */
     iparm[34] = 1;	/* PARDISO use C-style indexing for ia and ja arrays */
     maxfct = 1;	/* Maximum number of numerical factorizations. */
     mnum = 1;		/* Which factorization to use. */
     msglvl = 1;	/* Print statistical information in file */
     error = 0;		/* Initialize error flag */

     /* -------------------------------------------------------------------- */
     /* .. Initialize the internal solver memory pointer. This is only */
     /* necessary for the FIRST call of the PARDISO solver. */
     /* -------------------------------------------------------------------- */
     for (i = 0; i < 64; i++)
     {
	  pt[i] = 0;
     }
     /* -------------------------------------------------------------------- */
     /* .. Reordering and Symbolic Factorization. This step also allocates */
     /* all memory that is necessary for the factorization. */
     /* -------------------------------------------------------------------- */
     phase = 11;
     PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
	      &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
     if (error != 0)
     {
	  printf ("\nERROR during symbolic factorization: %d", error);
	  exit (1);
     }
     printf ("\nReordering completed ... ");
     printf ("\nNumber of nonzeros in factors = %d", iparm[17]);
     printf ("\nNumber of factorization MFLOPS = %d", iparm[18]);
/* -------------------------------------------------------------------- */
/* .. Numerical factorization. */
/* -------------------------------------------------------------------- */
     phase = 22;
     PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
	      &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
     if (error != 0)
     {
	  printf ("\nERROR during numerical factorization: %d", error);
	  exit (2);
     }
     printf ("\nFactorization completed ... ");
/* -------------------------------------------------------------------- */
/* .. Back substitution and iterative refinement. */
/* -------------------------------------------------------------------- */
     phase = 33;
     iparm[7] = 2;			/* Max numbers of iterative refinement steps. */
     /* Set right hand side to one. */
     for (i = 0; i < n; i++)
     {
	  b[i] = 1;
     }
     PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
	      &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &error);
     if (error != 0)
     {
	  printf ("\nERROR during solution: %d", error);
	  exit (3);
     }
     printf ("\nSolve completed ... ");
     printf ("\nThe solution of the system is: ");
     for (i = 0; i < n; i++)
     {
	  printf ("\n x [%d] = % f", i, x[i]);
     }
     printf ("\n");
/* -------------------------------------------------------------------- */
/* .. Termination and release of memory. */
/* -------------------------------------------------------------------- */
     phase = -1;			/* Release internal memory. */
     PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
	      &n, &ddum, ia, ja, &idum, &nrhs,
	      iparm, &msglvl, &ddum, &ddum, &error);
     return 0;
}

     

void BuildMatrixXc( const MatrixVector *const Mat, float *MatCouple, const Coupling_Node *const CNodes )
{

     int icoup;    /* Counter for the coupling nodes */
     int jcoup;    /* Another counter */
     
     for( icoup = 0; icoup < CNodes->Order; icoup++ ){

	  for (jcoup = icoup; jcoup < CNodes->Order; jcoup++){
	       
	       MatCouple[icoup*CNodes->Order + jcoup] = Mat->Array[(CNodes->Array[icoup]-1)*Mat->Cols + CNodes->Array[jcoup] - 1];
	       /* Now add the elements belonging to the same row as the current coupling
		* node but also belonging to the same column as the rest of the coupling
		* nodes */
	       MatCouple[jcoup*CNodes->Order + icoup] = MatCouple[icoup*CNodes->Order + jcoup];
	  }
     }
}

void BuildMatrixXcm( const MatrixVector *const Mat, MatrixVector *const VecXcm, const Coupling_Node *const CNodes )
{

     int Length;
     int icoup;      /* Counter for the coupling nodes */
     int jcoup;
     int incx, incy;
     int PosXcm, Acumulated_Length;

     incx = Mat->Rows;
     incy = 1;          /* The values in the Xm matrix are stored in columns. Therefore the stride
			 * has to be 1 because of the way FORTRAN handles arrays (major column ordering) */

     /* Since the matrix is symmetric and only a part of it is stored, this routine has to be splitted into two
      * parts. The first will copy the elements above the coupling node, while the second will focus on the
      * other part */

     /* Copy until the first coupling node */
     PosXcm = 0;
     Length = CNodes->Array[0] - 1;
     for ( jcoup = 0; jcoup < CNodes->Order; jcoup++ ){
	  PosXcm = jcoup*VecXcm->Rows;
	  scopy_( &Length, &Mat->Array[CNodes->Array[jcoup] - 1], &incx, &VecXcm->Array[PosXcm], &incy );
     }

     /* Copy until the last coupling node */
     Acumulated_Length = Length;
     for( icoup = 1; icoup < CNodes->Order; icoup++ ){
	     
	  Length = CNodes->Array[icoup] - CNodes->Array[icoup-1] - 1;
	  for ( jcoup = icoup; jcoup < CNodes->Order; jcoup++ ){
	 
	       PosXcm = jcoup*VecXcm->Rows + Acumulated_Length;
	       
	       scopy_( &Length, &Mat->Array[(CNodes->Array[jcoup] - 1) + (CNodes->Array[icoup-1])*Mat->Rows], &incx, &VecXcm->Array[PosXcm], &incy );
	  }
	  Acumulated_Length = Acumulated_Length + Length;
     }

     /* Do the same but starting from the opposite side */
     incx = 1;
     Length = Mat->Rows - CNodes->Array[CNodes->Order -1];
     for ( jcoup = CNodes->Order - 1; jcoup >= 0; jcoup = jcoup - 1 ){
	  PosXcm = VecXcm->Rows*(jcoup+1) - Length;
	  scopy_( &Length, &Mat->Array[(CNodes->Array[jcoup] - 1)*Mat->Rows + (CNodes->Array[CNodes->Order-1])], &incx, &VecXcm->Array[PosXcm], &incy );
     }
     Acumulated_Length = Length;
     incx = 1;
     for( icoup = CNodes->Order -2; icoup >= 0; icoup = icoup -1 ){
	  Length = CNodes->Array[icoup + 1] - CNodes->Array[icoup] - 1;
	  for( jcoup = icoup; jcoup >= 0; jcoup = jcoup - 1){
	       PosXcm = (jcoup+1)*VecXcm->Rows - Acumulated_Length - Length;
	       scopy_( &Length, &Mat->Array[(CNodes->Array[jcoup] - 1)*Mat->Rows + (CNodes->Array[icoup])], &incx, &VecXcm->Array[PosXcm], &incy );
	  }
	  Acumulated_Length = Acumulated_Length + Length;
     }
}
