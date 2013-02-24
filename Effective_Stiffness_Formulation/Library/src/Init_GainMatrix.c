/**
 * \file Init_GainMatrix.c
 * \author Ferran Obón Santacana
 * \version 1.0 PARDISO solver for matrix inversion is deprecated.
 * \date 9th of February 2013
 * 
 * \todo Add support for packed storage to reduce memory use.
 *
 * \brief Routines to compute the gain matrix.
 *
 * Routines for calculating the gain matrix. The routines make use of the BLAS and LAPACK libraries to perform the linear
 * algebra operations, including the matrix inversion in single and double precision. Sparse BLAS operations are supported
 * through the Intel MKL library, but since matrix inversion is a dense operations they still rely on LAPACK for this
 * operation. The PARDISO solver is no longer supported since it requires more memory and time than the LAPACK equivalent
 * routines.
 */
#include <stdlib.h>         /* For exit() */

#include "Print_Messages.h" /* For Print_Message() */
#include "Initiation.h"

#if _SPARSE_
#include <mkl.h>
#else 
#include "Netlib.h"
#endif

void IGainMatrix( MatrixVector_t *const Gain, const MatrixVector_t *const Mass, const MatrixVector_t *const Damp,
		  const MatrixVector_t *const Stiff, const Scalars Const )
{

     char uplo;    
     int lda, info; /* Leading dimension of the matrix, LAPACK error handling variable */

     int ione = 1;     /* Integer of value one */
     double one = 1.0; /* Double precision one for dlascl() parameter cfrom */

     uplo = 'L';  /* The lower part of the matrix will be used and the upper part will strictly
		   * not be referenced */
     lda = Max( 1, Gain->Rows );

     /* Gain = Const.Alpha*M + Const.Beta*C + Const.Gamma*K */
     Add3Mat( Gain, Mass, Damp, Stiff, Const );

     /* LAPACK: Compute the Cholesky factorization */
     dpotrf_( &uplo, &Gain->Rows, Gain->Array, &lda, &info );

     if ( info == 0 ){
	  Print_Message( SUCCESS, 2, STRING, "Cholesky factorization successfully completed." );
     }
     else if (info < 0){
	  Print_Message( ERROR, 3, STRING, "Cholesky factorization: the ", INT, -info, STRING, "th argument has an illegal value." );
	  exit( EXIT_FAILURE );
     } else if (info > 0){
	  Print_Message( ERROR, 3, STRING, "Cholesky factorization: the leading minor of order ", INT, info,
			 STRING, " is not positive definite, and the factorization could not be completed." );
     }

     /* LAPACK: Compute the inverse of the Gain matrix using the Cholesky factorization computed
      * by dpotrf_( ) */
     dpotri_( &uplo, &Gain->Rows, Gain->Array, &lda, &info );

     if ( info == 0 ){
	  Print_Message( SUCCESS, 1, STRING, "Matrix Inversion successfully completed." );
	  exit( EXIT_FAILURE );
     } else if (info < 0){
	  Print_Message( ERROR, 3, STRING, "Matrix inversion: the ", INT, -info, STRING, "th argument has an illegal value." );
     } else if (info > 0){
	  Print_Message( ERROR, 5, STRING, "Matrix Inversion: the (", INT, info, STRING, "," ,INT, info,
			 STRING, ") element of the factor U or L is zero, and the inverse could not be computed.\n" );
	  exit( EXIT_FAILURE );
     }

     dlascl_( &uplo, &ione, &ione, &one, Scalars.Lambda, &Gain->Rows, &Gain->Cols, Gain->Array, &lda, &info );
     if( info < 0 ){
	  Print_Message( ERROR, 3, STRING, "Gain matrix: the ", INT, -info, STRING, "th argument has an illegal value and cannot be scaled." );
	  exit( EXIT_FAILURE );
     } else {
	  Print_Message( SUCCESS, 1, STRING, "Gain matrix successfully calculated.\n" );
     }

}

#if _SPARSE_
void IGainMatrix_Sp( MatrixVector_t *const Gain, const MatrixVector_Sp_t *const Mass, const MatrixVector_Sp_t *const Damp,
		    const MatrixVector_Sp_t *const Stiff, const Scalars Const )
{

     char uplo;
     int lda, info;              /* Leading dimension of the gain matrix, LAPACK error handling variable */
     MatrixVector_Sp Sp_TempMat; /* Temporal sparse matrix */

     int ione = 1;     /* Integer of value one */
     double one = 1.0; /* Double precision one for dlascl() parameter cfrom */

     Init_MatrixVector_Sp( &Sp_TempMat, Damp->Rows, Damp->Cols, Damp->Num_Nonzero );
     /* Gain = Const.Alpha*M + Const.Beta*C + Const.Gamma*K */
     Add3Mat_Sp( &Sp_TempMat, Mass, Damp, Stiff, Const );

     /* Sparse to dense conversion. The Gain matrix will be symmetrical and only the upper part (lower part
      * in FORTRAN) will be referenced.
      */
     CSR_to_Dense( &Sp_TempMat, Gain, 0 );
     Destroy_MatrixVector_Sp( &Sp_TempMat );

     uplo = 'L';  /* The lower part of the matrix will be used and the upper part will strictly not be
		   * referenced */
     lda = Max( 1, Gain->Rows );

     /* LAPACK: Compute the Cholesky factorization */
     dpotrf_( &uplo, &Gain->Rows, Gain->Array, &lda, &info );

     if ( info == 0 ){
	  Print_Message( SUCCESS, 2, STRING, "Cholesky factorization successfully completed." );
     }
     else if (info < 0){
	  Print_Message( ERROR, 3, STRING, "Cholesky factorization: the ", INT, -info, STRING, "th argument has an illegal value." );
	  exit( EXIT_FAILURE );
     } else if (info > 0){
	  Print_Message( ERROR, 3, STRING, "Cholesky factorization: the leading minor of order ", INT, info,
			 STRING, " is not positive definite, and the factorization could not be completed." );
     }

     /* LAPACK: Compute the inverse of the Gain matrix using the Cholesky factorization computed
      * by dpotrf_( ) */
     dpotri_( &uplo, &Gain->Rows, Gain->Array, &lda, &info );

     if ( info == 0 ){
	  Print_Message( SUCCESS, 1, STRING, "Matrix Inversion successfully completed." );
	  exit( EXIT_FAILURE );
     } else if (info < 0){
	  Print_Message( ERROR, 3, STRING, "Matrix inversion: the ", INT, -info, STRING, "th argument has an illegal value." );
     } else if (info > 0){
	  Print_Message( ERROR, 5, STRING, "Matrix Inversion: the (", INT, info, STRING, "," ,INT, info,
			 STRING, ") element of the factor U or L is zero, and the inverse could not be computed.\n" );
	  exit( EXIT_FAILURE );
     }

     dlascl_( &uplo, &ione, &ione, &one, Scalars.Lambda, &Gain->Rows, &Gain->Cols, Gain->Array, &lda, &info );
     if( info < 0 ){
	  Print_Message( ERROR, 3, STRING, "Gain matrix: the ", INT, -info, STRING, "th argument has an illegal value and cannot be scaled." );
	  exit( EXIT_FAILURE );
     } else {
	  Print_Message( SUCCESS, 1, STRING, "Gain matrix successfully calculated.\n" );
     }

}

/* This function is deprecated */
void IGainMatrix_Pardiso( MatrixVector_t *const Gain, const MatrixVector_t *const Mass, const MatrixVector_t *const Damp,
			 const MatrixVector_t *const Stiff, const Scalars Const )
{

     int i;
     int iparm[64];
     void *pt[64];
     int maxfct, mnum, msglvl, error;
     int mtype = 2;                   /* Real symmetric matrix */
     int phase;
     MatrixVector IdentMatrix;
     MatrixVector TempMat;
     MatrixVector_Sp Sp_TempMat;
     double ddum;                       /* Dummy double */
     int idum;                          /* Dummy integer */
     char uplo;                         /* Used in dlascl() */
     int info, lda;                     /* Error handling and Leading dimension for dlascl() */

     Init_MatrixVector( &TempMat, Gain->Rows, Gain->Cols );

     Add3Mat( &TempMat, &(*Stiff), &(*Mass), &(*Damp), Const );

     Dense_to_CSR( &TempMat, &Sp_TempMat, 0 );   /* Transform into CSR format */
     Destroy_MatrixVector( &TempMat );           /* Destroy the dense matrix */

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
     iparm[11] = 2;
     iparm[12] = 1;	/* Maximum weighted matching algorithm is switched-off (default for symmetric).
			 * Try iparm[12] = 1 in case of inappropriate accuracy */
     iparm[13] = 0;	/* Output: Number of perturbed pivots */
     iparm[17] = -1;	/* Output: Number of nonzeros in the factor LU */
     iparm[18] = -1;	/* Output: Mflops for LU factorization */
     iparm[19] = 0;	/* Output: Numbers of CG Iterations */
     iparm[27] = 0;     /* Input/output matrices are double precision */
     iparm[34] = 0;	/* PARDISO uses 1 based indexing for ia and ja arrays */
     maxfct = 1;	/* Maximum number of numerical factorizations. */
     mnum = 1;		/* Which factorization to use. */
     msglvl = 1;	/* Print statistical information in file */
     error = 0;		/* Initialize error flag */

     /* --------------------------------------------------------------------
      * .. Initialize the internal solver memory pointer. This is only
      * necessary for the FIRST call of the PARDISO solver.
      * -------------------------------------------------------------------- */
     for (i = 0; i < 64; i++)
     {
	  pt[i] = 0;
     }

     /* --------------------------------------------------------------------
      * .. Reordering and Symbolic Factorization. This step also allocates
      * all memory that is necessary for the factorization.
      * -------------------------------------------------------------------- */
     phase = 11;
     PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &Sp_TempMat.Rows, Sp_TempMat.Values, Sp_TempMat.RowIndex,
	      Sp_TempMat.Columns, &idum, &Gain->Rows, iparm, &msglvl, &ddum, &ddum, &error);
     if (error != 0){
	  Print_Message( ERROR, 3, "PARDISO: Error during symbolic factorization: ", INT, error, STRING, "." );
	  exit (1);
     }
     Print_Message( SUCCESS, 1, STRING, "PARDISO: Reordering completed ... " );
     Print_Message( INFO, 3, STRING, "PARDISO: Number of nonzeros in factors = ", INT, iparm[17], STRING, "." ); 
     Print_Message( INFO, 3, STRING, "PARDISO: Number of factorization MFLOPS = ", INT, iparm[18], STRING, "." ); 

     /* -------------------------------------------------------------------- */
     /* .. Numerical factorization. */
     /* -------------------------------------------------------------------- */
     phase = 22;
     
     PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &Sp_TempMat.Rows, Sp_TempMat.Values, Sp_TempMat.RowIndex,
	      Sp_TempMat.Columns, &idum, &Gain->Rows, iparm, &msglvl, &ddum, &ddum, &error);
     if (error != 0){
	  Print_Message( ERROR, 3, STRING, "PARDISO: Error during numerical factorization: ", INT, error, STRING, "." );
	  exit (2);
     }
     Print_Message( SUCCESS, 1, STRING, "PARDISO: Factorization completed ... " );

     /* --------------------------------------------------------------------
      * .. Back substitution and iterative refinement.
      * -------------------------------------------------------------------- */
     phase = 33;
     iparm[7] = 10; /* Max numbers of iterative refinement steps. */

     /* Set right hand side to be the identity matrix. */     
     IdentMatrix = Generate_IdentityMatrix( Gain->Rows, Gain->Cols );

     PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &Sp_TempMat.Rows, Sp_TempMat.Values, Sp_TempMat.RowIndex,
	      Sp_TempMat.Columns, &idum, &Gain->Rows, iparm, &msglvl, IdentMatrix.Array, Gain->Array, &error);

     /* --------------------------------------------------------------------
      *	.. Termination and release of memory.
      *	-------------------------------------------------------------------- */
     phase = -1; /* Release internal memory. */
     PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &Sp_TempMat.Rows, &ddum, Sp_TempMat.RowIndex, Sp_TempMat.Columns,
	      &idum, &Gain->Rows, iparm, &msglvl, &ddum, &ddum, &error);

     Destroy_MatrixVector( &IdentMatrix );
     Destroy_MatrixVector_Sp( &Sp_TempMat );

     Print_Message( SUCCESS, 1, STRING, "PARDISO: Matrix Inversion successfully completed." );

     ddum = 1.0; idum = 1;
     lda = Max( 1, Gain->Rows );
     uplo = 'L';   /* The lower part of the matrix will be used (upper part in C) */
     dlascl_( &uplo, &idum, &idum, &ddum, Scalars.Lambda, &Gain->Rows, &Gain->Cols, Gain->Array, &lda, &info );
     if( info < 0 ){
	  Print_Message( ERROR, 3, STRING, "Gain matrix: the ", INT, -info, STRING, "th argument has an illegal value and cannot be scaled." );
	  exit( EXIT_FAILURE );
     } else {
	  Print_Message( SUCCESS, 1, STRING, "PARDISO: Gain matrix successfully calculated." );
     }

}

/* This function is deprecated */
void IGainMatrix_Pardiso_Sp( MatrixVector_t *const Gain, const MatrixVector_Sp_t *const Mass, const MatrixVector_Sp_t *const Damp,
		       const MatrixVector_Sp_t *const Stiff, const Scalars Const )
{
     int i;
     int iparm[64];
     void *pt[64];
     int maxfct, mnum, msglvl, error;
     int mtype = 2;                   /* Real symmetric matrix */
     int phase;
     MatrixVector IdentMatrix;
     MatrixVector_Sp Sp_TempMat;
     double ddum;                       /* Dummy double */
     int idum;                          /* Dummy integer */
     char uplo;                         /* Used in dlascl() */
     int info, lda;                     /* Error handling and Leading dimension for dlascl() */


     Init_MatrixVector_Sp( &Sp_TempMat, Damp->Rows, Damp->Cols, Damp->Num_Nonzero );

     Add3Mat_Sp( &Sp_TempMat, &(*Stiff), &(*Mass), &(*Damp), Const );

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
     iparm[11] = 2;
     iparm[12] = 1;	/* Maximum weighted matching algorithm is switched-off (default for symmetric).
			 *  Try iparm[12] = 1 in case of inappropriate accuracy */
     iparm[13] = 0;	/* Output: Number of perturbed pivots */
     iparm[17] = -1;	/* Output: Number of nonzeros in the factor LU */
     iparm[18] = -1;	/* Output: Mflops for LU factorization */
     iparm[19] = 0;	/* Output: Numbers of CG Iterations */
     iparm[27] = 0;     /* Input/output matrices are double precision */
     iparm[34] = 0;	/* PARDISO use 1 based indexing for ia and ja arrays */
     maxfct = 1;	/* Maximum number of numerical factorizations. */
     mnum = 1;		/* Which factorization to use. */
     msglvl = 1;	/* Print statistical information in file */
     error = 0;		/* Initialize error flag */

     /* --------------------------------------------------------------------
      *	.. Initialize the internal solver memory pointer. This is only
      * necessary for the FIRST call of the PARDISO solver.
      * -------------------------------------------------------------------- */
     for (i = 0; i < 64; i++)
     {
	  pt[i] = 0;
     }

     /* --------------------------------------------------------------------
      * .. Reordering and Symbolic Factorization. This step also allocates
      * all memory that is necessary for the factorization.
      * -------------------------------------------------------------------- */
     phase = 11;

     PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &Sp_TempMat.Rows, Sp_TempMat.Values, Sp_TempMat.RowIndex,
	      Sp_TempMat.Columns, &idum, &Gain->Rows, iparm, &msglvl, &ddum, &ddum, &error);

     if (error != 0) {
	  Print_Message( ERROR, 3, "PARDISO: Error during symbolic factorization: ", INT, error, STRING, "." );
	  exit (1);
     }
     Print_Message( SUCCESS, 1, STRING, "PARDISO: Reordering completed ... " );
     Print_Message( INFO, 3, STRING, "PARDISO: Number of nonzeros in factors = ", INT, iparm[17], STRING, "." ); 
     Print_Message( INFO, 3, STRING, "PARDISO: Number of factorization MFLOPS = ", INT, iparm[18], STRING, "." ); 

     /* --------------------------------------------------------------------
      * .. Numerical factorization.
      * -------------------------------------------------------------------- */
     phase = 22;
     PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &Sp_TempMat.Rows, Sp_TempMat.Values, Sp_TempMat.RowIndex,
	      Sp_TempMat.Columns, &idum, &Gain->Rows, iparm, &msglvl, &ddum, &ddum, &error);
     if (error != 0){
	  Print_Message( ERROR, 3, STRING, "PARDISO: Error during numerical factorization: ", INT, error, STRING, "." );
	  exit (2);
     }
     Print_Message( SUCCESS, 1, STRING, "PARDISO: Factorization completed ... " );

     /* --------------------------------------------------------------------
      * .. Back substitution and iterative refinement.
      * -------------------------------------------------------------------- */
     phase = 33;
     iparm[7] = 10;	/* Max numbers of iterative refinement steps. */

     /* Set right hand side to be the identity matrix. */
     IdentMatrix = Generate_IdentityMatrix( Gain->Rows, Gain->Cols );

     PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &Sp_TempMat.Rows, Sp_TempMat.Values, Sp_TempMat.RowIndex,
	      Sp_TempMat.Columns, &idum, &Gain->Rows, iparm, &msglvl, IdentMatrix.Array, Gain->Array, &error);

     /* --------------------------------------------------------------------
      * .. Termination and release of memory.
      * -------------------------------------------------------------------- */
     phase = -1;       /* Release internal memory. */
     PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
	      &Sp_TempMat.Rows, &ddum, Sp_TempMat.RowIndex, Sp_TempMat.Columns, &idum, &Gain->Rows,
	      iparm, &msglvl, &ddum, &ddum, &error);

     Destroy_MatrixVector( &IdentMatrix );
     Destroy_MatrixVector_Sp( &Sp_TempMat );

     Print_Message( SUCCESS, 1, STRING, "PARDISO: Matrix Inversion successfully completed." );

     ddum = 1.0; idum = 1;
     lda = Max( 1, Gain->Rows );
     uplo = 'L';   /* The lower part of the matrix will be used (upper part in C) */
     dlascl_( &uplo, &idum, &idum, &ddum, Scalars.Lambda, &Gain->Rows, &Gain->Cols, Gain->Array, &lda, &info );
     if( info < 0 ){
	  Print_Message( ERROR, 3, STRING, "Gain matrix: the ", INT, -info, STRING, "th argument has an illegal value and cannot be scaled." );
	  exit( EXIT_FAILURE );
     } else {
	  Print_Message( SUCCESS, 1, STRING, "PARDISO: Gain matrix successfully calculated." );
     }
}
#endif /* _SPARSE */


