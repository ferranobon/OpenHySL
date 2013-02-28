#include <stdlib.h>   /* For exit() */

#include "MatrixVector.h"
#include "MatrixVector_Sp.h"

#include "Initiation.h"
#include "Initiation_Sp.h"

#include "Auxiliary_Math.h"

#include "Print_Messages.h"

#include <mkl_lapack.h>
#include <mkl_pardiso.h>

void IGainMatrix_Sp( MatrixVector_t *const Gain, const MatrixVector_Sp_t *const Mass, const MatrixVector_Sp_t *const Damp,
		    const MatrixVector_Sp_t *const Stiff, const Scalars_t Const )
{

     char uplo;
     int lda, info;              /* Leading dimension of the gain matrix, LAPACK error handling variable */
     MatrixVector_Sp_t Sp_TempMat; /* Temporal sparse matrix */

     int ione = 1;     /* Integer of value one */
     double one = 1.0; /* Double precision one for dlascl() parameter cfrom */
     double Scalar;    /* Double precision scalar */

     MatrixVector_Create_Sp( Damp->Rows, Damp->Cols, Damp->Num_Nonzero, &Sp_TempMat );
     /* Gain = Const.Alpha*M + Const.Beta*C + Const.Gamma*K */
     MatrixVector_Add3Mat_Sp(Mass, Damp, Stiff, Const, &Sp_TempMat );

     /* Sparse to dense conversion. The Gain matrix will be symmetrical and only the upper part (lower part
      * in FORTRAN) will be referenced.
      */
     MatrixVector_CSR2Dense( &Sp_TempMat, 0, Gain );
     MatrixVector_Destroy_Sp( &Sp_TempMat );

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

     Scalar = Const.Lambda;
     dlascl_( &uplo, &ione, &ione, &one, &Scalar, &Gain->Rows, &Gain->Cols, Gain->Array, &lda, &info );
     if( info < 0 ){
	  Print_Message( ERROR, 3, STRING, "Gain matrix: the ", INT, -info, STRING, "th argument has an illegal value and cannot be scaled." );
	  exit( EXIT_FAILURE );
     } else {
	  Print_Message( SUCCESS, 1, STRING, "Gain matrix successfully calculated.\n" );
     }

}

/* This function is deprecated */
void IGainMatrix_Pardiso( MatrixVector_t *const Gain, const MatrixVector_t *const Mass, const MatrixVector_t *const Damp,
			 const MatrixVector_t *const Stiff, const Scalars_t Const )
{

     int i;
     int iparm[64];
     void *pt[64];
     int maxfct, mnum, msglvl, error;
     int mtype = 2;                   /* Real symmetric matrix */
     int phase;
     MatrixVector_t IdentMatrix;
     MatrixVector_t TempMat;
     MatrixVector_Sp_t Sp_TempMat;
     double ddum;                       /* Dummy double */
     int idum;                          /* Dummy integer */
     char uplo;                         /* Used in dlascl() */
     int info, lda;                     /* Error handling and Leading dimension for dlascl() */
     double Scalar;    /* Double precision scalar */

     MatrixVector_Create( Gain->Rows, Gain->Cols, &TempMat );

     MatrixVector_Add3Mat( Stiff, Mass, Damp, Const, &TempMat );

     MatrixVector_Dense2CSR( &TempMat, 0, &Sp_TempMat );   /* Transform into CSR format */
     MatrixVector_Destroy( &TempMat );           /* Destroy the dense matrix */

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
     PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &Sp_TempMat.Rows, Sp_TempMat.Values, Sp_TempMat.RowIndex,
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

     MatrixVector_Destroy( &IdentMatrix );
     MatrixVector_Destroy_Sp( &Sp_TempMat );

     Print_Message( SUCCESS, 1, STRING, "PARDISO: Matrix Inversion successfully completed." );

     ddum = 1.0; idum = 1;
     lda = Max( 1, Gain->Rows );
     uplo = 'L';   /* The lower part of the matrix will be used (upper part in C) */

     Scalar = Const.Lambda;
     dlascl_( &uplo, &idum, &idum, &ddum, &Scalar, &Gain->Rows, &Gain->Cols, Gain->Array, &lda, &info );
     if( info < 0 ){
	  Print_Message( ERROR, 3, STRING, "Gain matrix: the ", INT, -info, STRING, "th argument has an illegal value and cannot be scaled." );
	  exit( EXIT_FAILURE );
     } else {
	  Print_Message( SUCCESS, 1, STRING, "PARDISO: Gain matrix successfully calculated." );
     }

}

/* This function is deprecated */
void IGainMatrix_Pardiso_Sp( MatrixVector_t *const IGain, const MatrixVector_Sp_t *const Mass, const MatrixVector_Sp_t *const Damp,
			     const MatrixVector_Sp_t *const Stiff, const Scalars_t Const )
{
     int i;
     int iparm[64];
     void *pt[64];
     int maxfct, mnum, msglvl, error;
     int mtype = 2;                   /* Real symmetric matrix */
     int phase;
     MatrixVector_t IdentMatrix;
     MatrixVector_Sp_t Sp_TempMat;
     double ddum;                       /* Dummy double */
     int idum;                          /* Dummy integer */
     char uplo;                         /* Used in dlascl() */
     int info, lda;                     /* Error handling and Leading dimension for dlascl() */
     double Scalar;    /* Double precision scalar */

     MatrixVector_Create_Sp( Damp->Rows, Damp->Cols, Damp->Num_Nonzero, &Sp_TempMat );

     MatrixVector_Add3Mat_Sp( Stiff, Mass, Damp, Const, &Sp_TempMat );

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
	      Sp_TempMat.Columns, &idum, &IGain->Rows, iparm, &msglvl, &ddum, &ddum, &error);

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
	      Sp_TempMat.Columns, &idum, &IGain->Rows, iparm, &msglvl, &ddum, &ddum, &error);
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
     IdentMatrix = Generate_IdentityMatrix( IGain->Rows, IGain->Cols );

     PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &Sp_TempMat.Rows, Sp_TempMat.Values, Sp_TempMat.RowIndex,
	      Sp_TempMat.Columns, &idum, &IGain->Rows, iparm, &msglvl, IdentMatrix.Array, IGain->Array, &error);

     /* --------------------------------------------------------------------
      * .. Termination and release of memory.
      * -------------------------------------------------------------------- */
     phase = -1;       /* Release internal memory. */
     PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
	      &Sp_TempMat.Rows, &ddum, Sp_TempMat.RowIndex, Sp_TempMat.Columns, &idum, &IGain->Rows,
	      iparm, &msglvl, &ddum, &ddum, &error);

     MatrixVector_Destroy( &IdentMatrix );
     MatrixVector_Destroy_Sp( &Sp_TempMat );

     Print_Message( SUCCESS, 1, STRING, "PARDISO: Matrix Inversion successfully completed." );

     ddum = 1.0; idum = 1;
     lda = Max( 1, IGain->Rows );
     uplo = 'L';   /* The lower part of the matrix will be used (upper part in C) */
     Scalar = Const.Lambda;
     dlascl_( &uplo, &idum, &idum, &ddum, &Scalar, &IGain->Rows, &IGain->Cols, IGain->Array, &lda, &info );
     if( info < 0 ){
	  Print_Message( ERROR, 3, STRING, "Gain matrix: the ", INT, -info, STRING, "th argument has an illegal value and cannot be scaled." );
	  exit( EXIT_FAILURE );
     } else {
	  Print_Message( SUCCESS, 1, STRING, "PARDISO: Gain matrix successfully calculated." );
     }
}
