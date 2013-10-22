#include <stdio.h>          /* For printf(), fprintf() */
#include <stdlib.h>         /* For exit() */

#include "Auxiliary_Math.h"
#include "GainMatrix.h"
#include "MatrixVector.h"
#include "MatrixVector_Sp.h"
#include "Print_Messages.h"

#include "Definitions.h"

#include <mkl_blas.h>
#include <mkl_lapack.h>
#include <mkl_pardiso.h>

void IGainMatrix_Sp( MatrixVector_t *const IGain, const MatrixVector_Sp_t *const Mass, const MatrixVector_Sp_t *const Damp,
		     const MatrixVector_Sp_t *const Stiff, const Scalars_t Const )
{

     char uplo;
     int lda, info;                /* Leading dimension of the gain matrix, LAPACK error handling variable */
     MatrixVector_Sp_t Sp_TempMat; /* Temporal sparse matrix */

     int ione = 1;                 /* Integer of value one */
     HYSL_FLOAT one = 1.0;         /* HYSL_FLOAT precision one for dlascl() parameter cfrom */
     HYSL_FLOAT Scalar;            /* HYSL_FLOAT precision scalar */

     MatrixVector_Create_Sp( Damp->Rows, Damp->Cols, Damp->Num_Nonzero, &Sp_TempMat );
     /* Gain = Const.Alpha*M + Const.Beta*C + Const.Gamma*K */
     MatrixVector_Add3Mat_Sp( Mass, Damp, Stiff, Const, &Sp_TempMat );

     /* Sparse to dense conversion. The Gain matrix will be symmetrical and only the upper part (lower part in
      * FORTRAN) will be referenced.
      */
     MatrixVector_CSR2Dense( &Sp_TempMat, 0, IGain );
     MatrixVector_Destroy_Sp( &Sp_TempMat );

     uplo = 'L';  /* The lower part of the matrix will be used and the upper part will strictly not be
		   * referenced */
     lda = Max( 1, IGain->Rows );

     /* LAPACK: Compute the Cholesky factorization */
     hysl_potrf( &uplo, &IGain->Rows, IGain->Array, &lda, &info );

     if ( info == 0 ){
	  Print_Header( SUCCESS );
	  printf( "Cholesky factorization successfully completed.\n" );
     }
     else if (info < 0){
	  Print_Header( ERROR );
	  fprintf( stderr, "Cholesky factorization: the %d-th argument has an illegal value.\n", -info );
	  exit( EXIT_FAILURE );
     } else if (info > 0){
	  Print_Header( ERROR );
	  fprintf( stderr, "Cholesky factorization: the leading minor of order %d is not positive definite,", info );
	  fprintf( stderr, " and the factorization could not be completed.\n" );
	  exit( EXIT_FAILURE );
     }

     /* LAPACK: Compute the inverse of the Gain matrix using the Cholesky factorization computed by
      * dpotrf_() */
     hysl_potri( &uplo, &IGain->Rows, IGain->Array, &lda, &info );

     if ( info == 0 ){
	  Print_Header( SUCCESS );
	  printf( "Matrix Inversion successfully completed.\n" );
     } else if (info < 0){
	  Print_Header( ERROR );
	  fprintf( stderr, "Matrix inversion: the %d-th argument has an illegal value.\n", -info );
	  exit( EXIT_FAILURE );
     } else if (info > 0){
	  Print_Header( ERROR );
	  fprintf( stderr, "Matrix Inversion: the (%d,%d) element of the factor U or L is zero,", info, info );
	  fprintf( stderr, "and the inverse could not be computed.\n" );
	  exit( EXIT_FAILURE );
     }

     Scalar = Const.Lambda;
     hysl_lascl( &uplo, &ione, &ione, &one, &Scalar, &IGain->Rows, &IGain->Cols, IGain->Array, &lda, &info );

     if( info < 0 ){
	  Print_Header( ERROR );
	  fprintf( stderr, "Gain matrix: the %d-th argument has an illegal value and cannot be scaled.\n", -info );
	  exit( EXIT_FAILURE );
     } else {
	  Print_Header( SUCCESS );
	  printf( "Gain matrix successfully calculated.\n" );
     }
}

void IGainMatrix_Sp_PS( MatrixVector_t *const IGain, const MatrixVector_Sp_t *const Mass, const MatrixVector_Sp_t *const Damp,
		    const MatrixVector_Sp_t *const Stiff, const Scalars_t Const )
{

     char uplo;
     int info;                     /* LAPACK error handling variable */
     MatrixVector_Sp_t Sp_TempMat; /* Temporal sparse matrix */
     HYSL_FLOAT Scalar;                /* HYSL_FLOAT precision scalar */
     int incx, Length;

     uplo = 'L';  /* The lower part of the matrix will be used and the upper part will strictly not be
		   * referenced */
     incx = 1;

     MatrixVector_Create_Sp( Damp->Rows, Damp->Cols, Damp->Num_Nonzero, &Sp_TempMat );
     /* Gain = Const.Alpha*M + Const.Beta*C + Const.Gamma*K */
     MatrixVector_Add3Mat_Sp( Mass, Damp, Stiff, Const, &Sp_TempMat );

     /* Sparse to dense conversion. The Gain matrix will be symmetrical and only the upper part (lower part in
      * FORTRAN) will be referenced. */
     MatrixVector_CSR2Packed( &Sp_TempMat, IGain );
     MatrixVector_Destroy_Sp( &Sp_TempMat );

     /* LAPACK: Compute the Cholesky factorization */
     hysl_pptrf( &uplo, &IGain->Rows, IGain->Array, &info );

     if ( info == 0 ){
	  Print_Header( SUCCESS );
	  printf( "Cholesky factorization successfully completed.\n" );
     }
     else if (info < 0){
	  Print_Header( ERROR );
	  fprintf( stderr, "Cholesky factorization: the %d-th argument has an illegal value.\n", -info );
	  exit( EXIT_FAILURE );
     } else if (info > 0){
	  Print_Header( ERROR );
	  fprintf( stderr, "Cholesky factorization: the leading minor of order %d is not positive definite,", info );
	  fprintf( stderr, " and the factorization could not be completed.\n" );
	  exit( EXIT_FAILURE );
     }

     /* LAPACK: Compute the inverse of the Gain matrix using the Cholesky factorization computed by
      * dpptrf_() */
     hysl_pptri( &uplo, &IGain->Rows, IGain->Array, &info );

     if ( info == 0 ){
	  Print_Header( SUCCESS );
	  printf( "Matrix Inversion successfully completed.\n" );
     } else if (info < 0){
	  Print_Header( ERROR );
	  fprintf( stderr, "Matrix inversion: the %d-th argument has an illegal value.\n", -info );
	  exit( EXIT_FAILURE );
     } else if (info > 0){
	  Print_Header( ERROR );
	  fprintf( stderr, "Matrix Inversion: the (%d,%d) element of the factor U or L is zero,", info, info );
	  fprintf( stderr, "and the inverse could not be computed.\n" );
	  exit( EXIT_FAILURE );
     }

     Length = (IGain->Rows*IGain->Cols + IGain->Rows)/2;
     Scalar = Const.Lambda;
     hysl_scal( &Length, &Scalar, IGain->Array, &incx );
     
     Print_Header( SUCCESS );
     printf( "Gain matrix successfully calculated.\n" );
}

void IGainMatrix_Float2Double_Sp( MatrixVector_t *const IGain, const MatrixVector_Sp_t *const Mass,
				  const MatrixVector_Sp_t *const Damp, const MatrixVector_Sp_t *const Stiff,
				  const Scalars_t Const )
{

     char uplo;
     int lda, info;                 /* Leading dimension of the gain matrix, LAPACK error handling variable */
     MatrixVector_Sp_t Sp_TempMat;  /* Temporal sparse matrix */

     int ione = 1;                  /* Integer of value one */
     double one = 1.0;              /* double precision one for dlascl() parameter cfrom */
     double Scalar;                 /* double precision scalar */

     double *TempMatDouble = NULL;  /* Temporal double matrix for matrix inversion purposes. */
     size_t i;

     uplo = 'L';  /* The lower part of the matrix will be used and the upper part will strictly not be
		   * referenced */
     lda = Max( 1, IGain->Rows );

     MatrixVector_Create_Sp( Damp->Rows, Damp->Cols, Damp->Num_Nonzero, &Sp_TempMat );

     /* Gain = Const.Alpha*M + Const.Beta*C + Const.Gamma*K */
     MatrixVector_Add3Mat_Sp( Mass, Damp, Stiff, Const, &Sp_TempMat );

     /* Sparse to dense conversion. The Gain matrix will be symmetrical and only the upper part (lower part in
      * FORTRAN) will be referenced.
      */
     MatrixVector_CSR2Dense( &Sp_TempMat, 0, IGain );
     MatrixVector_Destroy_Sp( &Sp_TempMat );

     TempMatDouble = (double *) calloc ( (size_t) IGain->Rows*(size_t) IGain->Cols, sizeof(double) );

     if( TempMatDouble == NULL ){
	  Print_Header( ERROR );
	  fprintf( stderr, "IGainMatrix_Float2Double(): Out of memory.\n" );
	  exit( EXIT_FAILURE );
     }

#pragma omp parallel for
     for( i = 0; i < (size_t) IGain->Rows*(size_t) IGain->Cols; i++ ){
	  TempMatDouble[i] = (double) IGain->Array[i];
     }

     /* LAPACK: Compute the Cholesky factorization */
     dpotrf_( &uplo, &IGain->Rows, TempMatDouble, &lda, &info );

     if ( info == 0 ){
	  Print_Header( SUCCESS );
	  printf( "Cholesky factorization successfully completed.\n" );
     }
     else if (info < 0){
	  Print_Header( ERROR );
	  fprintf( stderr, "Cholesky factorization: the %d-th argument has an illegal value.\n", -info );
	  exit( EXIT_FAILURE );
     } else if (info > 0){
	  Print_Header( ERROR );
	  fprintf( stderr, "Cholesky factorization: the leading minor of order %d is not positive definite,", info );
	  fprintf( stderr, " and the factorization could not be completed.\n" );
	  exit( EXIT_FAILURE );
     }

     /* LAPACK: Compute the inverse of the Gain matrix using the Cholesky factorization computed by
      * dpotrf_() */
     dpotri_( &uplo, &IGain->Rows, TempMatDouble, &lda, &info );

     if ( info == 0 ){
	  Print_Header( SUCCESS );
	  printf( "Matrix Inversion successfully completed.\n" );
     } else if (info < 0){
	  Print_Header( ERROR );
	  fprintf( stderr, "Matrix inversion: the %d-th argument has an illegal value.\n", -info );
	  exit( EXIT_FAILURE );
     } else if (info > 0){
	  Print_Header( ERROR );
	  fprintf( stderr, "Matrix Inversion: the (%d,%d) element of the factor U or L is zero,", info, info );
	  fprintf( stderr, "and the inverse could not be computed.\n" );
	  exit( EXIT_FAILURE );
     }

     Scalar = (double) Const.Lambda;
     dlascl_( &uplo, &ione, &ione, &one, &Scalar, &IGain->Rows, &IGain->Cols, TempMatDouble, &lda, &info );

     if( info < 0 ){
	  Print_Header( ERROR );
	  fprintf( stderr, "Gain matrix: the %d-th argument has an illegal value and cannot be scaled.\n", -info );
	  exit( EXIT_FAILURE );
     } else {
	  Print_Header( SUCCESS );
	  printf( "Gain matrix successfully calculated.\n" );
     }


#pragma omp parallel for
     for( i = 0; i < (size_t) IGain->Rows*(size_t) IGain->Cols; i++ ){
	  IGain->Array[i] = (float) TempMatDouble[i];
     }

     /* Free the temporal HYSL_FLOAT matrix */
     free( TempMatDouble );
}

void IGainMatrix_Float2Double_Sp_PS( MatrixVector_t *const IGain, const MatrixVector_Sp_t *const Mass,
				     const MatrixVector_Sp_t *const Damp, const MatrixVector_Sp_t *const Stiff,
				     const Scalars_t Const )
{

     char uplo;
     int info;                      /* LAPACK error handling variable */
     MatrixVector_Sp_t Sp_TempMat;  /* Temporal sparse matrix */
     double Scalar;                 /* double precision scalar */
     int incx;

     double *TempMatDouble = NULL;  /* Temporal double matrix for matrix inversion purposes. */
     size_t Length, i;

     uplo = 'L';  /* The lower part of the matrix will be used and the upper part will strictly not be
		   * referenced */
     incx = 1;

     MatrixVector_Create_Sp( Damp->Rows, Damp->Cols, Damp->Num_Nonzero, &Sp_TempMat );
     /* Gain = Const.Alpha*M + Const.Beta*C + Const.Gamma*K */
     MatrixVector_Add3Mat_Sp( Mass, Damp, Stiff, Const, &Sp_TempMat );

     /* Sparse to dense conversion. The Gain matrix will be symmetrical and only the upper part (lower part in
      * FORTRAN) will be referenced. */
     MatrixVector_CSR2Packed( &Sp_TempMat, IGain );
     MatrixVector_Destroy_Sp( &Sp_TempMat );

     Length = ((size_t) IGain->Rows*(size_t) IGain->Cols + (size_t)IGain->Rows)/2;
     TempMatDouble = (double *) calloc ( Length, sizeof(double) );
     if( TempMatDouble == NULL ){
	  Print_Header( ERROR );
	  fprintf( stderr, "IGainMatrix_Float2Double(): Out of memory.\n" );
	  exit( EXIT_FAILURE );
     }

#pragma omp parallel for
     for( i = 0; i < Length; i++ ){
	  TempMatDouble[i] = (double) IGain->Array[i];
     }

     /* LAPACK: Compute the Cholesky factorization */
     dpptrf_( &uplo, &IGain->Rows, TempMatDouble, &info );

     if ( info == 0 ){
	  Print_Header( SUCCESS );
	  printf( "Cholesky factorization successfully completed.\n" );
     }
     else if (info < 0){
	  Print_Header( ERROR );
	  fprintf( stderr, "Cholesky factorization: the %d-th argument has an illegal value.\n", -info );
	  exit( EXIT_FAILURE );
     } else if (info > 0){
	  Print_Header( ERROR );
	  fprintf( stderr, "Cholesky factorization: the leading minor of order %d is not positive definite,", info );
	  fprintf( stderr, " and the factorization could not be completed.\n" );
	  exit( EXIT_FAILURE );
     }

     /* LAPACK: Compute the inverse of the Gain matrix using the Cholesky factorization computed by
      * dpptrf_() */
     dpptri( &uplo, &IGain->Rows, TempMatDouble, &info );

     if ( info == 0 ){
	  Print_Header( SUCCESS );
	  printf( "Matrix Inversion successfully completed.\n" );
     } else if (info < 0){
	  Print_Header( ERROR );
	  fprintf( stderr, "Matrix inversion: the %d-th argument has an illegal value.\n", -info );
	  exit( EXIT_FAILURE );
     } else if (info > 0){
	  Print_Header( ERROR );
	  fprintf( stderr, "Matrix Inversion: the (%d,%d) element of the factor U or L is zero,", info, info );
	  fprintf( stderr, "and the inverse could not be computed.\n" );
	  exit( EXIT_FAILURE );
     }

     Scalar = Const.Lambda;
     dscal( (int *) &Length, &Scalar, TempMatDouble, &incx );

#pragma omp parallel for
     for( i = 0; i < Length; i++ ){
	  IGain->Array[i] = (float) TempMatDouble[i];
     }

     /* Free the temporal HYSL_FLOAT matrix */
     free( TempMatDouble );
     
     Print_Header( SUCCESS );
     printf( "Gain matrix successfully calculated.\n" );
}

/* This function is deprecated */
void IGainMatrix_Pardiso( MatrixVector_t *const IGain, const MatrixVector_t *const Mass, const MatrixVector_t *const Damp,
			 const MatrixVector_t *const Stiff, const Scalars_t Const )
{

     int i;
     int iparm[64];
     void *pt[64];
     int maxfct, mnum, msglvl, error;
     int mtype = 2;                     /* Real symmetric matrix */
     int phase;
     MatrixVector_t IdentMatrix;
     MatrixVector_t TempMat;
     MatrixVector_Sp_t Sp_TempMat;
     HYSL_FLOAT ddum;                       /* Dummy HYSL_FLOAT */
     int idum;                          /* Dummy integer */
     char uplo;                         /* Used in dlascl() */
     int info, lda;                     /* Error handling and Leading dimension for dlascl() */
     HYSL_FLOAT Scalar;                     /* HYSL_FLOAT precision scalar */

     MatrixVector_Create( IGain->Rows, IGain->Cols, &TempMat );

     MatrixVector_Add3Mat( Mass, Damp, Stiff, Const, &TempMat );

     MatrixVector_Dense2CSR( &TempMat, 0, &Sp_TempMat );   /* Transform into CSR format */
     MatrixVector_Destroy( &TempMat );                     /* Destroy the dense matrix */

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
     iparm[12] = 1;	/* Maximum weighted matching algorithm is switched-off (default for symmetric).  Try
			 * iparm[12] = 1 in case of inappropriate accuracy */
     iparm[13] = 0;	/* Output: Number of perturbed pivots */
     iparm[17] = -1;	/* Output: Number of nonzeros in the factor LU */
     iparm[18] = -1;	/* Output: Mflops for LU factorization */
     iparm[19] = 0;	/* Output: Numbers of CG Iterations */
#if _FLOAT_
     iparm[27] = 1;     /* Input/output matrices are single precision */
#else
     iparm[27] = 0;     /* Input/output matrices are double precision */
#endif
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
	      Sp_TempMat.Columns, &idum, &IGain->Rows, iparm, &msglvl, &ddum, &ddum, &error);

     if (error != 0){
	  Print_Header( ERROR );
	  fprintf( stderr, "PARDISO: Error during symbolic factorization %d.\n" , error );
	  exit( EXIT_FAILURE );
     }

     Print_Header( SUCCESS );
     printf( "PARDISO: Reordering completed ...\n" );
     Print_Header( INFO );
     printf( "PARDISO: Number of nonzeros in factors = %d.\n", iparm[17] ); 
     Print_Header( INFO );
     printf( "PARDISO: Number of factorization MFLOPS = %d", iparm[18] ); 

     /* -------------------------------------------------------------------- */
     /* .. Numerical factorization. */
     /* -------------------------------------------------------------------- */
     phase = 22;
     
     PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &Sp_TempMat.Rows, Sp_TempMat.Values, Sp_TempMat.RowIndex,
	      Sp_TempMat.Columns, &idum, &IGain->Rows, iparm, &msglvl, &ddum, &ddum, &error);

     if (error != 0){
	  Print_Header( ERROR );
	  fprintf( stderr, "PARDISO: Error during numerical factorization %d.\n", error );
	  exit (2);
     }

     Print_Header( SUCCESS );
     printf( "PARDISO: Factorization completed ...\n" );

     /* --------------------------------------------------------------------
      * .. Back substitution and iterative refinement.
      * -------------------------------------------------------------------- */
     phase = 33;
     iparm[7] = 10; /* Max numbers of iterative refinement steps. */

     /* Set right hand side to be the identity matrix. */     
     IdentMatrix = Generate_IdentityMatrix( IGain->Rows, IGain->Cols );

     PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &Sp_TempMat.Rows, Sp_TempMat.Values, Sp_TempMat.RowIndex,
	      Sp_TempMat.Columns, &idum, &IGain->Rows, iparm, &msglvl, IdentMatrix.Array, IGain->Array, &error);

     /* --------------------------------------------------------------------
      *	.. Termination and release of memory.
      *	-------------------------------------------------------------------- */
     phase = -1; /* Release internal memory. */
     PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &Sp_TempMat.Rows, &ddum, Sp_TempMat.RowIndex, Sp_TempMat.Columns,
	      &idum, &IGain->Rows, iparm, &msglvl, &ddum, &ddum, &error);

     MatrixVector_Destroy( &IdentMatrix );
     MatrixVector_Destroy_Sp( &Sp_TempMat );

     Print_Header( SUCCESS );
     printf( "PARDISO: Matrix Inversion successfully completed.\n" );

     ddum = 1.0; idum = 1;
     lda = Max( 1, IGain->Rows );
     uplo = 'L';   /* The lower part of the matrix will be used (upper part in C) */

     Scalar = Const.Lambda;
     hysl_lascl( &uplo, &idum, &idum, &ddum, &Scalar, &IGain->Rows, &IGain->Cols, IGain->Array, &lda, &info );

     if( info < 0 ){
	  Print_Header( ERROR );
	  fprintf( stderr, "Gain matrix: the %d-th argument has an illegal value and cannot be scaled.\n", -info );
	  exit( EXIT_FAILURE );
     } else {
	  Print_Header( SUCCESS );
	  printf( "PARDISO: Gain matrix successfully calculated.\n" );
     }

}

/* This function is deprecated */
void IGainMatrix_Pardiso_Sp( MatrixVector_t *const IIGain, const MatrixVector_Sp_t *const Mass, const MatrixVector_Sp_t *const Damp,
			     const MatrixVector_Sp_t *const Stiff, const Scalars_t Const )
{
     int i;
     int iparm[64];
     void *pt[64];
     int maxfct, mnum, msglvl, error;
     int mtype = 2;                     /* Real symmetric matrix */
     int phase;
     MatrixVector_t IdentMatrix;
     MatrixVector_Sp_t Sp_TempMat;
     HYSL_FLOAT ddum;                       /* Dummy HYSL_FLOAT */
     int idum;                          /* Dummy integer */
     char uplo;                         /* Used in dlascl() */
     int info, lda;                     /* Error handling and Leading dimension for dlascl() */
     HYSL_FLOAT Scalar;                     /* HYSL_FLOAT precision scalar */

     MatrixVector_Create_Sp( Damp->Rows, Damp->Cols, Damp->Num_Nonzero, &Sp_TempMat );

     MatrixVector_Add3Mat_Sp( Mass, Damp, Stiff, Const, &Sp_TempMat );

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
     iparm[12] = 1;	/* Maximum weighted matching algorithm is switched-off (default for symmetric).  Try
			 *  iparm[12] = 1 in case of inappropriate accuracy */
     iparm[13] = 0;	/* Output: Number of perturbed pivots */
     iparm[17] = -1;	/* Output: Number of nonzeros in the factor LU */
     iparm[18] = -1;	/* Output: Mflops for LU factorization */
     iparm[19] = 0;	/* Output: Numbers of CG Iterations */
#if _FLOAT_
     iparm[27] = 1;     /* Input/output matrices are single precision */
#else
     iparm[27] = 0;     /* Input/output matrices are double precision */
#endif
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
	      Sp_TempMat.Columns, &idum, &IIGain->Rows, iparm, &msglvl, &ddum, &ddum, &error);

     if (error != 0){
	  Print_Header( ERROR );
	  fprintf( stderr, "PARDISO: Error during symbolic factorization %d.\n", error );
	  exit( EXIT_FAILURE );
     }

     Print_Header( SUCCESS );
     printf( "PARDISO: Reordering completed ...\n" );
     Print_Header( INFO );
     printf( "PARDISO: Number of nonzeros in factors = %d.\n", iparm[17] ); 
     Print_Header( INFO );
     printf( "PARDISO: Number of factorization MFLOPS = %d", iparm[18] );

     /* --------------------------------------------------------------------
      * .. Numerical factorization.
      * -------------------------------------------------------------------- */
     phase = 22;

     PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &Sp_TempMat.Rows, Sp_TempMat.Values, Sp_TempMat.RowIndex,
	      Sp_TempMat.Columns, &idum, &IIGain->Rows, iparm, &msglvl, &ddum, &ddum, &error);

     if (error != 0){
	  Print_Header( ERROR );
	  fprintf( stderr, "PARDISO: Error during numerical factorization %d.\n", error );
	  exit (2);
     }

     Print_Header( SUCCESS );
     printf( "PARDISO: Factorization completed ...\n" );

     /* --------------------------------------------------------------------
      * .. Back substitution and iterative refinement.
      * -------------------------------------------------------------------- */
     phase = 33;
     iparm[7] = 10;	/* Max numbers of iterative refinement steps. */

     /* Set right hand side to be the identity matrix. */
     IdentMatrix = Generate_IdentityMatrix( IIGain->Rows, IIGain->Cols );

     PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &Sp_TempMat.Rows, Sp_TempMat.Values, Sp_TempMat.RowIndex,
	      Sp_TempMat.Columns, &idum, &IIGain->Rows, iparm, &msglvl, IdentMatrix.Array, IIGain->Array, &error);

     /* --------------------------------------------------------------------
      * .. Termination and release of memory.
      * -------------------------------------------------------------------- */
     phase = -1;       /* Release internal memory. */
     PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
	      &Sp_TempMat.Rows, &ddum, Sp_TempMat.RowIndex, Sp_TempMat.Columns, &idum, &IIGain->Rows,
	      iparm, &msglvl, &ddum, &ddum, &error);

     MatrixVector_Destroy( &IdentMatrix );
     MatrixVector_Destroy_Sp( &Sp_TempMat );

     Print_Header( SUCCESS );
     printf( "PARDISO: Matrix Inversion successfully completed.\n" );

     ddum = 1.0; idum = 1;
     lda = Max( 1, IIGain->Rows );
     uplo = 'L';   /* The lower part of the matrix will be used (upper part in C) */
     Scalar = Const.Lambda;
     hysl_lascl( &uplo, &idum, &idum, &ddum, &Scalar, &IIGain->Rows, &IIGain->Cols, IIGain->Array, &lda, &info );

     if( info < 0 ){
	  Print_Header( ERROR );
	  fprintf( stderr, "Gain matrix: the %d-th argument has an illegal value and cannot be scaled.\n", -info );
	  exit( EXIT_FAILURE );
     } else {
	  Print_Header( SUCCESS );
	  printf( "PARDISO: Gain matrix successfully calculated.\n" );
     }
}
