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
#include "Substructure.h"
#include "Colors.h"

/* Matrix market format */
#include "mmio.h"

#if _SPARSE_
#include "mkl.h"
#include "mkl_pardiso.h"
#endif

void InitConstants( AlgConst *const InitConst, const char *FileName )
{

     ConfFile *Config;
     
     Config = ConfFile_Create( 70 );

     ConfFile_ReadFile( Config, FileName );

     /* Use Relative or absolute values */
     (*InitConst).Use_Absolute_Values = ConfFile_GetInt( Config, "General:Use_Absolute_Values" );
     if ( InitConst->Use_Absolute_Values != 0 && InitConst->Use_Absolute_Values != 1 ){
	  PrintErrorAndExit( "Invalid option for Use_Absolute_Values" );
     }

     InitConst->Read_Sparse = ConfFile_GetInt( Config, "General:Read_Sparse" );
     if ( InitConst->Read_Sparse != 0 && InitConst->Read_Sparse != 1 ){
	  PrintErrorAndExit( "Invalid option for Read_Sparse" );
     }

     InitConst->Use_Sparse = ConfFile_GetInt( Config, "General:Use_Sparse" );
     if ( InitConst->Use_Sparse != 0 && InitConst->Use_Sparse != 1 ){
	  PrintErrorAndExit( "Invalid option for Use_Sparse" );
     }

     InitConst->Use_Pardiso = ConfFile_GetInt( Config, "General:Use_Pardiso" );
     if ( InitConst->Use_Pardiso != 0 && InitConst->Use_Pardiso != 1 ){
	  PrintErrorAndExit( "Invalid option for Use_Pardiso" );
     }

     InitConst->Read_LVector = ConfFile_GetInt( Config, "General:Read_LVector" );
     if ( InitConst->Read_LVector != 0 && InitConst->Read_LVector != 1 ){
	  PrintErrorAndExit( "Invalid option for Read_LVector" );
     }

     /* Order of the matrices */
     (*InitConst).Order = ConfFile_GetInt( Config, "General:Order" );
     if ( InitConst->Order <= 0 ){
	  PrintErrorAndExit( "Invalid option for the order of the matrices" );
     }

     if( !InitConst->Read_LVector ){
	  InitConst->ExcitedDOF = Get_Excited_DOF( Config, "General:Excited_DOF" );
     } else {
	  InitConst->ExcitedDOF = NULL;
     }

     /* Number of steps and Time step */
     (*InitConst).Nstep = (unsigned int) ConfFile_GetInt( Config, "General:Num_Steps" );
     if ( InitConst->Nstep <= 0 ){
	  PrintErrorAndExit( "Invalid number of steps" );
     }

     (*InitConst).Delta_t = ConfFile_GetDouble( Config, "General:Delta" );
     if ( InitConst->Delta_t <= 0.0 ){
	  PrintErrorAndExit( "Invalid time step" );
     }

     /* Rayleigh values */
     (*InitConst).Rayleigh.Alpha = ConfFile_GetDouble( Config, "Rayleigh:Alpha" );
     (*InitConst).Rayleigh.Beta = ConfFile_GetDouble( Config, "Rayleigh:Beta" );

     /* Newmark integration constants */
     (*InitConst).Newmark.Gamma = ConfFile_GetDouble( Config, "Newmark:Gamma" );
     (*InitConst).Newmark.Beta = ConfFile_GetDouble( Config, "Newmark:Beta" );

     /* PID Constants */
     (*InitConst).PID.P = ConfFile_GetDouble( Config, "PID:P" );
     (*InitConst).PID.I = ConfFile_GetDouble( Config, "PID:I" );
     (*InitConst).PID.D = ConfFile_GetDouble( Config, "PID:D" );

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
     (*InitConst).FileM = strdup( ConfFile_GetString( Config, "FileNames:Mass_Matrix" ) );
     (*InitConst).FileK = strdup( ConfFile_GetString( Config, "FileNames:Stiffness_Matrix" ) );
     (*InitConst).FileC = strdup( ConfFile_GetString( Config, "FileNames:Damping_Matrix" ) );
     if( InitConst->Read_LVector ){
	  (*InitConst).FileLV = strdup( ConfFile_GetString( Config, "FileNames:Load_Vector" ) );

     } else {
	  InitConst->FileLV = NULL;
     }
     (*InitConst).FileCNodes = strdup( ConfFile_GetString( Config, "FileNames:Coupling_Nodes" ) );
     (*InitConst).FileData = strdup( ConfFile_GetString( Config, "FileNames:Ground_Motion" ) );
     (*InitConst).FileOutput = strdup( ConfFile_GetString( Config, "FileNames:OutputFile" ) );

     /* Get the communication protocol to be used */
     GetNetworkInformation( &InitConst->Remote, Config );

     /* Read the information regarding the numerical sub-structures */

     /* Number of substructures */
     (*InitConst).OrderSub = ConfFile_GetInt( Config, "Substructure:Order" );
     if ( InitConst->OrderSub < 0 ){
	  PrintErrorAndExit( "Invalid option for the number of sub-structuresr of the matrices" );
     }
     
     /* Number of substructures */
     InitConst->NSubstep = (unsigned int) ConfFile_GetInt( Config, "Substructure:Num_Substeps" );

     InitConst->DeltaT_Sub = InitConst->Delta_t/(double) InitConst->NSubstep;

     ConfFile_Free( Config );
}

int* Get_Excited_DOF( const ConfFile *const Config, const char *Expression )
{
     unsigned int i, j;
     int *DOF_Table;  
     char *FullString;
     char Temp[1];

     FullString = strdup( ConfFile_GetString( Config, Expression ) );

     /* The first position contains the number of degrees of Freedom per node present in the
      * structure */
     strncpy( Temp, &FullString[0], (size_t) 1 );

     DOF_Table = (int *) calloc( (size_t) (atoi(Temp)+1), sizeof(int) );
     DOF_Table[0] = atoi( Temp );

     j = 1;
     for( i = 1; i < strlen( FullString ); i++ ){
	  if ( FullString[i] != ' ' ){
	       strncpy( Temp, &FullString[i], (size_t) 1 );
	       DOF_Table[j] = atoi( Temp );
	       j = j + 1;
	  }
     }

     free( FullString );
     return DOF_Table;
}

void Delete_InitConstants( AlgConst *const InitConst )
{

     if( InitConst->ExcitedDOF != NULL ){
	  free( InitConst->ExcitedDOF );
     }
     free( InitConst->FileM );
     free( InitConst->FileK );

     if( InitConst->FileC != NULL ){
	  free( InitConst->FileC );
     }
     if( InitConst->FileLV != NULL ){
	  free( InitConst->FileLV );
     }
     free( InitConst->FileCNodes );
     free( InitConst->FileData );
     free( InitConst->FileOutput );

     Delete_NetworkInformation( &InitConst->Remote );
}

void Read_Coupling_Nodes( Coupling_Node *const CNodes, const int OrderSub, const double DeltaTSub, const char *Filename )
{
     FILE *InFile;
     int i, j;
     int itemp;
     double *ftemp;
     char *ctemp;
     EXP_Sub *Experimental;

     InFile = fopen( Filename, "r" );

     if ( InFile != NULL ){
	  /* The first value should be the number of Coupling nodes */
	  fscanf( InFile, "%i", &CNodes->Order );

	  if( CNodes->Order != OrderSub ){
	       fclose( InFile );
	       PrintErrorAndExit( "Invalid Number of Substructures" );
	  }
	  
	  /* Allocate the necessary memory */
	  CNodes->Array = (int *) calloc( (size_t) CNodes->Order, sizeof(int) );
	  CNodes->Sub = (Substructure *) malloc( (size_t) CNodes->Order*sizeof(Substructure) );
	  CNodes->u0c0 = (double *) calloc( (size_t) CNodes->Order, sizeof(double) );
	  /* Read the contents of the file */
	  for( i = 0; i < CNodes->Order; i++ ){
	       fscanf( InFile, "%i %i", &CNodes->Array[i], &CNodes->Sub[i].Type );
	       
	       switch (CNodes->Sub[i].Type) {
	       case USE_ADWIN:
		    CNodes->Sub[i].SimStruct = (void *) malloc( sizeof(EXP_Sub) );
		    Experimental = CNodes->Sub[i].SimStruct;
		    /* Dynamic string input. Reads everything between " " */
		    fgets( Experimental->Description, MAX_DESCRIPTION, InFile );

		    if( Experimental->Description[strlen(Experimental->Description) - 1] != '\n' && !feof(InFile) ){
			 fprintf( stderr, "The desciption of the substructure %i, was too long. Maximum 80 characters\n", i + 1 );
			 exit( EXIT_FAILURE );
		    }

		    break;
	       case USE_EXACT:
		    fscanf( InFile, "%i", &itemp );
		    if ( itemp != UHYDE_NUMPARAM_INIT ){
			 fprintf( stderr, "Wrong number of parameters for the substructue number %i of type TMD.\n", i );
			 fprintf( stderr, "The number of init parameters should be %i\n", EXACT_NUMPARAM_INIT );
			 exit( EXIT_FAILURE );
		    } else {
			 printf( "Simulating the sub-structure in the coupling node %d as an exact integration method.\n", CNodes->Array[i] );
			 CNodes->Sub[i].SimStruct = (void *) malloc( sizeof(TMD_Sim) );
			 ftemp = NULL;
			 ftemp = (double *) calloc( (size_t) EXACT_NUMPARAM_INIT, sizeof( double ) );

			 for( j = 0; j < EXACT_NUMPARAM_INIT; j++ ){
			      fscanf( InFile, "%lf", &ftemp[j] );
			 }
			 
			 ExactSolution_Init( ftemp[0], ftemp[1], ftemp[2], DeltaTSub, (TMD_Sim *) CNodes->Sub[i].SimStruct );
			 free( ftemp );
		    }
		    break;
	       case USE_UHYDE:
		    fscanf( InFile, "%i", &itemp );
		    if ( itemp != UHYDE_NUMPARAM_INIT ){
			 fprintf( stderr, "Wrong number of parameters for the substructue number %i of type UHYDE.\n", i );
			 fprintf( stderr, "The number of init parameters should be %i\n", UHYDE_NUMPARAM_INIT );
			 exit( EXIT_FAILURE );
		    } else {
			 printf( "Simulating the sub-structure in the coupling node %d as a UHYDE-fbr device.\n", CNodes->Array[i] );
			 CNodes->Sub[i].SimStruct = (void *) malloc( sizeof(UHYDE_Sim) );
			 ftemp = NULL;
			 ftemp = (double *) calloc( (size_t) UHYDE_NUMPARAM_INIT, sizeof( double ) );

			 for( j = 0; j < UHYDE_NUMPARAM_INIT; j++ ){
			      fscanf( InFile, "%lf", &ftemp[j] );
			 }
			 
			 Simulate_UHYDE_1D_Init( ftemp[0], ftemp[1], ftemp[2], (UHYDE_Sim *) CNodes->Sub[i].SimStruct );
			 free( ftemp );
		    }
		    break;
	       case USE_MEASURED:
		    CNodes->Sub[i].SimStruct = (void *) malloc( sizeof(EXP_Sub) );
		    Experimental = CNodes->Sub[i].SimStruct;
		    /* Dynamic string input. Reads everything between " " */
		    fgets( Experimental->Description, MAX_DESCRIPTION, InFile );

		    if( Experimental->Description[strlen(Experimental->Description) - 1] != '\n' && !feof(InFile) ){
			 fprintf( stderr, "The desciption of the substructure %i, was too long. Maximum 80 characters\n", i + 1 );
			 exit( EXIT_FAILURE );
		    }
	       }
	  }
	  /* Close the file */
	  fclose( InFile );
     } else ErrorFileAndExit( "It is not possible to read data because it was not possible to open: ", Filename );

}

void Delete_CouplingNodes( Coupling_Node *CNodes )
{
     int i;

     for( i = 0; i < CNodes->Order; i++ ){
	  free( CNodes->Sub[i].SimStruct );
     }

     free( CNodes->Sub );
     free( CNodes->Array );
     free( CNodes->u0c0 );
}

void CalculateMatrixC( const MatrixVector *const Mass, const MatrixVector *const Stif, MatrixVector *const Damp, const RayleighConst *const Rayleigh )
{
     char uplo;
     int ione;
     int incx, incy;
     int info, lda, ldb;
     double done;

     int Rows, Cols, Length;
     double alpha, beta;

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
     dlacpy_( &uplo, &Rows, &Cols, (*Mass).Array, &lda, (*Damp).Array, &ldb );

     /* LAPACK: C = Rayleigh.alpha*M = Rayleigh.alpha*C */
     dlascl_( &uplo, &ione, &ione, &done, &alpha, &(*Damp).Rows, &(*Damp).Cols, (*Damp).Array, &lda, &info );

     if ( info < 0 ){
	  LAPACKPErrorAndExit( "dlascl: The ", -info, "-th argument had an illegal value" );
     }

     /* BLAS: C = alpha*M + beta*K = C + beta*K. Only half of the matrix is calculated */
     for ( i = 0; i < (*Damp).Rows; i++ ){
	  Length = (*Damp).Rows - i;
	  daxpy_( &Length, &beta, &(*Stif).Array[i*(*Stif).Rows + i], &incx, &(*Damp).Array[i*(*Damp).Rows +i], &incy);
     }

     PrintSuccess( "Damping matrix successfully calculated.\n" );
}

#if _SPARSE_
void CalculateMatrixC_Sp( const Sp_MatrixVector *const Mass, const Sp_MatrixVector *const Stif, Sp_MatrixVector *const Damp, const RayleighConst *const Rayleigh )
{
     Sp_MatrixVector Temp;
     int Length, i;
     int incx, incy;
     double alpha, beta;
     char trans;
     int job, sort, info;

     alpha = Rayleigh->Alpha;
     beta = Rayleigh->Beta;

     Init_MatrixVector_Sp( &Temp, Mass->Rows, Mass->Cols, Mass->Num_Nonzero );

     incx = 1; incy = 1;
     Length = Temp.Num_Nonzero;
     dcopy_( &Length, Mass->Values, &incx, Temp.Values, &incy );

#pragma omp parallel for
     for (i = 0; i < Length; i++ ){
	  Temp.Columns[i] = Mass->Columns[i];
     }

     /* Scal the Values array */
     dscal_( &Length, &alpha, Temp.Values, &incx );

     /* Copy the RowIndex array */
     Length = Temp.Rows + 1;
#pragma omp parallel for
     for (i = 0; i < Length; i++ ){
	  Temp.RowIndex[i] = Mass->RowIndex[i];
     }

     trans = 'N';  /* The operation C = Temp + beta*B is performed */
     job = 0;      /* The routine computes the addition */
     sort = 0;     /* The routine does not perform any reordering */
     mkl_dcsradd( &trans, &job, &sort, &Temp.Rows, &Temp.Cols, Temp.Values, Temp.Columns, Temp.RowIndex,
		  &beta, Stif->Values, Stif->Columns, Stif->RowIndex,
		  Damp->Values, Damp->Columns, Damp->RowIndex, &Damp->Num_Nonzero, &info );

     /* Delete the previously allocated sparse matrix */
     Destroy_MatrixVector_Sparse( &Temp );

     if ( info > 0){
	  PrintErrorAndExit( "Number of elements exceeded while calculating the Damping matrix" );
     } else if ( info < 0 ){
	  PrintErrorAndExit( "I do not understand" );
     } else {
	  PrintSuccess( "Damping matrix successfully calculated.\n" );
     }
}
#endif

void CalculateMatrixKeinv( MatrixVector *const Keinv, const MatrixVector *const Mass, const MatrixVector *const Damp, const MatrixVector *const Stif, const Scalars Const )
{

     char uplo;
     int lda, info;

     uplo = 'L';  /* The lower part of the matrix will be used; the upper part will strictly not be referenced */
     lda = Max( 1, (*Keinv).Rows );

     /* Perform Meinv = [M + gamma*Delta_t*C + beta*Delta_t^2*K] */
     Add3Mat( &(*Keinv), &(*Stif), &(*Mass), &(*Damp), Const );

     /* LAPACK: Compute the Cholesky factorization of the symmetric positive definite matrix Meinv */
     dpotrf_( &uplo, &(*Keinv).Rows, Keinv->Array, &lda, &info );

     if ( info == 0 ){
	  PrintSuccess( "Cholesky factorization successfully completed.\n" );
     }
     else if (info < 0){
	  LAPACKPErrorAndExit( "Cholesky factorization: the ", -info, "th argument has an illegal value." );
     } else if (info > 0){
	  LAPACKPErrorAndExit("Cholesky factorization: the leading minor of order ", info, " is not positive definite, and the factorization could not be completed." );
     }

     /* LAPACK: Compute the inverse of Me using the Cholesky factorization computed by pdpotrf_( ) */
     dpotri_( &uplo, &(*Keinv).Rows, Keinv->Array, &lda, &info );

     if ( info == 0 ){
	  PrintSuccess( "Matrix Inversion successfully completed.\n" );
     } else if (info < 0){
	  LAPACKPErrorAndExit( "Matrix Inversion: the ", -info, "th argument has an illegal value." );
     } else if (info > 0){
	  fprintf( stderr, "Matrix Inversion: the (%d,%d) element of the factor U or L is zero, and the inverse could not be computed.\n", info, info );
	  fprintf( stderr, "Exiting program.\n" );
	  exit( EXIT_FAILURE );
     }
}

#if _SPARSE_
void CalculateMatrixKeinv_Pardiso( MatrixVector *const Keinv, const MatrixVector *const Mass, const MatrixVector *const Damp, const MatrixVector *const Stiff, const Scalars Const )
{
     int i;
     int iparm[64];
     void *pt[64];
     int maxfct, mnum, msglvl, error;
     int mtype = 2;                   /* Real symmetric matrix */
     int phase;
     MatrixVector IdentMatrix;
     MatrixVector TempMat;
     Sp_MatrixVector Sp_TempMat;
     double fdum;                       /* Dummy double */
     int idum;                          /* Dummy integer */

     Init_MatrixVector( &TempMat, Keinv->Rows, Keinv->Cols );

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
     iparm[12] = 1;	/* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
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
	      &Sp_TempMat.Rows, Sp_TempMat.Values, Sp_TempMat.RowIndex, Sp_TempMat.Columns, &idum, &Keinv->Rows, iparm, &msglvl, &fdum, &fdum, &error);
     if (error != 0)
     {
	  printf ("\nERROR during symbolic factorization: %d", error);
	  exit (1);
     }
     PrintSuccess( "\nReordering completed ... " );
     printf ("\nNumber of nonzeros in factors = %d", iparm[17]);
     printf ("\nNumber of factorization MFLOPS = %d", iparm[18]);
/* -------------------------------------------------------------------- */
/* .. Numerical factorization. */
/* -------------------------------------------------------------------- */
     phase = 22;
     
     PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
	      &Sp_TempMat.Rows, Sp_TempMat.Values, Sp_TempMat.RowIndex, Sp_TempMat.Columns, &idum, &Keinv->Rows, iparm, &msglvl, &fdum, &fdum, &error);
     if (error != 0)
     {
	  printf ("\nERROR during numerical factorization: %d", error);
	  exit (2);
     }
     PrintSuccess ( "\nFactorization completed ... " );
/* -------------------------------------------------------------------- */
/* .. Back substitution and iterative refinement. */
/* -------------------------------------------------------------------- */
     phase = 33;
     iparm[7] = 10;			/* Max numbers of iterative refinement steps. */

     /* Set right hand side to be the identity matrix. */     
     IdentMatrix = Generate_IdentityMatrix( Keinv->Rows, Keinv->Cols );

     PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
	      &Sp_TempMat.Rows, Sp_TempMat.Values, Sp_TempMat.RowIndex, Sp_TempMat.Columns, &idum, &Keinv->Rows, iparm, &msglvl, IdentMatrix.Array, Keinv->Array, &error);

/* -------------------------------------------------------------------- */
/* .. Termination and release of memory. */
/* -------------------------------------------------------------------- */
     phase = -1;			/* Release internal memory. */
     PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
	      &Sp_TempMat.Rows, &fdum, Sp_TempMat.RowIndex, Sp_TempMat.Columns, &idum, &Keinv->Rows,
	      iparm, &msglvl, &fdum, &fdum, &error);

     Destroy_MatrixVector( &IdentMatrix );
     Destroy_MatrixVector_Sparse( &Sp_TempMat );

     PrintSuccess( "Matrix Inversion successfully completed.\n" );
}

void CalculateMatrixKeinv_Pardiso_Sparse( MatrixVector *const Keinv, const Sp_MatrixVector *const Mass, const Sp_MatrixVector *const Damp, const Sp_MatrixVector *const Stiff, const Scalars Const )
{
     int i;
     int iparm[64];
     void *pt[64];
     int maxfct, mnum, msglvl, error;
     int mtype = 2;                   /* Real symmetric matrix */
     int phase;
     MatrixVector IdentMatrix;
     Sp_MatrixVector Sp_TempMat;
     double fdum;                       /* Dummy double */
     int idum;                          /* Dummy integer */

     Init_MatrixVector_Sp( &Sp_TempMat, Damp->Rows, Damp->Cols, Damp->Num_Nonzero );

     Add3Mat_Sparse( &Sp_TempMat, &(*Stiff), &(*Mass), &(*Damp), Const );
     MatrixVector_To_File_Sparse( &Sp_TempMat, "TempMatSp.txt" );

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
     iparm[12] = 1;	/* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
     iparm[13] = 0;	/* Output: Number of perturbed pivots */
     iparm[17] = -1;	/* Output: Number of nonzeros in the factor LU */
     iparm[18] = -1;	/* Output: Mflops for LU factorization */
     iparm[19] = 0;	/* Output: Numbers of CG Iterations */
     iparm[27] = 0;     /* Input/output matrices are single precision */
     iparm[34] = 0;	/* PARDISO use 1 based indexing for ia and ja arrays */
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
	      &Sp_TempMat.Rows, Sp_TempMat.Values, Sp_TempMat.RowIndex, Sp_TempMat.Columns, &idum, &Keinv->Rows, iparm, &msglvl, &fdum, &fdum, &error);
     if (error != 0)
     {
	  printf ("\nERROR during symbolic factorization: %d", error);
	  exit (1);
     }
     PrintSuccess ("\nReordering completed ... ");
     printf ("\nNumber of nonzeros in factors = %d", iparm[17]);
     printf ("\nNumber of factorization MFLOPS = %d", iparm[18]);
/* -------------------------------------------------------------------- */
/* .. Numerical factorization. */
/* -------------------------------------------------------------------- */
     phase = 22;
     PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
	      &Sp_TempMat.Rows, Sp_TempMat.Values, Sp_TempMat.RowIndex, Sp_TempMat.Columns, &idum, &Keinv->Rows, iparm, &msglvl, &fdum, &fdum, &error);
     if (error != 0)
     {
	  printf ("\nERROR during numerical factorization: %d", error);
	  exit (2);
     }
     PrintSuccess ("\nFactorization completed ... ");
/* -------------------------------------------------------------------- */
/* .. Back substitution and iterative refinement. */
/* -------------------------------------------------------------------- */
     phase = 33;
     iparm[7] = 10;			/* Max numbers of iterative refinement steps. */

     /* Set right hand side to be the identity matrix. */
     IdentMatrix = Generate_IdentityMatrix( Keinv->Rows, Keinv->Cols );

     PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
	      &Sp_TempMat.Rows, Sp_TempMat.Values, Sp_TempMat.RowIndex, Sp_TempMat.Columns, &idum, &Keinv->Rows, iparm, &msglvl, IdentMatrix.Array, Keinv->Array, &error);

/* -------------------------------------------------------------------- */
/* .. Termination and release of memory. */
/* -------------------------------------------------------------------- */
     phase = -1;			/* Release internal memory. */
     PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
	      &Sp_TempMat.Rows, &fdum, Sp_TempMat.RowIndex, Sp_TempMat.Columns, &idum, &Keinv->Rows,
	      iparm, &msglvl, &fdum, &fdum, &error);

     Destroy_MatrixVector( &IdentMatrix );
     Destroy_MatrixVector_Sparse( &Sp_TempMat );

     PrintSuccess( "Matrix Inversion successfully completed" );
}


#endif

MatrixVector Generate_IdentityMatrix( int Rows, int Cols )
{
     MatrixVector Identity;
     unsigned short int i;

     if( Rows != Cols ){
	  PrintErrorAndExit( "The number of rows and columns must be equal in order to generate an identity matrix" );
     }

     Init_MatrixVector( &Identity, Rows, Cols );
     
     for( i = 0; i < Rows; i++ ){
	  Identity.Array[i + Rows*i] = 1.0;
     }

     return Identity;
}


void BuildMatrixXc( const MatrixVector *const Mat, double *MatCouple, const Coupling_Node *const CNodes )
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

void BuildMatrixXcm( const MatrixVector *const Mat, MatrixVector *const MatXcm, const Coupling_Node *const CNodes )
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
	  PosXcm = jcoup*MatXcm->Rows;
	  dcopy_( &Length, &Mat->Array[CNodes->Array[jcoup] - 1], &incx, &MatXcm->Array[PosXcm], &incy );
     }

     /* Copy until the last coupling node */
     Acumulated_Length = Length;
     for( icoup = 1; icoup < CNodes->Order; icoup++ ){
	     
	  Length = CNodes->Array[icoup] - CNodes->Array[icoup-1] - 1;
	  for ( jcoup = icoup; jcoup < CNodes->Order; jcoup++ ){
	 
	       PosXcm = jcoup*MatXcm->Rows + Acumulated_Length;
	       
	       dcopy_( &Length, &Mat->Array[(CNodes->Array[jcoup] - 1) + (CNodes->Array[icoup-1])*Mat->Rows], &incx, &MatXcm->Array[PosXcm], &incy );
	  }
	  Acumulated_Length = Acumulated_Length + Length;
     }

     /* Do the same but starting from the opposite side */
     incx = 1;
     Length = Mat->Rows - CNodes->Array[CNodes->Order -1];
     for ( jcoup = CNodes->Order - 1; jcoup >= 0; jcoup = jcoup - 1 ){
	  PosXcm = MatXcm->Rows*(jcoup+1) - Length;
	  dcopy_( &Length, &Mat->Array[(CNodes->Array[jcoup] - 1)*Mat->Rows + (CNodes->Array[CNodes->Order-1])], &incx, &MatXcm->Array[PosXcm], &incy );
     }
     Acumulated_Length = Length;
     incx = 1;
     for( icoup = CNodes->Order -2; icoup >= 0; icoup = icoup -1 ){
	  Length = CNodes->Array[icoup + 1] - CNodes->Array[icoup] - 1;
	  for( jcoup = icoup; jcoup >= 0; jcoup = jcoup - 1){
	       PosXcm = (jcoup+1)*MatXcm->Rows - Acumulated_Length - Length;
	       dcopy_( &Length, &Mat->Array[(CNodes->Array[jcoup] - 1)*Mat->Rows + (CNodes->Array[icoup])], &incx, &MatXcm->Array[PosXcm], &incy );
	  }
	  Acumulated_Length = Acumulated_Length + Length;
     }
}

void Generate_LoadVectorForm( MatrixVector *const LoadVector, int *DOF )
{
     int i, j;

     i = 0;
     while( i < LoadVector->Rows ){
	  
	  for ( j = 1; j < DOF[0]; j++ ){
	       LoadVector->Array[i] = 1.0*(double) DOF[j];
	       i = i + 1;
	  }
     }
}
