static char help[] = "Reads in a Symmetric matrix in MatrixMarket format. Writes\n\
it using the PETSc sparse format. It also adds a Vector set to random values to the\n\
output file. Input parameters are:\n\
  -fin <filename> : input file\n\
  -fout <filename> : output file\n\n";

#include <petscmat.h>

int main( int argc, char**argv )
{

     Mat A;
     char InFile[PETSC_MAX_PATH_LEN], OutFile[PETSC_MAX_PATH_LEN];
     PetscInt i, m, n, nnz, col, row;
     PetscScalar val;
     PetscMPIInt size;
     PetscErrorCode ierr;
     FILE *TheFile;
     char temp;
     PetscViewer viewer;

     PetscInitialize( &argc, &argv, (char *) 0, help );

     ierr = MPI_Comm_size( PETSC_COMM_WORLD, &size); CHKERRQ( ierr );
     if ( size > 1 ){
	  SETERRQ( PETSC_COMM_WORLD, 1, "Uniprocessor routine only\n" );
     }

     /* Read the matrix */
     ierr = PetscOptionsGetString( PETSC_NULL, "-fin", InFile, PETSC_MAX_PATH_LEN, PETSC_NULL ); CHKERRQ( ierr );
     ierr = PetscFOpen( PETSC_COMM_SELF, InFile, "r", &TheFile ); CHKERRQ( ierr );


     /* The first line has the matrix dimensions */
     fscanf( TheFile, "%d %d %d", &m, &n, &nnz );
     printf ("m = %d, n = %d, nnz = %d\n", m, n, nnz);

     ierr = MatCreateSeqAIJ( PETSC_COMM_WORLD, m, n, nnz, 0, &A ); CHKERRQ( ierr );
     ierr = MatSetOption( A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE); CHKERRQ( ierr );

     for( i = 0; i < nnz; i++ ){
	  fscanf( TheFile, "%d%c%d%c%f", &row, &temp, &col, &temp, (float *) &val );
	  printf("%d %d %f\n", row, col, val );
	  ierr = MatSetValues( A, 1, &row, 1, &col, &val, INSERT_VALUES); CHKERRQ( ierr );
     }

     fclose( TheFile );

     ierr = MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY ); CHKERRQ( ierr );
     ierr = MatAssemblyEnd( A, MAT_FINAL_ASSEMBLY); CHKERRQ( ierr );

     ierr = PetscPrintf( PETSC_COMM_SELF,"Reading matrix complete.\n" ); CHKERRQ( ierr );
     ierr = PetscOptionsGetString( PETSC_NULL,"-fout",OutFile,PETSC_MAX_PATH_LEN,PETSC_NULL );CHKERRQ( ierr );
     ierr = PetscViewerBinaryOpen( PETSC_COMM_WORLD, OutFile, FILE_MODE_WRITE, &viewer ); CHKERRQ( ierr );
     ierr = MatView(A,viewer); CHKERRQ(ierr);

     ierr = MatDestroy( &A); CHKERRQ( ierr );

     ierr = PetscFinalize();

     return 0;
}

