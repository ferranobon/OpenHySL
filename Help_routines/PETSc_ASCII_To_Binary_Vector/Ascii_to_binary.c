static char help[] = "Reads in a Symmetric matrix in MatrixMarket format. Writes\n\
it using the PETSc sparse format. It also adds a Vector set to random values to the\n\
output file. Input parameters are:\n\
  -fin <filename> : input file\n\
  -fout <filename> : output file\n\n";

#include <petscmat.h>

int main( int argc, char**argv )
{

     Vec A;
     char InFile[PETSC_MAX_PATH_LEN], OutFile[PETSC_MAX_PATH_LEN];
     PetscInt i, row;
     PetscScalar val;
     PetscMPIInt size;
     PetscErrorCode ierr;
     FILE *TheFile;
     PetscViewer viewer;

     PetscInitialize( &argc, &argv, (char *) 0, help );

     ierr = MPI_Comm_size( PETSC_COMM_WORLD, &size); CHKERRQ( ierr );
     if ( size > 1 ){
	  SETERRQ( PETSC_COMM_WORLD, 1, "Uniprocessor routine only\n" );
     }

     /* Read the matrix */
     ierr = PetscOptionsGetString( PETSC_NULL, "-fin", InFile, PETSC_MAX_PATH_LEN, PETSC_NULL ); CHKERRQ( ierr );
     ierr = PetscFOpen( PETSC_COMM_SELF, InFile, "r", &TheFile ); CHKERRQ( ierr );


     /* The first line has the vector dimensions */
     fscanf( TheFile, "%d", &row );
     printf ("m = %d\n", row );

     ierr = VecCreateSeq( PETSC_COMM_SELF, row, &A ); CHKERRQ( ierr );

     for( i = 0; i < row; i++ ){
	  fscanf( TheFile, "%f", (float *) &val );
	  ierr = VecSetValue( A, i, val, INSERT_VALUES); CHKERRQ( ierr );
     }

     fclose( TheFile );

     ierr = VecAssemblyBegin( A ); CHKERRQ( ierr );
     ierr = VecAssemblyEnd( A); CHKERRQ( ierr );

     ierr = PetscPrintf( PETSC_COMM_SELF,"Reading matrix complete.\n" ); CHKERRQ( ierr );
     ierr = PetscOptionsGetString( PETSC_NULL,"-fout",OutFile,PETSC_MAX_PATH_LEN,PETSC_NULL );CHKERRQ( ierr );
     ierr = PetscViewerBinaryOpen( PETSC_COMM_WORLD, OutFile, FILE_MODE_WRITE, &viewer ); CHKERRQ( ierr );

     VecView( A, viewer);

     ierr = VecDestroy( &A); CHKERRQ( ierr );

     ierr = PetscFinalize();

     return 0;
}

