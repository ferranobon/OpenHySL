#include <petscsys.h>
#include <petscviewer.h>
#include <petscmat.h>

int PETSc_LoadMatrix_FromFile( MPI_Comm Comm, const char *FileName, Mat *Matrix, const PetscInt Rows, const PetscInt Cols )
{

     Mat Temp;
     int fd1;
     PetscInt header[4], GRows, GCols;
     PetscViewer viewer;

     /* Open the files */
     PetscViewerBinaryOpen( Comm, FileName, FILE_MODE_READ, &viewer );
     PetscViewerBinaryGetDescriptor( viewer, &fd1 );
     PetscBinaryRead( fd1, (char *) header, 4 , PETSC_INT);

     /* error checking on files */
     if ( header[0] != MAT_FILE_CLASSID){
	  SETERRQ( PETSC_COMM_WORLD, PETSC_ERR_FILE_UNEXPECTED, "The input file does not contain a valid PETSc matrix object");
     }

     /* Get the number of global rows and columns */
     GRows = header[1]; GCols = header[2];
     printf( "Global rows %d Global columns %d\n", GRows, GCols );
     if( GRows != Rows || GCols != Cols ){
	  SETERRQ( PETSC_COMM_WORLD, PETSC_ERR_FILE_UNEXPECTED, "The sizes of the input matrix are not the same as the specified ones.");
     }
     PetscViewerDestroy( &viewer );

     PetscViewerBinaryOpen( Comm, FileName, FILE_MODE_READ, &viewer );
     MatCreate( Comm, &Temp );
     MatSetSizes( Temp, PETSC_DECIDE, PETSC_DECIDE, GRows, GCols );
     MatSetType( Temp, MATMPIAIJ );
     MatLoad( Temp, viewer );

     PetscViewerDestroy( &viewer );
     *Matrix = Temp;
     return 0;
}

int PETSc_LoadMatrix_FromFile_Dense( MPI_Comm Comm, const char *FileName, Mat *Matrix, const PetscInt Rows, const PetscInt Cols )
{

     Mat Temp;
     int fd1;
     PetscInt header[4], GRows, GCols;
     PetscViewer viewer;

     /* Open the files */
     PetscViewerBinaryOpen( Comm, FileName, FILE_MODE_READ, &viewer );
     PetscViewerBinaryGetDescriptor( viewer, &fd1 );
     PetscBinaryRead( fd1, (char *) header, 4 , PETSC_INT);

     /* error checking on files */
     if ( header[0] != MAT_FILE_CLASSID){
	  SETERRQ( PETSC_COMM_WORLD, PETSC_ERR_FILE_UNEXPECTED, "The input file does not contain a valid PETSc matrix object");
     }

     /* Get the number of global rows and columns */
     GRows = header[1]; GCols = header[2];
     printf( "Global rows %d Global columns %d\n", GRows, GCols );
     if( GRows != Rows || GCols != Cols ){
	  SETERRQ( PETSC_COMM_WORLD, PETSC_ERR_FILE_UNEXPECTED, "The sizes of the input matrix are not the same as the specified ones.");
     }
     PetscViewerDestroy( &viewer );

     PetscViewerBinaryOpen( Comm, FileName, FILE_MODE_READ, &viewer );
     MatCreate( Comm, &Temp );
     MatSetSizes( Temp, PETSC_DECIDE, PETSC_DECIDE, GRows, GCols );
     MatSetType( Temp, MATMPIDENSE );
     MatLoad( Temp, viewer );

     PetscViewerDestroy( &viewer );
     *Matrix = Temp;
     return 0;
}

int PETSc_LoadVector_FromFile( MPI_Comm Comm, const char *FileName, Vec *Vector, const PetscInt Rows )
{
     Vec Temp;
     int fd1;
     PetscInt header[2], GRows;
     PetscViewer viewer;

     /* Open the files */
     PetscViewerBinaryOpen( Comm, FileName, FILE_MODE_READ, &viewer );
     PetscViewerBinaryGetDescriptor( viewer, &fd1 );
     PetscBinaryRead( fd1, (char *) header, 2 , PETSC_INT);

     /* error checking on files */
     if ( header[0] != VEC_FILE_CLASSID){
	  SETERRQ( PETSC_COMM_WORLD, PETSC_ERR_FILE_UNEXPECTED, "The input file does not contain a valid PETSc vector object");
     }

     /* Get the number of global rows and columns */
     GRows = header[1];

     if( GRows != Rows ){
	  SETERRQ( PETSC_COMM_WORLD, PETSC_ERR_FILE_UNEXPECTED, "The number or rows of the input vector is not the same as the specified one");
     }
     PetscViewerDestroy( &viewer );


     PetscViewerBinaryOpen( Comm, FileName, FILE_MODE_READ, &viewer );
     VecCreate( Comm, &Temp );
     VecSetType( Temp, VECMPI );
     VecSetUp( Temp );
     VecLoad( Temp, viewer );

     PetscViewerDestroy( &viewer );
     *Vector = Temp;
     return 0;
}
