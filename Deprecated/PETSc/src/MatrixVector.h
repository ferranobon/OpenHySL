#ifndef MATRIXVECTOR_H_
#define MATRIXVECTOR_H_

#include <petscmat.h>

/**
 * \brief A structure
 * Stores information of several constants that will be used in some matrix vector operations
 */
typedef struct Scal {
	  float Alpha;   /*!< \brief First Constant.*/
	  float Beta;    /*!< \brief Second Constant.*/
	  float Gamma;   /*!< \brief Third Constant.*/
} Scalars;


int PETSc_LoadMatrix_FromFile( MPI_Comm Comm, const char *FileName, Mat *Matrix, const PetscInt Rows, const PetscInt Cols );
int PETSc_LoadMatrix_FromFile_Dense( MPI_Comm Comm, const char *FileName, Mat *Matrix, const PetscInt Rows, const PetscInt Cols );
int PETSc_LoadVector_FromFile( MPI_Comm Comm, const char *FileName, Vec *Vector, const PetscInt Rows );

#endif /* MATRIXVECTOR_H_ */
