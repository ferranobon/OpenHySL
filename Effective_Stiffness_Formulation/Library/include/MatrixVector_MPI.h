/**
 * \file MatrixVector_MPI.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 18th of March 2013
 *
 * \brief PMatrixVector_t creation and manipulation prototypes.
 *
 * This file contains the prototypes of those functions involved in creating/destroying dense matrices and
 * vectors that are distributed across several processes using MPI.
 */

#ifndef MATRIXVECTOR_MPI_H_
#define MATRIXVECTOR_MPI_H_

#include "MatrixVector.h"

typedef struct DistInfo {
  int Row;
  int Col;
} DistInfo_t ;

typedef struct PMatrixVector {
	double *Array;
	int Desc[9];
	DistInfo_t GlobalSize; /* Stores the size of the global matrix: rows, columns */
	DistInfo_t LocalSize;  /* Stores the size of the local matrix: rows, columns */
	DistInfo_t BlockSize;  /* Block size (vertical, horizontal) */
} PMatrixVector_t;

void PMatrixVector_Create( int icntxt, const int NumRows, const int NumCols, const int BlRows, int const BlCols,
			   PMatrixVector_t *const Mat );
void PMatrixVector_Set2Value( const double Value, PMatrixVector_t *const Mat );
void PMatrixVector_ModifyElement( const int GRowIndex, const int GColIndex, const double Value,
				  const char *Operation, PMatrixVector_t *const Mat );
void PMatrixVector_Add3Mat( PMatrixVector_t *const MatA, PMatrixVector_t *const MatB,
			    PMatrixVector_t *const MatC, const Scalars_t Const, PMatrixVector_t *const MatY );
void PMatrixVector_FromFile( const char *Filename, PMatrixVector_t *const Mat );
void PMatrixVector_FromFile_MM( const char *Filename, PMatrixVector_t *const Mat );
void PMatrixVector_ToFile( PMatrixVector_t *const Mat, const char *Filename );
void PMatrixVector_Destroy( PMatrixVector_t *const Mat );

#endif /* MATRIXVECTOR_MPI_H_ */
