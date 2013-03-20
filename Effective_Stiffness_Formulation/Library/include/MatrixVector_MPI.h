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

void PMatrixVector_Create( int icntxt, PMatrixVector_t *const Mat, const int NumRows, const int NumCols, const int BlRows, int const BlCols );
void PMatrixVector_Set2Value( PMatrixVector_t *const Mat, const double Value );
void PMatrixVector_ModifyElement( PMatrixVector_t *const Mat, int GRowIndex, int GColIndex, const double Value, const char *Operation );
void PMatrixVector_Add3Mat( const PMatrixVector_t *const MatA, const PMatrixVector_t *const MatB, const PMatrixVector_t *const MatC,
			    PMatrixVector_t *const MatY, Scalars_t Const );
void PMatrixVector_FromFile( const char *Filename, PMatrixVector_t *const Mat );
void PMatrixVector_FromFile_MM( const char *Filename, PMatrixVector_t *const Mat );
void PMatrixVector_ToFile( PMatrixVector_t *const Mat, const char *Filename );
void PMatrixVector_Destroy( PMatrixVector_t *const Mat );

#endif /* MATRIXVECTOR_MPI_H_ */
