/*
 * PMatrixVector.h
 *
 *  Created on: 22/07/2011
 *      Author: ferran
 */

#ifndef PMATRIXVECTOR_H_
#define PMATRIXVECTOR_H_

typedef struct {
	float Alpha;
	float Beta;
	float Gamma;
} Scalars;

typedef struct {
  int Row;
  int Col;
} AdditionalInfo ;

typedef struct {
	float *Array;
	int Desc[9];
	AdditionalInfo GlobalSize; /* Stores the size of the global matrix: rows, columns */
	AdditionalInfo LocalSize;  /* Stores the size of the local matrix: rows, columns */
	AdditionalInfo BlockSize;  /* Block size (vertical, horizontal) */
} PMatrixVector;

void CreateDistMatrix( int icntxt, PMatrixVector *const Mat, const int NumRows, const int NumCols, const int BlRows, int const BlCols );
void Set2Value( PMatrixVector *const Mat, const float Value );
void ModifyElement( PMatrixVector *const Mat, int GRowIndex, int GColIndex, const float Value, const char *Operation );
void PAdd3Mat( PMatrixVector *const MatY, PMatrixVector *const MatA, PMatrixVector *const MatB, PMatrixVector *const MatC, Scalars Const );
void DistMatrixFromFile( PMatrixVector *const Mat, const char *Filename );
void DistMatrixToFile( PMatrixVector *const Mat, const char *Filename );
void DestroyDistMatrix( PMatrixVector *const Mat );
int Max ( const int a, const int b );

#endif /* PMATRIXVECTOR_H_ */
