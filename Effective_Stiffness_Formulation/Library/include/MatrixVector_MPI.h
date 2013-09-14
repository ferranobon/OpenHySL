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
	DistInfo_t GlobalSize; /*!< Stores the size of the global matrix: rows, columns */
	DistInfo_t LocalSize;  /*!< Stores the size of the local matrix: rows, columns */
	DistInfo_t BlockSize;  /*!< Block size (vertical, horizontal) */
} PMatrixVector_t;

/**
 * \brief Adds three distributed matrices of the same dimensions. MPI version
 * 
 * \warning The effect of using this routine on vectors is unknown.
 *
 * \pre
 * - All elements of type \c MatVec must be properly initialised through PMatrixVector_Create().
 * - The matrices must be symmetrical and only the upper part will be referenced (lower part in FORTRAN
 *   routines).
 * - \f$S(\mathcal Y) \geq max(S(\mathcal A),S(\mathcal B),S(\mathcal C))\f$ where \f$S(\mathcal X) =
 *   X.Rows*X.Cols\f$ is the size of the matrix.
 *
 * This routine adds three distributed matrices through the operation:
 *
 * \f[\mathcal Y = \alpha \mathcal A + \beta \mathcal B + \gamma \mathcal C\f]
 *
 * where:
 * - \f$\mathcal Y\f$, \f$\mathcal A\f$, \f$\mathcal B\f$, \f$\mathcal C\f$ are symmetric matrices,
 * - and \f$\alpha\f$, \f$\beta\f$, \f$\gamma\f$ are real scalars.
 *
 * It makes use of PBLAS and ScaLAPACK routines to perform the lineal algebra operations. For the non
 * distributed matrices the routines MatrixVector_Add3Mat(), MatrixVector_Add3Mat_PS() or
 * MatrixVector_Add3Mat_Sp() should be used depending on the used storage method.
 * 
 * \param[in]     MatA  (global) Symmetric distributed matrix \f$\mathcal A\f$ with only the upper part
 *                      referenced (lower part in FORTRAN routines).
 * \param[in]     MatB  (global) Symmetric distributed matrix \f$\mathcal B\f$ with only the upper part
 *                      referenced (lower part in FORTRAN routines).
 * \param[in]     MatC  (global) Symmetric distributed matrix \f$\mathcal C\f$ with only the upper part
 *                      referenced (lower part in FORTRAN routines).
 * \param[in]     Const (global) Scalars that multiply the matrices A (\f$\alpha\f$), \c Const.Aplha), B
 *                      (\f$\beta\f$), \c Const.Beta) and C (\f$\gamma\f$), \c Const.Gamma).
 * \param[in,out] MatY  (global) Symmetric distributed matrix \f$\mathcal Y\f$ with only the upper part
 *                      referenced (lower part in FORTRAN routines). On entry only the dimensions are
 *                      referenced.
 *
 * \post \c MatY is the result of the operation:
 * \f[\mathcal Y = \alpha \mathcal A + \beta \mathcal B + \gamma \mathcal C\f]
 *
 * \sa PMatrixVector_t and Scalars_t.
 */
void PMatrixVector_Add3Mat( PMatrixVector_t *const MatVecA, PMatrixVector_t *const MatVecB,
			    PMatrixVector_t *const MatVecC, const Scalars_t Const, PMatrixVector_t *const MatVecY );


void PMatrixVector_Create( int icntxt, const int NumRows, const int NumCols, const int BlRows, int const BlCols,
			   PMatrixVector_t *const MatVec );

/**
 * \brief The memory allocated in PMatrixVector_t is freed.
 *
 * \pre \c MatVec must be properly initialised through PMatrixVector_Create().
 * 
 * \param[out] MatVec (global) The distributed matrix or vector to be destroyed.
 *
 * \post
 * - \c MatVec.GlobalSize.Row and \c MatVec.GlobalSize.Col are set to 0.
 * - \c MatVec.LocalSize.Row and \c MatVec.LocalSize.Col are set to 0.
 * - \c MatVec.Array no longer points to allocated memory.
 *
 * \sa PMatrixVector_t.
 */
void PMatrixVector_Destroy( PMatrixVector_t *const MatVec );

/**
 * \brief Reads a matrix or a vector from an ASCII file in a dense format and distributes it to the processes
 * in the grid. MPI version.
 * 
 * \pre \c MatVec must be properly initialised through PMatrixVector_Create().
 *
 * The contents of a file are distributed across the processes of the grid which MatVec belongs to. To do so
 * it first reads all the contents into the process within the grid with coordinates (0,0) (creating a local
 * matrix with a size of \f$MxN\f$ where \em M and \em N are the rows and columns of the distributed matrix)
 * and then distributes all the contents with the rest. For non-distributed matrices or vectors the routine
 * MatrixVector_FromFile() should be used instead.

 * \param[in]     Filename (global) Name of the ASCII file to be opened.
 * \param[in,out] MatVec   On input only the number of rows and columns is referenced.
 *
 * \post \c MatVec.Array have the contents of the ASCII file. The contents vary depending on the process
 * coordinates.
 *
 * \sa PMatrixVector_t.
 */
void PMatrixVector_FromFile( const char *Filename, PMatrixVector_t *const MatVec );

/**
 * \brief Reads a matrix or a vector in the MatrixMarket format and distributes it to the processes in the
 * grid. MPI version.
 *
 * \warning This routine requires the MatrixMarket header files.
 *
 * \pre
 * - \c MatVec must be properly initialised through PMatrixVector_Create().
 * - \c Filename must be in MatrixMarket format and stored in a sparse way.
 *
 * This routine reads a matrix or a vector from a MatrixMarket (\cite MatrixMarket) formatted file. It can
 * handle only sparse formats but the output will always be a dense distributed matrix. The contents of a file
 * are distributed across the processes of the grid which MatVec belongs to. To do so it first reads all the
 * contents into the process within the grid with coordinates (0,0) (creating a local matrix with a size of
 * \f$MxN\f$ where \em M and \em N are the rows and columns of the distributed matrix) and then distributes
 * all the contents with the rest. For non-distributed matrices or vectors the routines
 * PMatrixVector_FromFile_MM(), PMatrixVector_FromFile_MM_Sp() or PMatrixVector_FromFile_MM_PS() should be
 * used instead.
 *
 * \param[in]     Filename (global) The file with a MatrixMarket format.
 * \param[in,out] MatVec   On input only the number of rows and columns is referenced.
 *
 * \post \c MatVec.Array have the contents of the file always in dense storage. The contents vary depending on
 * the process coordinates.
 *
 * \sa PMatrixVector_t.
 */
void PMatrixVector_FromFile_MM( const char *Filename, PMatrixVector_t *const MatVec );

/**
 * \brief Performs basic algebra operations on an element of the distributed matrix or vector. MPI version.
 * 
 * \pre
 * - \c MatVec must be properly initialised through PMatrixVector_Create().
 * - RowIndex and ColIndex must be in one based index.
 * - \c Operation must be \c Set, \c Add, \c Multiply or \c Divide.
 * 
 * A basic linear algebra operation is performed on one of the elements of the matrix or vector. The operation
 * performed is controlled through the variable \c Operation. Currently only four operations are supported.
 * 
 * - <tt>Operation = Set</tt>. \f$A(i,j) = \alpha\f$.
 * - <tt>Operation = Add</tt>. \f$A(i,j) = A(i,j) + \alpha\f$.
 * - <tt>Operation = Multiply</tt>. \f$A(i,j) = A(i,j)*\alpha\f$.
 * - <tt>Operation = Divide</tt>. \f$A(i,j) = frac{A(i,j)}{\alpha}\f$.
 *
 * If the operation is not supported, the routine calls <tt>exit( EXIT_FAILURE )</tt>.
 *
 * \param[in]     RowIndex  The row index \f$i\f$ (one based index).
 * \param[in]     ColIndex  The column index \f$j\f$ (one based index).
 * \param[in]     Alpha     The value to be set, added, multiplied and divided \f$\alpha\f$.
 * \param[in]     Operation Controls what operation is performed: \c Set, \c Add, \c Multiply or \c Divide.
 * \param[in,out] MatVec    Matrix or vector to be modified. On entry only the number of columns and the value
 *                          (except in the case when <tt>Operation = Set</tt>) are referenced. On output
 *                          \f$A(i,j)\f$, is modified accordingly.
 *
 * \post One of the supported operations is performed. If the operation is not supported, the routine calls
 * <tt>exit( EXIT_FAILURE )</tt>.
 *
 * \sa PMatrixVector_t
 */
void PMatrixVector_ModifyElement( int GRowIndex, int GColIndex, const double Value, const char *Operation,
				  PMatrixVector_t *const MatVec );

/**
 * \brief Sets all the members of a distributed matrix or vector to the specified value. MPI version.
 *
 * \pre \c MatVec must be properly initialised through PMatrixVector_Create().
 *
 * All the local elements in \c MatVec.Array are set to the specified value. Although performed on distributed
 * matrices or vectors, the operation is performed localy (BLAS routines).
 *
 * \param[in]     Value  (global) All elements of the matrix or vector will be set to this value.
 * \param[in,out] MatVec (local) Distributed matrix or vector. On input only the number of local rows and
 *                       columns is referenced. On output all its (local) elements are set to \c Value.
 *
 * \post All the elements in \c MatVec.Array are set to \c Value.
 */
void PMatrixVector_Set2Value( const double Value, PMatrixVector_t *const MatVec );

/**
 * \brief Writes a distributed matrix or a vector to an ASCII file in a dense format. MPI version.
 *
 * A distributed matrix is written into an ASCII file in a dense format. To do so a non-distributed version of
 * the matrix (size \f$MxN\f$ where \em M and \em N are the rows and columns of the distributed matrix) is
 * created in the process within the grid (which \c MatVec belongs to) with coordinates (0,0). The 'local'
 * content is then written into the specified file. For non-distributed matrices or vectors the routine
 * MatrixVector_FromFile() should be used instead.
 * 
 * \pre \c MatVec must be properly initialised through MatrixVector_Create().
 *
 * \param[in] MatVec   The matrix or vector to save to the file \c Filename.
 * \param[in] Filename (global) Name of the ASCII file to be opened.
 *
 * \post The ASCII file \c Filename has the contents of \c MatVec.Array in a dense format.
 *
 * \sa PMatrixVector_t.
 */
void PMatrixVector_ToFile( PMatrixVector_t *const MatVec, const char *Filename );

#endif /* MATRIXVECTOR_MPI_H_ */
