/**
 * \file MatrixVector_Sp.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 19th of February 2013
 *
 * \brief MatrixVector_Sp_t creation and manipulation prototypes.
 *
 * This file contains the prototypes of those functions involved in creating/destroying sparse matrices and
 * vectors.
 */

#ifndef MATRIXVECTOR_SP_H_
#define MATRIXVECTOR_SP_H_

#include "MatrixVector.h"  /* For Scalars_t */

/**
 * \brief Sparse matrix storage. MKL CSR-\em three \em array \em variation.
 *
 * Structure designed to represent sparse matrices and vectors within the library. The components of the
 * structure were selected in order to match those of Intel's MKL CSR-\em three \em array \em variation
 * format: \c Values, \c Columns and \c RowIndex in addition to the number of rows and columns. For dense
 * matrix representation, \c MatrixVector_t should be used instead.
 *
 * \sa MatrixVector_t.
 */
typedef struct MatVec_Sp {
     int Rows;        /*!< \brief Number of Rows of the matrix. */
     int Cols;        /*!< \brief Number of Columns of the matrix. */
     int Num_Nonzero; /*!< \brief Number of non-zero elements. */
     double *Values;  /*!< \brief A double precision array that contains the non-zero elements of a sparse
		       * matrix. The non-zero elements are mapped into the values array using the row-major
		       * upper triangular storage mapping. The lenght of the array is equal to the number of
		       * non-zero elements in the matrix. */
     int *Columns;    /*!< \brief Element \a i of the integer array columns is the number of the column that
		       * contains the i-th element in the values array. The lenght of the array is equal to
		       * the number of non-zero elements in the matrix. */
     int *RowIndex;   /*!< \brief Element \a j of the integer array rowIndex gives the index of the element in
		       * the values array that is first non-zero element in a row j. The length of the array
		       * is equal to the number of rows plus one. */
} MatrixVector_Sp_t;

/**
 * \brief Adds three matrices of the same dimensions. Sparse version
 *
 * \pre
 * - All elements of type \c MatrixVector_Sp_t must be properly intialised through the
 *   MatrixVector_Create_Sp() routine.
 * - The matrices must in Intel's MKL CSR-\em three \em array \em variation and in one based index.
 * - The matrices must be symmetrical and only the upper part will be referenced (lower part in FORTRAN
 *   routines).
 * - \f$S(\mathcal Y) \geq max(S(\mathcal A),S(\mathcal B),S(\mathcal C))\f$ where \f$S(\mathcal X) =
 *   X.Rows*X.Cols\f$ is the size of the matrix.
 * - The number of non-zero elements of \f$nnz(\mathcal Y) \geq nnz(\mathcal A + \mathcal B + \mathcal C)\f$.
 *
 * This routine adds three matrices through the operation:
 *
 * \f[\mathcal Y = \alpha \mathcal A + \beta \mathcal B + \gamma \mathcal C\f]
 *
 * where:
 * - \f$\mathcal Y\f$, \f$\mathcal A\f$, \f$\mathcal B\f$, \f$\mathcal C\f$ are symmetric matrices,
 * - and \f$\alpha\f$, \f$\beta\f$, \f$\gamma\f$ are real scalars.
 *
 * It makes use of BLAS and Sparse BLAS routines from the Intel Math Kernel Library (\cite MKL_2013) to
 * perform the linear algebra operations. For the full or packed storage cases use the MatrixVector_Add3Mat()
 * or MatrixVector_Add3Mat_PS() routines respectively.
 * 
 * \param[in]     MatA  Symmetric matrix \f$\mathcal A\f$ with only the upper part referenced (lower part in
 *                      FORTRAN routines).
 * \param[in]     MatB  Symmetric matrix \f$\mathcal B\f$ with only the upper part referenced (lower part in
 *                      FORTRAN routines).
 * \param[in]     MatC  Symmetric matrix \f$\mathcal C\f$ with only the upper part referenced (lower part in
 *                      FORTRAN routines).
 * \param[in]     Const Scalars that multiply the matrices A (\f$\alpha\f$), \c Const.Aplha), B (\f$\beta\f$),
 *                      \c Const.Beta) and C (\f$\gamma\f$), \c Const.Gamma).
 * \param[in,out] MatY  Symmetric matrix \f$\mathcal Y\f$ with only the upper part referenced (lower part in
 *                      FORTRAN routines). On entry only the dimensions and the number of non-zero elements
 *                      are referenced.
 *
 * \post \c MatY is the result of the operation:
 * \f[\mathcal Y = \alpha \mathcal A + \beta \mathcal B + \gamma \mathcal C\f]
 *
 * \sa MatrixVector_Sp_t Scalars_t.
 */
void MatrixVector_Add3Mat_Sp( const MatrixVector_Sp_t *const MatA, const MatrixVector_Sp_t *const MatB,
			      const MatrixVector_Sp_t *const MatC, const Scalars_t Const, MatrixVector_Sp_t *const MatY );
/**
 * \brief Allocates the memory for the \c MatrixVector_Sp_t type.
 *
 * \pre MatVec_Sp should have been pre-initialised through MatrixVector_SetRowsCols_Sp().
 *
 * Allocates the space for the MatrixVector_Sp_t according to the Intel's MKL CSR-\em three \em array \em
 * variation (\cite MKL_2013) format.
 * - \c MatVec_Sp.Values has a length equal to \c nnz.
 * - \c MatVec_Sp.Columns has a length equal to \c nnz.
 * - \c MatVec_Sp.RowIndex has a length equal to \c MatVec_Sp.Rows.
 *
 * All elements are initialised to 0.
 *
 * \param[in]     nnz       Number of non-zero elements of the array.
 * \param[in,out] MatVec_Sp Sparse matrix or vector to be initialised. On entry only the dimensions are
 *                          referenced.
 *
 * \post
 * - \c MatVec_Sp.Values has a length equal to \c nnz.
 * - \c MatVec_Sp.Columns has a length equal to \c nnz.
 * - \c MatVec_Sp.RowIndex has a length equal to \c MatVec_Sp.Rows.
 * - The memory should be deallocated through MatrixVector_Destroy_Sp().
 *
 * \sa MatrixVector_Sp_t.
 */
void MatrixVector_AllocateSpace_Sp( const int nnz, MatrixVector_Sp_t *const MatVec_Sp );

/**
 * \brief Counts the number of non-zero elements in a general matrix.
 *
 * \param[in] Matrix Single dimensional array of length \f$L = Rows*Cols\f$.
 * \param[in] Rows Number of rows.
 * \param[in] Cols Number of columns.
 * \returns Number of non-zero elements in the general matrix.
 */
int MatrixVector_CountNNZ_GE( const double *const Matrix, const int Rows, const int Cols );

/**
 * \brief Counts the number of non-zero elements in the upper triangular part of a symmetric matrix.
 *
 * \param[in] Sym_Matrix Single dimensional array of length \f$L = Rows^2\f$. Only the upper part will be
 *                       referenced.
 * \param[in] Rows       Number of rows.

 * \returns Number of non-zero elements in the upper triangular part of the symmetric matrix.
 */
int MatrixVector_CountNNZ_SY( const double *const Sym_Matrix, const int Rows );

/**
 * \brief Initialises a sparse matrix or vector.
 *
 * Initialises a sparse matrix or vector according to the Intel's MKL CSR-\em three \em array \em variation
 * (\cite MKL_2013) format.  For a full or packed storage representation the functions MatrixVector_Create()
 * or Matrixvector_Create_PS() should be used instead.
 *
 * \param[in]  Rows      The number of rows.
 * \param[in]  Cols      The number of columns.
 * \param[in]  nnz       Number of non-zero elements of the array.
 * \param[out] MatVec_Sp The matrix or vector to initialise.
 *
 * \post
 * - <tt>MatVec_Sp.Rows = Rows</tt> and <tt>MatVec_Sp.Cols = Cols</tt>.
 * - \c MatVec_Sp.Values has a length equal to \c nnz.
 * - \c MatVec_Sp.Columns has a length equal to \c nnz.
 * - \c MatVec_Sp.RowIndex has a length equal to \c MatVec_Sp.Rows.
 * - The memory should be allocated through MatrixVector_AllocateSpace_Sp().
 *
 * \sa MatrixVector_Sp_t.
 */
void MatrixVector_Create_Sp( const int Rows, const int Cols, const int nnz, MatrixVector_Sp_t *const MatVec_Sp );

/**
 * \brief Converts a sparse matrix/vector into a dense matrix/vector
 * 
 * \pre 
 * - \c MatVec must be properly initialised through the MatrixVector_Create() routine.
 * - \c MatVec_Sp must be properly intialised through the MatrixVector_Create_Sp() routine.
 * - The sparse matrix/vector must in Intel's MKL CSR-\em three \em array \em variation and in one based
 *   index.
 * - If the sparse matrix is symmetric only the upper triangular part must be present.
 * - The dimensions of the matrices/vectors must be the same.
 * - \c Operation must be either 0 (symmetric matrix) or 1 (dense matrix/vector).
 *
 * This routine converts a sparse matrix in Intel's MKL CSR-\em three \em array \em variation (\cite MKL_2013)
 * to a dense format. Currently two operations are available and controlled by \c MatVec_Type. If:
 * 
 * - \f$MatVec_Type = 1\f$ the input matrix is considered to be symmetric and only the upper will be
 *   considered to be present in the sparse matrix.
 * - \f$MatVec_Type = 0\f$ all the elements of the general matrix/vector are considered and transfered to the
 *   dense format.
 *
 * \param[in]     MatVec_Sp   Sparse matrix or vector.
 * \param[in]     MatVec_Type Type of the input sparse matrix: Symmetrical (\f$MatVec_Type=1\f$) or general
 *                            (\f$MatVec_Type = 1\f$).
 * \param[in,out] MatVec      Dense matrix or vector. On input only the dimensions are referenced.
 *
 * \post \c MatVec is a dense representation of the sparse matrix/vector in \c MatVec_Sp.
 *
 * \sa MatrixVector_t and MatrixVector_Sp_t. 
 */
void MatrixVector_CSR2Dense( const MatrixVector_Sp_t *const MatVec_Sp,  const int MatVec_Type, MatrixVector_t *const MatVec );

/**
 * \brief Converts a symmetric sparse matrix in CSR-\em three \em array \em variation into a dense packed
 * storage matrix.
 * 
 * \pre 
 * - \c MatVec_PS must be properly initialised through the MatrixVector_Create_PS() routine.
 * - \c MatVec_Sp must be properly intialised through the MatrixVector_Create_Sp() routine.
 * - The sparse matrix is symmetric, with only the upper triangular part present, in Intel's MKL CSR-\em three
 *   \em array \em variation and in one based index.
 * - The dimensions of the matrices/vectors must be the same.
 *
 * This routine converts a sparse matrix in Intel's MKL CSR-\em three \em array \em variation (\cite MKL_2013)
 * to a dense packed storage format.
 *
 * \param[in]     MatVec_Sp Symmetric sparse matrix.
 * \param[in,out] MatVec_PS Dense matrix in packed storage On input only the dimensions are referenced.
 *
 * \post \c MatVec_PS is a dense representation in packed storage of the sparse matrix/vector in \c
 * MatVec_Sp. Only the upper part is referenced
 *
 * \sa MatrixVector_t and MatrixVector_Sp_t. 
 */
void MatrixVector_CSR2Packed( const MatrixVector_Sp_t *const MatVec_Sp,  MatrixVector_t *const MatVec_PS );


/**
 * \brief Converts a dense matrix/vector into a sparse matrix/vector.
 * 
 * \pre 
 * - \c MatVec must be properly initialised through the MatrixVector_Create() routine.
 * - \c MatVec_Sp must not be intialised.
 * - If the dense matrix is symmetric only the upper triangular part will be referenced.
 * - The dimensions of the matrices/vectors must be the same.
 * - \c Operation must be either 0 (symmetric matrix) or 1 (dense matrix/vector).
 *
 * This routine converts a dense matrix/vector into a sparse matrix/vector in Intel's MKL CSR-\em three \em
 * array \em variation (\cite MKL_2013). Currently two operations are available and controlled by \c
 * MatVec_Type. If:
 * 
 * - \f$MatVec_Type = 1\f$ the input matrix is considered to be symmetric and only the upper will referenced.
 * - \f$MatVec_Type = 0\f$ all the elements of the dense matrix/vector are considered and transfered to the
 *   sparse format.
 *
 * \param[in]  MatVec      Dense matrix or vector.
 * \param[in]  MatVec_Type Type of the input sparse matrix: Symmetrical (\f$MatVec_Type=1\f$) or general
 *                         (\f$MatVec_Type = 1\f$).
 * \param[out] MatVec_Sp   Sparse matrix or vector.
 *
 * \post
 * - \c MatVec_Sp contains only the non-zero elements of the dense matrix/vector \c MatVec.
 * - If the dense matrix is symmetric, only the upper triangular part will be explicitlytransfered to the
 *   sparse matrix, ignoring the lower part.
 * - The sparse matrix/vector is in Intel's MKL CSR-\em three \em array \em variation and in one based index.
 *
 * \sa MatrixVector_t and MatrixVector_Sp_t. 
 */
void MatrixVector_Dense2CSR( const MatrixVector_t *const MatVec, const int MatVec_Type, MatrixVector_Sp_t *const MatVec_Sp );

/**
 * \brief The memory allocated in MatrixVector_Sp_t is freed.
 *
 * \pre \c MatVec_Sp must be properly initialised through MatrixVector_Create_Sp().
 * 
 * \param[out] MatVec_Sp The matrix or vector to destroy.
 *
 * \post
 * - \c MatVec_Sp.Rows, \c MatVec_Sp.Cols and \c MatVec_Sp.Num_Nonzero are set to 0.
 * - \c MatVec_Sp.Array, \c MatVec_Sp.Columns and \c MatVec_Sp.RowIndex no longer point to allocated memory.
 *
 * \sa MatrixVector_t.
 */
void MatrixVector_Destroy_Sp( MatrixVector_Sp_t *const MatVec_Sp );

/**
 * \brief Reads a matrix or a vector in the MatrixMarket format.
 *
 * \warning This routine requires the MatrixMarket header files.
 *
 * \pre
 * - \c MatVec_Sp must have the dimensions specified through MatrixVector_SetRowsCols_Sp().
 * - \c Filename must be in MatrixMarket format and stored in a sparse way.
 *
 * This routine reads a matrix or a vector from a MatrixMarket (\cite MatrixMarket) formatted file. It can
 * handle only sparse formats. If a dense matrix is desired the routine MatrixVector_FromFile_MM() should be
 * used instead.
 *
 * \param[in]     Filename  The file with a MatrixMarket format.
 * \param[in,out] MatVec_Sp On input only the number of rows and columns is referenced.
 *
 * \post \c MatVec_Sp.Array have the contents of the file always in dense storage.
 */
void MatrixVector_FromFile_MM_Sp( const char *Filename, MatrixVector_Sp_t *const MatVec_Sp );

/**
 * \brief Sets the number of rows, columns and non-zero elements of a sparse matrix.
 *
 * \pre \c MatVec_Sp should not be initialised.
 *
 * \param[in]  Rows      The number of rows.
 * \param[in]  Cols      The number of columns.
 * \param[out] MatVec_Sp The matrix or vector to initialise.
 *
 * \post
 * - <tt>MatVec.Rows = Rows</tt> and <tt>MatVec.Cols = Cols</tt>.
 * - The memory should be allocated through MatrixVector_AllocateSpace_Sp().
 *
 * \sa MatrixVector_Sp_t.
 */
void MatrixVector_SetRowsCols_Sp( const int Rows, const int Cols, MatrixVector_Sp_t *const MatVec_Sp );

/**
 * \brief Writes a sparse matrix vector into an ASCII file in Intel's MKL CSR-\em three \em array \em
 * variation format.
 *
 * \pre
 * - \c MatVec_Sp is a sparse matrix in Intel's MKL CSR-\em three \em array \em variation format.
 * - \c MatVec_Sp should be properly initialised through MatrixVector_Create_Sp().
 *
 * \param[in] MatVec_Sp A sparse matrix.
 * \param[in] Filename  Name of the desired ASCII file.
 *
 * \post Filename is an ASCII file with the sparse matrix represented in the Intel's MKL CSR-\em three \em
 * array \em variation format.
 */
void MatrixVector_ToFile_Sp( const MatrixVector_Sp_t *const MatVec_Sp, const char *Filename );

#endif /* MATRIXVECTOR_SP_H_ */
