/**
 * \file MatrixVector_PS.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 11th of March 2013
 *
 * \brief Creation and manipulation of matrices that use packed storage.
 *
 * This file contains the prototypes of those functions involved in creating/destroying dense matrices in
 * packed storage.
 */

#ifndef MATRIXVECTOR_PS_H_
#define MATRIXVECTOR_PS_H_

#include "MatrixVector.h"
#include "Definitions.h"

/**
 * \brief Adds three matrices of the same dimensions.
 * 
 * \warning The effect of using this routine on vectors is unknown.
 *
 * This routine adds three matrices through the operation:
 *
 * \f[\mathcal Y = \alpha \mathcal A + \beta \mathcal B + \gamma \mathcal C\f]
 *
 * where:
 * - \f$\mathcal Y\f$, \f$\mathcal A\f$, \f$\mathcal B\f$, \f$\mathcal C\f$ are symmetric matrices,
 * - and \f$\alpha\f$, \f$\beta\f$, \f$\gamma\f$ are real scalars.
 *
 * It makes use of BLAS routines to perform the lineal algebra operations. For the full storage and sparse
 * version use the MatrixVector_Add3Mat() and MatrixVector_Add3Mat_Sp() routines respectively.
 *
 * \pre
 * - All elements of type \c Matrix must be properly initialised through MatrixVector_Create_PS().
 * - The matrices must be symmetrical and in packed storage with the upper part referenced (lower part in
 *   FORTRAN routines).
 * - \f$S(\mathcal Y) \geq max(S(\mathcal A),S(\mathcal B),S(\mathcal C))\f$ where \f$S(\mathcal X) =
 *   \frac{(X.Rows*X.Cols + X.Rows)}{2}\f$ is the size of the matrix. 
 * 
 * \param[in]     MatA  Symmetric matrix \f$\mathcal A\f$ in packed storage with the upper part referenced
 *                      (lower part in FORTRAN routines).
 * \param[in]     MatB  Symmetric matrix \f$\mathcal B\f$ in packed storage with the upper part referenced
 *                      (lower part in FORTRAN routines).
 * \param[in]     MatC  Symmetric matrix \f$\mathcal C\f$ in packed storage with the upper part referenced
 *                      (lower part in FORTRAN routines).
 * \param[in]     Const Scalars that multiply the matrices A (\f$\alpha\f$), \c Const.Aplha), B (\f$\beta\f$),
 *                      \c Const.Beta) and C (\f$\gamma\f$), \c Const.Gamma).
 * \param[in,out] MatY  Symmetric matrix \f$\mathcal Y\f$ in packed storage with the upper part referenced
 *                      (lower part in FORTRAN routines). On entry only the dimensions are referenced.
 *
 * \post \c MatY is the result of the operation:
 * \f[\mathcal Y = \alpha \mathcal A + \beta \mathcal B + \gamma \mathcal C\f]
 *
 * \sa MatrixVector_t and Scalars_t.
 */
void MatrixVector_Add3Mat_PS( const MatrixVector_t *const MatA, const MatrixVector_t *const MatB,
			   const MatrixVector_t *const MatC, const Scalars_t Const,
			   MatrixVector_t *const MatY );
/**
 * \brief Allocates the memory for the \c MatrixVector_t type.
 *
 * A \c MatrixVector_t type is initialised in packed storage. The routine allocates an amount of memory as a
 * single dimenson array with length \f$L = (Rows*Cols + Rows)/2\f$. The number of rows and columns is also
 * stored and all elements of the array are initialised and set to 0.0. For matrices in full storage or sparse
 * matrices the routines MatrixVector_Create() or MatrixVector_Create_Sp() should be used respectively.
 *
 * \pre \f$ Rows \geq 0\f$ and \f$Cols \geq 0 \f$.
 *
 * \param[in]  Rows   The number of rows.
 * \param[in]  Cols   The number of columns.
 * \param[out] Matrix The matrix or vector to initialise.
 *
 * \post
 * - <tt>Matrix.Rows = Rows</tt> and <tt>Matrix.Cols = Cols</tt>.
 * - The length of the allocated double/single precision array is set to \f$L = Rows*Cols\f$ and all its
 * values initialised to 0.0.
 * - The memory should be deallocated through MatrixVector_Destroy_PS().
 *
 * \sa MatrixVector_t.
 */
void MatrixVector_Create_PS( const int Rows, const int Cols, MatrixVector_t *const Matrix );

/**
 * \brief The memory allocated in MatrixVector_t is freed.
 *
 * \pre \c Matrix must be properly initialised through MatrixVector_Create_PS().
 * 
 * \param[out] Matrix The matrix in packed storage to be destroyed.
 *
 * \post
 * - \c Matrix.Rows and \c Matrix.Cols are set to 0.
 * - \c Matrix.Array no longer points to allocated memory.
 *
 * \sa MatrixVector_t.
 */
void MatrixVector_Destroy_PS( MatrixVector_t *const Matrix );

/**
 * \brief Reads a symmetric matrix in full storage from an ASCII file and stores the contents in packed
 * storage format.
 * 
 * \pre
 * - \c Matrix must be properly initialised through MatrixVector_Create_PS().
 * - \c The file must contain a symmetric matrix with at least the upper part referenced. The the lower part
 *   must be present but its values will not be stored.
 *
 * \param[in]     Filename Name of the ASCII file to be opened.
 * \param[in,out] Matrix   On input only the number of rows and columns is referenced.
 *
 * \post \c Matrix.Array have the contents of the ASCII file in packed storage.
 */
void MatrixVector_FromFile_GE2PS( const char *Filename, MatrixVector_t *const Matrix );

/**
 * \brief Reads a matrix or a vector in the MatrixMarket format.
 *
 * \warning This routine requires the MatrixMarket header files.
 *
 * This routine reads a symmetric matrix from a MatrixMarket (\cite MatrixMarket) formatted file and stores
 * its contents in a packed storage format. It can handle only sparse formats but the output will always be in
 * packed storage matrix. If a sparse or a full storage format is desired the routines
 * MatrixVector_FromFile_MM_Sp() or MatrixVector_FromFile_MM() should be used instead.
 *
 * \pre
 * - \c Matrix must be properly initialised through MatrixVector_Create_PS().
 * - \c Filename must be in MatrixMarket format and stored in a sparse way.
 *
 *
 * \param[in]     Filename The file with a MatrixMarket format.
 * \param[in,out] Matrix   On input only the number of rows and columns is referenced.
 *
 * \post \c Matrix.Array have the contents of the file always in packed storage.
 */
void MatrixVector_FromFile_MM_PS( const char *Filename, MatrixVector_t *const Matrix );

/**
 * \brief Performs basic algebra operations on an element of a matrix in packed storage.
 * 
 * A basic linear algebra operation is performed on one of the elements of the matrix in packed storage
 * vector. The operation performed is controlled through the variable \c Operation. Currently only four
 * operations are supported.
 * 
 * - <tt>Operation = Set</tt>. \f$A(i,j) = \alpha\f$.
 * - <tt>Operation = Add</tt>. \f$A(i,j) = A(i,j) + \alpha\f$.
 * - <tt>Operation = Multiply</tt>. \f$A(i,j) = A(i,j)*\alpha\f$.
 * - <tt>Operation = Divide</tt>. \f$A(i,j) = frac{A(i,j)}{\alpha}\f$.
 *
 * If the operation is not supported, the routine calls <tt>exit( EXIT_FAILURE )</tt>.
 *
 * \pre
 * - \c Matrix must be properly initialised through MatrixVector_Create_PS().
 * - \c Matrix contains only the upper part of the whole matrix and it is in packed storage.
 * - RowIndex and ColIndex must be in one based index.
 * - \c Operation must be \c Set, \c Add, \c Multiply or \c Divide.
 *
 * \param[in]     RowIndex  The row index \f$i\f$ (one based index).
 * \param[in]     ColIndex  The column index \f$j\f$ (one based index).
 * \param[in]     Alpha     The value to be set, added, multiplied and divided \f$\alpha\f$.
 * \param[in]     Operation Controls what operation is performed: \c Set, \c Add, \c Multiply or \c Divide.
 * \param[in,out] Matrix    Matrix or vector to be modified. On entry only the number of columns and the value
 *                          (except in the case when <tt>Operation = Set</tt>) are referenced. On output
 *                          \f$A(i,j)\f$, is modified accordingly.
 *
 * \post One of the supported operations is performed. If the operation is not supported, the routine calls
 * <tt>exit( EXIT_FAILURE )</tt>.
 */
void MatrixVector_ModifyElement_PS( const int RowIndex, const int ColIndex, const hysl_float_t Alpha,
				 const char *Operation, MatrixVector_t *const Matrix );

/**
 * \brief Sets all the members to the specified value.
 *
 * All the elements in \c Matrix.Array are set to the specified value. For a full representation of the matrix
 * the routine MatrixVector_Set2Value_PS() should be used instead. It makes use of BLAS routines to perform
 * the lineal algebra operations.
 *
 * \pre \c Matrix must be properly initialised through MatrixVector_Create_PS() and be in packed storage.
 *
 * \param[in]     Value  All elements of the matrix or vector will be set to this value.
 * \param[in,out] Matrix Matrix in packed storage. On input only the number of rows and columns is
 *                       referenced. On output all its elements are set to \c Value.
 *
 * \post All the elements in \c Matrix.Array are set to \c Value.
 */
void MatrixVector_Set2Value_PS( const hysl_float_t Value, MatrixVector_t *const Matrix );

/**
 * \brief Writes a matrix in packed storage and with the upper part referenced to an ASCII file in a full
 * format.
 * 
 * \pre
 * - \c Matrix must be properly initialised through MatrixVector_Create_PS().
 * - \c Matrix must be symmetric and in packed storage with the upper part referenced.
 *
 * \param[in] Matrix   The symmetric matrix in packed storage to save to the file \c Filename.
 * \param[in] Filename Name of the ASCII file to be opened.
 *
 * \post The ASCII file \c Filename has the contents of \c Matrix.Array in a full format.
 */
void MatrixVector_ToFile_PS2Full( const MatrixVector_t *const Matrix, const char *Filename );

#endif /* MATRIXVECTOR_H_ */
