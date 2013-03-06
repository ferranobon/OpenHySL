/**
 * \file MatrixVector.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 19th of February 2013
 * \todo Add support for packaged storages to decrease the memory use
 *
 * \brief MatrixVector_t creation and manipulation prototypes.
 *
 * This file contains the prototypes of those functions involved in creating/destroying dense matrices and vectors.
 */

#ifndef MATRIXVECTOR_H_
#define MATRIXVECTOR_H_

/**
 * \brief General matrix storage.
 *
 * Structure designed to represent matrix and vector types within the library. It stores the number of rows and columns
 * as well as the values within the matrix. It uses a 1-dimensional array for this purpose. This is intended to be used
 * only for dense matrices and vectors. If a sparse format is desired, \c MatrixVector_Sp should be used instead. 
 *
 * \sa MatrixVector_Sp_t.
 */
typedef struct MatVect {
     int Rows;       /*!< Number of Rows of the matrix */
     int Cols;       /*!< \brief Number of Columns of the matrix */
     double *Array;  /*!< \brief Array of size \f$Size = Rows*Cols\f$ */
} MatrixVector_t;

/**
 * \brief Structure to handle constants to be used in some matrix/vector operations.
 */
typedef struct Scal{
     double Alpha;   /*!< \brief First constant.*/
     double Beta;    /*!< \brief Second constant.*/
     double Gamma;   /*!< \brief Third constant.*/
     double Lambda;  /*!< \brief Fourth constant.*/
} Scalars_t;

/**
 * \brief Adds three matrices of the same dimensions.
 * 
 * \warning The effect of using this routine on vectors is unknown.
 *
 * \pre
 * - All elements of type \c MatVec must be properly initialised through MatrixVector_Create().
 * - The matrices must be symmetrical and only the upper part will be referenced (lower part in FORTRAN routines).
 * - \f$S(\mathcal Y) \geq max(S(\mathcal A),S(\mathcal B),S(\mathcal C))\f$ where \f$S(\mathcal X) = X.Rows*X.Cols\f$ is
 * the size of the matrix.
 *
 * This routine adds three matrices through the operation:
 *
 * \f[\mathcal Y = \alpha \mathcal A + \beta \mathcal B + \gamma \mathcal C\f]
 *
 * where:
 * - \f$\mathcal Y\f$, \f$\mathcal A\f$, \f$\mathcal B\f$, \f$\mathcal C\f$ are symmetric matrices,
 * - and \f$\alpha\f$, \f$\beta\f$, \f$\gamma\f$ are real scalars.
 *
 * It makes use of BLAS and LAPACK routines to perform the lineal algebra operations. For the dense version use the
 *  MatrixVector_Add3Mat() routine.
 * 
 * \param[in] MatA Symmetric matrix \f$\mathcal A\f$ with only the upper part referenced (lower part in FORTRAN routines).
 * \param[in] MatB Symmetric matrix \f$\mathcal B\f$ with only the upper part referenced (lower part in FORTRAN routines).
 * \param[in] MatC Symmetric matrix \f$\mathcal C\f$ with only the upper part referenced (lower part in FORTRAN routines).
 * \param[in] Const Scalars that multiply the matrices A (\f$\alpha\f$), \c Const.Aplha), B (\f$\beta\f$), \c Const.Beta)
 * and C (\f$\gamma\f$), \c Const.Gamma).
 * \param[in] MatY Symmetric matrix \f$\mathcal Y\f$ with only the upper part referenced (lower part in FORTRAN routines). On entry
 * only the dimensions are referenced.
 *
 * \post \c MatY is the result of the operation:
 * \f[\mathcal Y = \alpha \mathcal A + \beta \mathcal B + \gamma \mathcal C\f]
 *
 * \sa MatrixVector_t Scalars_t.
 */
void MatrixVector_Add3Mat( const MatrixVector_t *const MatA, const MatrixVector_t *const MatB, const MatrixVector_t *const MatC,
			   const Scalars_t Const, MatrixVector_t *const MatY );
/**
 * \brief Allocates the memory for the \c MatrixVector_t type.
 *
 * \pre \f$ Rows \geq 0\f$ and \f$Cols \geq 0 \f$.
 *
 * A \c MatrixVector_t type is initialised. The routine allocates an amount of memory as a single dimenson array with length
 * \f$L = Rows*Cols\f$. The number of rows and columns is also stored and all elements of the array are initialised and
 * set to 0.0.
 *
 * \param[in] Rows The number of rows.
 * \param[in] Cols The number of columns.
 * \param[out] MatVec The matrix or vector to initialise.
 *
 * \post
 * - <tt>MatVec.Rows = Rows</tt> and <tt>MatVec.Cols = Cols</tt>.
 * - The length of the allocated double array is set to \f$L = Rows*Cols\f$ and all its values initialised to 0.0.
 * - The memory should be deallocated through MatrixVector_Destroy().
 *
 * \sa MatrixVector_t.
 */
void MatrixVector_Create( const int Rows, const int Cols, MatrixVector_t *const MatVec );

/**
 * \brief The memory allocated in MatrixVector_t is freed.
 *
 * \pre \c MatVec must be properly initialised through MatrixVector_Create().
 * 
 * \param[out] MatVec The matrix or vector to destroy.
 *
 * \post
 * - \c MatVec.Rows and \c MatVec.Cols are set to 0.
 * - \c MatVec.Array no longer points to allocated memory.
 *
 * \sa MatrixVector_t.
 */
void MatrixVector_Destroy( MatrixVector_t *const MatVec );

/**
 * \brief Reads a matrix or a vector from an ASCII file in a dense format.
 * 
 * \pre \c MatVec must be properly initialised through MatrixVector_Create().
 *
 * \param[in] Filename Name of the ASCII file to be opened.
 * \param[in,out] MatVec On input only the number of rows and columns is referenced.
 *
 * \post \c MatVec.Array have the contents of the ASCII file.
 */
void MatrixVector_FromFile( const char *Filename, MatrixVector_t *const MatVec );

/**
 * \brief Reads a matrix or a vector in the MatrixMarket format.
 *
 * \warning This routine requires the MatrixMarket header files.
 *
 * \pre
 * - \c MatVec must be properly initialised through MatrixVector_Create().
 * - \c Filename must be in MatrixMarket format and stored in a sparse way.
 *
 * This routine reads a matrix or a vector from a MatrixMarket (\cite MatrixMarket) formatted file. It can
 * handle only sparse formats but the output will always be a dense matrix. If a sparse matrix
 * is desired the routine MatrixVector_FromFile_MM_Sp() should be used instead.
 *
 * \param[in] Filename The file with a MatrixMarket format.
 * \param[in,out] MatVec On input only the number of rows and columns is referenced.
 *
 * \post \c MatVec.Array have the contents of the file always in dense storage.
 */
void MatrixVector_FromFile_MM( const char *Filename, MatrixVector_t *const MatVec );

/**
 * \brief Performs basic algebra operations on an element of the matrix or vector.
 * 
 * \pre
 * - \c MatVec must be properly initialised through MatrixVector_Create().
 * - RowIndex and ColIndex must be in one based index.
 * - \c Operation must be \c Set, \c Add, \c Multiply or \c Divide.
 * 
 * A basic linear algebra operation is performed on one of the elements of the matrix or vector. The operation
 * performed is controlled through the variable \c Operation. Currently only four operations are performed.
 * 
 * - <tt>Operation = Set</tt>. \f$A(i,j) = \alpha\f$.
 * - <tt>Operation = Add</tt>. \f$A(i,j) = A(i,j) + \alpha\f$.
 * - <tt>Operation = Multiply</tt>. \f$A(i,j) = A(i,j)*\alpha\f$.
 * - <tt>Operation = Divide</tt>. \f$A(i,j) = frac{A(i,j)}{\alpha}\f$.
 *
 * If the operation is not supported, the routine calls <tt>exit( EXIT_FAILURE )</tt>.
 *
 * \param[in] RowIndex The row index \f$i\f$ (one based index).
 * \param[in] ColIndex The column index \f$j\f$ (one based index).
 * \param[in] Alpha The value to be set, added, multiplied and divided \f$\alpha\f$.
 * \param[in] Alpha The value to be set, added, multiplied and divided \f$\alpha\f$.
 * \param[in] Operation Controls what operation is performed: \c Set, \c Add, \c Multiply or \c Divide.
 * \param[in,out] MatVec Matrix or vector to be modified. On entry only the number of columns and the value (except in
 * the case when <tt>Operation = Set</tt>) are referenced. On output \f$A(i,j)\f$, is modified accordingly.
 *
 * \post One of the supported operations is performed. If the operation is not supported, the routine calls
 * <tt>exit( EXIT_FAILURE )</tt>.
 */
void MatrixVector_ModifyElement( const int RowIndex, const int ColIndex, const double Alpha,
				 const char *Operation, MatrixVector_t *const MatVec );

/**
 * \brief Sets all the members to the specified value.
 *
 * \pre \c MatVec must be properly initialised through MatrixVector_Create().
 *
 * All the elements in \c MatVec.Array are set to the specified value.It makes use of BLAS routines to perform the lineal algebra operations.
 *
 * \param[in] Value All elements of the matrix or vector will be set to this value.
 * \param[in,out] MatVec Matrix or vector. On input only the number of rows and columns is referenced. On
 * output all its elements are set to \c Value.
 *
 * \post All the elements in \c MatVec.Array are set to \c Value.
 */
void MatrixVector_Set2Value( const double Value, MatrixVector_t *const MatVec );

/**
 * \brief Writes a matrix or a vector to an ASCII file in a dense format.
 * 
 * \pre \c MatVec must be properly initialised through MatrixVector_Create().
 *
 * \param[in] MatVec The matrix or vector to save to the file \c Filename.
 * \param[in] Filename Name of the ASCII file to be opened.
 *
 * \post The ASCII file \c Filename has the contents of \c MatVec.Array in a dense format.
 */
void MatrixVector_ToFile( const MatrixVector_t *const MatVec, const char *Filename );

#endif /* MATRIXVECTOR_H_ */
