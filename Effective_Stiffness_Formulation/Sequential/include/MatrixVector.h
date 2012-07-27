/**
 * \file MatrixVector.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 10th of August 2011
 * \todo Add support for packaged storages to decrease the memory use
 *
 * \brief MatrixVector creation and manipulation prototypes.
 *
 * This file contains the prototypes of those functions involved in creating/destroying matrices and vectors. Essential
 * matrix/vector manipulations are also contemplated.
 */

#ifndef MATRIXVECTOR_H_
#define MATRIXVECTOR_H_

/**
 * \brief General Matrix Storage
 *
 * The matrix structure contains information to store a general matrix. The elements are
 * stored in row major order.
 */
typedef struct MatVec {
	int Rows; /*!< \brief Number of Rows of the matrix */
	int Cols; /*!< \brief Number of Columns of the matrix */
	float *Array; /*!< \brief Array of size \f$Size = Rows*Cols\f$ */
} MatrixVector;

/**
 * \brief A structure
 * Stores information of several constants that will be used in some matrix vector operations
 */
typedef struct Scal {
	  float Alpha;   /*!< \brief First Constant.*/
	  float Beta;    /*!< \brief Second Constant.*/
	  float Gamma;   /*!< \brief Third Constant.*/
} Scalars;

/**
 * \brief Creates and initialises to zero the elements in the structure MatrixVector
 *
 * Function used to create and initialise the structure MatrixVector. It makes use of the calloc() function to allocate
 * the memory dynamically and initialise all the elements to zero. The elements are located using a 1-dimensional array
 * for both matrices and vectors, being the size of it the product of \f$Rows*Columns\f$.
 *
 * \param[out] Mat The matrix or vector to create
 * \param[in] Rows The number of rows of the matrix/vector
 * \param[in] Cols The number of columns of the matrix/vector
 *
 * \sa MatrixVector.
 */
void Init_MatrixVector( MatrixVector *const Mat, const int Rows, const int Cols );

/**
 * \brief Reads a Matrix or a vector from a file
 *
 * The contents of the desired text file, will be read and stored into the Matrix or Vector.
 *
 * \pre The file is supposed to be an ASCII file and the elements need to be in the proper order (row or major or the
 * appropiate packing storage) since the routine stores the elements in the memory sequentially.
 * Also, the data structure MatrixVector should be properly initialised through the Init_MatrixVector() routine
 *
 * \param[in,out] Mat The matrix or vector where the content of the file will be stored
 * \param[in] Filename The name of the file to be opened.
 *
 * \sa MatrixVector.
 */
void MatrixVector_From_File( MatrixVector *const Mat, const char *Filename );

/**
 * \brief Sets all the elements to the desired value
 *
 * All the elements are set to the desired value. This function makes use of the BLAS routine \c dcopy_( ) with
 * \f$\textrm{incx} = 0\f$ and \f$\textrm{incy} = 1\f$
 *
 * \pre The data structure MatrixVector should be properly initialised through the Init_MatrixVector() routine.
 *
 * \param[in,out] Mat The matrix or vector to be modified
 * \param[in] Value Contains the value that will be assigned to all the elements of \c Mat.
 *
 * \sa MatrixVector.
 */
void Set2Value( MatrixVector *const Mat, const float Value );

/**
 * \brief Modify the element given by the coordinates of the matrix or vector by a given operation
 *
 * The element defined by its coordinates, \f$(RowIndex, ColIndex)\f$ is modified according to the basic operations: Set, Add, Multiply
 * and Divide. The operation is selected during function call through the argument operation. The Matrix must contain elements of type float.
 *
 * \pre The data structure MatrixVector should be properly initialised through the Init_MatrixVector() routine.
 *
 * \param[in,out] Mat The matrix or vector whose element will be modified.
 * \param[in] RowIndex stores the row index of the element.
 * \param[in] ColIndex stores the column index of the element.
 * \param[in] Value stores the value that will be used in the selected operations.
 * \param[in] Operation It can have four possible values: 'Set', 'Add', 'Multiply' and 'Divide'. Assuming (i,j) to be the
 * desired element in the Matrix Mat:
 * \li If \c operation = 'Set', then \f$M(i,j) = Value\f$.
 * \li If \c operation = 'Add', then \f$M(i,j) = M(i,j) + Value\f$. Note that a negative Value will result in a substraction.
 * \li If \c operation = 'Multiply', then  \f$M(i,j) = M(i,j)*Value\f$.
 * \li If \c operation = 'Divide', then \f$M(i,j) = M(i,j)/Value\f$.
 *
 * \sa MatrixVector.
 */
void Modify_Element( MatrixVector *const Mat, const int RowIndex, const int ColIndex, const float Value, const char *Operation );

/**
 * \brief Adds three symmetric matrices
 *
 * This routine adds three symmetric matrices into a new matrix Y through the formula \f$ [Y] = \alpha [A] + \beta [B] + \gamma [C]\f$
 * Only the upper part is referenced (Lower part in FORTRAN routines). It makes use of the dlacpy_( ), dlascl_( ) and daxpy_( ) routines.
 *
 * \pre The data structures of type MatrixVector should be properly initialised through the Init_MatrixVector() routine.
 * Furthermore, the matrices should be symmetric and the elements of the upper part (Lower part in FORTRAN routines) must be
 * stored in general storage. Furthermore, \f$S(Y) \geq max(S(A),S(B),S(C))\f$ where \f$S(X) = X.Rows*X.Cols\f$ is the size of
 * the matrix.
 *
 * \param[in,out] MatY will store the upper part of the operation.
 * \param[in] MatA matrix that will be multiplied by the first constant.
 * \param[in] MatB matrix that will be multiplied by the second constant.
 * \param[in] MatC matrix that will be multiplied by the third constant.
 * \param[in] Const an argument of type Scalars. On entry, \c Const.Alpha \f$= \alpha\f$, \c Const.Beta \f$= \beta\f$ and \c Const.Gamma \f$= \gamma\f$.
 *
 * \sa MatrixVector, Scalars.
 */
void Add3Mat( MatrixVector *const MatY, const MatrixVector *const MatA, const MatrixVector *const MatB, const MatrixVector *const MatC, const Scalars Const );

/**
 * \brief Writes a matrix or a vector to a file
 *
 * The contents of the matrix or vector are saved into a file in row-major order. If the file cannot be opened, the program
 * exits abnormally.
 *
 * \pre The data structure must be properly initialised and their elements must in column-major order.
 *
 * \param[in] Mat The matrix or vector whose array must be stored into a file
 * \param[in] Filename The name of the file where the contents of \c Mat will be stored.
 *
 * \sa MatrixVector.
 */
void MatrixVector_To_File( const MatrixVector *const Mat, const char *Filename );

/**
 * \brief Deallocates the memory of the desired MatrixVector.
 *
 * The memory is deallocated using free() and the number of rows and columns is set to 0. This routine should be used as soon as a
 * matrix or vector becomes useless to free the memory or before changing size of the matrix or vector.
 *
 * \pre The data structure MatrixVector should be properly initialised.
 *
 * \param[out] Mat The matrix or vector to be destroyed.
 *
 * \sa MatrixVector.
 */
void Destroy_MatrixVector( MatrixVector *const Mat);

/**
 * \brief Returns the maximum of two values
 *
 * The maximum of two float values is returned at the end of the function.
 *
 * \param a First value
 * \param b Second value
 * \return max(a,b)
 */
int Max ( const int a, const int b );

#endif /* MATRIXVECTOR_H_ */
