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
	double *Array; /*!< \brief Array of size \f$Size = Rows*Cols\f$ */
} MatrixVector;

/**
 * \brief Sparse Matrix Storage. MKL CSR-three array variation.
 *
 * This structure is used in order to store sparse matrices. The Intel MKL CSR (Compressed Sparse Row) three array
 * variation [REFERENCE MKL needed] is used: \c Values, \c Columns and \c RowIndex.
 *
 */
typedef struct SpMatVec{
     int Rows;        /*!< \brief Number of Rows of the matrix. */
     int Cols;        /*!< \brief Number of Columns of the matrix. */
     int Num_Nonzero; /*!< \brief Number of non-zero elements. */
     double *Values;   /*!< \brief A real or complex array that contains the non-zero elements of a sparse matrix. The non-zero
		       * elements are mapped into the values array using the row-major upper triangular storage mapping. The
		       * lenght of the array is equal to the number of non-zero elements in the matrix. */
     int *Columns;    /*!< \brief Element \a i of the integer array columns is the number of the column that contains the
		       * i-th element in the values array. The lenght of the array is equal to the number of non-zero elements
		       * in the matrix. */
     int *RowIndex;   /*!< \brief Element \a j of the integer array rowIndex gives the index of the element in the values
		       * array that is first non-zero element in a row j. The length of the array is equal to the number of
		       * rows plus one. */
} Sp_MatrixVector;

/**
 * \brief A structure
 * Stores information of several constants that will be used in some matrix vector operations
 */
typedef struct Scal {
	  double Alpha;   /*!< \brief First Constant.*/
	  double Beta;    /*!< \brief Second Constant.*/
	  double Gamma;   /*!< \brief Third Constant.*/
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

void Set_RowsCols_Sp( Sp_MatrixVector *const Mat, const int Rows, const int Cols );
void AllocateSpace_Sp( Sp_MatrixVector *const Mat, const int nnz );
void Init_MatrixVector_Sp( Sp_MatrixVector *const Mat, const int Rows, const int Cols, const int nnz );

#if _SPARSE_
/**
 * \brief Converts a dense matrix into CSR-three array variation format.
 *
 * The dense matrix is converted into the CSR-three array variation format of the Intel MKL library. It first counts the
 * number of non-zero elements and afterwards it allocates the necessary memory for the Values and Columns arrays
 * (\sa Sp_MatVec). The CSR matrix will be zero-based indexing and will contain the upper triangular part of the dense matrix.
 *
 * \param[in] Mat Dense matrix.
 * \param[out] Sp_Mat Sparse matrix stored in CSR-three variation array. It contains the upper triangular part of the dense
 * matrix and follows a one-based indexing.
 * \param[in] Operation If
 * - \f$Operation = 0\f$ the matrix is considered to be symmetric and only the upper triangular
 * part is considered.
 * - \f$Operation = 1\f$ the matrix is considered to be general.
 *
 */
void Dense_to_CSR( const MatrixVector *const Mat, Sp_MatrixVector *const Sp_Mat, const int Operation );

void CSR_to_Dense( const Sp_MatrixVector *const Sp_Mat, MatrixVector *const Mat, const int Operation );

void MatrixVector_From_File_Sp( Sp_MatrixVector *const Mat, const char *Filename );
#endif

/**
 * \brief Counts the non-zero elements in a Symmetric Matrix.
 *
 * The non-zero elements of a Symmetric Matrix are counted. Only the upper triangular part is explicitly used (Lower diagonal
 * part in Fortran).
 *
 * \param[in] Sym_Matrix The considered symmetrical matrix.
 * \param[in] Rows The number of rows of \c Sym_Matrix.
 *
 * \return Number of non-zero elements in the upper triangular part of the matrix.
 *
 */
int Count_Nonzero_Elements_SY( const double *const Sym_Matrix, const int Rows );

/**
 * \brief Counts the non-zero elements in a General Matrix.
 *
 * The non-zero elements of a general Matrix are counted.
 *
 * \param[in] Matrix The considered general matrix.
 * \param[in] Rows The number of rows of \c Matrix.
 * \param[in] Columns The number of columns of \c Matrix.
 *
 * \return Number of non-zero elements in the general matrix.
 *
 */
int Count_Nonzero_Elements_GE( const double *const Sym_Matrix, const int Rows, const int Cols );

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
 * \brief Reads a sparse matrix/vector from a file and stores it in dense format.
 *
 * A sparse matrix/vector will be read from a file and it will be stored in dense
 * format. The file must be formatted as follows, and the row and column number must
 * begin at zero:
 *
 * 0:0:45.2700
 * 1:1:67.8991
 * 5:0:51.2414
 * ...
 *
 * Note that the order of the elements does not matter in this case. This routine
 * will read the contents within the file until it reaches the end, making use of
 * the feof() function. Therefore it is important to make sure that the amount of 
 * dynamically allocated memory is enough to handle properly the contents of the
 * file.
 *
 * \pre
 * - The file is supposed to be an ASCII file with the format mentioned above.
 * - The data structure MatrixVector should be properly initialised through the Init_MatrixVector() routine.
 *
 * \param[in,out] Mat The matrix or vector where the content of the file will be
 * stored in a dense format.
 * \param[in] Filename The name of the file to be opened.
 *
 * \sa MatrixVector.
 */
void MatrixVector_From_File_Sp2Dense( MatrixVector *const Mat, const char *Filename );

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
void Set2Value( MatrixVector *const Mat, const double Value );

/**
 * \brief Modify the element given by the coordinates of the matrix or vector by a given operation
 *
 * The element defined by its coordinates, \f$(RowIndex, ColIndex)\f$ is modified according to the basic operations: Set, Add, Multiply
 * and Divide. The operation is selected during function call through the argument operation. The Matrix must contain elements of type double.
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
void Modify_Element( MatrixVector *const Mat, const int RowIndex, const int ColIndex, const double Value, const char *Operation );

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

void Add3Mat_Sparse( Sp_MatrixVector *const MatY, const Sp_MatrixVector *const MatA, const Sp_MatrixVector *const MatB, const Sp_MatrixVector *const MatC, const Scalars Const );

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
 * \brief Writes a sparse matrix or a vector to a file.
 *
 * The contents of the sparse matrix or vector are saved into a file in
 * row-major order and in CSR representation. If the file cannot be opened, the
 * program exits abnormally.
 *
 * \pre The data structure must be properly initialised and their elements must be in
 * CSR-three array variation and zero-based index.
 *
 * \param[in] Sp_Mat The sparse matrix or vector whose array must be stored into a file.
 * \param[in] Filename The name of the file where the contents of \c Sp_Mat will be stored.
 *
 * \sa Sp_MatrixVector.
 */
void MatrixVector_To_File_Sparse( const Sp_MatrixVector *const Sp_Mat, const char *Filename );

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
 * \brief Deallocates the memory of the desired Sp_MatrixVector in CSR-three array variation format.
 *
 * The memory is deallocated using free() and the number of rows, columns and non-zero elements is set to 0.
 * This routine should be used as soon as a matrix or vector becomes useless to free the memory or before
 * changing size of the matrix or vector.
 *
 * \pre - The data structure Sp_MatrixVector should be properly initialised.
 *      - The format of the Sparse matrix should be compliant with the three array variation of Intel MKL libraries.
 *
 * \param[out] Sp_Mat The sparse matrix or vector to be destroyed.
 *
 * \sa Sp_MatrixVector.
 */
void Destroy_MatrixVector_Sparse( Sp_MatrixVector *const Sp_Mat );

/**
 * \brief Returns the maximum of two values
 *
 * The maximum of two double values is returned at the end of the function.
 *
 * \param a First value
 * \param b Second value
 * \return max(a,b)
 */
int Max ( const int a, const int b );

#endif /* MATRIXVECTOR_H_ */