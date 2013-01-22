/**
 * \file Initiation.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 10th of August 2011
 * \todo Add support for packaged storages to decrease the memory use.
 *
 * \brief Prototypes of the functions used in the Initiation phase.
 *
 * This file defines the routines used during the initiation phase of the substructure algorithm: construction of the Proportional Viscous
 * Damping Matrix, the Effective Mass Matrix and the Gain Matrix and their non-coupling and coupling part.
 */

#ifndef INITIATION_H_
#define INITIATION_H_

#include "MatrixVector.h"
#include "Conf_Parser.h"
#include "Send_Receive_Data.h"
#include "Substructure.h"

/**
 * \brief Structure to handle the PID error compensator.
 *
 * This structure handles the PID error compensator type, by storing the three parameters
 * that characterise a PID: the proportional, the integral and derivative part.
 */
typedef struct {
     double P; /*!< \brief Proportional constant */
     double I; /*!< \brief Integral constant */
     double D; /*!< \brief Derivative constant */
} PIDValues;

/**
 * \brief Structure to store the Rayleigh Viscous Damping.
 *
 * This struct consists of two variables to store the parameters required to build
 * the Viscous Damping Matrix using Rayleigh, that is \f$\mathcal{C} = \alpha\mathcal{M} + \beta\mathcal{K}\f$, where:
 * \li \f$\mathcal{C}\f$ is the Viscous Damping Matrix,
 * \li \f$\mathcal{M}\f$ is the Mass Matrix,
 * \li \f$\mathcal{K}\f$ is the Stiffness Matrix,
 * \li and \f$\alpha\f$ and \f$beta\f$ are the parameters that multiply the mass and stiffness matrix respectively.
 */
typedef struct {
     double Alpha; /*!< \brief First constant.*/
     double Beta;  /*!< \brief Second constant.*/
} RayleighConst;

/**
 * \brief Structure to handle the Newmark Constants.
 *
 * This structure handles the Newmark constants, by storing the two parameters required by the Newmark-\f$beta\f$
 * method: \f$\gamma\f$ and \f$\beta\f$.
 */
typedef struct {
     double Gamma; /*!< \brief First constant.*/
     double Beta;  /*!< \brief Second constant.*/
} NewmarkConst;

/**
 * \brief Structure to wrap constants and filenames.
 *
 * This structure stores several constants that will be used in different parts of the substructure
 * algorithm, such as the order of the matrices, the position of the coupling nodes, number of steps
 * filenames and other contants. It act as a configuration file and it has to be modified by hand
 * before compiling.
 *
 * \sa RayleighConst, NewmarkConst, PIDValues and InitConstants()
 *
 */
typedef struct {

     int Order;               /*!< \brief Order of the matrices */
     int OrderSub;   /*!< \brief Number of substructures */

     unsigned int Nstep;      /*!< \brief Number of steps */
     unsigned int NSubstep;   /*!< \brief Number of sub-steps */

     int Use_Absolute_Values; /*!< \brief Variable to control whether to use absolute values in the equation of motion or relative values. Affects how the input load is calculated */

     int Read_Sparse;         /*!< \brief The matrices will be read using a sparse format. \sa MatrixVector_From_File_Sp2Dense() and MatrixVector_From_File_Sp(). */
     int Use_Sparse;          /*!< \brief Sparse routines will be used instead of dense ones. Only available with the support of the Intel MKL libraries.
				\sa MatrixVector_From_File_Sp(), CalculateMatrixC_Sp(), CalculateMatrixKeinv_Pardiso_Sparse(), Calc_Input_Load_AbsValues_Sparse(),
				Calc_Input_Load_RelValues_Sparse() and EffK_Calc_Effective_Force_Sparse(). */
     int Use_Pardiso;         /*!< \brief Uses PARDISO solver instead of LAPACK to compute the matrix inversion. It is only available throught the MKL libraries.
				\sa CaclulateMatrixKeinv_Pardiso() and CalculateMatrixKeinv_Pardiso_Sparse(). */
     int Read_LVector;        /*!< \brief Read the load vector instead of generating it. */
     int *ExcitedDOF;

     double Delta_t;           /*!< \brief Time increment \f$\Delta t\f$ */
     double DeltaT_Sub;        /*!< \brief Time increment for the sub-stepping process */

     RayleighConst Rayleigh;  /*!< \brief Stores Rayleigh Constants alpha (\c Rayleigh.Alpha or \f$\alpha_R\f$) and beta (\c Rayleigh.Beta or \f$\beta_R\f$) */
     NewmarkConst Newmark;    /*!< \brief Stores Newmark Constants gamma (\c Newmark.Gamma or \f$\gamma_N\f$) and (\c Newmark.Beta or \f$\beta_N\f$) */
     PIDValues PID;           /*!< \brief Stores proportional (\c PID.P), integral (\c PID.I) and derivative (\c PID.D) part. */

     /* Constants for Calculations */
     double Const1;           /*!< \brief \f$Const_1=\beta_N\Delta t^2\f$ */
     double Const2;           /*!< \brief \f$Const_2 = (0.5 - 2*\beta_N + \gamma_N)*\Delta t^2\f$ */
     double Const3;           /*!< \brief \f$Const_3 = (0.5 + \beta_N - \gamma_N)*\Delta t^2\f$ */

     /* Constants for Step ending */
     double a0;               /*!< \brief \f$a_0 = \frac{1}{\beta_N\Delta t^2}\f$ */
     double a1;               /*!< \brief \f$a_1 = \frac{\gamma_N}{\beta_N\Delta t}\f$ */
     double a2;               /*!< \brief \f$a_2 = \frac{1}{\beta_N\Delta t}\f$ */
     double a3;               /*!< \brief \f$a_3 = \frac{1}{2\beta_N\Delta t} - 1\f$ */
     double a4;               /*!< \brief \f$a_4 = \frac{\gamma_N}{\beta_N} - 1\f$ */
     double a5;               /*!< \brief \f$a_5 = \Delta_t\biggl(\frac{\gamma_N}{2\beta_N} - 1\biggr)\f$ */
     double a6;               /*!< \brief \f$a_6= \frac{1 - \gamma_N}{\Delta t}\f$ */
     double a7;               /*!< \brief \f$a_7 = \gamma_N\Delta t\f$ */

     /* Files where data are located */
     char* FileM;            /*!< \brief Stores the name of the file that contains the Mass Matrix */
     char* FileK;            /*!< \brief Stores the name of the file that contains the Stiffness Matrix */
     char* FileC;            /*!< \brief Stores the name of the file that contains the Damping Matrix */
     char* FileLV;            /*!< \brief Stores the name of the file that contains the Load Vector */
     char* FileCNodes;       /*!< \brief Stores the name of the file that contains the vector of coupling nodes. */
     char* FileData;         /*!< \brief Stores the name of the file that contains displacement, velocity and acceleration */
     char* FileOutput;       /*!< \brief Name of the file to store the output values of the process */

     /* Information regarding the type of communication */
     Remote_Machine_Info Remote; /*!< \brief Stores the data for the Remote site. \sa Remote_Machine_Info */
} AlgConst;

/**
 * \brief Definition of constant values and filenames that will be used during the Algorithm.
 *
 * This routine reads the values specified in a configuration file. In it, several constants like the order of the matrices,
 * the number of steps and file names where the Mass and Stiffness matrices are stored (the damping matrix is optional).
 *
 * \param[out] InitConst A structure that comprises of several constants.
 * \param[in] FileName Name of the configuration file.
 *
 * \post
 * - The size of the matrices will determine the memory that will be allocated when defining a MatrixVector type and also how
 * many elements will be read/written from/to the files.
 * - The number of steps must be equal to the number of rows of the file "DataFile".
 * - The values of the Newmark integration, PID and Rayleigh constants must be coherent/feasible. The algorithm will not perform checks
 * on them.
 * \sa RayleighConst, NewmarkConst and PIDValues.
 *
 */
void InitConstants( AlgConst *const InitConst, const char* FileName );

/**
 * \brief Frees the memory allocated during the InitConstants() routine.
 *
 * The memory allocated in the InitConstants() by the strdup() function is deallocated.
 * This includes basically the filenames and the IP, Port, Login and Password variables.
 * 
 * \param[out] InitConst Structure containing the data to be deallocated.
 *
 * \sa InitConstants Delete_ServerInformation.
 */
void Delete_InitConstants( AlgConst *const InitConst );

/**
 * \brief Reads the coupling nodes from a file.
 *
 * The coupling nodes are read from a file and stored sequentially in a dynamically
 * allocated array. The first number of the file must be always the number of 
 * coupling nodes to be readen.
 * 
 * \pre
 * - The file must be an ASCII file with the first value meaning the number
 * of nodes to be read.
 * - The datastructure Coupling_Nodes should not be initialised, since this is done
 * in this routine.
 *
 * \param[out] CNodes Data structure to store both: the number of coupling nodes and
 * a list of them.
 * \param[in] Filename The name of the file to be opened.
 *
 * \post CNodes must contain a list of the coupling nodes in increasing row order and the number of them.
 */
void Read_Coupling_Nodes( Coupling_Node *const CNodes, const int OrderSub, const double DeltaTSub, const char *Filename );

/**
 * \brief Construction of Proportional Viscous Damping Matrix using Rayleigh Damping.
 *
 * This routine calculates the Proportional Viscous Damping Matrix using Rayleigh Damping through the equation \f$\mathcal{C} = \alpha\mathcal{M} \cdot \beta \mathcal{K}\f$
 * (see Dynamics of Structures p. 234 ). The coefficients \f$\alpha\f$ and \f$\beta\f$ have to be provided. It makes use of the dlacpy_( ), dlascl_( ) and daxpy_( ) routines.
 *
 * \pre The matrices must be symmetrical and only the upper part of it will be referenced (lower part in FORTRAN routines)
 *
 * \param[in] Mass The Mass matrix.
 * \param[in] Stif The Stiffness matrix.
 * \param[in,out] Damp The Proportional Viscous Damping Matrix. As an input, only the size of the matrix is referenced, not its elements.
 * \param[in] Rayleigh It contains the values of the proportional coefficients \f$\alpha\f$ (\c Rayleigh.Alpha) and \f$\beta\f$ (\c Rayleigh.Beta).
 *
 * \post \c Damp is a symmetric matrix in general storage with only the upper part referenced (Lower part in FORTRAN routines).
 * It contains the result of \f$\mathcal{C} = \alpha \mathcal{M} \cdot \beta \mathcal{K}\f$.
 *
 * \sa MatrixVector and RayleighConst.
 */
void CalculateMatrixC( const MatrixVector *const Mass, const MatrixVector *const Stif, MatrixVector *const Damp, const RayleighConst *const Rayleigh );

/**
 * \brief Construction of Proportional Viscous Damping Matrix using Rayleigh Damping and sparse matrices in Intel's MKL CSR-\em three \em array \em variation format.
 *
 * This is the sparse version of the routine to calculate the Proportional Viscous Damping Matrix using Rayleigh Damping through the equation
 * \f$\mathcal{C} = \alpha \mathcal{M} \cdot \beta \mathcal{K}\f$ (see Dynamics of Structures p. 234 ). The coefficients \f$\alpha\f$ and \f$\beta\f$ have to be provided.
 * It makes use of the external functions scopy_(), sscal_() and mkl_scsradd(). The format is assumed to be the CSR-\em three \em array \em variation.
 *
 * \pre The matrices must be in Intel's MKL CSR-\em three \em array \em variation format in one based index. Due to the symmetry, only the upper part
 * (lower part in FORTRAN) is referenced.
 *
 * \param[in] Mass The Mass matrix.
 * \param[in] Stif The Stiffness matrix.
 * \param[in,out] Damp The Proportional Viscous Damping Matrix. As an input, only the size of the matrix is referenced, not its elements.
 * \param[in] Rayleigh It contains the values of the proportional coefficients \f$\alpha\f$ (Rayleigh.Alpha) and \f$\beta\f$ (Rayleigh.Beta).
 *
 * \post \c Damp is a symmetric matrix in Intel's MKL CSR-\em three \em array \em variation format.
 * It contains the result of \f$\mathcal{C} = \alpha \mathcal{M} \cdot \beta \mathcal{K}\f$.
 *
 * \sa MatrixVector and RayleighConst.
 */
void CalculateMatrixC_Sp( const Sp_MatrixVector *const Mass, const Sp_MatrixVector *const Stif, Sp_MatrixVector *const Damp, const RayleighConst *const Rayleigh );

/**
 * \brief Construction of the inverse of the Effective Stiffness Matrix using LAPACK routines.
 *
 * This routine calculates the so called Effective Stiffness Matrix through the equation \f$\mathcal{K}^{1}_e = \biggl[\mathcal{K} + \frac{1}{\beta\Delta\ t^2}\mathcal{M} + \frac{\gamma}{\beta\Delta t} \mathcal{C}\biggr]^{-1}\f$.
 * It makes use of the Add3Mat() to add the three matrices and the LAPACK routines dpotrf_() to calculate the Cholesky factorisation and dpotri_() to calculate the inverse. Note that the
 * matrix inversion is performed in double precision and afterwards it is converted to single precision. This is done in order to avoid numerical problems that arised when inverting the
 * matrix in single precision.
 *
 * \pre The matrices must be symmetrical, in general storage and only the upper part of it will be referenced (lower part in FORTRAN routines)
 *
 * \param[in,out] Keinv The inverse of the effective stiffness matrix. As an input, only the size of the matrix is referenced, not its elements.
 * \param[in] Mass The Mass matrix.
 * \param[in] Stif The Stiffness matrix.
 * \param[in] Damp The Viscous Damping matrix.
 * \param[in] Const The values are:
 * - \c Const.Alpha \f$ = 1\f$.
 * - \c Const.Beta \f$= \frac{1}{\beta\Delta\ t^2}\f$.
 * - \c Const.Gamma \f$= \frac{\gamma}{\beta\Delta t}\f$.
 *
 * \post \c Keinv is a symmetric matrix in general storage with only the upper part referenced (Lower part in FORTRAN routines). It contains the result of
 * \f$\mathcal{K}^{1}_e = \biggl[\mathcal{K} + \frac{1}{\beta\Delta\ t^2}\mathcal{M} + \frac{\gamma}{\beta\Delta t} \mathcal{C}\biggr]^{-1}\f$.
 *
 * \sa MatrixVector, Scalars and Add3Mat().
 */
void CalculateMatrixKeinv( MatrixVector *const Keinv, const MatrixVector *const Mass, const MatrixVector *const Damp, const MatrixVector *const Stif, const Scalars Const );

/**
 * \brief Construction of the inverse of the Effective Stiffness Matrix using PARDISO solver and dense matrices as input. It requires Intel MKL.
 *
 * This routine calculates the so called Effective Stiffness Matrix through the equation \f$\mathcal{K}^{1}_e = \biggl[\mathcal{K} + \frac{1}{\beta\Delta\ t^2}\mathcal{M} + \frac{\gamma}{\beta\Delta t}\mathcal{C}\biggr]^{-1}\f$.
 * It makes use of the Pardiso solver, converting the result of Add3Mat() to Intel's MKL CSR-\em three \em array \em variation format in order to compute the matrix inversion. Note that the
 * matrix inversion is performed in double precision and afterwards it is converted to single precision. This is done in order to avoid numerical problems that arised when inverting the
 * matrix in single precision. The Intel MKL libraries are required since it makes use of the mkl_sdnscsr() (in Dense_to_CSR()) routine and the Pardiso solver.
 *
 * \pre The matrices must be symmetrical, in general storage and only the upper part of it will be referenced (lower part in FORTRAN routines)
 *
 * \param[in,out] Keinv The inverse of the effective stiffness matrix. As an input, only the size of the matrix is referenced, not its elements.
 * \param[in] Mass The Mass matrix.
 * \param[in] Stif The Stiffness matrix.
 * \param[in] Damp The Viscous Damping matrix.
 * \param[in] Const The values are:
 * - \c Const.Alpha \f$ = 1\f$.
 * - \c Const.Beta \f$= \frac{1}{\beta\Delta\ t^2}\f$.
 * - \c Const.Gamma \f$= \frac{\gamma}{\beta\Delta t}\f$.
 *
 * \post \c Keinv is a symmetric matrix in general storage with only the upper part referenced (Lower part in FORTRAN routines). It contains the result of
 * \f$\mathcal{K}^{1}_e = \biggl[\mathcal{K} + \frac{1}{\beta\Delta\ t^2}\mathcal{M} + \frac{\gamma}{\beta\Delta t}\mathcal{C}\biggr]^{-1}\f$.
 *
 * \sa MatrixVector, Scalars, Add3Mat(), Dense_to_CSR() and Generate_IdentityMatrix().
 */
void CalculateMatrixKeinv_Pardiso( MatrixVector *const Keinv, const MatrixVector *const Mass, const MatrixVector *const Damp, const MatrixVector *const Stiff, const Scalars Const );

/**
 * \brief Construction of the inverse of the Effective Stiffness Matrix using PARDISO solver and sparse matrices in Intel's MKL CSR-\em three \em array \em variation format as input.
 * It requires Intel MKL.
 *
 * This routine calculates the so called Effective Stiffness Matrix through the equation \f$\mathcal{K}^{1}_e = \biggl[\mathcal{K} + \frac{1}{\beta\Delta\ t^2}\mathcal{M} + \frac{\gamma}{\beta\Delta t}\mathcal{C}\biggr]^{-1}\f$.
 * It makes use of the Pardiso solver and the matrices in Intel's MKL CSR-\em three \em array \em variation format in order to compute the matrix inversion. Note that the
 * matrix inversion is performed in double precision and afterwards it is converted to single precision. This is done in order to avoid numerical problems that arised when inverting the
 * matrix in single precision. The Intel MKL libraries are required since it makes use of the mkl_sdnscsr() (in Dense_to_CSR()) routine and the Pardiso solver.
 *
 * \pre
 * - The \c Keinv matrix must be symmetrical and only the upper part of it will be referenced (lower part in FORTRAN routines). It has to be initialised with the
 * Init_MatrixVector() routine.
 * - The Mass, Viscous and Stiffness matrices must be in Intel's MKL CSR-\em three \em array \em variation format in one based index. Due to their symmetry, only the upper part
 * (lower part in FORTRAN) is referenced.
 *
 * \param[in,out] Keinv The inverse of the effective stiffness matrix. As an input, only the size of the matrix is referenced, not its elements.
 * \param[in] Mass The Mass matrix.
 * \param[in] Stif The Stiffness matrix.
 * \param[in] Damp The Viscous Damping matrix.
 * \param[in] Const The values are:
 * - \c Const.Alpha \f$ = 1\f$.
 * - \c Const.Beta \f$= \frac{1}{\beta\Delta\ t^2}\f$.
 * - \c Const.Gamma \f$= \frac{\gamma}{\beta\Delta t}\f$.
 *
 * \post \c Keinv is a symmetric matrix in Intel's MKL CSR-\em three \em array \em variation format with only the upper part referenced (Lower part in FORTRAN routines).
 * It contains the result of \f$\mathcal{K}^{1}_e = \biggl[\mathcal{K} + \frac{1}{\beta\Delta\ t^2}\mathcal{M} + \frac{\gamma}{\beta\Delta t}\mathcal{C}\biggr]^{-1}\f$.
 *
 * \sa Sp_MatrixVector, Scalars, Add3Mat_Sparse(), Dense_to_CSR() and Generate_IdentityMatrix().
 */
void CalculateMatrixKeinv_Pardiso_Sparse( MatrixVector *const Keinv, const Sp_MatrixVector *const Mass, const Sp_MatrixVector *const Damp, const Sp_MatrixVector *const Stiff, const Scalars Const );

/**
 * \brief Generation of a Identity Matrix.
 *
 * The identity matrix with sizes equal to \c Rows and \c Cols (number of rows and columns respectively) is generated. The output format is in general storage.
 * 
 * \ore The number of rows must be equal to the number of columns, otherwise the routine will fail
 *
 * \param[in] Rows The number of rows.
 * \param[in] Cols The number of columns.
 *
 * \return A identity matrix of \f$Size= Rows*Cols\f$ in general storage.
 */
MatrixVector Generate_IdentityMatrix( int Rows, int Cols );

/**
 * \brief Construction of the coupling matrix.
 *
 * This routine copies the values in the coupling positions of a symmetric matrix constructing a sub-matrix of \f$ Size = Number~of~coupling~nodes^2\f$. Since this matrix has
 * to be sent to the experimental facility in order to perform the sub-stepping process, it is treated as a full matrix (although it is symmetrical). Therefore, all elements
 * are copied if \f$Number~of~coupling~nodes \geq 1\f$. For example, given the symmetric matrix \f$\mathcal{A}\f$ and coupling positions in 2 and 5, the resulting \f$\mathcal{A}_{Couple}\f$ would
 * be as follows:
 *
 * \f[ \mathcal{A} = \begin{pmatrix}
 *   1 & -1 & 3  & 4  & 5\\
 *   * & \mathbf{5}  & 4  & -3 & \mathbf{2}\\
 *   * & *  & 3  & 6  & 7\\
 *   * & *  & *  & 11 & 4\\
 *   * & *  & *  & *  & \mathbf{8}\\
 * \end{pmatrix}
 * \Longrightarrow \mathcal{A}_{Couple} = \begin{pmatrix}
 * 5 & 2\\
 * 2 & 8\\
 * \end{pmatrix}\f]
 *
 * \pre
 * - The input matrix has to be a symmetrical matrix containing at least the upper part (lower part in FORTRAN routines) in
 * general storage.
 * - The coupling nodes are assumed to be in increasing order of rows and in one based index.
 * - MatCouple must be properly initialised and its size should be \f$Size \geq Number~of~coupling~nodes^2\f$.
 *
 * \param[in] Mat The matrix that will be decoupled.
 * \param[out] MatCouple The matrix where the coupling nodes are saved.
 * \param[in] CNodes Structure containing the coupling nodes in increasing order of rows.
 *
 * \post \c MatCouple is a symmetrical matrix \f$Size = Number~of~coupling~nodes^2\f$ in general storage that contains the values in the coupling position of the specified matrix.
 *
 * \sa MatrixVector and CNodes.
 *
 */
void BuildMatrixXc( const MatrixVector *const Mat, double *MatCouple, const Coupling_Node *const CNodes );

/**
 * \brief Construction of the non-coupling part of a given matrix.
 *
 * This routine copies the non-coupling values of a column with coupling degrees of freedom of the symmetric matrix Mat,
 * constructing a matrix of \f$Size = (Order - Order_C)\cdot Order_C\f$. Where:
 *
 * - \f$Order\f$ is the number of rows and columns of the input matrix.
 * - \f$Order_C\f$ is the number of coupling degrees of freedom.
 *
 * For example, given the given the symmetric matrix \f$\mathcal{A}\f$ and coupling positions in 2 and 5, the resulting \f$\mathcal{A}_{cm}\f$ would
 * be as follows:
 * 
 * \f[ \mathcal{A} = \begin{pmatrix}
 *   1 & \mathbf{-1} & 3  & 4  & \mathbf{5}\\
 *   * & 5  & \mathbf{4}  & \mathbf{-3} & 2\\
 *   * & *  & 3  & 6  & \mathbf{7}\\
 *   * & *  & *  & 11 & \mathbf{4}\\
 *   * & *  & *  & *  & 8\\
 * \end{pmatrix}
 * \Longrightarrow \mathcal{A}_{cm} = \begin{pmatrix}
 * -1 & 5\\
 *  4 & 7\\
 * -3 & 4\\
 * \end{pmatrix}\f]
 *
 * \pre
 * - The input matrix has to be a symmetrical matrix containing at least the upper part (lower part in FORTRAN routines) in
 * general storage.
 * - \f$Order > Order_C\f$.
 * - \c Matcm must be of \f$Size = (Order - Order_C)\cdot Order_C\f$.
 * - The coupling nodes are assumed to be in increasing order of rows and in one based index.
 *
 * \param[in] Mat The matrix that will be decoupled.
 * \param[in,out] Matcm The matrix where the non-coupling elemets of a column with a coupling node are stored.
 * \param[in] CNodes Structure containing the coupling nodes in increasing order of rows.
 *
 * \post \c Matcm is a general matrix of \f$Size = (Order - Order_C)\cdot Order_C\f$ with the non-coupling elements of the columns with coupling nodes.
 *
 * \sa MatrixVector.
 */
void BuildMatrixXcm( const MatrixVector *const Mat, MatrixVector *const Matcm,  const Coupling_Node *const CNodes );

void Generate_LoadVectorForm( MatrixVector *const LoadVector, int *DOF );
int* Get_Excited_DOF( const ConfFile *const Config, const char *Expression );
#endif /* INITIATION_H_ */
