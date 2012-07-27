/**
 * \file Initiation.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 10th of August 2011
 * \todo Add support for packaged storages to decrease the memory use.
 * \todo Make the routines BuildMatrixXc() and BuildMatrixXcm() to be able to handle non consecutive coupling nodes.
 *
 * \brief Prototypes of the functions used in the Initiation phase.
 *
 * This file defines the routines used during the initiation phase of the substructure algorithm: construction of the Proportional Viscous
 * Damping Matrix, the Effective Mass Matrix and the Gain Matrix and their non-coupling and coupling part.
 */

#ifndef INITIATION_H_
#define INITIATION_H_

#include "MatrixVector.h"
/**
 * \brief Structure to handle the PID error compensator.
 *
 * This structure handles the PID error compensator type, by storing the three parameters
 * that characterise a PID: the proportional, the integral and derivative part.
 */
typedef struct {
     float P; /*!< \brief Proportional constant */
     float I; /*!< \brief Integral constant */
     float D; /*!< \brief Derivative constant */
} PIDValues;

/**
 * \brief Structure to store the Rayleigh Viscous Damping.
 *
 * This struct consists of two variables to store the parameters required to build
 * the Viscous Damping Matrix using Rayleigh, that is \f$[C] = \alpha[M] + \beta[K]\f$, where:
 * \li \f$[C]\f$ is the Viscous Damping Matrix,
 * \li \f$[M]\f$ is the Mass Matrix,
 * \li \f$[K]\f$ is the Stiffness Matrix,
 * \li and \f$\alpha\f$ and \f$beta\f$ are the parameters that multiply the mass and stiffness matrix respectively.
 */
typedef struct {
     float Alpha; /*!< \brief First constant.*/
     float Beta;  /*!< \brief Second constant.*/
} RayleighConst;

/**
 * \brief Structure to handle the Newmark Constants.
 *
 * This structure handles the Newmark constants, by storing the two parameters required by the Newmark-\f$beta\f$
 * method: \f$\gamma\f$ and \f$\beta\f$.
 */
typedef struct {
     float Gamma; /*!< \brief First constant.*/
     float Beta;  /*!< \brief Second constant.*/
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

     int OrderC;              /*!< \brief Number of coupling DOF */
     int PosCouple;           /*!< \brief Position of the coupling DOF */

     unsigned int Nstep;               /*!< \brief Number of steps */
     int Use_Absolute_Values; /*!< \brief Variable to control whether to use absolute values in the equation of motion or relative values. Affects how the input load is calculated */
     float Delta_t;           /*!< \brief Time increment \f$\Delta t\f$ */

     RayleighConst Rayleigh;  /*!< \brief Stores Rayleigh Constants alpha (\c Rayleigh.Alpha or \f$\alpha_R\f$) and beta (\c Rayleigh.Beta or \f$\beta_R\f$) */
     NewmarkConst Newmark;    /*!< \brief Stores Newmark Constants gamma (\c Newmark.Gamma or \f$\gamma_N\f$) and (\c Newmark.Beta or \f$\beta_N\f$) */
     PIDValues PID;           /*!< \brief Stores proportional (\c PID.P), integral (\c PID.I) and derivative (\c PID.D) part. */

     /* Constants for Calculations */
     float Const1;           /*!< \brief \f$Const_1=\beta_N\Delta t^2\f$ */
     float Const2;           /*!< \brief \f$Const_2 = (0.5 - 2*\beta_N + \gamma_N)*\Delta t^2\f$ */
     float Const3;           /*!< \brief \f$Const_3 = (0.5 + \beta_N - \gamma_N)*\Delta t^2\f$ */

     /* Constants for Step ending */
     float a0;               /*!< \brief \f$a_0 = \frac{1}{\beta_N\Delta t^2}\f$ */
     float a1;               /*!< \brief \f$a_1 = \frac{\gamma_N}{\beta_N\Delta t}\f$ */
     float a2;               /*!< \brief \f$a_1 = \frac{1}{\beta_N\Delta t}\f$ */
     float a3;               /*!< \brief \f$a_3 = \frac{1}{2\beta_N\Delta t} - 1\f$ */
     float a4;               /*!< \brief \f$a_4 = \frac{\gamma_N}{\beta_N} - 1\f$ */
     float a5;               /*!< \brief \f$a_5 = \Delta_t\biggl(\frac{\gamma_N}{2\beta_N} - 1\biggr)\f$ */
     float a6;               /*!< \brief \f$a_6= \frac{1 - \gamma_N}{\Delta t}\f$ */
     float a7;               /*!< \brief \f$a_7 = \gamma_N\Delta t\f$ */

     /* Files where data are located */
     const char* FileM;       /*!< \brief Stores the name of the file that contains the Mass Matrix */
     const char* FileK;       /*!< \brief Stores the name of the file that contains the Stiffness Matrix */
     const char* FileC;       /*!< \brief Stores the name of the file that contains the Damping Matrix */
     const char* FileLVector; /*!< \brief Stores the name of the file that contains the vector used for the load. This vector usually contains 1 and 0 */
     const char* FileData;    /*!< \brief Stores the name of the file that contains displacement, velocity and acceleration */

     /* Information regarding the type of communication */
     int Type_Protocol;           /*!< \brief Identifies the protocol to be used. 1 = Custom, 2 = NSEP, 3 = OpenFresco */
} AlgConst;

/**
 * \brief Definition of constant values and filenames that will be used during the Algorithm.
 *
 * This routine is intended to be used as a configuration file. In it, several constants like the order of the matrices,
 * the number of steps and file names where the Mass and Stiffness matrices are stored (the damping matrix is optional). Although
 * this function could be implemented reading this values from a configuration file, thus skiping the need to compile it each time one
 * of these values is changed, the actual approach has some benefits from the performance point of view. This is specially important in
 * the case of loops, since the compiler knows before hand the size of them.
 *
 * \param[out] AConst A structure that comprises of several constants.
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
void InitConstants( AlgConst *const AConst );

/**
 * \brief Identify the communication protocol to be used
 *
 * The communication protocol to be used is identified, and a proper return value is given. It makes use of the function
 * Get_Server_Information().
 *
 * \pre The first line of the file \c Connection.txt must contain the desired protocol type
 *
 * \return 
 * - 0 if the desired protocol is of type \c Custom.
 * - 1 if the desired protocol is of type \c PNSE.
 * - 2 if the desired protocol is of type \c OpenFresco.
 * - -1 if the desired protocol is not recognised.
 */
int Get_Type_Protocol( void );

/**
 * \brief Construction of Proportional Viscous Damping Matrix using Rayleigh Damping.
 *
 * This routine calculates the Proportional Viscous Damping Matrix using Rayleigh Damping through the equation \f$[C] = \alpha [M] \cdot \beta [K]\f$
 * (see Dynamics of Structures p. 234 ). The coefficients \f$\alpha\f$ and \f$\beta\f$ have to be provided. It makes use of the dlacpy_( ), dlascl_( ) and daxpy_( ) routines.
 *
 * \pre The matrices must be symmetrical and only the upper part of it will be referenced (lower part in FORTRAN routines)
 *
 * \param[in] Mass The Mass matrix.
 * \param[in] Stif The Stiffness matrix.
 * \param[in,out] Damp The Proportional Viscous Damping Matrix. As an input, only the size of the matrix is referenced, not its elements.
 * \param[in] Rayleigh It contains the values of the proportional coefficients \f$\alpha\f$ (Rayleigh.Alpha) and \f$\beta\f$ (Rayleigh.Beta).
 *
 * \post \c Damp is a symmetric matrix in general storage with only the upper part referenced (Lower part in FORTRAN routines).
 * It contains the result of \f$[C] = \alpha [M] \cdot \beta [K]\f$.
 *
 * \sa MatrixVector and RayleighConst.
 */
void CalculateMatrixC( const MatrixVector *const Mass, const MatrixVector *const Stif, MatrixVector *const Damp, const RayleighConst *const Rayleigh );

/**
 * \brief Construction of the inverse of the Effective Mass Matrix.
 *
 * This routine calculates the so called Effective Mass Matrix through the equation \f$M_{e,inv} = [M + \gamma\Delta tC + \beta\Delta t^2 K]^{-1}\f$. It makes use
 * of the Add3Mat( ) to add the three matrices and the LAPACK routines dpotrf_( ) to calculate the Cholesky factorisation and pdpotri_( ) to calculate the inverse.
 *
 * \pre The matrices must be symmetrical and only the upper part of it will be referenced (lower part in FORTRAN routines)
 *
 * \param[in,out] Meinv The inverse of the effective mass matrix. As an input, only the size of the matrix is referenced, not its elements.
 * stiffness, mass and damping matrices). It will store the upper part of the Rayleigh Damping Matrix in general storage as an output.
 * \param[in] Mass The Mass matrix.
 * \param[in] Stif The Stiffness matrix.
 * \param[in] Damp The Viscous Damping matrix.
 * \param[in] Const The values are:
 * - \c Const.Alpha \f$ = 1\f$.
 * - \c Const.Beta \f$= \gamma\Delta t\f$.
 * - \c Const.Gamma \f$ = \beta\Delta t^2\f$.
 *
 * \post \c Meinv is a symmetric matrix in general storage with only the upper part referenced (Lower part in FORTRAN routines). It contains the result of
 * \f$M_{e,inv} = [M + \gamma\Delta tC + \beta\Delta t^2 K]^{-1}\f$.
 *
 * \sa MatrixVector, Scalars and Add3Mat( ).
 */
void CalculateMatrixKeinv( MatrixVector *const Meinv, const MatrixVector *const Mass, const MatrixVector *const Damp, const MatrixVector *const Stif, const Scalars Const );

/**
 * \brief Construction of the Gain Matrix
 *
 * This routine calculates the Gain Matrix \f$G = \beta\Delta t^2 M_{e,inv}\f$. It makes use of the dlacpy_( ) and dlascl_( ) LAPACK routines.
 *
 * \pre The matrices must be symmetrical and only the upper part of it will be referenced (lower part in FORTRAN routines)
 *
 * \param[in,out] Gain The Gain matrix. As an input, only the size of the matrix is referenced, not its elements.
 * \param[in] Meinv The inverse of the effective mass matrix.
 * \param[in] Const On entry \f$Const = \beta\Delta t^2\f$.
 *
 * \post \c Gain is a symmetric matrix in general storage with only the upper part referenced (Lower part in FORTRAN routines). It contains the result of \f$G = \beta\Delta t^2 M_{e,inv}\f$.
 *
 * \sa MatrixVector.
 */
void CalculateMatrixG( MatrixVector *const Gain, const MatrixVector *const Meinv, float Const );

/**
 * \brief Construction of the coupling nodes
 *
 * This routine copies the coupling nodes starting at the Coupling Position (PosCpl) constructing a matrix of orderc. Since this matrix has to be sent
 * to the experimental facility in order to perform the sub-stepping process, it is treated as a full matrix (although it is symmetrical). Therefore, all elements are copied
 * if \f$orderc \geq 1\f$.
 *
 * \pre
 * - The input matrix has to be a symmetrical matrix containing at least the upper part (Lower part in FORTRAN routines) in
 * general storage.
 * - The coupling nodes are assumed to be in consecutive rows (the first row is given by PosCpl).
 * - MatCouple must be properly initialised and its size should be \f$size \geq OrderC\cdot OrderC\f$.
 *
 * \param[in] Mat The matrix that will be decoupled.
 * \param[out] MatCouple The matrix where the coupling nodes are saved.
 * \param[in] PosCpl The position of the first coupling node.
 * \param[in] OrderC The number of coupling nodes. In case that \f$OrderC > 1\f$, the routine assumes that they are consecutive.
 *
 * \post \c MatCouple is a symmetrical matrix \f$OrderC\cdot OrderC\f$ in general storage that contains the coupling nodes.
 *
 * \sa MatrixVector.
 *
 */
void BuildMatrixXc( const MatrixVector *const Mat, float *MatCouple, const int PosCpl, const int OrderC );

/**
 * \brief Construction of the non-coupling part of the row where the Coupling node is located.
 *
 * This routine copies the non-coupling part of the matrix Mat, constructing a matrix of (Order - Order)xOrderC. This means that
 * the rows were the coupling nodes are located, starting at PosClp and ending at PosClp+OrderC, are copied to the new matrix.
 *
 * \pre The input matrix has to be a symmetrical matrix containing at least the upper part (Lower part in FORTRAN routines) in
 * general storage. The coupling nodes are assumed to be in consecutive rows (the first row is given by PosCpl).
 * VecXcm must be of \f$size = (Order - OrderC)\cdot OrderC\f$, where:
 * - \e Order is the number of rows and columns of the input matrix.
 * - \e OrderC is the number of coupling degrees of freedom
 *
 * \param[in] Mat The matrix that will be decoupled.
 * \param[in,out] VecXcm The matrix where the non-coupling elemets of a row with a couping node are stored.
 * \param[in] PosCpl The position of the first coupling node.
 * \param[in] OrderC The number of coupling nodes. In case that \f$OrderC > 1\f$, the routine assumes that they are consecutive.
 *
 * \post \c VecXcm is a general matrix of size \f$size = (Order - OrderC)\cdot OrderC\f$ with the non-coupling elements of the row with a coupling node.
 *
 * \sa MatrixVector.
 */
void BuildMatrixXcm( const MatrixVector *const Mat, MatrixVector *const VecXcm, const int PosCpl, const int OrderC );


#endif /* INITIATION_H_ */
