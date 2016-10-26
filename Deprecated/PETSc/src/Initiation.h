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

#include <petscmat.h>

#include "MatrixVector.h"
#include "Conf_Parser.h"
#include "Send_Receive_Data.h"

/**
 * \brief Structure to handle the PID error compensator.
 *
 * This structure handles the PID error compensator type, by storing the three parameters
 * that characterise a PID: the proportional, the integral and derivative part.
 */
typedef struct {
     PetscScalar P; /*!< \brief Proportional constant */
     PetscScalar I; /*!< \brief Integral constant */
     PetscScalar D; /*!< \brief Derivative constant */
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
     PetscScalar Alpha; /*!< \brief First constant.*/
     PetscScalar Beta;  /*!< \brief Second constant.*/
} RayleighConst;

/**
 * \brief Structure to handle the Newmark Constants.
 *
 * This structure handles the Newmark constants, by storing the two parameters required by the Newmark-\f$beta\f$
 * method: \f$\gamma\f$ and \f$\beta\f$.
 */
typedef struct {
     PetscScalar Gamma; /*!< \brief First constant.*/
     PetscScalar Beta;  /*!< \brief Second constant.*/
} NewmarkConst;

/**
 * \brief Structure to store the coupling nodes.
 * 
 * This structure is used in order to store the coupling nodes that will be used
 * during a test. The nodes are stored sequentially and in increasing order.
 */
typedef struct {
     PetscInt *Array;  /*!< \brief Array containing the coupling nodes */
     PetscInt Order;   /*!< \brief Number of coupling nodes */
} Coupling_Node;

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

     PetscInt Order;               /*!< \brief Order of the matrices */

     PetscInt Nstep;      /*!< \brief Number of steps */

     PetscInt Use_Absolute_Values; /*!< \brief Variable to control whether to use absolute values in the equation of motion or relative values. Affects how the input load is calculated */

     PetscScalar Delta_t;           /*!< \brief Time increment \f$\Delta t\f$ */

     RayleighConst Rayleigh;  /*!< \brief Stores Rayleigh Constants alpha (\c Rayleigh.Alpha or \f$\alpha_R\f$) and beta (\c Rayleigh.Beta or \f$\beta_R\f$) */
     NewmarkConst Newmark;    /*!< \brief Stores Newmark Constants gamma (\c Newmark.Gamma or \f$\gamma_N\f$) and (\c Newmark.Beta or \f$\beta_N\f$) */
     PIDValues PID;           /*!< \brief Stores proportional (\c PID.P), integral (\c PID.I) and derivative (\c PID.D) part. */

     /* Constants for Calculations */
     PetscScalar Const1;           /*!< \brief \f$Const_1=\beta_N\Delta t^2\f$ */
     PetscScalar Const2;           /*!< \brief \f$Const_2 = (0.5 - 2*\beta_N + \gamma_N)*\Delta t^2\f$ */
     PetscScalar Const3;           /*!< \brief \f$Const_3 = (0.5 + \beta_N - \gamma_N)*\Delta t^2\f$ */

     /* Constants for Step ending */
     PetscScalar a0;               /*!< \brief \f$a_0 = \frac{1}{\beta_N\Delta t^2}\f$ */
     PetscScalar a1;               /*!< \brief \f$a_1 = \frac{\gamma_N}{\beta_N\Delta t}\f$ */
     PetscScalar a2;               /*!< \brief \f$a_1 = \frac{1}{\beta_N\Delta t}\f$ */
     PetscScalar a3;               /*!< \brief \f$a_3 = \frac{1}{2\beta_N\Delta t} - 1\f$ */
     PetscScalar a4;               /*!< \brief \f$a_4 = \frac{\gamma_N}{\beta_N} - 1\f$ */
     PetscScalar a5;               /*!< \brief \f$a_5 = \Delta_t\biggl(\frac{\gamma_N}{2\beta_N} - 1\biggr)\f$ */
     PetscScalar a6;               /*!< \brief \f$a_6= \frac{1 - \gamma_N}{\Delta t}\f$ */
     PetscScalar a7;               /*!< \brief \f$a_7 = \gamma_N\Delta t\f$ */

     /* Files where data are located */
     char* FileM;       /*!< \brief Stores the name of the file that contains the Mass Matrix */
     char* FileK;       /*!< \brief Stores the name of the file that contains the Stiffness Matrix */
     char* FileC;       /*!< \brief Stores the name of the file that contains the Damping Matrix */
     char* FileLVector; /*!< \brief Stores the name of the file that contains the vector used for the load. This vector usually contains 1 and 0 */
     char* FileCNodes;  /*!< \brief Stores the name of the file that contains the vector of coupling nodes. */
     char* FileData;    /*!< \brief Stores the name of the file that contains displacement, velocity and acceleration */
     char* FileOutput;    /*!< \brief Name of the file to store the output values of the process */

     /* Information regarding the type of communication */
     Remote_Machine_Info Remote;
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

void BroadcastConfFile( AlgConst *const InitConst );

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

void Read_Coupling_Nodes( Coupling_Node *const CNodes, const char *Filename );
void BroadCast_Coupling_Nodes( Coupling_Node *const CNodes );

void CalculateMatrixC( const Mat *const Mass, const Mat *const Stiff, Mat *const Damp, const RayleighConst *const Rayleigh );

void EffK_Caclulate_Keinv( Mat Mass, Mat Stiff, Mat Damp, Mat Keinv, const Scalars Const );

void BuildMatrixXc( MPI_Comm Comm, Mat Matrix, PetscScalar *MatCouple, const Coupling_Node *const CNodes );
void BuildMatrixXcm( MPI_Comm Comm, Mat Matrix, Mat MatXcm, const Coupling_Node *const CNodes );
PetscScalar* GetMat_Value( MPI_Comm Comm, Mat Matrix, PetscInt NumRows, PetscInt *Rows, PetscInt NumCols, PetscInt *Cols );
PetscInt GetOwner_Position( MPI_Comm Comm, Mat Matrix, PetscInt Row );

#endif
