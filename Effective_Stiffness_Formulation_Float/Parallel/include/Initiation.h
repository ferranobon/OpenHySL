/*
 * Initiation.h
 *
 *  Created on: 22/07/2011
 *      Author: ferran
 */

#ifndef INITIATION_H_
#define INITIATION_H_

#include "Netlib.h"
#include "PMatrixVector.h"
#include "Send_Receive_Data.h"

/**
 * \brief Structure to handle the PID error compensator
 *
 * This structure handles the PID error compensator type, by storing the three parameters
 * that characterise a PID: the proportional, the integral and derivative part
 */
typedef struct {
     float P; /*!< Proportional constant */
     float I; /*!< Integral constant */
     float D; /*!< Derivative constant */
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
     float Alpha;
     float Beta;
} RayleighConst;

/**
 * \brief Structure to handle the Newmark Constants.
 *
 * This structure handles the Newmark constants, by storing the two parameters required by the Newmark-\f$beta\f$
 * method: \f$\gamma\f$ and \f$\beta\f$.
 */
typedef struct {
     float Gamma;
     float Beta;
} NewmarkConst;

typedef struct {
     int Rows;
     int Cols;
} GridInfo;

/**
 * \brief Structure to store the coupling nodes.
 * 
 * This structure is used in order to store the coupling nodes that will be used
 * during a test. The nodes are stored sequentially and in increasing order.
 */
typedef struct {
     int *Array;  /*!< \brief Array containing the coupling nodes */
     int Order;   /*!< \brief Number of coupling nodes */
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

     GridInfo ProcessGrid;    /* Stores the information regarding the desired process grid */
     GridInfo BlockSize;      /* Stores the information regarding the desired block size */

     int Order;               /*!< \brief Order of the matrices */

     int Nstep;               /*!< \brief Number of steps */

     int Use_Absolute_Values; /*!< \brief Variable to control whether to use absolute values in the equation of motion or relative values. Affects how the input load is calculated */

     float Delta_t;          /*!< \brief Time increment \f$\Delta t\f$ */

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

void InitConstants( AlgConst *const AConst, const int size );
void BroadcastConfFile( AlgConst *const InitConst );
void Delete_InitConstants( AlgConst *const InitConst );
void Read_Coupling_Nodes( Coupling_Node *const CNodes, const char *Filename );
void BroadCast_Coupling_Nodes( Coupling_Node *const CNodes );
void CalculateMatrixC( PMatrixVector *const Mass, PMatrixVector *const Stif, PMatrixVector *const Damp, RayleighConst Rayleigh );
void CalculateMatrixKeinv( PMatrixVector *const Keinv, PMatrixVector *const Mass, PMatrixVector *const Damp, PMatrixVector *const Stif, Scalars Const );
void CalculateMatrixG( PMatrixVector *const Gain, PMatrixVector *const Keinv, float Const );
void BuildMatrixXc( MPI_Comm Comm, PMatrixVector *const Mat, float *MatCouple, const Coupling_Node *const CNodes );
void BuildMatrixXcm( MPI_Comm Comm, PMatrixVector *const Mat, PMatrixVector *const VecXcm, const Coupling_Node *const CNodes );

#endif /* INITIATION_H_ */