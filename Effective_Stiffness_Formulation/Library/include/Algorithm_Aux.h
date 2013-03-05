#ifndef _ALGORITHM_AUX_H_
#define _ALGORITHM_AUX_H_

#include "Initiation.h"         /* For Rayleigh_t */
#include "Error_Compensation.h" /* For PID_t */
#include "Conf_Parser.h"        /* For ConfFile_t */

/**
 * \brief Structure to handle the time integration constants.
 *
 * This structure handles the Newmark constants, by storing the two parameters required by the Newmark-\f$beta\f$
 * method: \f$\gamma\f$ and \f$\beta\f$.
 */
typedef struct TIntegration{
     double Gamma; /*!< \brief First constant.*/
     double Beta;  /*!< \brief Second constant.*/
} TIntegration_t;

/**
 * \brief Structure to wrap constants and filenames.
 *
 * This structure stores several constants that will be used in different parts of the substructure algorithm, such as
 * the order of the matrices, the position of the coupling nodes, number of steps filenames and other contants. It act
 * as a configuration file and it has to be modified by hand before compiling.
 *
 * \sa Rayleigh_t, Newmark_t, PID_t and Algorithm_Init().
 */
typedef struct AlgConst{

     int Order;               /*!< \brief Order of the matrices */
     int OrderSub;            /*!< \brief Number of substructures */

     unsigned int NStep;      /*!< \brief Number of steps */
     unsigned int NSubstep;   /*!< \brief Number of sub-steps */

     int Use_Absolute_Values; /*!< \brief Variable to control whether to use absolute values in the equation of motion
			       * or relative values. Affects how the input load is calculated
			       */

     int Read_Sparse;         /*!< \brief The matrices will be read using a sparse format. \sa
			       * MatrixVector_From_File_Sp2Dense() and MatrixVector_From_File_Sp().
			       */
     int Use_Sparse;          /*!< \brief Sparse routines will be used instead of dense ones. Only available with the
			       * support of the Intel MKL libraries.
			       * \sa MatrixVector_From_File_Sp(), CalculateMatrixC_Sp(),
			       * CalculateMatrixKeinv_Pardiso_Sparse(), Calc_Input_Load_AbsValues_Sparse(),
			       * Calc_Input_Load_RelValues_Sparse() and EffK_Calc_Effective_Force_Sparse().
			       */
     int Read_LVector;        /*!< \brief Read the load vector instead of generating it. */
     int *ExcitedDOF;         /*!< \brief Integer string containing wich elements of the external load will be applied.
			       * \sa AlgorithM_GetExcitedDOF().
			       */
     double Delta_t;          /*!< \brief Time increment \f$\Delta t\f$ */
     double DeltaT_Sub;       /*!< \brief Time increment for the sub-stepping process */

     Rayleigh_t Rayleigh;     /*!< \brief Stores Rayleigh Constants alpha (\c Rayleigh.Alpha or \f$\alpha_R\f$) and beta
			       * (\c Rayleigh.Beta or \f$\beta_R\f$)
			       */
     TIntegration_t Newmark;  /*!< \brief Stores Newmark Constants gamma (\c Newmark.Gamma or \f$\gamma_N\f$) and (\c
			       * Newmark.Beta or \f$\beta_N\f$)
			       */
     PID_t PID;               /*!< \brief Stores proportional (\c PID.P), integral (\c PID.I) and derivative (\c PID.D)
			       * part.
			       */

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
     char* FileLV;           /*!< \brief Stores the name of the file that contains the Load Vector */
     char* FileCNodes;       /*!< \brief Stores the name of the file that contains the vector of coupling nodes. */
     char* FileData;         /*!< \brief Stores the name of the file that contains displacement, velocity and
			      * acceleration */
     char* FileOutput;       /*!< \brief Name of the file to store the output values of the process */
} AlgConst_t;


/**
 * \brief Definition of constant values and filenames that will be used during the Algorithm.
 *
 * This routine reads the values specified in a configuration file. In it, several constants like the order of the
 * matrices, the number of steps and file names where the Mass and Stiffness matrices are stored (the damping matrix is
 * optional).
 *
 *
 * \pre The values of the Newmark integration, PID and Rayleigh constants must be coherent/feasible. The algorithm will
 * not perform checks on them.
 *
 * \param[in] FileName Name of the configuration file.
 * \param[out] InitConst A structure that comprises of several constants.
 *
 * \post
 * - The size of the matrices will determine the memory that will be allocated when defining a \c MatrixVector_t type
 * and also also how many elements will be read/written from/to the files.
 *- The number of steps must be equal to the number of rows of the file "DataFile".
 *
 * \sa AlgoConst_t, PID_t, Rayleigh_t and TIntegration_t.
 */
void Algorithm_Init( const char *FileName, AlgConst_t *const InitConst );

/**
 * \brief Frees the memory allocated during the Algorithm_Init() routine.
 *
 * \pre InitConst must be properly initialised through Algorithm_Init().
 * 
 * \param[out] InitConst Structure containing the data to be deallocated.
 *
 * \sa AlgConst_t.
 */
void Algorithm_Destroy( AlgConst_t *const InitConst );

/**
 * \brief Transforms a string with the desired DOFs to be excited into an integer array.
 *
 * \pre
 * - The first element of the list should contain the number of DOFs per node.
 * - The next \f$n\f$ entries are the degrees of freedom that are being excited (1) or not (0) by the external load.
 * - The order of which DOFs will have an external load applied must remain constant and keep the same patern.
 * - The rest of \f$\frac{Order}{n}\f$ must be equal to 0, where \f$Order\f$ is the length of the external load vector.
 * - \c Config must be properly initialised through ConfFile_Create() and have entries stored through ConfFile_ReadFile().
 * - \c Expression must be a valid entry in \c Config.
 *
 * Transforms a string with the desired sequence of DOFs that will have an external load applied. The string consist of
 * a first element that determines how many values will follow, \f$n\f$. The next \f$n\f$ values are the pattern and
 * they identify if a degree of freedom is excited (1) or not (0). The pattern will be applied as many times as \f$N =
 * \frac{Order}{n}\f$, where \f$Order\f$ is the length of the external load vector.
 *
 * \b Example \b 1:
 *
 * <tt>6 1 1 1 0 0 0 </tt>
 *
 * The first three elements of the external load vector will have a value \f$u_g*1.0\f$ while the other three will be
 * set to 0.0. The pattern will repeat for every 6 elements in the load vector.
 *
 * \b Example \b 2:
 *
 * <tt>4 1 0 1 1</tt>
 *
 * Only the second element of every four positions in the external load vector will be set to 0.0.
 *
 * \param[in] Config Structure containing the entries of the configuration file.
 * \param[in] Expression String with the desired DOFs to be excited.
 */
int* Algorithm_GetExcitedDOF( const ConfFile_t *const Config, const char *Expression );

#endif /* _ALGORITHM_AUX_H_ */