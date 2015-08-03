/**
 * \file Algorithm_Aux.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 16th of March 2013
 *
 * \brief Routines that depend on the implementation of the algorithm.
 *
 * This file includes routines that depend and vary with the implementation or purpose of the algorithm, like
 * reading earthquake records, configuration files, command line options, etc. Since these routines are
 * implementation dependant, they are not included in the main libraries and should be added to the compile
 * line.
 */
#ifndef _ALGORITHM_AUX_H_
#define _ALGORITHM_AUX_H_

#include <stdbool.h>
#include <sys/time.h>           /* For struct timeval */

#include "Rayleigh.h"           /* For Rayleigh_t */
#include "Error_Compensation.h" /* For PID_t */
#include "Conf_Parser.h"        /* For ConfFile_t */
#include "Definitions.h"

#define MPI_TIME_SLENGTH 30 /*!< \brief Length of the array \c Date_time in MPI since HDF5 does not support
			     * variable types in mpi mode. */

/**
 * \brief Structure to handle the elapsed time.
 *
 * This structure is used to store the starting time, ending time and elapsed time of the integration
 * process. It is used in non-mpi versions of the program and their values are stored in the HDF5 file.
 *
 * \sa SaveTime_MPI_t.
 */
typedef struct SaveTime {
     time_t Date_start;     /*!< \brief Date of the test. */
     char *Date_time;       /*!< \brief Stores the conversion in string format of the value stored in
			      \c Date_start. */
     double Elapsed_time;   /*!< \brief Elapsed time between the start of the test and the end of it in
			      millisencods. */
     struct timeval Start;  /*!< \brief Clock start time. Resolution of nanoseconds. */
     struct timeval End;    /*!< \brief Clock end time. Resolution of nanoseconds. */
} SaveTime_t;

/**
 * \brief Structure to handle the elapsed time.
 *
 * This structure is used to store the starting time, ending time and elapsed time of the integration
 * process. It is used in the mpi versions of the program and their values are stored in the HDF5 file.
 *
 * \sa SaveTime_t.
 */
typedef struct SaveTime_MPI {
     time_t Date_start;                /*!< \brief Date of the test. */
     char Date_time[MPI_TIME_SLENGTH]; /*!< \brief Fixed length because HDF5 does not support variable types
					* in mpi mode. */
     double Elapsed_time;              /*!< \brief Elapsed time between the start of the test and the end of
					* it in millisencods. */
     double Start;                     /*!< \brief MPI Wall start time. Resolution of milliconds. */
     double End;                       /*!< \brief MPI Wall end time. Resolution of miliseconds. */
} SaveTime_MPI_t;

/**
 * \brief Structure to handle the time integration constants.
 *
 * This structure handles the Newmark constants, by storing the two parameters required by the
 * Newmark-\f$beta\f$ method: \f$\gamma\f$ and \f$\beta\f$.
 */
typedef struct TIntegration{
     HYSL_FLOAT Gamma; /*!< \brief First constant.*/
     HYSL_FLOAT Beta;  /*!< \brief Second constant.*/
} TIntegration_t;

/**
 * \brief Structure to handle grid and block sizes in the MPI version of the Library.
 */
typedef struct MPIInfo {
     int Rows;     /*!< \brief Number of rows.*/
     int Cols;     /*!< \brief Number of columns.*/
} MPIInfo_t;

/**
 * \brief Structure to wrap constants and filenames.
 *
 * This structure stores several constants that will be used in different parts of the substructure algorithm,
 * such as the order of the matrices, the position of the coupling nodes, number of steps filenames and other
 * contants.
 *
 * \sa Rayleigh_t, Newmark_t, PID_t and Algorithm_Init().
 */
typedef struct AlgConst{

     MPIInfo_t ProcessGrid;   /*!< \brief Stores information regarding the process grid that will be used by
			       * BLACS/ScaLAPACK/PBLAS routines. Only referenced if MPI is used.*/
     MPIInfo_t BlockSize;     /*!< \brief Sizes of the blocks that will be used in order to distribute a
			       * matrix or a vector across the diferent processes within a grid. Only
			       * referenced if MPI is used.*/
     int Order;               /*!< \brief Order of the matrices */
     int OrderSub;            /*!< \brief Number of substructures */

     unsigned int NStep;      /*!< \brief Number of steps */
     unsigned int NSubstep;   /*!< \brief Number of sub-steps */

     int Use_Absolute_Values; /*!< \brief Variable to control whether to use absolute
			       * values in the equation of motion or relative
			       * values. Affects how the input load is calculated
			       */

     int Read_Sparse;         /*!< \brief The matrices will be read using a sparse format. \sa
			       * MatrixVector_From_File_Sp2Dense() and MatrixVector_From_File_Sp().
			       */
     int Use_Sparse;          /*!< \brief Sparse routines will be used instead of dense ones. Only available
			       * with the support of the Intel MKL libraries.  \sa
			       * MatrixVector_From_File_Sp(), CalculateMatrixC_Sp(),
			       * CalculateMatrixKeinv_Pardiso_Sparse(), Calc_Input_Load_AbsValues_Sparse(),
			       * Calc_Input_Load_RelValues_Sparse() and EffK_Calc_Effective_Force_Sparse().
			       */
     int Use_Packed;          /*!< \brief Matrices will use packed storage if possible. */
     int Read_LVector;        /*!< \brief Read the load vector instead of generating it. */
     int *ExcitedDOF;         /*!< \brief Integer string containing wich elements of the external load will be
			       * applied.  \sa AlgorithM_GetExcitedDOF().
			       */
     int Read_CMatrix;        /*!< \brief Read the damping matrix from a file instead of using Rayleigh Damping. */
     HYSL_FLOAT Delta_t;      /*!< \brief Time increment \f$\Delta t\f$ */
     HYSL_FLOAT DeltaT_Sub;   /*!< \brief Time increment for the sub-stepping process */

     HYSL_FLOAT Scale_Factor; /*!< \brief Scale factor for the time history of the input load */

     Rayleigh_t Rayleigh;     /*!< \brief Stores Rayleigh Constants alpha (\c Rayleigh.Alpha or
			       * \f$\alpha_R\f$) and beta (\c Rayleigh.Beta or \f$\beta_R\f$)
			       */
     TIntegration_t Newmark;  /*!< \brief Stores Newmark Constants gamma (\c Newmark.Gamma or \f$\gamma_N\f$)
			       * and (\c Newmark.Beta or \f$\beta_N\f$)
			       */
     PID_t PID;               /*!< \brief Stores proportional (\c PID.P), integral (\c PID.I) and derivative
			       * (\c PID.D) part.
			       */

     /* Constants for Step ending */
     HYSL_FLOAT a0;           /*!< \brief \f$a_0 = \frac{1}{\beta_N\Delta t^2}\f$ */
     HYSL_FLOAT a1;           /*!< \brief \f$a_1 = \frac{\gamma_N}{\beta_N\Delta t}\f$ */
     HYSL_FLOAT a2;           /*!< \brief \f$a_2 = \frac{1}{\beta_N\Delta t}\f$ */
     HYSL_FLOAT a3;           /*!< \brief \f$a_3 = \frac{1}{2\beta_N\Delta t} - 1\f$ */
     HYSL_FLOAT a4;           /*!< \brief \f$a_4 = \frac{\gamma_N}{\beta_N} - 1\f$ */
     HYSL_FLOAT a5;           /*!< \brief \f$a_5 = \Delta_t\biggl(\frac{\gamma_N}{2\beta_N} - 1\biggr)\f$ */
     HYSL_FLOAT a6;           /*!< \brief \f$a_6= \frac{1 - \gamma_N}{\Delta t}\f$ */
     HYSL_FLOAT a7;           /*!< \brief \f$a_7 = \gamma_N\Delta t\f$ */
     HYSL_FLOAT a8;           /*!< \brief \f$a_8 = \beta_N\Delta t^2\f$ */
     HYSL_FLOAT a9;           /*!< \brief \f$a_9 = \Delta t\f$ */
     HYSL_FLOAT a10;          /*!< \brief \f$a_{10} = \biggl(\frac{1}{2} - \beta_N\biggr)\Delta t^2\f$ */
     HYSL_FLOAT a16;          /*!< \brief \f$a_{16} = (1-2\gamma_N)\Delta t\f$ */
     HYSL_FLOAT a17;          /*!< \brief \f$a_{17} = \left(\frac{1}{2} - 2\beta_N + \gamma_N\right)\Delta t^2\f$ */
     HYSL_FLOAT a18;          /*!< \brief \f$a_{18} = \left(\frac{1}{2} + \beta_N - \gamma_N\right)\Delta t^2\f$ */

     /* Files where data are located */
     char* FileM;            /*!< \brief Stores the name of the file that contains the Mass Matrix */
     char* FileK;            /*!< \brief Stores the name of the file that contains the Stiffness Matrix */
     char* FileC;            /*!< \brief Stores the name of the file that contains the Damping Matrix */
     char* FileLV;           /*!< \brief Stores the name of the file that contains the Load Vector */
     char* FileCNodes;       /*!< \brief Stores the name of the file that contains the vector of coupling
			      * nodes. */
     char* FileData;         /*!< \brief Stores the name of the file that contains displacement, velocity and
			      * acceleration */
     char* FileOutput;       /*!< \brief Name of the file to store the output values of the process */
} AlgConst_t;


/**
 * \brief Definition of constant values and filenames that will be used during the Algorithm.
 *
 * This routine reads the values specified in a configuration file, such as the order of the matrices, the
 * number of steps and file names of the Mass and Stiffness matrices (the damping matrix is optional), and are
 * stored for further reference.
 *
 * \pre The values of the Newmark integration, PID and Rayleigh constants must be coherent/feasible. The
 * algorithm will not perform checks on them.
 *
 * \param[in]  FileName  Name of the configuration file.
 * \param[out] InitConst A structure that comprises of several constants.
 *
 * \post
 * - The size of the matrices will determine the memory that will be allocated when defining a \c
 *   MatrixVector_t type and also also how many elements will be read/written from/to the files.
 * - The number of steps must be equal to the number of rows of the file "DataFile".
 * - The memory should be deallocated through the Algorithm_Destroy() routine.
 *
 * \sa AlgoConst_t, PID_t, Rayleigh_t and TIntegration_t.
 */
void Algorithm_Init( const char *FileName, AlgConst_t *const InitConst );

/**
 * \brief Checks if the given value is correct.
 * 
 * \param[in] Value Integer to be checked.
 *
 * \return \c true if \c Value is 0 or 1, \c false otherwise.
 */
bool Valid_Value( const int Value );

/**
 * \brief Checks if a file exists and can be opened.
 * 
 * \param[in] Filename Name of the file to be opened.
 * 
 * \return \c true if the file can be opened in read mode, \c false otherwise.
 */
bool Valid_File( const char *Filename );

/**
 * \brief Changes the name of the file so that it does not overwrite an existing file.
 *
 * \param[in,out] Name On input it contains the conflicting name file. On output, a non-conflicting file name
 *                     is given.
 */
char* Change_Filename( const char *const Name );


/**
 * \brief Definition of constant values and filenames that will be used during the Algorithm. MPI version.
 *
 * \warning This routine should only be called by a single MPI process (usually process 0) and then
 * broadcasted to the rest through the Algorithm_BroadcastConfFile() routine.
 *
 * This routine reads the values specified in a configuration file, such as the order of the matrices, the
 * number of steps and file names of the Mass and Stiffness matrices (the damping matrix is optional), and are
 * stored for further reference. This is the MPI version.
 *
 * \pre The values of the Newmark integration, PID and Rayleigh constants must be coherent/feasible. The
 * algorithm will not perform checks on them.
 *
 * \param[in]  FileName  Name of the configuration file.
 * \param[out] InitConst A structure that comprises of several constants.
 *
 * \post
 * - The size of the matrices will determine the memory that will be allocated when defining a \c
 *   MatrixVector_t type and also also how many elements will be read/written from/to the files.
 * - The number of steps must be equal to the number of rows of the file "DataFile".
 * - The memory should be deallocated through the Algorithm_Destroy() routine.
 *
 * \sa AlgoConst_t, PID_t, Rayleigh_t and TIntegration_t, MPIInfo_t.
 */
void Algorithm_Init_MPI( const char *FileName, AlgConst_t *const InitConst );

/**
 * \brief Broadcasts the AlgConst_t struct to the rest of the MPI processes.
 *
 * \warning This routine requires MPI.
 *
 * \pre \c InitConst must be properly initialised through Algorithm_Init() in MPI process 0.
 *
 * \param[in] InitConst Struct to be broadcasted.
 *
 * \sa AlgConst_t.
 */
void Algorithm_BroadcastConfFile( AlgConst_t *const InitConst );

/**
 * \brief Broadcasts a string to the rest of the MPI processes and allocates the necessary memory for the
 * string to be stored in oder processes.
 *
 * \warning This routine requires MPI.
 *
 * \pre MPI must be properly initialised and the calling process must be the root process. The string must
 * not be allocated in the other processes.
 *
 * \param[in] String String to be broadcasted.
 *
 * \return A pointer to the allocated memory where the string \c String is stored.
 */
char* Algorithm_BroadcastString( char *String );

/**
 * \brief Broadcasts the array containing which are the DOFs that get excited by the ground acceleration.
 *
 * \warning This routine requires MPI.
 *
 * \pre \c InitConst must be properly initialised through Algorithm_Init() in MPI process 0.
 *
 * \param[in] IntArray The pointer to the integer array that has to be passed on.
 *
 * \return A pointer to the allocated memory where the integer array IntArray is stored.
 */
int* Algorithm_BroadcastExcitedDOF( int *IntArray );

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
 * Transforms a string with the desired sequence of DOFs that will have an external load applied. The string
 * consist of a first element that determines how many values will follow, \f$n\f$. The next \f$n\f$ values
 * are the pattern and they identify if a degree of freedom is excited (1) or not (0). The pattern will be
 * applied as many times as \f$N = \frac{Order}{n}\f$, where \f$Order\f$ is the length of the external load
 * vector.
 *
 * \b Example \b 1:
 *
 * <tt>6 1 1 1 0 0 0 </tt>
 *
 * The first three elements of the external load vector will have a value \f$u_g*1.0\f$ while the other three
 * will be set to 0.0. The pattern will repeat for every 6 elements in the load vector.
 *
 * \b Example \b 2:
 *
 * <tt>4 1 0 1 1</tt>
 *
 * Only the second element of every four positions in the external load vector will be set to 0.0.
 *
 * \pre
 * - The first element of the list should contain the number of DOFs per node.
 * - The next \f$n\f$ entries are the degrees of freedom that are being excited (1) or not (0) by the external
 *   load.
 * - The order of which DOFs will have an external load applied must remain constant and keep the same patern.
 * - The rest of \f$\frac{Order}{n}\f$ must be equal to 0, where \f$Order\f$ is the length of the external
 *   load vector.
 * - \c Config must be properly initialised through ConfFile_Create() and have entries stored through
 *   ConfFile_ReadFile().
 * - \c Expression must be a valid entry in \c Config.
 *
 * \param[in] Config Structure containing the entries of the configuration file.
 * \param[in] Expression String with the desired DOFs to be excited.
 */
int* Algorithm_GetExcitedDOF( const ConfFile_t *const Config, const char *Expression );

/**
 * \brief Reads the acceleration, velocity and displacement of an earthquake from a file.
 *
 * The data from an earthquake is readen from the specified file and its components are stored in different
 * variables, namely: acceleration, velocity and displacement. If the specified file does not exist, an
 * appropriate error is displayed and the program exists with \c EXIT_FAILURE.
 *
 * \pre
 * - The input text file should be in ASCII format and it must contain no header and four columns in the
 *   following order:
 *    -# It can contain anything, since it will not be stored, it is recommended that it has the step number
 *       to detect any possible errors (see the next precondition).
 *    -# Acceleration in \f$[mm/s^2]\f$.
 *    -# Velocity in \f$[mm/s]\f$.
 *    -# Displacement in \f$[mm]\f$.
 * - Note also that the number of steps in the file should be equal to the number of steps defined in
 *   Algorithm_Init()InitConstants(). No check is performed.
 * - The arrays must be properly initialised and its length should be \f$l = Number~of~Steps\f$.
 *
 * \param[in]  NumSteps     The number of steps in which the earthquake is divided.
 * \param[in]  Filename     The name of the file where the earthquake is recorded.
 * \param[in]  Scale_Factor Factor that will be multiplied to the readen values.
 * \param[out] Acceleration Array that will store the second column of the data file, that is, the
 *                          acceleration of the recorder earthquake converted to \f$[mm/s^2]\f$.
 * \param[out] Velocity     Array that will store the third column of the data file, that is, the velocity of
 *                          the recorded earthquake converted to \f$[m/s]\f$.
 * \param[out] Displacement Array that will store the fourth column of the data file, that is, the
 *                          displacement of the recorded earthquake converted to \f$[m]\f$.
 *
 * \post \c Acceleration, \c Velocity and \c Displacement will store the acceleration, velocity and
 * displacement of the recorded earthquake in \f$[m/s^2]\f$, \f$[m/s]\f$ and \f$[m]\f$ respectively.
 *
 * \sa Algorithm_Initnit().
 */
void Algorithm_ReadDataEarthquake( const unsigned int NumSteps, const char *Filename,
				   const HYSL_FLOAT Scale_Factor, HYSL_FLOAT *const Acceleration,
				   HYSL_FLOAT *const Velocity, HYSL_FLOAT *const Displacement );

/**
 * \brief Prints a help text with the different options available when launching the program.
 *
 * \param[in] Program_Name Name of the executable.
 */
void Algorithm_PrintHelp( const char *Program_Name );

/**
 * \brief Concatenates different strings into one.
 *
 * This function concatenates different strings into one. It takes any number of arguments but the first
 * argument must be always the number of strings that are going to be passed. The strings are concatenated in
 * the same order that is given in the arguments. For example:
 * 
 * <tt>StrConc = Concatenate_Strings( 2, "This is a ", "concatenated string." );</tt>
 *
 * will result in <tt>StrConc = "This is a concatenated string."</tt>
 *
 * \pre \c count is the number of strings that are going to be concatenated. All strings must be correctly
 * initialised except for the newly created string that is given as a return.
 *
 * \param[in] count Number of strings to concatenate.
 *
 * \return A pointer to the allocated memory containing the concatenated strings.
 *
 * \post The initial strings are still accessible while a new string, with all the concatenated strings in it
 * is given as return.
 */
char* Concatenate_Strings( int count, ... );

#endif /* _ALGORITHM_AUX_H_ */
