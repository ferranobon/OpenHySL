/**
 * \file Substructure_Exact.h
 * \author Ferran Obón Santacana
 * \version 1.0 
 * \date 13th of September 2013
 *
 * \brief Routines to deal with substructures that can be solved through exact integration routines.
 *
 * These routines deal with the creation, destruction and simulation of those substructures that can be solved
 * using exact integration methods. Therefore, they are used mainly for vibrating substructures like Tunned
 * Mass Dampers (TMDs) or other structural components (like in the Arianne IV tests)
 */
#ifndef SUBSTRUCTURE_EXACT_H_
#define SUBSTRUCTURE_EXACT_H_

#include "Substructure.h"
#include "MatrixVector.h"
#include "Rayleigh.h"
#include "Definitions.h"

#define EXACTMDOF_NUMPARAM_INIT 3  /*!< \brief Number of required parameters in order to initialise a
				    * substructure of type Exact fbr */
#define EXACTSDOF_NUMPARAM_INIT 3  /*!< \brief Number of required parameters in order to initialise a
				    * substructure of type Exact fbr */
typedef struct ExactSim {
     MatrixVector_t Mass, Stiff, Damp;
     MatrixVector_t EValues, EVectors;
     MatrixVector_t Damping_Ratios;
     MatrixVector_t Load;
     MatrixVector_t Init_Disp, Init_Vel;
     MatrixVector_t End_Disp, End_Vel, End_Acc;

     HYSL_FLOAT Disp0, DispT;
     HYSL_FLOAT Vel0, VelT, VelTdT;
     HYSL_FLOAT Acc0, AccT, AccTdT;

     char *Description;  /*!< \brief Optional description of the substructure. */

     Rayleigh_t Ray_Sub, Ray_Main;
     HYSL_FLOAT a0, a2, a3, a6, a7;
} ExactSim_t;

typedef struct ExactSimESP {
     HYSL_FLOAT Mass, Stiff, Damp, xi, omega;
     HYSL_FLOAT Init_Disp, Init_Vel;
     HYSL_FLOAT End_Disp, End_Vel, End_Acc;
     HYSL_FLOAT Force_0, Force_1;
     HYSL_FLOAT Disp0, DispT;
     HYSL_FLOAT Vel0, VelT, VelTdT;
     HYSL_FLOAT Acc0, AccT, AccTdT;

     HYSL_FLOAT A, B, C, D, E, F, G, H;
     HYSL_FLOAT a0, a2, a3, a6, a7;
     char *Description;  /*!< \brief Optional description of the substructure. */
} ExactSimESP_t;

/**
 * \param[in]  Mass
 * \param[in]  Stiff
 * \param[in]  NDOF
 * \param[in]  RayM_Alpha
 * \param[in]  RayM_Beta
 * \param[in]  RayS_Alpha
 * \param[in]  RayS_Beta
 * \param[in]  a0
 * \param[in]  a2
 * \param[in]  a3
 * \param[in]  a6
 * \param[in]  a7
 * \param[in]  Description
 * \param[out] Sub
 */
void Substructure_ExactSolutionMDOF_Init( const HYSL_FLOAT *const Mass, const HYSL_FLOAT *const Stiff, const int NDOF,
					  const HYSL_FLOAT RayM_Alpha, const HYSL_FLOAT RayM_Beta,
					  const HYSL_FLOAT RayS_Alpha, const HYSL_FLOAT RayS_Beta,
					  const HYSL_FLOAT a0, const HYSL_FLOAT a2, const HYSL_FLOAT a3,
					  const HYSL_FLOAT a6, const HYSL_FLOAT a7, const char *Description,
					  ExactSim_t *const Sub );
/**
 * \param[in]  Mass
 * \param[in]  Damp
 * \param[in]  Stiff
 * \param[in]  a0
 * \param[in]  a2
 * \param[in]  a3
 * \param[in]  a6
 * \param[in]  a7
 * \param[in]  Description
 * \param[out] Sub
 */
void Substructure_ExactSolutionSDOF_Init( const HYSL_FLOAT Mass, const HYSL_FLOAT Damp, const HYSL_FLOAT Stiff,
					  const HYSL_FLOAT a0, const HYSL_FLOAT a2, const HYSL_FLOAT a3,
					  const HYSL_FLOAT a6, const HYSL_FLOAT a7, const char *Description,
					  ExactSim_t *const Sub );
/**
 *
 * \param[in]  a0
 * \param[in]  a1
 * \param[in]  Num_DOF
 * \param[in]  Eigen_Values
 * \param[out] Damping_Ratios
 */
void Compute_DampingRatios_Rayleigh( const HYSL_FLOAT a0, const HYSL_FLOAT a1, const int Num_DOF,
				     const HYSL_FLOAT *const Eigen_Values, HYSL_FLOAT *const Damping_Ratios );

/**
 * \param[in]     DispTdT
 * \param[in]     ramp
 * \param[in]     GAcc
 * \param[in]     DeltaT
 * \param[in,out] Sub
 * \param[out]    fc
 */
void Substructure_ExactSolutionMDOF( const HYSL_FLOAT DispTdT, const HYSL_FLOAT ramp, const HYSL_FLOAT GAcc,
				      const HYSL_FLOAT DeltaT, ExactSim_t *const Sub, HYSL_FLOAT *const fc );

/**
 * \param[in]     DispTdT
 * \param[in]     ramp
 * \param[in]     GAcc
 * \param[in]     DeltaT
 * \param[in,out] Sub
 * \param[out]    fc
 */
void Substructure_ExactSolutionSDOF( const HYSL_FLOAT DispTdT, const HYSL_FLOAT ramp, const HYSL_FLOAT GAcc,
				      const HYSL_FLOAT DeltaT, ExactSim_t *const Sub, HYSL_FLOAT *const fc );

/**
 * \param[in]  Mass
 * \param[in]  Eigen_Values
 * \param[in]  Eigen_Vectors
 * \param[in]  Damping_Ratios
 * \param[in]  Init_Disp
 * \param[in]  Init_Velo
 * \param[in]  Load
 * \param[out] End_Disp
 * \param[out] End_Velo
 * \param[out] End_Accel
 * \param[in]  Delta_T
 */
void Duhamel_Integral( const MatrixVector_t *const Mass, const MatrixVector_t *const Eigen_Values,
		       const MatrixVector_t *const Eigen_Vectors, const MatrixVector_t *const Damping_Ratios,
		       const MatrixVector_t *const Init_Disp, const MatrixVector_t *const Init_Velo, 
		       const MatrixVector_t *const Load, MatrixVector_t *const End_Disp,
		       MatrixVector_t *const End_Velo, MatrixVector_t *const End_Accel, const HYSL_FLOAT Delta_T );

//void Substructure_ExactSolution_SDOF( const HYSL_FLOAT u0c, const HYSL_FLOAT DeltaT, ExactSim_t *const Sub, HYSL_FLOAT
//*const fc );

/**
 * \param[in,out] Sub
 */
void Substructure_ExactSolutionMDOF_Destroy( ExactSim_t *const Sub );

/**
 * \param[in,out] Sub
 */
void Substructure_ExactSolutionSDOF_Destroy( ExactSim_t *const Sub );

/**
 * \param[in]  Mass
 * \param[in]  Damp
 * \param[in]  Stiff
 * \param[in]  DeltaT
 * \param[in]  Description
 * \param[out] Sub
 */
void Substructure_ExactSolutionESP_Init( const HYSL_FLOAT Mass, const HYSL_FLOAT Damp, const HYSL_FLOAT Stiff, const HYSL_FLOAT DeltaT, const char *Description, ExactSimESP_t *const Sub );

/**
 * \param[in]     u0c
 * \param[in]     DeltaT
 * \param[in,out] Sub
 * \param[out]    fc
 */
void Substructure_ExactSolutionESP_SDOF( const HYSL_FLOAT u0c, const HYSL_FLOAT ramp, const HYSL_FLOAT GAcc, const HYSL_FLOAT DeltaT, ExactSimESP_t *const Sub, HYSL_FLOAT *const fc );

/**
 * \param[out] Sub
 */
void Substructure_ExactSolutionESP_Destroy( ExactSimESP_t *const Sub );

#endif /* SUBSTRUCTURE_EXACT_H_ */
