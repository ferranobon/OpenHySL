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

#define EXACTMDOF_NUMPARAM_INIT 3  /*!< Subber of required parameters in order to initialise a substructure of
				    * type Exact fbr */
#define EXACTSDOF_NUMPARAM_INIT 3  /*!< Subber of required parameters in order to initialise a substructure of
				    * type Exact fbr */
typedef struct ExactSim {
     MatrixVector_t Mass, Stiff, Damp;
     MatrixVector_t EValues, EVectors;
     MatrixVector_t Damping_Ratios;
     MatrixVector_t Load;
     MatrixVector_t Init_Disp, Init_Vel;
     MatrixVector_t End_Disp, End_Vel, End_Acc;

     double Disp0, DispT;
     double Vel0, VelT, VelTdT;
     double Acc0, AccT, AccTdT;

     char *Description;  /*!< Optional description of the substructure. */

     Rayleigh_t Ray_Sub, Ray_Main;
     double a0, a2, a3, a6, a7;
} ExactSim_t;

typedef struct ExactSimESP {
     double Mass, Stiff, Damp;
     double Disp0, Vel0;
     double Disp, Vel;
     double Force_0, Force_1;
     double u0c_old;
     double v0c;

     double A, B, C, D, E, F, G, H;

     char *Description;  /*!< Optional description of the substructure. */
} ExactSimESP_t;

void Substructure_ExactSolutionMDOF_Init( const double *const Mass, const double *const Stiff, const int NDOF,
					  const double RayM_Alpha, const double RayM_Beta,
					  const double RayS_Alpha, const double RayS_Beta,
					  const double a0, const double a2, const double a3,
					  const double a6, const double a7, const char *Description,
					  ExactSim_t *const Sub );
void Substructure_ExactSolutionSDOF_Init( const double Mass, const double Damp, const double Stiff,
					  const double a0, const double a2, const double a3,
					  const double a6, const double a7, const char *Description,
					  ExactSim_t *const Sub );
void Compute_DampingRatios_Rayleigh( const double a0, const double a1, const int Num_DOF,
				     const double *const Eigen_Values, double *const Damping_Ratios );
void Substructure_ExactSolutionMDOF( const double DispTdT, const double ramp, const double GAcc,
				      const double DeltaT, ExactSim_t *const Sub, double *const fc );
void Substructure_ExactSolutionSDOF( const double DispTdT, const double ramp, const double GAcc,
				      const double DeltaT, ExactSim_t *const Sub, double *const fc );
void Duhamel_Integral( const MatrixVector_t *const Mass, const MatrixVector_t *const Eigen_Values,
		       const MatrixVector_t *const Eigen_Vectors, const MatrixVector_t *const Damping_Ratios,
		       const MatrixVector_t *const Init_Disp, const MatrixVector_t *const Init_Velo, 
		       const MatrixVector_t *const Load, MatrixVector_t *const End_Disp,
		       MatrixVector_t *const End_Velo, MatrixVector_t *const End_Accel, const double Delta_T );

//void Substructure_ExactSolution_SDOF( const double u0c, const double DeltaT, ExactSim_t *const Sub, double *const fc );
void Substructure_ExactSolutionMDOF_Destroy( ExactSim_t *const Sub );
void Substructure_ExactSolutionSDOF_Destroy( ExactSim_t *const Sub );


void Substructure_ExactSolutionESP_Init( const double Mass, const double Damp, const double Stiff, const double DeltaT, const char *Description, ExactSimESP_t *const Sub );
void Substructure_ExactSolutionESP_SDOF( const double u0c, const double DeltaT, ExactSimESP_t *const Sub, double *const fc );
void Substructure_ExactSolutionESP_Destroy( ExactSimESP_t *const Sub );
#endif /* SUBSTRUCTURE_EXACT_H_ */
