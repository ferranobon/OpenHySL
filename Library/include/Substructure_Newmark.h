/**
 * \file Substructure_Newmark.h
 * \author Ferran Obón Santacana
 * \version 1.0 
 * \date 4th of July 2013
 *
 * \brief Routines to deal with substructures that can be solved through Newmark integration routines.
 *
 * These routines deal with the creation, destruction and simulation of those substructures that can be solved
 * using Newmark. Therefore, they are used mainly for vibrating substructures like Tunned Mass Dampers (TMDs)
 * or other structural components (like in the Arianne IV tests)
 */
#ifndef SUBSTRUCTURE_NEWMARK_H_
#define SUBSTRUCTURE_NEWMARK_H_

#include "Substructure.h"
#include "MatrixVector.h"
#include "Rayleigh.h"
#include "Definitions.h"

#define NEWMARK_NUMPARAM_INIT 3  /*!< \brief Number of required parameters in order to initialise a
				    * substructure of type Newmark */

typedef struct NewmarkSim {

     HYSL_FLOAT Beta, Gamma;

     HYSL_FLOAT Mass, Stiff, Damp, G;
     HYSL_FLOAT dT, vT, aT, l;
     HYSL_FLOAT dTdT, vTdT, aTdT;

     HYSL_FLOAT Disp0, DispT;
     HYSL_FLOAT Vel0, VelT, VelTdT;
     HYSL_FLOAT Acc0, AccT, AccTdT;

     HYSL_FLOAT a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10;
     HYSL_FLOAT A0, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10;
     char *Description;  /*!< \brief Optional description of the substructure. */
} NewmarkSim_t;

/**
 * \param[in]  Mass
 * \param[in]  Damp
 * \param[in]  Stiff
 * \param[in]  DeltaT
 * \param[in]  Description
 * \param[out] Sub
 */
void Substructure_Newmark_Init( const HYSL_FLOAT Mass, const HYSL_FLOAT Damp, const HYSL_FLOAT Stiff, const HYSL_FLOAT DeltaTSub, const HYSL_FLOAT DeltaT, const HYSL_FLOAT Beta,
				const HYSL_FLOAT Gamma, const char *Description, NewmarkSim_t *const Sub );

/**
 * \param[in]     u0c
 * \param[in]     DeltaT
 * \param[in,out] Sub
 * \param[out]    fc
 */
void Substructure_Newmark_SDOF( const HYSL_FLOAT DispTdT, const HYSL_FLOAT ramp, const HYSL_FLOAT GAcc, NewmarkSim_t *const Sub, HYSL_FLOAT *const fc );

/**
 * \param[out] Sub
 */
void Substructure_Newmark_Destroy( NewmarkSim_t *const Sub );
#endif /* SUBSTRUCTURE_NEWMARK_H_ */
