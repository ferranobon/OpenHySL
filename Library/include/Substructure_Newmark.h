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

     hysl_float_t Beta, Gamma;

     hysl_float_t Mass, Stiff, Damp, G;
     hysl_float_t dT, vT, aT, l;
     hysl_float_t dTdT, vTdT, aTdT;

     hysl_float_t Disp0, DispT;
     hysl_float_t Vel0, VelT, VelTdT;
     hysl_float_t Acc0, AccT, AccTdT;

     hysl_float_t a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10;
     hysl_float_t A0, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10;
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
void Substructure_Newmark_Init( const hysl_float_t Mass, const hysl_float_t Damp, const hysl_float_t Stiff, const hysl_float_t DeltaTSub, const hysl_float_t DeltaT, const hysl_float_t Beta,
				const hysl_float_t Gamma, const char *Description, NewmarkSim_t *const Sub );

/**
 * \param[in]     u0c
 * \param[in]     DeltaT
 * \param[in,out] Sub
 * \param[out]    fc
 */
void Substructure_Newmark_SDOF( const hysl_float_t DispTdT, const hysl_float_t ramp, const hysl_float_t GAcc, NewmarkSim_t *const Sub, hysl_float_t *const fc );

/**
 * \param[out] Sub
 */
void Substructure_Newmark_Destroy( NewmarkSim_t *const Sub );
#endif /* SUBSTRUCTURE_NEWMARK_H_ */
