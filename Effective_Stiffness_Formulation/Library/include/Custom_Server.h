/**
 * \file Custom_Server.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 19th of June 2013
 *
 * \brief Routines that depend on the implementation of the custom server.
 *
 */
#ifndef _CUSTOM_SEVER_H_
#define _CUSTOM_SEVER_H_

#include "Algorithm_Aux.h" /* For TIntegration_t */

typedef struct CustomSever{

     int Order_Couple;        /*!< \brief Order of the coupling degrees of freedom */

     unsigned int NSubstep;   /*!< \brief Number of sub-steps */

     double Delta_t;          /*!< \brief Time increment \f$\Delta t\f$ */
     double DeltaT_Sub;       /*!< \brief Time increment for the sub-stepping process */

     TIntegration_t Newmark;  /*!< \brief Stores Newmark Constants gamma (\c Newmark.Gamma or \f$\gamma_N\f$)
			       * and (\c Newmark.Beta or \f$\beta_N\f$)
			       */

     /* Constants for Step ending */
     double a0;               /*!< \brief \f$a_0 = \frac{1}{\beta_N\Delta t^2}\f$ */
     double a1;               /*!< \brief \f$a_1 = \frac{\gamma_N}{\beta_N\Delta t}\f$ */
     double a2;               /*!< \brief \f$a_2 = \frac{1}{\beta_N\Delta t}\f$ */
     double a3;               /*!< \brief \f$a_3 = \frac{1}{2\beta_N\Delta t} - 1\f$ */
     double a4;               /*!< \brief \f$a_4 = \frac{\gamma_N}{\beta_N} - 1\f$ */
     double a5;               /*!< \brief \f$a_5 =
			       * \Delta_t\biggl(\frac{\gamma_N}{2\beta_N} - 1\biggr)\f$ */
     double a6;               /*!< \brief \f$a_6= \frac{1 - \gamma_N}{\Delta t}\f$ */
     double a7;               /*!< \brief \f$a_7 = \gamma_N\Delta t\f$ */

     /* Files where data are located */
     char* FileCNodes;       /*!< \brief Stores the name of the file that contains the vector of coupling
			      * nodes. */
} CustomServer_t;

void CustomServer_Init( const char *FileName, AlgConst_t *const InitConst );
void CustomServer_Destroy( AlgConst_t *const InitConst );
void CustomServer_PrintHelp( const char *Program_Name );

#endif /* _CUSTOM_SEVER_H_ */
