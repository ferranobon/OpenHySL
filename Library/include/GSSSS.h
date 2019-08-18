/**
 * \file GSSSS.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 26th of June 2014
 *
 * \brief Implementation of the GSSSS framework.
 *
 * This file contains the routines that implement the GSSSS framework into the Subfeed algorithm in the
 * A-Form, V-Form and D-Form. The in-depth description is given in \cite[Tamma_2012].
 * line.
 */
#ifndef _GSSSS_H_
#define _GSSSS_H_

#include "MatrixVector.h"
#include "MatrixVector_Sp.h"
#include "Error_Compensation.h"
#include "Definitions.h"

/**
 * \brief Structure to handle the time integration parameters in the GSSSS framework.
 *
 * This structure handles the integration parameters in the GSSSS algorithms. See \cite[Tamma_2004].
 */
typedef struct TIntegration_GSSSS{
     hysl_float_t w1;      /*!< \brief 1st parameter of the degenerated scalar polynomial function. See Eq. 46
		          *          in \cite[Tamma_2004].*/
     hysl_float_t w2;      /*!< \brief 2nd parameter of the degenerated scalar polynomial function. See Eq. 46
		          *          in \cite[Tamma_2004].*/
     hysl_float_t w3;      /*!< \brief 3rd parameter of the degenerated scalar polynomial function. See Eq. 46
		          *          in \cite[Tamma_2004].*/
     hysl_float_t A1W1;    /*!< \brief Parameter \f$\Lambda_1 W_1\f$. See Eq. 53 in \cite[Tamma_2004].*/
     hysl_float_t A2W2;    /*!< \brief Parameter \f$\Lambda_2 W_2\f$. See Eq. 53 in \cite[Tamma_2004].*/
     hysl_float_t A3W3;    /*!< \brief Parameter \f$\Lambda_3 W_3\f$. See Eq. 53 in \cite[Tamma_2004].*/
     hysl_float_t A4W1;    /*!< \brief Parameter \f$\Lambda_4 W_1\f$. See Eq. 53 in \cite[Tamma_2004].*/
     hysl_float_t A5W2;    /*!< \brief Parameter \f$\Lambda_5 W_2\f$. See Eq. 53 in \cite[Tamma_2004].*/
     hysl_float_t A6W1;    /*!< \brief Parameter \f$\Lambda_6 W_1\f$. See Eq. 53 in \cite[Tamma_2004].*/
     hysl_float_t W1;      /*!< \brief Weight value (Forces).*/
     hysl_float_t l1;      /*!< \brief Integration constant \f$\lambda_1\f$.*/
     hysl_float_t l2;      /*!< \brief Integration constant \f$\lambda_2\f$.*/
     hysl_float_t l3;      /*!< \brief Integration constant \f$\lambda_3\f$.*/
     hysl_float_t l4;      /*!< \brief Integration constant \f$\lambda_4\f$.*/
     hysl_float_t l5;      /*!< \brief Integration constant \f$\lambda_5\f$.*/
     hysl_float_t rho_inf; /*!< \brief Minimum absolute value of the eigenvalues of the amplification matrix
			  *          \f$\rho_\infty\f$. See Eq. 61 in \cite[Tamma_2004].*/
} TIntegration_GSSSS_t;

void GSSSS_Init( const hysl_float_t Rho1, const hysl_float_t Rho2, const hysl_float_t Rho3, TIntegration_GSSSS_t *const GSSSS );

void GSSSS_EffectiveForce_AForm( const MatrixVector_t *const Mass, const MatrixVector_t *const Damp,
				 const MatrixVector_t *const Stiff, const MatrixVector_t *const DispT,
				 const MatrixVector_t *const VelT, const MatrixVector_t *const AccT,
				 MatrixVector_t *const Tempvec, const TIntegration_GSSSS_t *const GSSSS,
				 const hysl_float_t DeltaT, MatrixVector_t *const Eff_ForceT );

void GSSSS_EffectiveForce_AForm_PS( const MatrixVector_t *const Mass, const MatrixVector_t *const Damp,
				    const MatrixVector_t *const Stiff, const MatrixVector_t *const DispT,
				    const MatrixVector_t *const VelT, const MatrixVector_t *const AccT,
				    MatrixVector_t *const Tempvec, const TIntegration_GSSSS_t *const GSSSS,
				    const hysl_float_t DeltaT, MatrixVector_t *const Eff_ForceT );

void GSSSS_EffectiveForce_AForm_Sp( const MatrixVector_Sp_t *const Mass, const MatrixVector_Sp_t *const Damp,
				    const MatrixVector_Sp_t *const Stiff, const MatrixVector_t *const DispT,
				    const MatrixVector_t *const VelT, const MatrixVector_t *const AccT,
				    MatrixVector_t *const Tempvec, const TIntegration_GSSSS_t *const GSSSS,
				    const hysl_float_t DeltaT, MatrixVector_t *const Eff_ForceT );

void GSSSS_Compute_NewState( const MatrixVector_t *const IGain, const MatrixVector_t *const Eff_ForceT,
			     const MatrixVector_t *const LoadTdT, const MatrixVector_t *const Err_ForceT,
			     const MatrixVector_t *const LoadT, const TIntegration_GSSSS_t *const GSSSS,
			     MatrixVector_t *const Tempvec, MatrixVector_t *const VecTdT_0 );

void GSSSS_Compute_NewState_PS( const MatrixVector_t *const IGain, const MatrixVector_t *const Eff_ForceT,
				const MatrixVector_t *const LoadTdT, const MatrixVector_t *const Err_ForceT,
				const MatrixVector_t *const LoadT, const TIntegration_GSSSS_t *const GSSSS,
				MatrixVector_t *const Tempvec, MatrixVector_t *const VecTdT_0 );

void GSSSS_ComputeDisplacement_AForm( const MatrixVector_t *const DispT, const MatrixVector_t *const VelT,
				      const MatrixVector_t *const AccT, const MatrixVector_t *const AccTdT,
				      const TIntegration_GSSSS_t *const GSSSS, const hysl_float_t DeltaT,
				      MatrixVector_t *const DispTdT );

void GSSSS_ComputeVelocity_AForm( const MatrixVector_t *const VelT, const MatrixVector_t *const AccT,
				  const MatrixVector_t *const AccTdT, const TIntegration_GSSSS_t *const GSSSS,
				  const hysl_float_t DeltaT, MatrixVector_t *const VelTdT );

void GSSSS_ErrorForce_PID( const MatrixVector_t *const Mass, const MatrixVector_t *const Damp, const MatrixVector_t *Stiff,
			   const MatrixVector_t *const AccTdT, const MatrixVector_t *const AccT, const MatrixVector_t *const VelT,
			   const MatrixVector_t *const DispT, const MatrixVector_t *const fc, const MatrixVector_t *const LoadTdT,
			   const MatrixVector_t *const LoadT, const TIntegration_GSSSS_t *const GSSSS, const hysl_float_t DeltaT,
			   const MatrixVector_t *const Tempvec, const PID_t *const PID, MatrixVector_t *const fe );

void GSSSS_ErrorForce_PID_PS( const MatrixVector_t *const Mass, const MatrixVector_t *const Damp, const MatrixVector_t *Stiff,
			      const MatrixVector_t *const AccTdT, const MatrixVector_t *const AccT, const MatrixVector_t *const VelT,
			      const MatrixVector_t *const DispT, const MatrixVector_t *const fc, const MatrixVector_t *const LoadTdT,
			      const MatrixVector_t *const LoadT, const TIntegration_GSSSS_t *const GSSSS, const hysl_float_t DeltaT,
			      const MatrixVector_t *const Tempvec, const PID_t *const PID, MatrixVector_t *const fe );

void GSSSS_ErrorForce_PID_Sp( const MatrixVector_Sp_t *const Mass, const MatrixVector_Sp_t *const Damp,
			      const MatrixVector_Sp_t *Stiff, const MatrixVector_t *const AccTdT,
			      const MatrixVector_t *const AccT, const MatrixVector_t *const VelT,
			      const MatrixVector_t *const DispT, const MatrixVector_t *const fc,
			      const MatrixVector_t *const LoadTdT, const MatrixVector_t *const LoadT,
			      const TIntegration_GSSSS_t *const GSSSS, const hysl_float_t DeltaT,
			      const MatrixVector_t *const Tempvec, const PID_t *const PID,
			      MatrixVector_t *const fe );

#endif
