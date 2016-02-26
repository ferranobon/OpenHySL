#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Substructure_StoneDrums.h"
#include "Auxiliary_Math.h"
#include "Print_Messages.h"
#include "Definitions.h"

void Substructure_StoneDrums_Init (const int PrevDOF, const HYSL_FLOAT alpha, const HYSL_FLOAT ko,
				   const HYSL_FLOAT beta, const HYSL_FLOAT gamma, const HYSL_FLOAT n,
				   const HYSL_FLOAT A0, const HYSL_FLOAT deltaA, const HYSL_FLOAT nu0,
				   const HYSL_FLOAT deltaNu, const HYSL_FLOAT eta0, const HYSL_FLOAT deltaEta,
				   const int BoucWen_Type, const char *Description, StoneDrums_t *const Sub )
{
     Sub->Description = strdup( Description );
  
     Sub->PrevDOF = PrevDOF;

     Sub->BoucWen_Type = BoucWen_Type;

     /* Basic Bouc-Wen parameters */
     Sub->alpha = alpha;
     Sub->ko = ko;
     Sub->beta = beta;
     Sub->gamma = gamma;
     Sub->n = n;

     /* Parameters for material degradation */
     Sub->A0 = A0;
     Sub->deltaA = deltaA;
     Sub->nu0 = nu0;
     Sub->deltaNu = deltaNu;
     Sub->eta0 = eta0;
     Sub->deltaEta = deltaEta;

     /* Initialisation of Newton-Raphson variables */
     Sub->maxIter = 100;
#if _FLOAT_
     Sub->tolerance = 1E-3f;
     Sub->startPoint = 0.01f;
     Sub->z_old = 0.0f;
     Sub->e_old = 0.0f;
     Sub->DispT = 0.0f;
     Sub->Result_Fc = 0.0f;
#else 
     Sub->tolerance = 1E-3;
     Sub->startPoint = 0.01;
     Sub->z_old = 0.0;
     Sub->e_old = 0.0;
     Sub->DispT = 0.0;
     Sub->Result_Fc = 0.0;
#endif

}

void Substructure_StoneDrums ( HYSL_FLOAT DispTdT, StoneDrums_t *const Sub, HYSL_FLOAT *force )
{
     int count;
     HYSL_FLOAT sign;
     HYSL_FLOAT e_new, A_new, nu_new, eta_new, Theta, Phi;     
     HYSL_FLOAT e_new_, A_new_, nu_new_, eta_new_, Phi_;
     HYSL_FLOAT vs_1, vs_2, zu, hze_exp, hze;
     HYSL_FLOAT hze_;
     HYSL_FLOAT fz_new, fz_new_;
     HYSL_FLOAT z_new, z_new_p, z_eval;
     HYSL_FLOAT deltaDisp;
     
     count = 0;
#if _FLOAT_
     z_new = 1.0f;
#else
     z_new = 1.0;
#endif
     z_new_p = Sub->startPoint; z_eval = Sub->startPoint;
     deltaDisp = DispTdT - Sub->DispT;

     while ( fabs(z_new_p - z_new) > Sub->tolerance && count < Sub->maxIter ){
	  e_new = Sub->e_old + (1.0 - Sub->alpha)*Sub->ko*deltaDisp*z_eval;
	  A_new = Sub->A0 - Sub->deltaA*e_new;
	  nu_new = Sub->nu0 + Sub->deltaNu*e_new;
	  eta_new = Sub->eta0 + Sub->deltaEta*e_new;
	  
	  sign = signum(deltaDisp*z_eval);
	  Theta = Sub->gamma + Sub->beta*sign;
	  Phi = A_new - pow(fabs(z_eval),Sub->n)*Theta*nu_new;
	  
	  vs_1 = Sub->vs0*(1.0 - exp(-Sub->p*e_new));
	  vs_2 = (Sub->psi0 + Sub->deltaPsi*e_new)*(Sub->lambda + vs_1);
	  zu = pow(1.0/((Sub->nu0 + Sub->deltaNu*e_new)*(Sub->beta + Sub->gamma)), 1.0/Sub->n);
	  hze_exp = exp(-pow(z_eval*signum(deltaDisp) -Sub->q*zu, 2.0)/pow(vs_2, 2.0));
	  hze = 1.0 - vs_1*hze_exp;
	       
	  fz_new = z_eval - Sub->z_old - hze*Phi/eta_new*deltaDisp;
	  
	  /* Evaluate function derivatives with respect to z_eval for the Newton-Rhapson scheme */
	  e_new_ = (1.0 - Sub->alpha)*Sub->ko*deltaDisp;
	  A_new_ = -Sub->deltaA*e_new_;
	  nu_new_ = Sub->deltaNu*e_new_;
	  eta_new_ = Sub->deltaEta*e_new_;
	  sign = signum(z_eval);
	  Phi_ = A_new_ - Sub->n*pow(fabs(z_eval), Sub->n - 1.0)*sign*Theta*nu_new - pow(fabs(z_eval), Sub->n)*Theta*nu_new_;
	  
	  hze_ = deltaDisp*Sub->ko*Sub->p*Sub->vs0*exp(-pow(z_eval*signum(deltaDisp) - Sub->q*pow(1.0/((Sub->beta + Sub->gamma)*(Sub->nu0 + Sub->deltaNu*(Sub->e_old - deltaDisp*Sub->ko*z_eval*(Sub->alpha - 1.0)))), 1.0/Sub->n), 2.0)/(pow(Sub->lambda - Sub->vs0*(exp(-Sub->p*(Sub->e_old - deltaDisp*Sub->ko*z_eval*(Sub->alpha - 1.0))) - 1.0), 2.0)*pow(Sub->psi0 + Sub->deltaPsi*(Sub->e_old - deltaDisp*Sub->ko*z_eval*(Sub->alpha - 1.0)), 2.0)))*exp(-Sub->p*(Sub->e_old - deltaDisp*Sub->ko*z_eval*(Sub->alpha - 1.0)))*(Sub->alpha - 1.0) - Sub->vs0*exp(-pow(z_eval*signum(deltaDisp) - Sub->q*pow(1.0/((Sub->beta + Sub->gamma)*(Sub->nu0 + Sub->deltaNu*(Sub->e_old - deltaDisp*Sub->ko*z_eval*(Sub->alpha - 1.0)))),1.0/Sub->n), 2.0)/(pow(Sub->lambda - Sub->vs0*(exp(-Sub->p*(Sub->e_old - deltaDisp*Sub->ko*z_eval*(Sub->alpha - 1.0))) - 1.0), 2.0)*pow(Sub->psi0 + Sub->deltaPsi*(Sub->e_old - deltaDisp*Sub->ko*z_eval*(Sub->alpha - 1.0)), 2.0)))*(exp(-Sub->p*(Sub->e_old - deltaDisp*Sub->ko*z_eval*(Sub->alpha - 1.0))) - 1.0)*((2*(signum(deltaDisp) - (Sub->deltaNu*deltaDisp*Sub->ko*Sub->q*(Sub->alpha - 1.0)*pow(1.0/((Sub->beta + Sub->gamma)*(Sub->nu0 + Sub->deltaNu*(Sub->e_old - deltaDisp*Sub->ko*z_eval*(Sub->alpha - 1.0)))),1.0/Sub->n - 1.0))/(Sub->n*(Sub->beta + Sub->gamma)*pow(Sub->nu0 + Sub->deltaNu*(Sub->e_old - deltaDisp*Sub->ko*z_eval*(Sub->alpha - 1.0)), 2.0)))*(z_eval*signum(deltaDisp) - Sub->q*pow(1.0/((Sub->beta + Sub->gamma)*(Sub->nu0 + Sub->deltaNu*(Sub->e_old - deltaDisp*Sub->ko*z_eval*(Sub->alpha - 1.0)))),1.0/Sub->n)))/(pow(Sub->lambda - Sub->vs0*(exp(-Sub->p*(Sub->e_old - deltaDisp*Sub->ko*z_eval*(Sub->alpha - 1.0))) - 1.0), 2.0)*pow(Sub->psi0 + Sub->deltaPsi*(Sub->e_old - deltaDisp*Sub->ko*z_eval*(Sub->alpha - 1.0)), 2.0)) + (2*Sub->deltaPsi*deltaDisp*Sub->ko*(Sub->alpha - 1.0)*pow(z_eval*signum(deltaDisp) - Sub->q*pow(1.0/((Sub->beta + Sub->gamma)*(Sub->nu0 + Sub->deltaNu*(Sub->e_old - deltaDisp*Sub->ko*z_eval*(Sub->alpha - 1.0)))),1.0/Sub->n), 2.0))/(pow(Sub->lambda - Sub->vs0*(exp(-Sub->p*(Sub->e_old - deltaDisp*Sub->ko*z_eval*(Sub->alpha - 1.0))) - 1.0), 2.0)*pow(Sub->psi0 + Sub->deltaPsi*(Sub->e_old - deltaDisp*Sub->ko*z_eval*(Sub->alpha - 1.0)),3)) + (2*deltaDisp*Sub->ko*Sub->p*Sub->vs0*exp(-Sub->p*(Sub->e_old - deltaDisp*Sub->ko*z_eval*(Sub->alpha - 1.0)))*(Sub->alpha - 1.0)*pow(z_eval*signum(deltaDisp) - Sub->q*pow(1.0/((Sub->beta + Sub->gamma)*(Sub->nu0 + Sub->deltaNu*(Sub->e_old - deltaDisp*Sub->ko*z_eval*(Sub->alpha - 1.0)))),1.0/Sub->n), 2.0))/(pow(Sub->lambda - Sub->vs0*(exp(-Sub->p*(Sub->e_old - deltaDisp*Sub->ko*z_eval*(Sub->alpha - 1.0))) - 1.0),3.0)*pow(Sub->psi0 + Sub->deltaPsi*(Sub->e_old - deltaDisp*Sub->ko*z_eval*(Sub->alpha - 1.0)), 2.0)));
	      
	  fz_new_ = 1.0 - (hze_*Phi*eta_new + hze*Phi_*eta_new - hze*Phi*eta_new_)/pow(eta_new, 2.0)*deltaDisp;

	  /* Perform a new step */
	  z_new = z_eval - fz_new/fz_new_;
	  
	  /* Update the root */
	  z_new_p = z_eval;
	  z_eval = z_new;
	  
	  count = count + 1;
	  
	  /* Warning if there is no convergence */
	  if (count == Sub->maxIter) {
	       Print_Header( WARNING );
	       fprintf( stderr, "Substructure_BoucWen(): could not find the root z_{i+1}, after %i iterations and norm: %lE\n", Sub->maxIter, fabs(z_new_p - z_new) );
	  }
	  
	  // Compute restoring force.
	  *force = Sub->alpha*Sub->ko*DispTdT + (1.0 - Sub->alpha)*Sub->ko*z_eval;
	  
	  // Compute material degradation parameters.
	  e_new = Sub->e_old + (1.0 - Sub->alpha)*Sub->ko*deltaDisp*z_eval;
	  A_new = Sub->A0 - Sub->deltaA*e_new;
	  nu_new = Sub->nu0 + Sub->deltaNu*e_new;
	  eta_new = Sub->eta0 + Sub->deltaEta*e_new;
     }
     Sub->DispT = DispTdT;
     Sub->e_old = e_new;
     Sub->z_old = z_eval;

     Sub->Result_Fc = *force;
}

#if _1_
void Substructure_StoneDrums ( const HYSL_FLOAT DispTdT, const HYSL_FLOAT ramp, StoneDrums_t * const Sub, HYSL_FLOAT *const fc )
{
     HYSL_FLOAT dStrain, Tstrain, strain;


     /* Newton-Rhapson */
     
}
#endif
     

void Substructure_StoneDrums_Destroy( StoneDrums_t *const Sub )
{
     free( Sub->Description );
}
