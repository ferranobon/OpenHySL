#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>

#include "Substructure_BoucWen.h"
#include "Auxiliary_Math.h"
#include "Print_Messages.h"
#include "Definitions.h"

void Substructure_BoucWen_Init (const HYSL_FLOAT alpha, const HYSL_FLOAT ko, const HYSL_FLOAT Fy,
				const HYSL_FLOAT beta, const HYSL_FLOAT gamma, const HYSL_FLOAT n,
				const HYSL_FLOAT A0, const HYSL_FLOAT deltaA, const HYSL_FLOAT nu0,
				const HYSL_FLOAT deltaNu, const HYSL_FLOAT eta0, const HYSL_FLOAT deltaEta,
				const HYSL_FLOAT vs0, const HYSL_FLOAT p, const HYSL_FLOAT q,
			        const HYSL_FLOAT lambda, const HYSL_FLOAT psi0, const HYSL_FLOAT deltaPsi,
				const int BoucWen_Type, const char *Description, BoucWen_t *const Sub )
{
     Sub->Description = strdup( Description );
  
     Sub->BoucWen_Type = BoucWen_Type;

     /* Basic Bouc-Wen parameters */
     Sub->alpha = alpha;
     Sub->ko = ko;
     Sub->Fy = Fy;
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

     /* Parameters for pinching modelling */
     Sub->vs0 = vs0;
     Sub->p = p;
     Sub->q = q;
     Sub->psi0 = psi0;
     Sub->deltaPsi = deltaPsi;
     Sub->lambda = lambda;
     
     /* Initialisation of Newton-Raphson variables */
     Sub->maxIter = 100;
#if _FLOAT_
     Sub->tolerance = 1E-12f;
     Sub->startPoint = 0.01f;
     Sub->z_old = 0.0f;
     Sub->e_old = 0.0f;
     Sub->DispT = 0.0f;
#else 
     Sub->tolerance = 1E-12;
     Sub->startPoint = 0.01;
     Sub->z_old = 0.0;
     Sub->e_old = 0.0;
     Sub->DispT = 0.0;
#endif

}

void Substructure_BoucWenSurface_Init (const HYSL_FLOAT alpha, const HYSL_FLOAT ko,
				       const HYSL_FLOAT Fy, const HYSL_FLOAT beta, const HYSL_FLOAT gamma, const HYSL_FLOAT n,
				       const HYSL_FLOAT A0, const HYSL_FLOAT deltaA, const HYSL_FLOAT nu0,
				       const HYSL_FLOAT deltaNu, const HYSL_FLOAT eta0, const HYSL_FLOAT deltaEta,
				       const int BoucWen_Type, const char *Description, BoucWenSurface_t *const Sub )
{
     int i;
     
     Sub->Description = strdup( Description );
  
     Sub->BoucWen_Type = BoucWen_Type;

     /* Basic Bouc-Wen parameters */
     Sub->alpha = alpha;
     Sub->ko = ko;
     Sub->Fy = Fy;
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

     Sub->vs0 = 0.0;
     Sub->p = 0.0;
     Sub->q = 0.0;
     Sub->psi0 = 0.0;
     Sub->deltaPsi = 0.0;
     Sub->lambda = 0.0;
     
     /* Initialisation of Newton-Raphson variables */
     Sub->maxIter = 100;
#if _FLOAT_
     Sub->tolerance = 1E-12f;
     for ( i = 0; i < 2; i++){
	  Sub->z_old[i] = 0.0f;
	  Sub->z[i] = 0.0f;
	  Sub->DispT[i] = 0.0f;
     }
#else 
     Sub->tolerance = 1E-12;
     for ( i = 0; i < 2; i++){
	  Sub->z_old[i] = 0.0;
	  Sub->z[i] = 0.0;
	  Sub->DispT[i] = 0.0;
     }
#endif

}

void Substructure_BoucWen ( const HYSL_FLOAT DispTdT, BoucWen_t *const Sub, HYSL_FLOAT *const force )
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
	  e_new = Sub->e_old + (1.0 - Sub->alpha)*Sub->ko/Sub->Fy*deltaDisp*z_eval;
	  A_new = Sub->A0 - Sub->deltaA*e_new;
	  nu_new = Sub->nu0 + Sub->deltaNu*e_new;
	  eta_new = Sub->eta0 + Sub->deltaEta*e_new;
	  
	  sign = signum(deltaDisp*z_eval);
	  Theta = Sub->gamma + Sub->beta*sign;
	  Phi = A_new - pow(fabs(z_eval),Sub->n)*Theta*nu_new;

	  if (Sub->BoucWen_Type == BOUC_WEN_BABER_NOORI){
	       vs_1 = Sub->vs0*(1.0 - exp(-Sub->p*e_new));
	       vs_2 = (Sub->psi0 + Sub->deltaPsi*e_new)*(Sub->lambda + vs_1);
	       zu = pow(1.0/((Sub->nu0 + Sub->deltaNu*e_new)*(Sub->beta + Sub->gamma)), 1.0/Sub->n);
	       hze_exp = exp(-pow(z_eval*signum(deltaDisp) -Sub->q*zu, 2.0)/pow(vs_2, 2.0));
	       hze = 1.0 - vs_1*hze_exp;

	       fz_new = z_eval - Sub->z_old - hze*Phi/eta_new*deltaDisp;
	  } else {
	       fz_new = z_eval - Sub->z_old - Phi/eta_new*deltaDisp*Sub->ko/Sub->Fy;
	  }
	  
	  /* Evaluate function derivatives with respect to z_eval for the Newton-Rhapson scheme */
	  e_new_ = (1.0 - Sub->alpha)*Sub->ko/Sub->Fy*deltaDisp;
	  A_new_ = -Sub->deltaA*e_new_;
	  nu_new_ = Sub->deltaNu*e_new_;
	  eta_new_ = Sub->deltaEta*e_new_;
	  sign = signum(z_eval);
	  Phi_ = A_new_ - Sub->n*pow(fabs(z_eval), Sub->n - 1.0)*sign*Theta*nu_new - pow(fabs(z_eval), Sub->n)*Theta*nu_new_;

	  if (Sub->BoucWen_Type == BOUC_WEN_BABER_NOORI){
	       hze_ = deltaDisp*Sub->ko/Sub->Fy*Sub->p*Sub->vs0*exp(-pow(z_eval*signum(deltaDisp) - Sub->q*pow(1.0/((Sub->beta + Sub->gamma)*(Sub->nu0 + Sub->deltaNu*(Sub->e_old - deltaDisp*Sub->ko/Sub->Fy*z_eval*(Sub->alpha - 1.0)))), 1.0/Sub->n), 2.0)/(pow(Sub->lambda - Sub->vs0*(exp(-Sub->p*(Sub->e_old - deltaDisp*Sub->ko/Sub->Fy*z_eval*(Sub->alpha - 1.0))) - 1.0), 2.0)*pow(Sub->psi0 + Sub->deltaPsi*(Sub->e_old - deltaDisp*Sub->ko/Sub->Fy*z_eval*(Sub->alpha - 1.0)), 2.0)))*exp(-Sub->p*(Sub->e_old - deltaDisp*Sub->ko/Sub->Fy*z_eval*(Sub->alpha - 1.0)))*(Sub->alpha - 1.0) - Sub->vs0*exp(-pow(z_eval*signum(deltaDisp) - Sub->q*pow(1.0/((Sub->beta + Sub->gamma)*(Sub->nu0 + Sub->deltaNu*(Sub->e_old - deltaDisp*Sub->ko/Sub->Fy*z_eval*(Sub->alpha - 1.0)))),1.0/Sub->n), 2.0)/(pow(Sub->lambda - Sub->vs0*(exp(-Sub->p*(Sub->e_old - deltaDisp*Sub->ko/Sub->Fy*z_eval*(Sub->alpha - 1.0))) - 1.0), 2.0)*pow(Sub->psi0 + Sub->deltaPsi*(Sub->e_old - deltaDisp*Sub->ko/Sub->Fy*z_eval*(Sub->alpha - 1.0)), 2.0)))*(exp(-Sub->p*(Sub->e_old - deltaDisp*Sub->ko/Sub->Fy*z_eval*(Sub->alpha - 1.0))) - 1.0)*((2*(signum(deltaDisp) - (Sub->deltaNu*deltaDisp*Sub->ko/Sub->Fy*Sub->q*(Sub->alpha - 1.0)*pow(1.0/((Sub->beta + Sub->gamma)*(Sub->nu0 + Sub->deltaNu*(Sub->e_old - deltaDisp*Sub->ko/Sub->Fy*z_eval*(Sub->alpha - 1.0)))),1.0/Sub->n - 1.0))/(Sub->n*(Sub->beta + Sub->gamma)*pow(Sub->nu0 + Sub->deltaNu*(Sub->e_old - deltaDisp*Sub->ko/Sub->Fy*z_eval*(Sub->alpha - 1.0)), 2.0)))*(z_eval*signum(deltaDisp) - Sub->q*pow(1.0/((Sub->beta + Sub->gamma)*(Sub->nu0 + Sub->deltaNu*(Sub->e_old - deltaDisp*Sub->ko/Sub->Fy*z_eval*(Sub->alpha - 1.0)))),1.0/Sub->n)))/(pow(Sub->lambda - Sub->vs0*(exp(-Sub->p*(Sub->e_old - deltaDisp*Sub->ko/Sub->Fy*z_eval*(Sub->alpha - 1.0))) - 1.0), 2.0)*pow(Sub->psi0 + Sub->deltaPsi*(Sub->e_old - deltaDisp*Sub->ko/Sub->Fy*z_eval*(Sub->alpha - 1.0)), 2.0)) + (2*Sub->deltaPsi*deltaDisp*Sub->ko/Sub->Fy*(Sub->alpha - 1.0)*pow(z_eval*signum(deltaDisp) - Sub->q*pow(1.0/((Sub->beta + Sub->gamma)*(Sub->nu0 + Sub->deltaNu*(Sub->e_old - deltaDisp*Sub->ko/Sub->Fy*z_eval*(Sub->alpha - 1.0)))),1.0/Sub->n), 2.0))/(pow(Sub->lambda - Sub->vs0*(exp(-Sub->p*(Sub->e_old - deltaDisp*Sub->ko/Sub->Fy*z_eval*(Sub->alpha - 1.0))) - 1.0), 2.0)*pow(Sub->psi0 + Sub->deltaPsi*(Sub->e_old - deltaDisp*Sub->ko/Sub->Fy*z_eval*(Sub->alpha - 1.0)),3)) + (2*deltaDisp*Sub->ko/Sub->Fy*Sub->p*Sub->vs0*exp(-Sub->p*(Sub->e_old - deltaDisp*Sub->ko/Sub->Fy*z_eval*(Sub->alpha - 1.0)))*(Sub->alpha - 1.0)*pow(z_eval*signum(deltaDisp) - Sub->q*pow(1.0/((Sub->beta + Sub->gamma)*(Sub->nu0 + Sub->deltaNu*(Sub->e_old - deltaDisp*Sub->ko/Sub->Fy*z_eval*(Sub->alpha - 1.0)))),1.0/Sub->n), 2.0))/(pow(Sub->lambda - Sub->vs0*(exp(-Sub->p*(Sub->e_old - deltaDisp*Sub->ko/Sub->Fy*z_eval*(Sub->alpha - 1.0))) - 1.0),3.0)*pow(Sub->psi0 + Sub->deltaPsi*(Sub->e_old - deltaDisp*Sub->ko/Sub->Fy*z_eval*(Sub->alpha - 1.0)), 2.0)));

	       fz_new_ = 1.0 - (hze_*Phi*eta_new + hze*Phi_*eta_new - hze*Phi*eta_new_)/pow(eta_new, 2.0)*deltaDisp;	       
	  } else {
	       fz_new_ = 1.0 - (Phi_*eta_new - Phi*eta_new_)/pow(eta_new, 2.0)*deltaDisp*Sub->ko/Sub->Fy;
	  }

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
	  *force = -Sub->alpha*Sub->ko*DispTdT - (1.0 - Sub->alpha)*Sub->Fy*z_eval;
	  
	  // Compute material degradation parameters.
	  e_new = Sub->e_old + (1.0 - Sub->alpha)*Sub->ko/Sub->Fy*deltaDisp*z_eval;
	  A_new = Sub->A0 - Sub->deltaA*e_new;
	  nu_new = Sub->nu0 + Sub->deltaNu*e_new;
	  eta_new = Sub->eta0 + Sub->deltaEta*e_new;
     }

     Sub->DispT = DispTdT;
     Sub->e_old = e_new;
     Sub->z_old = z_eval;
     Sub->BW_Force = -(*force);
}

void Substructure_BoucWenSurface ( const HYSL_FLOAT DispTdT1, const HYSL_FLOAT DispTdT2, BoucWenSurface_t *const Sub, HYSL_FLOAT *const force1, HYSL_FLOAT *const force2 )
{
     int count, i;
     HYSL_FLOAT fz_new[2], fz_new_[4];
     HYSL_FLOAT delta_z[2];
     HYSL_FLOAT deltaDisp[2];
     HYSL_FLOAT z_norm;
     HYSL_FLOAT tmp1, tmp2, tmp3, tmp4;
     
     count = 0;
     
     deltaDisp[0] = DispTdT1 - Sub->DispT[0];
     deltaDisp[1] = DispTdT2 - Sub->DispT[1];

     do {
	  z_norm = norm( 2, Sub->z );
	  if( z_norm == 0.0 ){
	       z_norm = DBL_EPSILON;
	  }
	  
	  /* Use Park's formulation */
	  tmp1 = Sub->z[0]*deltaDisp[0] + Sub->z[1]*deltaDisp[1];
	  tmp2 = Sub->gamma + Sub->beta*signum(tmp1);
	  tmp3 = pow(z_norm, Sub->n - 2.0)*tmp1*tmp2;
	  tmp4 = pow(z_norm, Sub->n - 4.0)*tmp2;
	  
	  // Function f(z)
	  fz_new[0] = Sub->z[0] - Sub->z_old[0] - (Sub->A0*deltaDisp[0] - Sub->z[0]*tmp3)*Sub->ko/Sub->Fy;
	  fz_new[1] = Sub->z[1] - Sub->z_old[1] - (Sub->A0*deltaDisp[1] - Sub->z[1]*tmp3)*Sub->ko/Sub->Fy;
	  
	  fz_new_[0] = 1.0 + tmp4*Sub->ko/Sub->Fy*(pow(Sub->z[1], 3.0)*deltaDisp[1] + 2.0*pow(Sub->z[1], 2.0)*Sub->z[0]*deltaDisp[0]
						   + Sub->z[1]*pow(Sub->z[0], 2.0)*deltaDisp[1]*(Sub->n - 1.0) + pow(Sub->z[0], 3.0)*deltaDisp[0]*Sub->n);
	  fz_new_[1] = Sub->z[1]*tmp4*Sub->ko/Sub->Fy*(pow(Sub->z[1], 2.0)*deltaDisp[0] + Sub->z[0]*Sub->z[1]*deltaDisp[1]*(Sub->n - 2.0)
						  + pow(Sub->z[0], 2.0)*deltaDisp[0]*(Sub->n - 1.0));
	  fz_new_[2] = Sub->z[0]*tmp4*Sub->ko/Sub->Fy*(pow(Sub->z[0],2.0)*deltaDisp[1] + Sub->z[0]*Sub->z[1]*deltaDisp[0]*(Sub->n - 2.0)
						  + pow(Sub->z[1], 2.0)*deltaDisp[1]*(Sub->n - 1.0));
	  fz_new_[3] = 1.0 + tmp4*Sub->ko/Sub->Fy*(pow(Sub->z[0], 3.0)*deltaDisp[0] + 2.0*pow(Sub->z[0],2.0)*Sub->z[1]*deltaDisp[1]
						   + Sub->z[0]*pow(Sub->z[1], 2.0)*deltaDisp[0]*(Sub->n - 1.0) + pow(Sub->z[1], 3.0)*deltaDisp[1]*Sub->n);

	  if ((fabs(fz_new_[0]) <= DBL_EPSILON) || (fabs(fz_new_[3]) <= DBL_EPSILON)){
	       Print_Header( WARNING );
	       fprintf( stderr, "WARNING: zero Jacobian in Newton-Raphson scheme for hysteretic evolution parameter z.\n");
         
	  }
	  
	  right_matrix_division( fz_new, fz_new_, delta_z );
	  
	  Sub->z[0] = Sub->z[0] - delta_z[0];
	  Sub->z[1] = Sub->z[1] - delta_z[1];

	  count = count + 1;
	       
     } while( (norm(2, delta_z) >= Sub->tolerance) && (count < Sub->maxIter) ); 

     Sub->z_old[0] = Sub->z[0]; Sub->z_old[1] = Sub->z[1];
     Sub->DispT[0] = DispTdT1; Sub->DispT[1] = DispTdT2;

     *force1 = -Sub->alpha*Sub->ko*DispTdT1 - (1.0 - Sub->alpha)*Sub->Fy*Sub->z[0];
     *force2 = -Sub->alpha*Sub->ko*DispTdT2 - (1.0 - Sub->alpha)*Sub->Fy*Sub->z[1];

     Sub->BW_Force[0] = -(*force1);
     Sub->BW_Force[1] = -(*force2);
}

void right_matrix_division (const HYSL_FLOAT *const vector, const HYSL_FLOAT *const matrix, HYSL_FLOAT *const output)
{
     int i;
     HYSL_FLOAT inv[4];
     HYSL_FLOAT detM;

     /* Calculate the determinant of a 2x2 matrix */
     detM = matrix[0]*matrix[3] - matrix[1]*matrix[2];

     /* Calculate the inverse of a 2x2 matrix */
     inv[0] = matrix[3];
     inv[1] = -matrix[1];
     inv[2] = -matrix[2];
     inv[3] = matrix[0];

     for (i = 0; i < 4; i++){
	  inv[i] = inv[i]/detM;
     }

     /* Perform the right matrix division */
     output[0] = vector[0]*inv[0] + vector[1]*inv[2];
     output[1] = vector[0]*inv[1] + vector[1]*inv[3];
}

void Substructure_BoucWen_Destroy( BoucWen_t *const Sub )
{
     free( Sub->Description );
}

void Substructure_BoucWenSurface_Destroy( BoucWenSurface_t *const Sub )
{
     free( Sub->Description );
}
