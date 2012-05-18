#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Substructure.h"

void Init_Constants_Substructure( ConstSub *const Constants )
{
     (*Constants).Order_Couple = 1;
     (*Constants).Num_Steps = 4096;
     (*Constants).Num_Sub = 4;
     (*Constants).DeltaT = 0.01;
     (*Constants).DeltaT_Sub = (*Constants).DeltaT/(float)(*Constants).Num_Sub;
}

void Simulate_Substructure( const float *const u0c, float *const uc, float *const fcprev, float *const fc, const int OrderC, const int NSub, const float Deltat_Sub )
{

     fc[0] = 0.0;
     fcprev[0] = 0.0;

     uc[0] = u0c[0];
}


void ExactSolution_Init( const float Mass, const float Damp, const float Stiff, const float DeltaT, TMD_Sim *const Num )
{
     float omega, omegaD, omega2;
     float xi;
     float expdt, sinDdt, cosDdt, omegadt, omegaDdt;

     /* Init the variables */
     Num->Disp0 = 0.0;
     Num->Vel0 = 0.0;

     Num->Disp = 0.0;
     Num->Vel = 0.0;

     Num->Force_0 = 0.0;
     Num->Force_1 = 0.0;

     Num->u0c_old = 0.0;
     Num->v0c = 0.0;

     /* Calculate the free vibration frequency */
     omega2 = Stiff/Mass;

     /* Check the value of omega */
     if ( omega2 > 0.0 ){
	  omega = sqrtf( omega2 );

	  /* Calculate the damping ratio */
	  xi = Damp/(2*Mass*omega);
	  
	  /* Check the value of xi */
	  if ( xi >= 0.0 && xi < 1 ){
	       /* Calculate the damping vibration frequency */
	       omegaD = omega*sqrtf(1-xi*xi);

	       /* Calculate several constants */
	       expdt = exp(-xi*omega*DeltaT);
	       sinDdt = sin(omegaD*DeltaT);
	       cosDdt = cos(omegaD*DeltaT);
	       omegadt = omega*DeltaT;
	       omegaDdt = omegaD*DeltaT;

	       Num->A = expdt*(xi*omega*sinDdt/omega + cosDdt);
	       Num->B = expdt*sinDdt/omegaD;
	       Num->C = (expdt*(((1.0-2.0*xi*xi)/omegaDdt-xi*omega/omegaD)*sinDdt-(1.0+2.0*xi/omegadt)*cosDdt)+2.0*xi/omegadt)/Stiff;
	       Num->D = (expdt*((2.0*xi*xi-1.0)*sinDdt/omegaDdt + 2.0*xi*cosDdt/omegadt)+(1.0-2*xi/omegadt))/Stiff;
	       Num->E = -expdt * omega * omega * sinDdt / omegaD;
	       Num->F = expdt * (cosDdt - xi*omega*sinDdt/omegaD);
	       Num->G = (expdt*((omega*omega/omegaD + xi*omega/omegaDdt)*sinDdt + cosDdt/DeltaT)-1.0/DeltaT)/Stiff;
	       Num->H = (-expdt*(xi*omega*sinDdt/omegaD + cosDdt) + 1.0)/(Stiff*DeltaT);

	   } else {
	       fprintf( stderr, "Exact Solution: Negative or overcritical damping.\n" );
	       exit( EXIT_FAILURE );
	  }
     } else {
	  fprintf( stderr, "Exact Solution: Negative vibration frerquency.\n" );
	  exit( EXIT_FAILURE );
     }
}



void ExactSolution_SDOF( const float Mass, const float Damp, const float Stiff, const float u0c, const float DeltaT, TMD_Sim *const Num, float *const fc )
{

     /* Compute initial velocity */
     Num->v0c = (u0c - Num->u0c_old)/DeltaT;

     /* Backup and computer new forcce */
     Num->Force_0 = Num->Force_1;
     Num->Force_1 = Stiff*u0c + Damp*Num->v0c;

     /* Backup the displacements */
     Num->Disp0 = Num->Disp;
     Num->Vel0 = Num->Vel;

     /* Compute the new state */
     Num->Disp = Num->A*Num->Disp0 + Num->B*Num->Vel0 + Num->C*Num->Force_0 + Num->D*Num->Force_1;
     Num->Vel = Num->E*Num->Disp0 + Num->F*Num->Vel0 + Num->G*Num->Force_0 + Num->H*Num->Force_1;

     /* Compute the coupling force */
     *fc = Stiff*(Num->Disp - u0c) + Damp*(Num->Vel - Num->v0c);

     /* Backup u0c vector */
     Num->u0c_old = u0c;
	 
}

void Simulate_UHYDE( const float u0c, const float Deltat_Sub, float *const Friction_Force ){
  
}
