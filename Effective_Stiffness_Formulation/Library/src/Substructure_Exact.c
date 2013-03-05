#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "Print_Messages.h"
#include "Substructure_Exact.h"

void Substructure_ExactSolution_Init( const double Mass, const double Damp, const double Stiff, const double DeltaT, const char *Description, ExactSim_t *const Sub )
{
     double omega, omegaD, omega2;
     double xi;
     double expdt, sinDdt, cosDdt, omegadt, omegaDdt;

     Sub->Description = strdup( Description );

     /* Init the variables */
     Sub->Mass = Mass;
     Sub->Damp = Damp;
     Sub->Stiff = Stiff;

     Sub->Disp0 = 0.0;
     Sub->Vel0 = 0.0;

     Sub->Disp = 0.0;
     Sub->Vel = 0.0;

     Sub->Force_0 = 0.0;
     Sub->Force_1 = 0.0;

     Sub->u0c_old = 0.0;
     Sub->v0c = 0.0;

     /* Calculate the free vibration frequency */
     omega2 = Stiff/Mass;

     /* Check the value of omega */
     if ( omega2 > 0.0 ){
	  omega = sqrt( omega2 );

	  /* Calculate the damping ratio */
	  xi = Damp/(2*Mass*omega);
	  
	  /* Check the value of xi */
	  if ( xi >= 0.0 && xi < 1.0 ){
	       /* Calculate the damping vibration frequency */
	       omegaD = omega*sqrt(1.0-xi*xi);

	       /* Calculate several constants */
	       expdt = exp(-xi*omega*DeltaT);
	       sinDdt = sin(omegaD*DeltaT);
	       cosDdt = cos(omegaD*DeltaT);
	       omegadt = omega*DeltaT;
	       omegaDdt = omegaD*DeltaT;

	       Sub->A = expdt*(xi*omega*sinDdt/omega + cosDdt);
	       Sub->B = expdt*sinDdt/omegaD;
	       Sub->C = (expdt*(((1.0-2.0*xi*xi)/omegaDdt-xi*omega/omegaD)*sinDdt-(1.0+2.0*xi/omegadt)*cosDdt)+2.0*xi/omegadt)/Stiff;
	       Sub->D = (expdt*((2.0*xi*xi-1.0)*sinDdt/omegaDdt + 2.0*xi*cosDdt/omegadt)+(1.0-2.0*xi/omegadt))/Stiff;
	       Sub->E = -expdt * omega * omega * sinDdt / omegaD;
	       Sub->F = expdt * (cosDdt - xi*omega*sinDdt/omegaD);
	       Sub->G = (expdt*((omega*omega/omegaD + xi*omega/omegaDdt)*sinDdt + cosDdt/DeltaT)-1.0/DeltaT)/Stiff;
	       Sub->H = (-expdt*(xi*omega*sinDdt/omegaD + cosDdt) + 1.0)/(Stiff*DeltaT);

	   } else {
	       Print_Header( ERROR );
	       fprintf( stderr, "Exact Solution: Negative or overcritical damping.\n" );
	       exit( EXIT_FAILURE );
	  }
     } else {
	  Print_Header( ERROR );
	  fprintf( stderr, "Exact Solution: Negative vibration frerquency.\n" );
	  exit( EXIT_FAILURE );
     }
}

void Substructure_ExactSolution_SDOF( const double u0c, const double DeltaT, ExactSim_t *const Sub, double *const fc )
{

     /* Compute initial velocity */
     Sub->v0c = (u0c - Sub->u0c_old)/DeltaT;

     /* Backup and computer new forcce */
     Sub->Force_0 = Sub->Force_1;
     Sub->Force_1 = Sub->Stiff*u0c + Sub->Damp*Sub->v0c;

     /* Backup the displacements */
     Sub->Disp0 = Sub->Disp;
     Sub->Vel0 = Sub->Vel;

     /* Compute the new state */
     Sub->Disp = Sub->A*Sub->Disp0 + Sub->B*Sub->Vel0 + Sub->C*Sub->Force_0 + Sub->D*Sub->Force_1;
     Sub->Vel = Sub->E*Sub->Disp0 + Sub->F*Sub->Vel0 + Sub->G*Sub->Force_0 + Sub->H*Sub->Force_1;

     /* Compute the coupling force */
     *fc = Sub->Stiff*(Sub->Disp - u0c) + Sub->Damp*(Sub->Vel - Sub->v0c);

     /* Backup u0c vector */
     Sub->u0c_old = u0c;
	 
}


void Substructure_ExactSolution_Destroy( ExactSim_t *const Sub )
{
     free( Sub->Description );
}
