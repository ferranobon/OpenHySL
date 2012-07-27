#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "Substructure.h"
#include "ErrorHandling.h"

void Init_Constants_Substructure( ConstSub *const Constants )
{
     (*Constants).Order_Couple = 1;
     (*Constants).Num_Steps = 4096;
     (*Constants).Num_Sub = 4;
     (*Constants).DeltaT = 0.01f;
     (*Constants).DeltaT_Sub = (*Constants).DeltaT/(float)(*Constants).Num_Sub;
}

void Simulate_Substructure_Measured_Values( const char *FileName, const float *const Keinv, const float *const u0c, float *const uc, float *const fcprev, float *const fc, const unsigned int OrderC, const unsigned int NSub )
{
     unsigned int i, Substep;
     float ramp0, ramp;

     static float u0c0 = 0.0;
     static int Is_First_Time = 1;
     static FILE *Fc_File;

     if( Is_First_Time ){
	  Fc_File = fopen( FileName, "r" );
   
	  if ( Fc_File == NULL ){
	       ErrorFileAndExit( "Cannot open the file: ", FileName );
	  } else {
	       Is_First_Time = 0;
	  }

     }

     for ( Substep = 1; Substep <= NSub; Substep++ ){

	  for ( i = 0; i < OrderC; i++ ){
	       /* Backup data so that fcprev contains always the last coupling force */
	       fcprev[i] = fc[i];
	  }

	  ramp = (float) Substep / (float) NSub;
	  ramp0 = 1 - ramp;   

	  uc[0] = ramp0*u0c0 + ramp*u0c[0] + Keinv[0]*fc[0];

	  /* Read the new value of fc */
	  fscanf( Fc_File, "%f", &fc[0] );
     }
}

void Simulate_Substructure( void *const Num, const int Mode, const float *const Keinv, const float *const u0c, float *const uc, float *const fcprev, float *const fc, const unsigned int OrderC, const unsigned int NSub, const float DeltaT_Sub )
{

     unsigned int i, Substep;
     float ramp0, ramp;
     static float u0c0 = 0.0;
     
     for ( Substep = 1; Substep <= NSub; Substep++ ){

	  for ( i = 0; i < OrderC; i++ ){
	       /* Backup data so that fcprev contains always the last coupling force */
	       fcprev[i] = fc[i];
	  }

	  ramp = (float) Substep / (float) NSub;
	  ramp0 = 1.0f - ramp;   

	  uc[0] = ramp0*u0c0 + ramp*u0c[0] + Keinv[0]*fc[0];

	  /* Compute the new fc */
	  if( Mode == USE_EXACT ){
	       ExactSolution_SDOF( u0c[0], DeltaT_Sub, (TMD_Sim *) Num, &fc[0] );
	  } else if ( Mode == USE_UHYDE ){
	       Simulate_UHYDE_1D( u0c[0], DeltaT_Sub, (UHYDE_Sim *) Num, &fc[0] );
	  } else assert( Mode < USE_EXACT || Mode > USE_UHYDE );
     }

     /* Backup u0c */
     u0c0 = u0c[0];
}


void ExactSolution_Init( const float Mass, const float Damp, const float Stiff, const float DeltaT, TMD_Sim *const Num )
{
     float omega, omegaD, omega2;
     float xi;
     float expdt, sinDdt, cosDdt, omegadt, omegaDdt;

     /* Init the variables */
     Num->Mass = Mass;
     Num->Damp = Damp;
     Num->Stiff = Stiff;

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
	       expdt = expf(-xi*omega*DeltaT);
	       sinDdt = sinf(omegaD*DeltaT);
	       cosDdt = cosf(omegaD*DeltaT);
	       omegadt = omega*DeltaT;
	       omegaDdt = omegaD*DeltaT;

	       Num->A = expdt*(xi*omega*sinDdt/omega + cosDdt);
	       Num->B = expdt*sinDdt/omegaD;
	       Num->C = (expdt*(((1.0f-2.0f*xi*xi)/omegaDdt-xi*omega/omegaD)*sinDdt-(1.0f+2.0f*xi/omegadt)*cosDdt)+2.0f*xi/omegadt)/Stiff;
	       Num->D = (expdt*((2.0f*xi*xi-1.0f)*sinDdt/omegaDdt + 2.0f*xi*cosDdt/omegadt)+(1.0f-2.0f*xi/omegadt))/Stiff;
	       Num->E = -expdt * omega * omega * sinDdt / omegaD;
	       Num->F = expdt * (cosDdt - xi*omega*sinDdt/omegaD);
	       Num->G = (expdt*((omega*omega/omegaD + xi*omega/omegaDdt)*sinDdt + cosDdt/DeltaT)-1.0f/DeltaT)/Stiff;
	       Num->H = (-expdt*(xi*omega*sinDdt/omegaD + cosDdt) + 1.0f)/(Stiff*DeltaT);

	   } else {
	       fprintf( stderr, "Exact Solution: Negative or overcritical damping.\n" );
	       exit( EXIT_FAILURE );
	  }
     } else {
	  fprintf( stderr, "Exact Solution: Negative vibration frerquency.\n" );
	  exit( EXIT_FAILURE );
     }
}



void ExactSolution_SDOF( const float u0c, const float DeltaT, TMD_Sim *const Num, float *const fc )
{

     /* Compute initial velocity */
     Num->v0c = (u0c - Num->u0c_old)/DeltaT;

     /* Backup and computer new forcce */
     Num->Force_0 = Num->Force_1;
     Num->Force_1 = Num->Stiff*u0c + Num->Damp*Num->v0c;

     /* Backup the displacements */
     Num->Disp0 = Num->Disp;
     Num->Vel0 = Num->Vel;

     /* Compute the new state */
     Num->Disp = Num->A*Num->Disp0 + Num->B*Num->Vel0 + Num->C*Num->Force_0 + Num->D*Num->Force_1;
     Num->Vel = Num->E*Num->Disp0 + Num->F*Num->Vel0 + Num->G*Num->Force_0 + Num->H*Num->Force_1;

     /* Compute the coupling force */
     *fc = Num->Stiff*(Num->Disp - u0c) + Num->Damp*(Num->Vel - Num->v0c);

     /* Backup u0c vector */
     Num->u0c_old = u0c;
	 
}

void Simulate_UHYDE_1D_Init( const float qyield, const float yield_factor, const float Friction, UHYDE_Sim *const Num )
{

     Num->u0c_old = 0.0;
     Num->q = 0.0;

     Num->qyield = qyield;
     Num->qplastic = qyield/yield_factor;
     Num->k = Friction/Num->qplastic;

}

void Simulate_UHYDE_1D( const float u0c, const float DeltaT, UHYDE_Sim *const Num, float *const Friction_Force )
{

     float v;
     float hq;
  
     /* The notation here follows the one presented in Cascade Report No. 1 Seismic qualification
      * of passive mitigation devices page 40.
      */

     v = (u0c - Num->u0c_old)/DeltaT;

     if ( ( fabsf(Num->q) <= Num->qyield) || Num->q*v <= 0.0 ){
	  hq = 1.0;
     } else if ( Num->qyield < fabsf(Num->q) || Num->q*v > 0.0 ){
	  hq = (Num->qplastic - Num->q)/(Num->qplastic - Num->qyield);
     } else {
	  assert( 0 );
     }

     Num->q = Num->q + hq*(u0c - Num->u0c_old);

     if( Num->q > Num->qplastic ){
	  Num->q = Num->qplastic;
     } else if ( Num->q < -Num->qplastic ){
	  Num->q = -Num->qplastic;
     }

     *Friction_Force = -Num->k*Num->q;
     
     /* Backup the displacement */
     Num->u0c_old = u0c;

}
