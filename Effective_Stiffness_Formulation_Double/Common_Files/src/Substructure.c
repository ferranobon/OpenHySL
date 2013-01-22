#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "Substructure.h"
#include "ErrorHandling.h"
#include "Netlib.h"
#include "Conf_Parser.h"

void Init_Constants_Substructure( ConstSub *const Constants, const char *Filename )
{
     ConfFile *Config;

     Config = ConfFile_Create( 45 );

     ConfFile_ReadFile( Config, Filename );

     Constants->Order_Couple = (unsigned int) ConfFile_GetInt( Config, "Substructure:Order" );
     Constants->Num_Steps = (unsigned int) ConfFile_GetInt( Config, "General:Num_Steps" );
     Constants->Num_Sub = (unsigned int) ConfFile_GetInt( Config, "Substructure:Num_Substeps");
     Constants->DeltaT = ConfFile_GetDouble( Config, "General:Delta" );
     Constants->DeltaT_Sub = (*Constants).DeltaT/(double)(*Constants).Num_Sub;
     
     ConfFile_Free( Config );
}

void Simulate_Substructure_Measured_Values( const char *FileName, const double *const Keinv, const double *const u0c, double *const uc, double *const fcprev, double *const fc, const unsigned int OrderC, const unsigned int NSub )
{
     unsigned int i, Substep;
     double ramp0, ramp;

     static double u0c0 = 0.0;
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

	  ramp = (double) Substep / (double) NSub;
	  ramp0 = 1 - ramp;   

	  uc[0] = ramp0*u0c0 + ramp*u0c[0] + Keinv[0]*fc[0];

	  /* Read the new value of fc */
	  fscanf( Fc_File, "%lf", &fc[0] );
     }
     u0c0 = u0c[0];
}

void Simulate_Substructures( Coupling_Node *const CNodes, double *Keinv, double *const u0c, double *const uc, double *const fcprev, double *const fc, const unsigned int NSubstep, const double DeltaT_Sub )
{

     unsigned int i, Substep;
     double ramp0, ramp;
     int incx = 1, incy = 1;
     double One;
     int Length;
     char uplo = 'L';
     TMD_Sim *TMD;
     UHYDE_Sim *UHYDE;  

     Length = (int) CNodes->Order;
     One = 1.0;

     for ( Substep = 1; Substep <= NSubstep; Substep++ ){

	  /* Backup data so that fcprev contains always the last coupling force */
	  dcopy_( &Length, fc, &incx, fcprev, &incy );
	       
	  ramp = (double) Substep / (double) NSubstep;

	  ramp0 = 1.0 - ramp;   

	  if ( CNodes->Order > 1 ){
	       dcopy_( &Length, CNodes->u0c0, &incx, uc, &incy );
	       dscal_( &Length, &ramp0, uc, &incx );
	       daxpy_( &Length, &ramp, u0c, &incx, uc, &incy );
	       dsymv_( &uplo, &Length, &One, Keinv, &Length, fc, &incx, &One, uc, &incy ); 
	  } else {
	       uc[0] = ramp0*CNodes->u0c0[0] + ramp*u0c[0] + Keinv[0]*fc[0];
	  }
	  
	  /* Compute the new fc */
	  for( i = 0; i < CNodes->Order; i ++ ){
	       if( CNodes->Sub[i].Type == USE_EXACT ){
		    TMD = (TMD_Sim *) CNodes->Sub[i].SimStruct;
		    ExactSolution_SDOF( u0c[i], DeltaT_Sub, &TMD[i], &fc[i] );
	       } else if ( CNodes->Sub[i].Type == USE_UHYDE ){
		    UHYDE = (UHYDE_Sim *) CNodes->Sub[i].SimStruct;
		    Simulate_UHYDE_1D( u0c[i], DeltaT_Sub, &UHYDE[i], &fc[i] );
	       } else assert( CNodes->Sub[i].Type < USE_EXACT || CNodes->Sub[i].Type > USE_UHYDE );
	  }
	  
     }

     /* Backup u0c */
     dcopy_( &Length, u0c, &incx, CNodes->u0c0, &incy );
}

void Simulate_Substructure( void *const Num, const int Mode, double *const Keinv, double *const u0c0, double *const u0c, double *const uc, double *const fcprev, double *const fc, const unsigned int OrderC, const unsigned int NSub, const double DeltaT_Sub )
{

     unsigned int i, Substep;
     double ramp0, ramp;
     int incx = 1, incy = 1;
     double One;
     int Length;
     char uplo = 'L';
     TMD_Sim *TMD;
     UHYDE_Sim *UHYDE;  

     Length = (int) OrderC;
     One = 1.0;

     for ( Substep = 1; Substep <= NSub; Substep++ ){


	  /* Backup data so that fcprev contains always the last coupling force */
	  dcopy_( &Length, fc, &incx, fcprev, &incy );
	  
	  ramp = (double) Substep / (double) NSub;

	  ramp0 = 1.0 - ramp;   

	  if ( OrderC > 1 ){
	       dcopy_( &Length, u0c0, &incx, uc, &incy );
	       dscal_( &Length, &ramp0, uc, &incx );
	       daxpy_( &Length, &ramp, u0c, &incx, uc, &incy );
	       dsymv_( &uplo, &Length, &One, Keinv, &Length, fc, &incx, &One, uc, &incy ); 
	  } else {
	       uc[0] = ramp0*u0c0[0] + ramp*u0c[0] + Keinv[0]*fc[0];
	  }
	  
	  /* Compute the new fc */
	  for( i = 0; i < OrderC; i ++ ){
	       if( Mode == USE_EXACT ){
		    TMD = (TMD_Sim *) Num;
		    ExactSolution_SDOF( u0c[i], DeltaT_Sub, &TMD[i], &fc[i] );
	       } else if ( Mode == USE_UHYDE ){
		    UHYDE = (UHYDE_Sim *) Num;
		    Simulate_UHYDE_1D( u0c[i], DeltaT_Sub, &UHYDE[i], &fc[i] );
	       } else assert( Mode < USE_EXACT || Mode > USE_UHYDE );
	  }
	  
     }

     /* Backup u0c */
     dcopy_( &Length, u0c, &incx, u0c0, &incy );
}


void ExactSolution_Init( const double Mass, const double Damp, const double Stiff, const double DeltaT, TMD_Sim *const Num )
{
     double omega, omegaD, omega2;
     double xi;
     double expdt, sinDdt, cosDdt, omegadt, omegaDdt;

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

	       Num->A = expdt*(xi*omega*sinDdt/omega + cosDdt);
	       Num->B = expdt*sinDdt/omegaD;
	       Num->C = (expdt*(((1.0-2.0*xi*xi)/omegaDdt-xi*omega/omegaD)*sinDdt-(1.0+2.0*xi/omegadt)*cosDdt)+2.0*xi/omegadt)/Stiff;
	       Num->D = (expdt*((2.0*xi*xi-1.0)*sinDdt/omegaDdt + 2.0*xi*cosDdt/omegadt)+(1.0-2.0*xi/omegadt))/Stiff;
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



void ExactSolution_SDOF( const double u0c, const double DeltaT, TMD_Sim *const Num, double *const fc )
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

void Simulate_UHYDE_1D_Init( const double qyield, const double yield_factor, const double Friction, UHYDE_Sim *const Num )
{

     Num->u0c_old = 0.0;
     Num->q = 0.0;

     Num->qyield = qyield;
     Num->qplastic = qyield/yield_factor;
     Num->k = Friction/Num->qplastic;

}

void Simulate_UHYDE_1D( const double u0c, const double DeltaT, UHYDE_Sim *const Num, double *const Friction_Force )
{

     double v;
     double hq;
  
     /* The notation here follows the one presented in Cascade Report No. 1 Seismic qualification
      * of passive mitigation devices page 40.
      */

     v = (u0c - Num->u0c_old)/DeltaT;

     if ( ( fabs(Num->q) <= Num->qyield) || Num->q*v <= 0.0 ){
	  hq = 1.0;
     } else if ( Num->qyield < fabs(Num->q) || Num->q*v > 0.0 ){
	  hq = (Num->qplastic - Num->q)/(Num->qplastic - Num->qyield);
     } else {
	  printf("qyield %e, abs(q) %e, q*v %e\n", Num->qyield, fabs(Num->q), Num->q*v );
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
