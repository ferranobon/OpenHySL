#include <stdlib.h>

#include <math.h>
#include <assert.h>
#include <string.h>

#include "Substructure_UHYDEfbr.h"
#include "Definitions.h"

void Substructure_SimUHYDE_1D_Init( const hysl_float_t qyield, const hysl_float_t yield_factor, const hysl_float_t Friction, const char *Description, UHYDEfbrSim_t *const Num )
{

     Num->Description = strdup( Description );

     Num->u0c_old = 0.0;
     Num->q = 0.0;

     Num->qyield = qyield;
     Num->qplastic = qyield/yield_factor;
     Num->k = Friction/Num->qplastic;

}

void Substructure_SimUHYDE_1D( const hysl_float_t u0c, const hysl_float_t DeltaT, UHYDEfbrSim_t *const Num, hysl_float_t *const Friction_Force )
{

     hysl_float_t v;
     hysl_float_t hq;
  
     /* The notation here follows the one presented in Cascade Report No. 1 Seismic qualification
      * of passive mitigation devices page 40.
      */

     v = (u0c - Num->u0c_old)/DeltaT;

     if ( ( fabs(Num->q) <= Num->qyield) || Num->q*v <= 0.0 ){
	  hq = 1.0;
     } else if ( Num->qyield < fabs(Num->q) || Num->q*v > 0.0 ){
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

void Substructure_SimUHYDE_Destroy( UHYDEfbrSim_t *const Sub )
{
     free( Sub->Description );
}
