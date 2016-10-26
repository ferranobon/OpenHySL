#include <stdlib.h>

#include <math.h>
#include <assert.h>
#include <string.h>

#include "Substructure_UHYDEfbr.h"
#include "Definitions.h"

void Substructure_SimUHYDE_1D_Init( const HYSL_FLOAT qyield, const HYSL_FLOAT yield_factor, const HYSL_FLOAT Friction, const char *Description, UHYDEfbrSim_t *const Num )
{

     Num->Description = strdup( Description );

     Num->u0c_old = 0.0;
     Num->q = 0.0;

     Num->qyield = qyield;
     Num->qplastic = qyield/yield_factor;
     Num->k = Friction/Num->qplastic;

}

void Substructure_SimUHYDE_1D( const HYSL_FLOAT u0c, const HYSL_FLOAT DeltaT, UHYDEfbrSim_t *const Num, HYSL_FLOAT *const Friction_Force )
{

     HYSL_FLOAT v;
     HYSL_FLOAT hq;
  
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
