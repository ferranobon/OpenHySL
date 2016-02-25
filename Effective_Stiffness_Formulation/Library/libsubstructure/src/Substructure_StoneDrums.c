#include <stdlib.h>

#include "Substructure_StoneDrums.h"

void Substructure_StoneDrums_Init (const int PrevDOF, const char *Description, StoneDrums_t *const Sub )
{
     Sub->Description = strdup( Description );
  
     Sub->PrevDOF = PrevDOF;

#if _FLOAT_
     Sub->a0 = 4E4f;
     Sub->a2 = 4E2f;
     Sub->a3 = 1.0f;
     Sub->a6 = 5E-3f;
     Sub->a7 = 5E-3f;
#else
     Sub->a0 = 4E4;
     Sub->a2 = 4E2;
     Sub->a3 = 1.0;
     Sub->a6 = 5E-3;
     Sub->a7 = 5E-3;
#endif

     Sub->Result_Fc = 0.0;
     Sub->Disp0 = 0.0; Sub->DispT = 0.0;
     Sub->Vel0 = 0.0; Sub->VelT = 0.0; Sub->VelTdT = 0.0;
     Sub->Acc0 = 0.0; Sub->AccT = 0.0; Sub->AccTdT = 0.0;

}

void Substructure_StoneDrums ( const HYSL_FLOAT DispTdT, const HYSL_FLOAT ramp, StoneDrums_t * const Sub, HYSL_FLOAT *const fc )
{

#if _FLOAT_
     Sub->AccTdT = (1.0f - ramp)*(-Sub->a0*Sub->Disp0 -Sub->a2*Sub->Vel0 -Sub->a3*Sub->Acc0)
	  + ramp*(-Sub->a0*Sub->DispT -Sub->a2*Sub->VelT -Sub->a3*Sub->AccT) + Sub->a0*DispTdT;

     Sub->VelTdT = (1.0f - ramp)*(Sub->Vel0 + Sub->a6*Sub->Acc0) + ramp*(Sub->VelT + Sub->a6*Sub->AccT) + Sub->a7*Sub->AccTdT;
#else
     Sub->AccTdT = (1.0 - ramp)*(-Sub->a0*Sub->Disp0 -Sub->a2*Sub->Vel0 -Sub->a3*Sub->Acc0)
	  + ramp*(-Sub->a0*Sub->DispT -Sub->a2*Sub->VelT -Sub->a3*Sub->AccT) + Sub->a0*DispTdT;

     Sub->VelTdT = (1.0 - ramp)*(Sub->Vel0 + Sub->a6*Sub->Acc0) + ramp*(Sub->VelT + Sub->a6*Sub->AccT) + Sub->a7*Sub->AccTdT;
#endif

     if ( Sub->VelTdT > 0.0 ){
	  *fc = -35000;
     } else if (Sub->VelTdT < 0.0){
	  *fc = 35000;
     } else {
	  *fc = 0.0;
     }

     Sub->Result_Fc = *fc;
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
