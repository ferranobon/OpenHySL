#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "Print_Messages.h"
#include "Substructure_Newmark.h"
#include "MatrixVector.h"
#include "Auxiliary_Math.h"
#include "Definitions.h"

#if _MKL_
#include "mkl_blas.h"
#else
#include "Netlib.h"
#endif

void Substructure_Newmark_Init( const HYSL_FLOAT Mass, const HYSL_FLOAT Damp, const HYSL_FLOAT Stiff, const HYSL_FLOAT DeltaTSub, HYSL_FLOAT DeltaT,
				const HYSL_FLOAT Beta, const HYSL_FLOAT Gamma, const char *Description, NewmarkSim_t *const Sub )
{

     Sub->Description = strdup( Description );

     /* Constants for Ending Step */
#if _FLOAT_
     Sub->a0 = 1.0f/(Beta*DeltaTSub*DeltaTSub);
     Sub->a2 = 1.0f/(Beta*DeltaTSub);
     Sub->a3 = 1.0f/(2.0f*Beta) - 1.0f;
     Sub->a4 = Gamma/Beta - 1.0f;
     Sub->a5 = (DeltaTSub/2.0f)*(Gamma/Beta - 2.0f);
     Sub->a6 = (1.0f - Gamma)*DeltaTSub;
     Sub->a10 = (0.5f - Beta)*DeltaTSub*DeltaTSub;
#else
     Sub->a0 = 1.0/(Beta*DeltaTSub*DeltaTSub);
     Sub->a2 = 1.0/(Beta*DeltaTSub);
     Sub->a3 = 1.0/(2.0*Beta) - 1.0;
     Sub->a4 = Gamma/Beta - 1.0;
     Sub->a5 = (DeltaTSub/2.0)*(Gamma/Beta - 2.0);
     Sub->a6 = (1.0 - Gamma)*DeltaTSub;
     Sub->a10 = (0.5 - Beta)*DeltaTSub*DeltaTSub;
#endif

     Sub->a1 = Gamma/(Beta*DeltaTSub);
     Sub->a7 = Gamma*DeltaTSub;
     Sub->a8 = Beta*DeltaTSub*DeltaTSub;
     Sub->a9 = DeltaTSub;


#if _FLOAT_
     Sub->A0 = 1.0f/(Beta*DeltaT*DeltaT);
     Sub->A2 = 1.0f/(Beta*DeltaT);
     Sub->A3 = 1.0f/(2.0f*Beta) - 1.0f;
     Sub->A4 = Gamma/Beta - 1.0f;
     Sub->A5 = (DeltaT/2.0f)*(Gamma/Beta - 2.0f);
     Sub->A6 = (1.0f - Gamma)*DeltaT;
     Sub->A10 = (0.5f - Beta)*DeltaT*DeltaT;
#else
     Sub->A0 = 1.0/(Beta*DeltaT*DeltaT);
     Sub->A2 = 1.0/(Beta*DeltaT);
     Sub->A3 = 1.0/(2.0*Beta) - 1.0;
     Sub->A4 = Gamma/Beta - 1.0;
     Sub->A5 = (DeltaT/2.0)*(Gamma/Beta - 2.0);
     Sub->A6 = (1.0 - Gamma)*DeltaT;
     Sub->A10 = (0.5 - Beta)*DeltaT*DeltaT;
#endif

     Sub->A1 = Gamma/(Beta*DeltaT);
     Sub->A7 = Gamma*DeltaT;
     Sub->A8 = Beta*DeltaT*DeltaT;
     Sub->A9 = DeltaT;

     /* Init the variables */
     Sub->Mass = Mass;
     Sub->Damp = Damp;
     Sub->Stiff = Stiff;

     Sub->G = 1.0/(Sub->Stiff + Sub->a0*Sub->Mass + Sub->a1*Sub->Damp);

     Sub->dT = 0.0;
     Sub->vT = 0.0;
     Sub->aT = 0.0;

     Sub->dTdT = 0.0;
     Sub->vTdT = 0.0;
     Sub->aTdT = 0.0;

     Sub->l = 0.0;

     Sub->DispT = 0.0;
     Sub->VelTdT = 0.0;

     printf("A0 %lE A6 %lE a0 %lE a6 %lE\n", Sub->A0, Sub->A6, Sub->a0, Sub->a6 );
}

void Substructure_Newmark_SDOF( const HYSL_FLOAT DispTdT, const HYSL_FLOAT ramp, const HYSL_FLOAT GAcc, NewmarkSim_t *const Sub, HYSL_FLOAT *const fc )
{

     /* Compute initial velocity */
#if _FLOAT_
     Sub->AccTdT = (1.0f - ramp)*(-Sub->A0*Sub->Disp0 -Sub->A2*Sub->Vel0 -Sub->A3*Sub->Acc0)
	  + ramp*(-Sub->A0*Sub->DispT -Sub->A2*Sub->VelT -Sub->A3*Sub->AccT) + Sub->A0*DispTdT;

     Sub->VelTdT = (1.0f - ramp)*(Sub->Vel0 + Sub->A6*Sub->Acc0) + ramp*(Sub->VelT + Sub->A6*Sub->AccT) + Sub->A7*Sub->AccTdT;
#else
     Sub->AccTdT = (1.0 - ramp)*(-Sub->A0*Sub->Disp0 -Sub->A2*Sub->Vel0 -Sub->A3*Sub->Acc0)
	  + ramp*(-Sub->A0*Sub->DispT -Sub->A2*Sub->VelT -Sub->A3*Sub->AccT) + Sub->A0*DispTdT;

     Sub->VelTdT = (1.0 - ramp)*(Sub->Vel0 + Sub->A6*Sub->Acc0) + ramp*(Sub->VelT + Sub->A6*Sub->AccT) + Sub->A7*Sub->AccTdT;
#endif

     /* Load */
     Sub->l = -Sub->Mass*(GAcc + Sub->AccTdT);

     /* New Displacements */
     Sub->dTdT = Sub->G*(Sub->l + Sub->Mass*(Sub->a0*Sub->dT + Sub->a2*Sub->vT + Sub->a3*Sub->aT) + Sub->Damp*(Sub->a1*Sub->dT + Sub->a4*Sub->vT + Sub->a5*Sub->aT));

     /* Accelerations */
     Sub->aTdT = Sub->a0*(Sub->dTdT - Sub->dT) - Sub->a2*Sub->vT - Sub->a3*Sub->aT;

     /* Velocity */
     Sub->vTdT = Sub->vT + Sub->a6*Sub->aT + Sub->a7*Sub->aTdT;

     /* Coupling force */
     *fc = Sub->Stiff*Sub->dTdT + Sub->Damp*Sub->vTdT;

     /* Backup */
     Sub->dT = Sub->dTdT;
     Sub->vT = Sub->vTdT;
     Sub->aT = Sub->aTdT;

}


void Substructure_Newmark_Destroy( NewmarkSim_t *const Sub )
{
     free( Sub->Description );
}
