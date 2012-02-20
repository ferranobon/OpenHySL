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
