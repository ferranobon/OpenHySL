#include "Substructure.h"

void Init_Constants_Substructure( ConstSub *const Constants )
{
     (*Constants).Order_Couple = 1;
     (*Constants).Num_Steps = 4096;
     (*Constants).Num_Sub = 4;
     (*Constants).DeltaT = 0.01;
     (*Constants).DeltaT_Sub = (*Constants).DeltaT/(float)(*Constants).Num_Sub;
}
