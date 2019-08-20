#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include "Substructure_BoucWen.h"
#include "Substructure_StoneDrums.h"

#include "Auxiliary_Math.h"
#include "Print_Messages.h"
#include "Definitions.h"

void Substructure_StoneDrums_Init (const int JointID, const int PrevJointID, const int NextJointID, const int Dir, const hysl_float_t *const BW, const int *const BoucWen_Type, const char *Description,
        StoneDrums_t *const Sub) {
    int Pos;

    Sub->Description = strdup(Description);

    Sub->JointID = JointID;
    Sub->PrevJointID = PrevJointID;
    Sub->NextJointID = NextJointID;
    Sub->Dir = Dir;

    Pos = 0;
    Substructure_BoucWenSurface_Init(BW[Pos], BW[Pos + 1], BW[Pos + 2], BW[Pos + 3], BW[Pos + 4], BW[Pos + 5], BW[Pos + 6], BW[Pos + 7], BW[Pos + 8], BW[Pos + 9], BW[Pos + 10], BW[Pos + 11],
            BoucWen_Type[0], Sub->Description, &Sub->BoucWen_Sliding);
    Pos = Pos + BOUCWENDEG_NUMPARAM_INIT;
    Substructure_BoucWen_Init(BW[Pos], BW[Pos + 1], BW[Pos + 2], BW[Pos + 3], BW[Pos + 4], BW[Pos + 5], BW[Pos + 6], BW[Pos + 7], BW[Pos + 8], BW[Pos + 9], BW[Pos + 10], BW[Pos + 11], 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, BoucWen_Type[1], Sub->Description, &Sub->BoucWen_Torsion);
    Pos = Pos + BOUCWENDEG_NUMPARAM_INIT;
    Substructure_BoucWenSurface_Init(BW[Pos], BW[Pos + 1], BW[Pos + 2], BW[Pos + 3], BW[Pos + 4], BW[Pos + 5], BW[Pos + 6], BW[Pos + 7], BW[Pos + 8], BW[Pos + 9], BW[Pos + 10], BW[Pos + 11],
            BoucWen_Type[2], Sub->Description, &Sub->BoucWen_Bending);

}

void Substructure_StoneDrums (const hysl_float_t *const DispTdT, const int Order, StoneDrums_t *const Sub, hysl_float_t *force) {
    hysl_float_t alpha;
    hysl_float_t Rel_DispXY, Rel_Torsion, Rel_Bending;
    hysl_float_t Force_XY, Moment_Torsion, Moment_Bending;

    alpha = atan((DispTdT[0] - DispTdT[1]) / (DispTdT[0] - DispTdT[1]));

    Rel_DispXY = sqrt(pow(DispTdT[0] - DispTdT[1], 2.0) + pow(DispTdT[0] - DispTdT[1], 2.0));
    Rel_Torsion = DispTdT[0];
    Rel_Bending = DispTdT[1];
    printf("DO NOT EXECUTE\n");
//     Substructure_BoucWen ( Rel_DispXY, &Sub->BoucWen_Sliding, &Force_XY );
    //   Substructure_BoucWen ( Rel_Torsion, &Sub->BoucWen_Torsion, &Moment_Torsion );
    // Substructure_BoucWen ( Rel_Bending, &Sub->BoucWen_Bending, &Moment_Bending );

    force[0] = Force_XY * cos(alpha);
    force[1] = Force_XY * sin(alpha);
    force[2] = Moment_Torsion;
    force[3] = Moment_Bending;

}

void Substructure_StoneDrums_Destroy (StoneDrums_t *const Sub) {
    free(Sub->Description);
}
