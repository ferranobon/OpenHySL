#ifndef SUBSTRUCTURE_H
#define SUBSTRUCTURE_H

#define NUM_CHANNELS 22  /* Number of channels to save from ADwin */

typedef struct{
     int Order_Couple;   /* Order of the coupling nodes. */
     int Num_Sub;        /* Number of sub-steps. */
     int Num_Steps;      /* Number of steps. */
     float DeltaT;       /* Time step. */
     float DeltaT_Sub;   /* Time between substeps. */
} ConstSub;

void Init_Constants_Substructure( ConstSub *const Constants );

#endif /* SUBSTRUCTURE_H */
