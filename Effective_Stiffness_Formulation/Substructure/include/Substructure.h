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

typedef struct{
     float Disp0, Vel0;
     float Disp, Vel;
     float Force_0, Force_1;
     float u0c_old;      /* Backup of the displacement at the coupling node */
     float v0c;          /* Velocity at the coupling node */

     float A,B,C,D,E,F,G,H;  /* Several constants */
} TMD_Sim;


void Init_Constants_Substructure( ConstSub *const Constants );


void ExactSolution_Init( const float Mass, const float Damp, const float Stiff, const float DeltaT, TMD_Sim *const Num );
void ExactSolution_SDOF( const float Mass, const float Damp, const float Stif, const float u0c, const float DeltaT, TMD_Sim *const Num, float *const fc );

#endif /* SUBSTRUCTURE_H */
