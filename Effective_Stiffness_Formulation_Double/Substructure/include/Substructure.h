#ifndef SUBSTRUCTURE_H
#define SUBSTRUCTURE_H

#define USE_ADWIN    0    /* Run using ADwin */
#define USE_EXACT    1    /* Simulate the substructure using the exact solution */
#define USE_UHYDE    2    /* Simulate the substructure using the exact solution */
#define USE_MEASURED 3    /* Simulate the substructure using measured values */

typedef struct{
     unsigned int Order_Couple;   /* Order of the coupling nodes. */
     unsigned int Num_Sub;        /* Number of sub-steps. */
     unsigned int Num_Steps;      /* Number of steps. */
     double DeltaT;       /* Time step. */
     double DeltaT_Sub;   /* Time between substeps. */
} ConstSub;

typedef struct{
     double Mass, Damp, Stiff;
     double Disp0, Vel0;
     double Disp, Vel;
     double Force_0, Force_1;
     double u0c_old;      /* Backup of the displacement at the coupling node */
     double v0c;          /* Velocity at the coupling node */

     double A,B,C,D,E,F,G,H;  /* Several constants */
} TMD_Sim;

typedef struct{
     double u0c_old;          /* Old displacement (the one from the previous sub-step */
     double q;                /* Displacement in the device */
     double qyield;           /* Yield displacement */
     double qplastic;         /* Plastic displacement */
     double k;                /* Initial stiffness */
} UHYDE_Sim;

void Init_Constants_Substructure( ConstSub *const Constants, const char *Filename );

void Simulate_Substructure_Measured_Values( const char *FileName, const double *const Keinv, const double *const u0c, double *const uc, double *const fcprev, double *const fc, const unsigned int OrderC, const unsigned int NSub );

void Simulate_Substructure( void *const Num, const int Mode, double *const Keinv, double *const u0c0, double *const u0c, double *const uc, double *const fcprev, double *const fc, const unsigned int OrderC, const unsigned int NSub, const double Deltat_Sub );

void ExactSolution_Init( const double Mass, const double Damp, const double Stiff, const double DeltaT, TMD_Sim *const Num );
void ExactSolution_SDOF( const double u0c, const double DeltaT, TMD_Sim *const Num, double *const fc );

void Simulate_UHYDE_1D_Init( const double qyield, const double yield_factor, const double Friction, UHYDE_Sim *const Num );

void Simulate_UHYDE_1D( const double u0c, const double DeltaT, UHYDE_Sim *const Num, double *const Friction_Force );

#endif /* SUBSTRUCTURE_H */
