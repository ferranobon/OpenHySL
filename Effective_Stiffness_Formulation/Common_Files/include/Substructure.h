#ifndef SUBSTRUCTURE_H
#define SUBSTRUCTURE_H

#define USE_ADWIN 0    /* Run using ADwin */
#define USE_EXACT 1    /* Simulate the substructure using the exact solution */
#define USE_UHYDE 2    /* Simulate the substructure using the UHYDE-fbr device */
#define USE_MEASURED 3    /* Simulate the substructure using measured values */

#define UHYDE_NUMPARAM_INIT 3
#define EXACT_NUMPARAM_INIT 3

typedef struct {
     void *SimStruct;
     int Type;
     float (*fInit)( int, float*, ... );
     void (*fCalc)( const float, const float, void *const, float *const );
} Substructure;

/**
 * \brief Structure to store the coupling nodes.
 * 
 * This structure is used in order to store the coupling nodes that will be used
 * during a test. The nodes are stored sequentially and in increasing order.
 */
typedef struct {
     int *Array;  /*!< \brief Array containing the coupling nodes */
     int Order;   /*!< \brief Number of coupling nodes */
     float *u0c0;
     Substructure *Sub;
} Coupling_Node;

typedef struct{
     unsigned int Order_Couple;   /* Order of the coupling nodes. */
     unsigned int Num_Sub;        /* Number of sub-steps. */
     unsigned int Num_Steps;      /* Number of steps. */
     float DeltaT;       /* Time step. */
     float DeltaT_Sub;   /* Time between substeps. */
} ConstSub;

typedef struct{
     float Mass, Damp, Stiff;
     float Disp0, Vel0;
     float Disp, Vel;
     float Force_0, Force_1;
     float u0c_old;      /* Backup of the displacement at the coupling node */
     float v0c;          /* Velocity at the coupling node */

     float A,B,C,D,E,F,G,H;  /* Several constants */
} TMD_Sim;

typedef struct{
     float u0c_old;          /* Old displacement (the one from the previous sub-step */
     float q;                /* Displacement in the device */
     float qyield;           /* Yield displacement */
     float qplastic;         /* Plastic displacement */
     float k;                /* Initial stiffness */
} UHYDE_Sim;

void Init_Constants_Substructure( ConstSub *const Constants, const char *Filename );

void Simulate_Substructure_Measured_Values( const char *FileName, const float *const Keinv, const float *const u0c, float *const uc, float *const fcprev, float *const fc, const unsigned int OrderC, const unsigned int NSub );

void Simulate_Substructures( Coupling_Node *const CNodes, float *Keinv, float *const u0c, float *const uc, float *const fcprev, float *const fc, const unsigned int NSubstep, const float DeltaT_Sub );

void Simulate_Substructure( void *const Num, const int Mode, float *const Keinv, float *const u0c0, float *const u0c, float *const uc, float *const fcprev, float *const fc, const unsigned int OrderC, const unsigned int NSub, const float Deltat_Sub );

void ExactSolution_Init( const float Mass, const float Damp, const float Stiff, const float DeltaT, TMD_Sim *const Num );
void ExactSolution_SDOF( const float u0c, const float DeltaT, TMD_Sim *const Num, float *const fc );

void Simulate_UHYDE_1D_Init( const float qyield, const float yield_factor, const float Friction, UHYDE_Sim *const Num );

void Simulate_UHYDE_1D( const float u0c, const float DeltaT, UHYDE_Sim *const Num, float *const Friction_Force );

#endif /* SUBSTRUCTURE_H */
