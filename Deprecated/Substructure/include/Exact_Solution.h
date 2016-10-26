#ifndef EXACT_SOLUTION_H
#define EXACT_SOLUTION_H

typedef struct{
     float *Mass, *Damp, *Stiff;
     float Num_DOF;
     float Disp0, Vel0;
     float Disp, Vel;
     float Force_0, Force_1;
     float u0c_old;      /* Backup of the displacement at the coupling node */
     float v0c;          /* Velocity at the coupling node */

     float A,B,C,D,E,F,G,H;  /* Several constants */
} TMD_Sim;


void Compute_EigenValues_EigenVectors ( const int Num_DOF, float *const MatrixA, float *const MatrixB, float *const Eigen_Values, float *const Eigen_Vectors );
void Compute_Damp_Matrix( const float a0, const float a1, const int Num_DOF, const float *const Mass, const float *const Stiff, float *const Damp );
void Compute_Damping_Ratios_and_Matrix ( const float a0, const float a1, const int Num_DOF, const float *const Eigen_Values, float *const Damping_Ratios );
void ExactSolution_Init( const int Num_DOF, const float a0, const float a1, const char* MassFile, const char *StiffFile, TMD_Sim *const Num );
void Exact_Solution ( const float *const u0c, const float DeltaT, const int Num_DOF, float *const fc );
