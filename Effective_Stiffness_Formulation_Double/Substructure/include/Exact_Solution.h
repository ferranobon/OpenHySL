#ifndef EXACT_SOLUTION_H
#define EXACT_SOLUTION_H

typedef struct{
     double *Mass, *Damp, *Stiff;
     double Num_DOF;
     double Disp0, Vel0;
     double Disp, Vel;
     double Force_0, Force_1;
     double u0c_old;      /* Backup of the displacement at the coupling node */
     double v0c;          /* Velocity at the coupling node */

     double A,B,C,D,E,F,G,H;  /* Several constants */
} TMD_Sim;


void Compute_EigenValues_EigenVectors ( const int Num_DOF, double *const MatrixA, double *const MatrixB, double *const Eigen_Values, double *const Eigen_Vectors );
void Compute_Damp_Matrix( const double a0, const double a1, const int Num_DOF, const double *const Mass, const double *const Stiff, double *const Damp );
void Compute_Damping_Ratios_and_Matrix ( const double a0, const double a1, const int Num_DOF, const double *const Eigen_Values, double *const Damping_Ratios );
void ExactSolution_Init( const int Num_DOF, const double a0, const double a1, const char* MassFile, const char *StiffFile, TMD_Sim *const Num );
void Exact_Solution ( const double *const u0c, const double DeltaT, const int Num_DOF, double *const fc );
