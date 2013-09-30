#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "Print_Messages.h"
#include "Substructure_Exact.h"
#include "MatrixVector.h"
#include "Auxiliary_Math.h"

#if _MKL_
#include "mkl_blas.h"
#else
#include "Netlib.h"
#endif

void Substructure_ExactSolutionMDOF_Init( const double *const Mass, const double *const Stiff, const int NDOF,
					  const double RayM_Alpha, const double RayM_Beta,
					  const double RayS_Alpha, const double RayS_Beta,
					  const double a0, const double a2, const double a3,
					  const double a6, const double a7, const char *Description,
					  ExactSim_t *const Sub )
{
     int i;

     Sub->Description = strdup( Description );

     MatrixVector_Create( NDOF, NDOF, &Sub->Mass );
     MatrixVector_Create( NDOF, NDOF, &Sub->Stiff );
     MatrixVector_Create( NDOF, NDOF, &Sub->EVectors );

     MatrixVector_Create( NDOF, 1, &Sub->EValues );
     MatrixVector_Create( NDOF, 1, &Sub->Damping_Ratios );
     MatrixVector_Create( NDOF, 1, &Sub->Init_Disp );
     MatrixVector_Create( NDOF, 1, &Sub->Init_Vel );
     MatrixVector_Create( NDOF, 1, &Sub->Load );
     MatrixVector_Create( NDOF, 1, &Sub->End_Disp );
     MatrixVector_Create( NDOF, 1, &Sub->End_Vel );
     MatrixVector_Create( NDOF, 1, &Sub->End_Acc );

     /* Init the variables */
     for( i = 0; i < NDOF*NDOF; i++ ){
	  Sub->Mass.Array[i] = Mass[i];
	  Sub->Stiff.Array[i] = Stiff[i];
     }

     Compute_Eigenvalues_Eigenvectors( &Sub->Stiff, &Sub->Mass, &Sub->EValues, &Sub->EVectors );
     Compute_DampingRatios_Rayleigh( RayS_Alpha, RayS_Beta, Sub->EValues.Rows, Sub->EValues.Array,
				     Sub->Damping_Ratios.Array );

     Sub->Disp0 = 0.0; Sub->DispT = 0.0;
     Sub->Vel0 = 0.0; Sub->VelT = 0.0; Sub->VelTdT = 0.0;
     Sub->Acc0 = 0.0; Sub->AccT = 0.0; Sub->AccTdT = 0.0;

     Sub->a0 = a0; Sub->a2 = a2; Sub->a3 = a3;
     Sub->a6 = a6; Sub->a7 = a7;

     Sub->Ray_Main.Alpha = RayM_Alpha; Sub->Ray_Main.Beta = RayM_Beta; 
     Sub->Ray_Sub.Alpha = RayS_Alpha; Sub->Ray_Sub.Beta = RayS_Beta;
}

void Substructure_ExactSolutionSDOF_Init( const double Mass, const double Damp, const double Stiff,
					  const double a0, const double a2, const double a3,
					  const double a6, const double a7, const char *Description,
					  ExactSim_t *const Sub )
{
     Sub->Description = strdup( Description );

     MatrixVector_Create( 1, 1, &Sub->Mass );
     MatrixVector_Create( 1, 1, &Sub->Stiff );
     MatrixVector_Create( 1, 1, &Sub->Damp );
     MatrixVector_Create( 1, 1, &Sub->EValues );
     MatrixVector_Create( 1, 1, &Sub->EVectors );

     MatrixVector_Create( 1, 1, &Sub->Damping_Ratios );
     MatrixVector_Create( 1, 1, &Sub->Init_Disp );
     MatrixVector_Create( 1, 1, &Sub->Init_Vel );
     MatrixVector_Create( 1, 1, &Sub->Load );
     MatrixVector_Create( 1, 1, &Sub->End_Disp );
     MatrixVector_Create( 1, 1, &Sub->End_Vel );
     MatrixVector_Create( 1, 1, &Sub->End_Acc );

     /* Init the variables */
     Sub->Mass.Array[0] = Mass;
     Sub->Stiff.Array[0] = Stiff;
     Sub->Damp.Array[0] = Damp;
     
     Compute_Eigenvalues_Eigenvectors( &Sub->Stiff, &Sub->Mass, &Sub->EValues, &Sub->EVectors );
     Sub->Damping_Ratios.Array[0] = Sub->Damp.Array[0]/(2.0*sqrt(Sub->EValues.Array[0])*Sub->Mass.Array[0]);

     Sub->Disp0 = 0.0; Sub->DispT = 0.0;
     Sub->Vel0 = 0.0; Sub->VelT = 0.0; Sub->VelTdT = 0.0;
     Sub->Acc0 = 0.0; Sub->AccT = 0.0; Sub->AccTdT = 0.0;

     Sub->a0 = a0; Sub->a2 = a2; Sub->a3 = a3;
     Sub->a6 = a6; Sub->a7 = a7;
}

void Compute_DampingRatios_Rayleigh( const double Ray_Alpha, const double Ray_Beta, const int Num_DOF,
				     const double *const Eigen_Values, double *const Damping_Ratios )
{
     int i; /* A counter */
     double omega;
     
#pragma omp parallel for
     for ( i = 0; i < Num_DOF; i++ ){
	  omega = sqrt( Eigen_Values[i] );
	  Damping_Ratios[i] = Ray_Alpha/(2.0*omega) + Ray_Beta*omega/2.0;
     }
}

void Substructure_ExactSolutionMDOF( const double DispTdT, const double ramp, const double GAcc, const double DeltaT,
				     ExactSim_t *const Sub, double *const fc )
{
     int i;
     int Length = Sub->Mass.Rows;

     /* Calculate the new acceleration and velocity for the coupling node */
     Sub->AccTdT = (1.0 - ramp)*(-Sub->a0*Sub->Disp0 -Sub->a2*Sub->Vel0 -Sub->a3*Sub->Acc0)
	  + ramp*(-Sub->a0*Sub->DispT -Sub->a2*Sub->VelT -Sub->a3*Sub->AccT) + Sub->a0*DispTdT;

     Sub->VelTdT = (1.0 - ramp)*(Sub->Vel0 + Sub->a6*Sub->Acc0) + ramp*(Sub->VelT + Sub->a6*Sub->AccT) + Sub->a7*Sub->AccTdT;

     /* Calculate the force acting on the substructure */
     for( i = 0; i < Length; i++ ){
	  Sub->Load.Array[i] = -(GAcc + Sub->AccTdT + Sub->Ray_Main.Alpha*Sub->VelTdT)*Sub->Mass.Array[i*Length + i];
     }

     /* Calculate the behaviour of the substructure */
     Duhamel_Integral( &Sub->Mass, &Sub->EValues, &Sub->EVectors, &Sub->Damping_Ratios, &Sub->Init_Disp,
		       &Sub->Init_Vel, &Sub->Load, &Sub->End_Disp, &Sub->End_Vel, &Sub->End_Acc, DeltaT );

     /* Calculate the new coupling force */
     if( Sub->Stiff.Rows > 1 ){
	  *fc = (Sub->Stiff.Array[0] - Sub->Stiff.Array[Sub->Stiff.Cols + 1])*(Sub->Ray_Sub.Beta*Sub->End_Vel.Array[0] + Sub->End_Disp.Array[0]);
     } else if (Sub->Stiff.Rows == 1) {
	  *fc = Sub->Stiff.Array[0]*(Sub->Ray_Sub.Beta*Sub->End_Vel.Array[0] + Sub->End_Disp.Array[0]);
     } else assert(0);

     /* Backup substructure values */
     for(i = 0; i < Length; i++){
	  Sub->Init_Disp.Array[i] = Sub->End_Disp.Array[i];
	  Sub->Init_Vel.Array[i] = Sub->End_Vel.Array[i];
     }
}

void Substructure_ExactSolutionSDOF( const double DispTdT, const double ramp, const double GAcc,
				      const double DeltaT, ExactSim_t *const Sub, double *const fc )
{

     /* Calculate the new acceleration and velocity for the coupling node */
     Sub->AccTdT = (1.0 - ramp)*(-Sub->a0*Sub->Disp0 -Sub->a2*Sub->Vel0 -Sub->a3*Sub->Acc0)
	  + ramp*(-Sub->a0*Sub->DispT -Sub->a2*Sub->VelT -Sub->a3*Sub->AccT) + Sub->a0*DispTdT;

     /* Calculate the force acting on the substructure */
     Sub->Load.Array[0] = -(GAcc + Sub->AccTdT)*Sub->Mass.Array[0];

     /* Calculate the behaviour of the substructure */
     Duhamel_Integral( &Sub->Mass, &Sub->EValues, &Sub->EVectors, &Sub->Damping_Ratios, &Sub->Init_Disp,
		       &Sub->Init_Vel, &Sub->Load, &Sub->End_Disp, &Sub->End_Vel, &Sub->End_Acc, DeltaT );

     /* Calculate the new coupling force */
     *fc = Sub->Stiff.Array[0]*Sub->End_Disp.Array[0] + Sub->Damp.Array[0]*Sub->End_Vel.Array[0];
//     *fc = Sub->Load.Array[0] - Sub->Mass.Array[0]*Sub->End_Acc.Array[0];
     //printf( "fc = %le, fc1 = %le\n", fc[0], fc1 );

     /* Backup substructure values */
     Sub->Init_Disp.Array[0] = Sub->End_Disp.Array[0];
     Sub->Init_Vel.Array[0] = Sub->End_Vel.Array[0];
}

void Duhamel_Integral( const MatrixVector_t *const Mass, const MatrixVector_t *const Eigen_Values,
		       const MatrixVector_t *const Eigen_Vectors, const MatrixVector_t *const Damping_Ratios,
		       const MatrixVector_t *const Init_Disp, const MatrixVector_t *const Init_Vel, 
		       const MatrixVector_t *const Load, MatrixVector_t *const End_Disp,
		       MatrixVector_t *const End_Vel, MatrixVector_t *const End_Acc, const double Delta_T )
{

     int i;                        /* Counter */
     int incx, incy;               /* Stride in the vectors */   
     char uplo = 'L';              /* The lower part (upper part in C) will be used and the upper part (lower
				    * part in C) will strictly not be referenced */
     char trans;
     double Alpha, Beta;           /* Constants for the BLAS routines */
     double om2, om, omp, zeta=0;
     double expo, sino, coso;
     double x_part, x_hom, a, b;    

     MatrixVector_t modal_force;
     MatrixVector_t modal_disp, modal_velo;
     MatrixVector_t temp1, temp2;

     MatrixVector_Create( Mass->Rows, 1, &modal_force );
     MatrixVector_Create( Mass->Rows, 1, &modal_disp );
     MatrixVector_Create( Mass->Rows, 1, &modal_velo );

     /* Transform into modal coordinates */
     MatrixVector_Create( Mass->Rows, 1, &temp1 );
     MatrixVector_Create( Mass->Rows, 1, &temp2 );

     Alpha = 1.0; Beta = 0.0;
     incx = 1; incy = 1;
     dsymv( &uplo, &temp1.Rows, &Alpha, Mass->Array, &temp1.Rows, Init_Disp->Array, &incx, &Beta, temp1.Array, &incy );
     dsymv( &uplo, &temp2.Rows, &Alpha, Mass->Array, &temp2.Rows, Init_Vel->Array, &incx, &Beta, temp2.Array, &incy );

     trans = 'T';
     dgemv( &trans, &temp1.Rows, &temp1.Rows, &Alpha, Eigen_Vectors->Array, &temp1.Rows, temp1.Array, &incx, &Beta,
	     modal_disp.Array, &incy );
     dgemv( &trans, &temp2.Rows, &temp1.Rows, &Alpha, Eigen_Vectors->Array, &temp2.Rows, temp2.Array, &incx, &Beta,
	     modal_velo.Array, &incy );
     dgemv( &trans, &temp1.Rows, &temp1.Rows, &Alpha, Eigen_Vectors->Array, &temp1.Rows, Load->Array, &incx, &Beta,
	     modal_force.Array, &incy );

     MatrixVector_Destroy( &temp1 );
     MatrixVector_Destroy( &temp2 );

     for (i = 0; i < Mass->Rows; i++) {

	  om2	= Eigen_Values->Array[i];

	  if (Damping_Ratios != NULL){
	       zeta	= Damping_Ratios->Array[i];
	  }

	  if (om2 > 0.0){
	       om = sqrt(om2);

	       if (zeta >= 0.0 && zeta < 1.0) {
		    omp	= om*sqrt(1.0-zeta*zeta);
		    x_part	= modal_force.Array[i]/om2;
		    a	= -x_part + modal_disp.Array[i];
		    b	= (zeta*om*a + modal_velo.Array[i])/omp;
		    coso	= cos(omp*Delta_T);
		    sino	= sin(omp*Delta_T);
		    expo	= exp(-zeta*om*Delta_T);
		    x_hom	= expo*(a*coso + b*sino);
		    modal_velo.Array[i]	= -zeta*om*x_hom + expo*omp*(-a*sino + b*coso);
		    modal_disp.Array[i]	= x_hom + x_part;
		    /* use modal_force for accelerations*/
		    modal_force.Array[i]	-= om2*modal_disp.Array[i];
		    modal_force.Array[i]	-= 2.0*om*zeta*modal_velo.Array[i];
	       } else {
		    Print_Header( ERROR );
		    fprintf( stderr, "Duhamel Integration: negative or overcritical damping: %le.\n", zeta );
		    exit( EXIT_FAILURE );
	       }
	  } else {
	       Print_Header( ERROR );
	       fprintf( stderr, "Duhamel Integration: non-positive eigenvalue.\n");
	       exit( EXIT_FAILURE );
	  }
     }

     /* Transform into physical coordinates */
     trans = 'N';
     dgemv( &trans, &modal_disp.Rows, &modal_disp.Rows, &Alpha, Eigen_Vectors->Array, &modal_disp.Rows,
	     modal_disp.Array, &incx, &Beta, End_Disp->Array, &incy );
     dgemv( &trans, &modal_velo.Rows, &modal_velo.Rows, &Alpha, Eigen_Vectors->Array, &modal_velo.Rows,
	     modal_velo.Array, &incx, &Beta, End_Vel->Array, &incy );
     dgemv( &trans, &modal_force.Rows, &modal_force.Rows, &Alpha, Eigen_Vectors->Array, &modal_force.Rows,
	     modal_force.Array, &incx, &Beta, End_Acc->Array, &incy );

     /* Free dynamically allocated memory */
     MatrixVector_Destroy( &modal_force );
     MatrixVector_Destroy( &modal_disp );
     MatrixVector_Destroy( &modal_velo );
}


void Substructure_ExactSolutionMDOF_Destroy( ExactSim_t *const Sub )
{
     free( Sub->Description );

     MatrixVector_Destroy( &Sub->Mass );
     MatrixVector_Destroy( &Sub->Stiff );
     MatrixVector_Destroy( &Sub->EValues );
     MatrixVector_Destroy( &Sub->EVectors );
     MatrixVector_Destroy( &Sub->Damping_Ratios );
     MatrixVector_Destroy( &Sub->Load );
     MatrixVector_Destroy( &Sub->Init_Disp );
     MatrixVector_Destroy( &Sub->Init_Vel );
     MatrixVector_Destroy( &Sub->End_Disp );
     MatrixVector_Destroy( &Sub->End_Vel );
     MatrixVector_Destroy( &Sub->End_Acc );
}

void Substructure_ExactSolutionSDOF_Destroy( ExactSim_t *const Sub )
{
     free( Sub->Description );

     MatrixVector_Destroy( &Sub->Mass );
     MatrixVector_Destroy( &Sub->Stiff );
     MatrixVector_Destroy( &Sub->Damp );
     MatrixVector_Destroy( &Sub->EValues );
     MatrixVector_Destroy( &Sub->EVectors );
     MatrixVector_Destroy( &Sub->Damping_Ratios );
     MatrixVector_Destroy( &Sub->Load );
     MatrixVector_Destroy( &Sub->Init_Disp );
     MatrixVector_Destroy( &Sub->Init_Vel );
     MatrixVector_Destroy( &Sub->End_Disp );
     MatrixVector_Destroy( &Sub->End_Vel );
     MatrixVector_Destroy( &Sub->End_Acc );
}

void Substructure_ExactSolutionESP_Init( const double Mass, const double Damp, const double Stiff, const double DeltaT, const char *Description, ExactSimESP_t *const Sub )
{
     double omega, omegaD, omega2;
     double xi;
     double expdt, sinDdt, cosDdt, omegadt, omegaDdt;

     Sub->Description = strdup( Description );

     /* Init the variables */
     Sub->Mass = Mass;
     Sub->Damp = Damp;
     Sub->Stiff = Stiff;

     Sub->Disp0 = 0.0;
     Sub->Vel0 = 0.0;

     Sub->Disp = 0.0;
     Sub->Vel = 0.0;

     Sub->Force_0 = 0.0;
     Sub->Force_1 = 0.0;

     Sub->u0c_old = 0.0;
     Sub->v0c = 0.0;

     /* Calculate the free vibration frequency */
     omega2 = Stiff/Mass;

     /* Check the value of omega */
     if ( omega2 > 0.0 ){
	  omega = sqrt( omega2 );

	  /* Calculate the damping ratio */
	  xi = Damp/(2.0*Mass*omega);
	  
	  /* Check the value of xi */
	  if ( xi >= 0.0 && xi < 1.0 ){
	       /* Calculate the damping vibration frequency */
	       omegaD = omega*sqrt(1.0-xi*xi);

	       /* Calculate several constants */
	       expdt = exp(-xi*omega*DeltaT);
	       sinDdt = sin(omegaD*DeltaT);
	       cosDdt = cos(omegaD*DeltaT);
	       omegadt = omega*DeltaT;
	       omegaDdt = omegaD*DeltaT;

	       Sub->A = expdt*(xi*omega*sinDdt/omega + cosDdt);
	       Sub->B = expdt*sinDdt/omegaD;
	       Sub->C = (expdt*(((1.0-2.0*xi*xi)/omegaDdt-xi*omega/omegaD)*sinDdt-(1.0+2.0*xi/omegadt)*cosDdt)+2.0*xi/omegadt)/Stiff;
	       Sub->D = (expdt*((2.0*xi*xi-1.0)*sinDdt/omegaDdt + 2.0*xi*cosDdt/omegadt)+(1.0-2.0*xi/omegadt))/Stiff;
	       Sub->E = -expdt * omega * omega * sinDdt / omegaD;
	       Sub->F = expdt * (cosDdt - xi*omega*sinDdt/omegaD);
	       Sub->G = (expdt*((omega*omega/omegaD + xi*omega/omegaDdt)*sinDdt + cosDdt/DeltaT)-1.0/DeltaT)/Stiff;
	       Sub->H = (-expdt*(xi*omega*sinDdt/omegaD + cosDdt) + 1.0)/(Stiff*DeltaT);

	   } else {
	       Print_Header( ERROR );
	       fprintf( stderr, "Exact Solution: Negative or overcritical damping.\n" );
	       exit( EXIT_FAILURE );
	  }
     } else {
	  Print_Header( ERROR );
	  fprintf( stderr, "Exact Solution: Negative vibration frequency.\n" );
	  exit( EXIT_FAILURE );
     }
}

void Substructure_ExactSolutionESP_SDOF( const double u0c, const double DeltaT, ExactSimESP_t *const Sub, double *const fc )
{

     /* Compute initial velocity */
     Sub->v0c = (u0c - Sub->u0c_old)/DeltaT;

     /* Backup and computer new forcce */
     Sub->Force_0 = Sub->Force_1;
     Sub->Force_1 = Sub->Stiff*u0c + Sub->Damp*Sub->v0c;

     /* Backup the displacements */
     Sub->Disp0 = Sub->Disp;
     Sub->Vel0 = Sub->Vel;

     /* Compute the new state */
     Sub->Disp = Sub->A*Sub->Disp0 + Sub->B*Sub->Vel0 + Sub->C*Sub->Force_0 + Sub->D*Sub->Force_1;
     Sub->Vel = Sub->E*Sub->Disp0 + Sub->F*Sub->Vel0 + Sub->G*Sub->Force_0 + Sub->H*Sub->Force_1;

     /* Compute the coupling force */
     *fc = Sub->Stiff*(Sub->Disp - u0c) + Sub->Damp*(Sub->Vel - Sub->v0c);

     /* Backup u0c vector */
     Sub->u0c_old = u0c;
	 
}


void Substructure_ExactSolutionESP_Destroy( ExactSimESP_t *const Sub )
{
     free( Sub->Description );
}
