#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "Print_Messages.h"
#include "Substructure_Exact.h"
#include "MatrixVector.h"
#include "Auxiliary_Math.h"
#include "Definitions.h"

#if _MKL_
#include "mkl_blas.h"
#else
#include "Netlib.h"
#endif

void Substructure_ExactSolutionMDOF_Init( const HYSL_FLOAT *const Mass, const HYSL_FLOAT *const Stiff, const int NDOF,
					  const HYSL_FLOAT RayM_Alpha, const HYSL_FLOAT RayM_Beta,
					  const HYSL_FLOAT RayS_Alpha, const HYSL_FLOAT RayS_Beta,
					  const HYSL_FLOAT a0, const HYSL_FLOAT a2, const HYSL_FLOAT a3,
					  const HYSL_FLOAT a6, const HYSL_FLOAT a7, const char *Description,
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

void Substructure_ExactSolutionSDOF_Init( const HYSL_FLOAT Mass, const HYSL_FLOAT Damp, const HYSL_FLOAT Stiff,
					  const HYSL_FLOAT a0, const HYSL_FLOAT a2, const HYSL_FLOAT a3,
					  const HYSL_FLOAT a6, const HYSL_FLOAT a7, const char *Description,
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
#if _FLOAT_
     Sub->Damping_Ratios.Array[0] = Sub->Damp.Array[0]/(2.0f*sqrtf(Sub->EValues.Array[0])*Sub->Mass.Array[0]);
#else
     Sub->Damping_Ratios.Array[0] = Sub->Damp.Array[0]/(2.0*sqrt(Sub->EValues.Array[0])*Sub->Mass.Array[0]);
#endif

     Sub->Disp0 = 0.0; Sub->DispT = 0.0;
     Sub->Vel0 = 0.0; Sub->VelT = 0.0; Sub->VelTdT = 0.0;
     Sub->Acc0 = 0.0; Sub->AccT = 0.0; Sub->AccTdT = 0.0;

     Sub->a0 = a0; Sub->a2 = a2; Sub->a3 = a3;
     Sub->a6 = a6; Sub->a7 = a7;

     printf("a0 %le a2 %le a3 %le a6%le a7 %le\n", a0, a2, a3, a6, a7 );
}

void Compute_DampingRatios_Rayleigh( const HYSL_FLOAT Ray_Alpha, const HYSL_FLOAT Ray_Beta, const int Num_DOF,
				     const HYSL_FLOAT *const Eigen_Values, HYSL_FLOAT *const Damping_Ratios )
{
     int i; /* A counter */
     HYSL_FLOAT omega;
     
     for ( i = 0; i < Num_DOF; i++ ){
#if _FLOAT_
	  omega = sqrtf( Eigen_Values[i] );
	  Damping_Ratios[i] = Ray_Alpha/(2.0f*omega) + Ray_Beta*omega/2.0f;
#else
	  omega = sqrt( Eigen_Values[i] );
	  Damping_Ratios[i] = Ray_Alpha/(2.0*omega) + Ray_Beta*omega/2.0;
#endif
     }
}

void Substructure_ExactSolutionMDOF( const HYSL_FLOAT DispTdT, const HYSL_FLOAT ramp, const HYSL_FLOAT GAcc, const HYSL_FLOAT DeltaT,
				     ExactSim_t *const Sub, HYSL_FLOAT *const fc )
{
     int i;
     int Length = Sub->Mass.Rows;

     /* Calculate the new acceleration and velocity for the coupling node */
#if _FLOAT_
     Sub->AccTdT = (1.0f - ramp)*(-Sub->a0*Sub->Disp0 -Sub->a2*Sub->Vel0 -Sub->a3*Sub->Acc0)
	  + ramp*(-Sub->a0*Sub->DispT -Sub->a2*Sub->VelT -Sub->a3*Sub->AccT) + Sub->a0*DispTdT;

     Sub->VelTdT = (1.0f - ramp)*(Sub->Vel0 + Sub->a6*Sub->Acc0) + ramp*(Sub->VelT + Sub->a6*Sub->AccT) + Sub->a7*Sub->AccTdT;
#else
     Sub->AccTdT = (1.0 - ramp)*(-Sub->a0*Sub->Disp0 -Sub->a2*Sub->Vel0 -Sub->a3*Sub->Acc0)
	  + ramp*(-Sub->a0*Sub->DispT -Sub->a2*Sub->VelT -Sub->a3*Sub->AccT) + Sub->a0*DispTdT;

     Sub->VelTdT = (1.0 - ramp)*(Sub->Vel0 + Sub->a6*Sub->Acc0) + ramp*(Sub->VelT + Sub->a6*Sub->AccT) + Sub->a7*Sub->AccTdT;
#endif

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

void Substructure_ExactSolutionSDOF( const HYSL_FLOAT DispTdT, const HYSL_FLOAT ramp, const HYSL_FLOAT GAcc,
				      const HYSL_FLOAT DeltaT, ExactSim_t *const Sub, HYSL_FLOAT *const fc )
{

     /* Calculate the new acceleration and velocity for the coupling node */
#if _FLOAT_
     Sub->AccTdT = (1.0f - ramp)*(-Sub->a0*Sub->Disp0 -Sub->a2*Sub->Vel0 -Sub->a3*Sub->Acc0)
	  + ramp*(-Sub->a0*Sub->DispT -Sub->a2*Sub->VelT -Sub->a3*Sub->AccT) + Sub->a0*DispTdT;
#else
     Sub->AccTdT = (1.0 - ramp)*(-Sub->a0*Sub->Disp0 -Sub->a2*Sub->Vel0 -Sub->a3*Sub->Acc0)
	  + ramp*(-Sub->a0*Sub->DispT -Sub->a2*Sub->VelT -Sub->a3*Sub->AccT) + Sub->a0*DispTdT;
#endif

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
		       MatrixVector_t *const End_Vel, MatrixVector_t *const End_Acc, const HYSL_FLOAT Delta_T )
{

     int i;                        /* Counter */
     int incx, incy;               /* Stride in the vectors */   
     char uplo = 'L';              /* The lower part (upper part in C) will be used and the upper part (lower
				    * part in C) will strictly not be referenced */
     char trans;
     HYSL_FLOAT Alpha, Beta;           /* Constants for the BLAS routines */
     HYSL_FLOAT om2, om, omp, zeta=0.0;
     HYSL_FLOAT expo, sino, coso;
     HYSL_FLOAT x_part, x_hom, a, b;    

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
     hysl_symv( &uplo, &temp1.Rows, &Alpha, Mass->Array, &temp1.Rows, Init_Disp->Array, &incx, &Beta, temp1.Array, &incy );
     hysl_symv( &uplo, &temp2.Rows, &Alpha, Mass->Array, &temp2.Rows, Init_Vel->Array, &incx, &Beta, temp2.Array, &incy );

     trans = 'T';
     hysl_gemv( &trans, &temp1.Rows, &temp1.Rows, &Alpha, Eigen_Vectors->Array, &temp1.Rows, temp1.Array, &incx, &Beta,
	     modal_disp.Array, &incy );
     hysl_gemv( &trans, &temp2.Rows, &temp1.Rows, &Alpha, Eigen_Vectors->Array, &temp2.Rows, temp2.Array, &incx, &Beta,
	     modal_velo.Array, &incy );
     hysl_gemv( &trans, &temp1.Rows, &temp1.Rows, &Alpha, Eigen_Vectors->Array, &temp1.Rows, Load->Array, &incx, &Beta,
	     modal_force.Array, &incy );

     MatrixVector_Destroy( &temp1 );
     MatrixVector_Destroy( &temp2 );

     for (i = 0; i < Mass->Rows; i++) {

	  om2	= Eigen_Values->Array[i];

	  if (Damping_Ratios != NULL){
	       zeta	= Damping_Ratios->Array[i];
	  }

	  if (om2 > 0.0){
#if _FLOAT_
	       om = sqrtf(om2);
#else
	       om = sqrt(om2);
#endif

	       if (zeta >= 0.0 && zeta < 1.0) {
#if _FLOAT_
		    omp	= om*sqrtf(1.0f-zeta*zeta);
#else
		    omp	= om*sqrt(1.0-zeta*zeta);
#endif
		    x_part	= modal_force.Array[i]/om2;
		    a	= -x_part + modal_disp.Array[i];
		    b	= (zeta*om*a + modal_velo.Array[i])/omp;
#if _FLOAT_
		    coso	= cosf(omp*Delta_T);
		    sino	= sinf(omp*Delta_T);
		    expo	= expf(-zeta*om*Delta_T);
#else
		    coso	= cos(omp*Delta_T);
		    sino	= sin(omp*Delta_T);
		    expo	= exp(-zeta*om*Delta_T);
#endif
		    x_hom	= expo*(a*coso + b*sino);
		    modal_velo.Array[i]	= -zeta*om*x_hom + expo*omp*(-a*sino + b*coso);
		    modal_disp.Array[i]	= x_hom + x_part;
		    /* use modal_force for accelerations*/
		    modal_force.Array[i]	-= om2*modal_disp.Array[i];
#if _FLOAT_
		    modal_force.Array[i]	-= 2.0f*om*zeta*modal_velo.Array[i];
#else
		    modal_force.Array[i]	-= 2.0*om*zeta*modal_velo.Array[i];
#endif
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
     hysl_gemv( &trans, &modal_disp.Rows, &modal_disp.Rows, &Alpha, Eigen_Vectors->Array, &modal_disp.Rows,
	     modal_disp.Array, &incx, &Beta, End_Disp->Array, &incy );
     hysl_gemv( &trans, &modal_velo.Rows, &modal_velo.Rows, &Alpha, Eigen_Vectors->Array, &modal_velo.Rows,
	     modal_velo.Array, &incx, &Beta, End_Vel->Array, &incy );
     hysl_gemv( &trans, &modal_force.Rows, &modal_force.Rows, &Alpha, Eigen_Vectors->Array, &modal_force.Rows,
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

void Substructure_ExactSolutionESP_Init( const HYSL_FLOAT Mass, const HYSL_FLOAT Damp, const HYSL_FLOAT Stiff, const HYSL_FLOAT DeltaT, const char *Description, ExactSimESP_t *const Sub )
{
     HYSL_FLOAT omega, omegaD, omega2;
     HYSL_FLOAT xi;
     HYSL_FLOAT expdt, sinDdt, cosDdt, omegadt, omegaDdt;

     Sub->Description = strdup( Description );

#if _FLOAT_
     Sub->a0 = 4E4f;
     Sub->a2 = 4E2f;
     Sub->a3 = 1.0f;
     Sub->a6 = 5E-3f;
     Sub->a7 = 5E-3f;
#else
     Sub->a0 = 4E4;
     Sub->a2 = 4E2;
     Sub->a3 = 1.0;
     Sub->a6 = 5E-3;
     Sub->a7 = 5E-3;
#endif

     /* Init the variables */
     Sub->Mass = Mass;
     Sub->Damp = Damp;
     Sub->Stiff = Stiff;

     Sub->Init_Disp = 0.0;
     Sub->Init_Vel = 0.0;

     Sub->End_Disp = 0.0;
     Sub->End_Vel = 0.0;

     Sub->Force_0 = 0.0;
     Sub->Force_1 = 0.0;

     Sub->DispT = 0.0;
     Sub->VelTdT = 0.0;

     /* Calculate the free vibration frequency */
     omega2 = Stiff/Mass;

     /* Check the value of omega */
     if ( omega2 > 0.0 ){
#if _FLOAT_
	  omega = sqrtf( omega2 );
#else
	  omega = sqrt( omega2 );
#endif

	  /* Calculate the damping ratio */
#if _FLOAT_
	  xi = Damp/(2.0f*Mass*omega);
#else
	  xi = Damp/(2.0*Mass*omega);
#endif
	  
	  /* Check the value of xi */
	  if ( xi >= 0.0 && xi < 1.0 ){
	       /* Calculate the damping vibration frequency */
#if _FLOAT_
	       omegaD = omega*sqrtf(1.0f-xi*xi);
#else
	       omegaD = omega*sqrt(1.0-xi*xi);
#endif

	       /* Calculate several constants */
#if _FLOAT_
	       expdt = expf(-xi*omega*DeltaT);
	       sinDdt = sinf(omegaD*DeltaT);
	       cosDdt = cosf(omegaD*DeltaT);
#else 
	       expdt = exp(-xi*omega*DeltaT);
	       sinDdt = sin(omegaD*DeltaT);
	       cosDdt = cos(omegaD*DeltaT);
#endif
	       omegadt = omega*DeltaT;
	       omegaDdt = omegaD*DeltaT;

	       Sub->A = expdt*(xi*omega*sinDdt/omegaD + cosDdt);
	       Sub->B = expdt*sinDdt/omegaD;
#if _FLOAT_
	       Sub->C = (expdt*(((1.0f-2.0f*xi*xi)/omegaDdt-xi*omega/omegaD)*sinDdt-(1.0f+2.0f*xi/omegadt)*cosDdt)+2.0f*xi/omegadt)/Stiff;
	       Sub->D = (expdt*((2.0f*xi*xi-1.0f)*sinDdt/omegaDdt + 2.0f*xi*cosDdt/omegadt)+(1.0f-2.0f*xi/omegadt))/Stiff;
#else
	       Sub->C = (expdt*(((1.0-2.0*xi*xi)/omegaDdt-xi*omega/omegaD)*sinDdt-(1.0+2.0*xi/omegadt)*cosDdt)+2.0*xi/omegadt)/Stiff;
	       Sub->D = (expdt*((2.0*xi*xi-1.0)*sinDdt/omegaDdt + 2.0*xi*cosDdt/omegadt)+(1.0-2.0*xi/omegadt))/Stiff;
#endif
	       Sub->E = -expdt * omega * omega * sinDdt / omegaD;
	       Sub->F = expdt * (cosDdt - xi*omega*sinDdt/omegaD);
#if _FLOAT_
	       Sub->G = (expdt*((omega*omega/omegaD + xi*omega/omegaDdt)*sinDdt + cosDdt/DeltaT)-1.0f/DeltaT)/Stiff;
	       Sub->H = (-expdt*(xi*omega*sinDdt/omegaD + cosDdt) + 1.0f)/(Stiff*DeltaT);
#else
	       Sub->G = (expdt*((omega*omega/omegaD + xi*omega/omegaDdt)*sinDdt + cosDdt/DeltaT)-1.0/DeltaT)/Stiff;
	       Sub->H = (-expdt*(xi*omega*sinDdt/omegaD + cosDdt) + 1.0)/(Stiff*DeltaT);
#endif


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

void Substructure_ExactSolutionESP_SDOF( const HYSL_FLOAT DispTdT, const HYSL_FLOAT ramp, const HYSL_FLOAT GAcc,const HYSL_FLOAT DeltaT, ExactSimESP_t *const Sub, HYSL_FLOAT *const fc )
{

     /* Compute initial velocity */
#if _FLOAT_
     Sub->AccTdT = (1.0f - ramp)*(-Sub->a0*Sub->Disp0 -Sub->a2*Sub->Vel0 -Sub->a3*Sub->Acc0)
	  + ramp*(-Sub->a0*Sub->DispT -Sub->a2*Sub->VelT -Sub->a3*Sub->AccT) + Sub->a0*DispTdT;

     Sub->VelTdT = (1.0f - ramp)*(Sub->Vel0 + Sub->a6*Sub->Acc0) + ramp*(Sub->VelT + Sub->a6*Sub->AccT) + Sub->a7*Sub->AccTdT;
#else
     Sub->AccTdT = (1.0 - ramp)*(-Sub->a0*Sub->Disp0 -Sub->a2*Sub->Vel0 -Sub->a3*Sub->Acc0)
	  + ramp*(-Sub->a0*Sub->DispT -Sub->a2*Sub->VelT -Sub->a3*Sub->AccT) + Sub->a0*DispTdT;

     Sub->VelTdT = (1.0 - ramp)*(Sub->Vel0 + Sub->a6*Sub->Acc0) + ramp*(Sub->VelT + Sub->a6*Sub->AccT) + Sub->a7*Sub->AccTdT;
#endif

     /* Backup and computer new forcce */
     Sub->Force_0 = Sub->Force_1;
     Sub->Force_1 = Sub->Stiff*DispTdT + Sub->Damp*Sub->VelTdT;
//     Sub->Force_1 = -Sub->Mass*(GAcc + Sub->AccTdT);

     /* Backup the displacements */
     Sub->Init_Disp = Sub->End_Disp;
     Sub->Init_Vel = Sub->End_Vel;

     /* Compute the new state */
     Sub->End_Disp = Sub->A*Sub->Init_Disp + Sub->B*Sub->Init_Vel + Sub->C*Sub->Force_0 + Sub->D*Sub->Force_1;
     Sub->End_Vel = Sub->E*Sub->Init_Disp + Sub->F*Sub->Init_Vel + Sub->G*Sub->Force_0 + Sub->H*Sub->Force_1;

     /* Compute the coupling force */
     *fc = Sub->Stiff*(Sub->End_Disp - DispTdT) + Sub->Damp*(Sub->End_Vel - Sub->VelTdT);
//     *fc = Sub->Stiff*Sub->End_Disp + Sub->Damp*Sub->End_Vel;

}


void Substructure_ExactSolutionESP_Destroy( ExactSimESP_t *const Sub )
{
     free( Sub->Description );
}
