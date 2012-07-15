1#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omph.h>

#include "Netlib.h"
#include "ErrorHandling.h"
#include "Exact_Solution.h"

void Compute_EigenValues_EigenVectors ( const int Num_DOF, float *const MatrixA, float *const MatrixB, float *const Eigen_Values, float *const Eigen_Vectors )
{
     int i, j;   /* Counters */
     int one = 1, Size;
     int lda, ldb, info;
     int lwork; /* Dimension of the array work */
     float *work, *TempMat, temp;

     Size = Num_DOF;
     lwork = 3*Size - 1;
     lda = max (1, Size);
     ldb = lda;

     TempMat = (float*) malloc(Num_DOF*Num_DOF*sizeof (float) );
     work = (float*) malloc( lwork*sizeof (float) );

     /* DSYGV_:On Entry Eigen_Vectors must contain the Matrix A */
     Size = Num_DOF * Num_DOF;
     scopy_( &Size, MatrixA, &one, Eigen_Vectors, &one );
     scopy_( &Size, MatrixB, &one, TempMat, &one );

     Size = Num_DOF;
     ssygv_( &one, "V", "L", &Size, Eigen_Vectors, &lda, TempMat, &ldb, Eigen_Values, work, &lwork, &info );

     for ( i = 0; i < Num_DOF - 1; i++){
	  
	  /* Order the Eigenvalues and eigenvectors in ascendent order */
	  if ( Eigen_Values[i] > Eigen_Values[i+1] ){
	       /* Swap Eigenvalues */
	       temp = Eigen_Values[i];
	       Eigen_Values[i] = Eigen_Values[i+1];
	       Eigen_Values[i+1] = temp;
	       /* Now Swap EigenVectors */
	       for( j = 0; j < Num_DOF; j++ ){
		    temp = Eigen_Vectors[Num_DOF*(i+1) + j];
		    Eigen_Vectors[Num_DOF*i + j] = Eigen_Vectors[Num_DOF*(i+1) + j];
		    Eigen_Vectors[Num_DOF*(i+1) + j] = temp;
	       }
	  }
     }

     /* Free the dynamically allocated memory */
     free( TempMat );
     free( work );
}

void Compute_Damp_Matrix( const float a0, const float a1, const int Num_DOF, const float *const Mass, const float *const Stiff, float *const Damp )
{

     int i;  /* A counter */

#pragma omp parallel for
     for ( i = 0; i < Num_DOF*Num_DOF; i++ ){
	  Damp[i] = Mass[i]*a0 + Stiff[i]*a1;
     }
}

void Compute_Damping_Ratios_and_Matrix ( const float a0, const float a1, const int Num_DOF, const float *const Eigen_Values, float *const Damping_Ratios )
{
     int i; /* A counter */
     float omega;
     
     for ( i = 0; i < Num_DOF; i++ ){
	  omega = sqrt( Eigen_Values[i] );
	  Damping_Ratios[i] = a0/(2.0*omega) + a1*omega/2.0;
     }
}

void ExactSolution_Init( const int Num_DOF, const float a0, const float a1, const char* MassFile, const char *StiffFile, TMD_Sim *const Num )
{
     fopen *InFile_Mass, *InFile_Stiff;
     float Eigen_Values;

     /* Init the variables */
     Num->Mass = (float *) calloc( Num_DOF*Num_DOF, sizeof(float) );
     Num->Stiff = (float *) calloc( Num_DOF*Num_DOF, sizeof(float) );
     Num->Damp = (float *) calloc( Num_DOF*Num_DOF, sizeof(float) );
     Num->Eigen_Vectors = (float *) calloc( Num_DOF, sizeof(float) );

     Num->Damping_Ratios = (float *) calloc( Num_DOF, sizeof(float) );

     Num->Omega = (float *) calloc( Num_DOF, sizeof(float) );

     Num->Disp0 = (float *) calloc( Num_DOF, sizeof(float) );
 
     Eigen_Values = (float *) calloc( Num_DOF, sizeof(float) );

     /* Read the matrices from a text file */
     InFile_Mass = fopen( MassFile, "r" );
     InFile_Stiff = fopen( StiffFile, "r" );

     if( InFile_Mass == NULL || InFile_Stiff == NULL ){
	  if( InFile_Mass == NULL ){
	       ErrorFileAndExit( "Could not open: ", "Mass.txt" );
	  } else {
	       ErrorFileAndExit( "Could not open: ", "Stiff.txt" );
	  }
     } else {
	  for( i = 0; i < Num_DOF*Num_DOF; i++ ){
	       fscanf( InFile_Mass, "%f", &Num->Mass[0] );
	       fscanf( InFile_Mass, "%f", &Num->Stiff[0] );
	  }
     }

     fclose( InFile_Mass );
     fclose( InFile_Stiff );

     /* Compute the eigen values and eigenvectors */
     Compute_EigenValues_EigenVectors ( Num_DOF, Num->Mass, Num->Stiff, Num->Eigen_Values,
					Num->Eigen_Vectors );

     /* Calculate the damping matrix */
     Compute_Damp_Matrix ( a0, a1, Num_DOF, Num->Mass, Num->Stiff, Num->Damp );

     /* Calculate the damping ratios */
     Compute_Damping_Ratios ( a0, a1, Num_DOF, Eigen_Values, Num->Damping_Ratios );


     /* Calculate omega = sqrt( eigen_values ) */
     for ( i = 0; i < Num_DOF; i++ ){
	  if ( Eigen_Values[i] > 0 ){
	       Num->Omega[i] = sqrtf( Num->Eigen_Values[i] );
	  } else {
	       PrintErrorAndExit( "Duhamel Integration: non-positive eigenvalue");
	  }
     }

     free( Eigen_Values );
}

void Exact_Solution ( const float *const u0c, const float DeltaT, const int Num_DOF, float *const fc )
{

     int i; /* Counter */

     /* Variables required by the BLAS routines */
     int Size;
     static int incx = 1, incy = 1;     /* Stride in the vectors */
     static char uplo = 'L';
     static char trans;
     static float Alpha, Beta;

     float om2, om, omp, zeta=0;
     float expo, sino, coso;
     float x_part, x_hom, a, b;    

     float *Init_Disp, *Init_Vel, *Force;
     float *End_Disp, *End_Vel, *End_Acc;
     float *Modal_Force, *Modal_Disp, *Modal_Vel, *temp;

     /* Calculate new velocity. v = (u0c - Disp0)/Delta_T */
     Size = Num_DOF;
     scopy_( &Size, u0c, &incx, &Init_Vel, &incy );
     Alpha = -1.0;
     saxpy_( &Size, &Alpha, Num->Disp0, &incx, Num->Init_Vel, &incy );
     Alpha = 1.0/DeltaT;
     sscal_( &Size, &Alpha, Num->Init_Vel, &incx );


     Init_Disp = (float *) calloc( Num_DOF, sizeof (float) );
     Init_Vel = (float *) calloc( Num_DOF, sizeof (float) );
     Force = (float *) calloc( Num_DOF, sizeof (float) );

     Modal_Force = (float *) calloc( Num_DOF, sizeof (float) );
     Modal_Disp = (float *) calloc( Num_DOF, sizeof (float) );
     Modal_Vel = (float *) calloc( Num_DOF, sizeof (float) );


     temp = (float *) calloc( Num_DOF, sizeof (float) );

     /* Compute new input force. Force = M*u0c + C*Init_Vel */
     Size = Num_DOF; Alpha = 1.0; Beta = 0.0;
     ssymv( &uplo, &Size, &Alpha, Num->Mass, u0c, &incx, &Beta, temp, &incy );
     Beta = 1.0;
     ssymv( &uplo, &Size, &Alpha, Num->Damp, Init_Vel, &incx, &Beta, Force, &incy );

     /* Transform into modal coordinates */
     Size = Num_DOF*Num_DOF;
     scopy_( &Size, Init_Disp, &incx, temp, &incy );
     Size = Num_DOF; Beta = 0.0;
     ssymv( &uplo, &Size, &Alpha, Num->Mass, temp, &incx, &Beta, Init_Disp );

     Size = Num_DOF*Num_DOF;
     scopy_( &Size, Num->Init_Vel, &incx, temp, &incy );
     Size = Num_DOF;
     ssymv( &uplo, &Size, &Alpha, Num->Mass, temp, &incx, &Beta, Init_Vel );
     
     free( temp );

     trans = 'T';
     sgemv( &trans, &Size, &Size, &Alpha, Num->Eigen_Vectors, &Size, Init_Disp, &incx, &Beta, Modal_Disp, &incy );
     sgemv( &trans, &Size, &Size, &Alpha, Num->Eigen_Vectors, &Size, Init_Vel, &incx, &Beta, Modal_Vel, &incy );
     sgemv( &trans, &Size, &Size, &Alpha, Num->Eigen_Vectors, &Size, Force, &incx, &Beta, Modal_Force, &incy );

     free( Init_Disp );
     free( Init_Vel );
     free( Force );

#pragma omp parallel for
     for (i = 0; i < Num_DOF; i++) {

	  om2	= Num->Eigen_Values[i];

	  if (Damping_Ratios != NULL){
	       zeta = Num->Damping_Ratios[i];
	  }

	  if (om2 > 0){
	       om = sqrt(om2);
	       if (zeta >= 0 && zeta < 1) {
		    omp	= om*sqrt(1-zeta*zeta);
		    x_part	= Modal_Force[i]/om2;
		    a	= -x_part + Modal_Disp[i];
		    b	= (zeta*om*a + Modal_Vel[i])/omp;
		    coso	= cos(omp*Delta_T);
		    sino	= sin(omp*Delta_T);
		    expo	= exp(-zeta*om*Delta_T);
		    x_hom	= expo*(a*coso + b*sino);
		    Modal_Vel[i]	= -zeta*om*x_hom
			 + expo*omp*(-a*sino + b*coso);
		    Modal_Disp[i]	= x_hom + x_part;
		    /* use Modal_Force for accelerations*/
		    Modal_Force[i]	-= om2*Modal_Disp[i];
		    Modal_Force[i]	-= 2*om*zeta*Modal_Vel[i];
	       } else {
		    fprintf( stderr, "Duhamel Integration: negative or overcritical damping.\n");
		    exit( EXIT_FAILURE );
	       }
	  } else {
	       fprintf( stderr, "Duhamel Integration: non-positive eigenvalue.\n");
	  }
     }

     End_Disp = (float *) calloc( Num_DOF, sizeof (float) );
     End_Vel = (float *) calloc( Num_DOF, sizeof (float) );
     End_Acc = (float *) calloc( Num_DOF, sizeof (float) );

     /* Transform into physical coordinates */
     Multiply_Matrix_Vector( Num_DOF, Num_DOF, "N", 1.0, Eigen_Vectors, Modal_Disp, 0.0, End_Disp );
     Multiply_Matrix_Vector( Num_DOF, Num_DOF, "N", 1.0, Eigen_Vectors, Modal_Vel, 0.0, End_Velo );
     Multiply_Matrix_Vector( Num_DOF, Num_DOF, "N", 1.0, Eigen_Vectors, Modal_Force, 0.0, End_Accel);

#pragma omp parallel for
     for( i = 0; i< Num_DOF; i++ ){
	  fc[i] = End_Accel[i] - Acc0[i];
	  End_Accel[i] = Acc0[i];
	  Disp0[i] = u0c[i]
     }   

     /* Free dynamically allocated memory */
     free( Modal_Force );
     free( Modal_Disp );
     free( Modal_Vel );

     free( End_Disp );
     free( End_Vel );
     free( End_Acc );
}
