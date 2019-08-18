#include "GSSSS.h"
#include "MatrixVector.h"
#include "Error_Compensation.h"
#include "Definitions.h"

#if _MKL_
#include <mkl_blas.h>
#else
#include "Netlib.h"
#endif

void GSSSS_Init( const hysl_float_t Rho1, const hysl_float_t Rho2, const hysl_float_t Rho3, TIntegration_GSSSS_t *const GSSSS )
{

     GSSSS->w1 = -15.0*(1.0 - 2.0*Rho3)/(1.0 - 4.0*Rho3);
     GSSSS->w2 = 15.0*(3.0 - 4.0*Rho3)/(1.0 - 4.0*Rho3);
     GSSSS->w3 = -35.0*(1.0 - Rho3)/(1.0 - 4.0*Rho3);


     GSSSS->W1 = (0.5 + GSSSS->w1/3.0 + GSSSS->w2/4.0 + GSSSS->w3/5.0)/(1.0 + GSSSS->w1/2.0 + GSSSS->w2/3.0 + GSSSS->w3/4.0);

     GSSSS->A1W1 = 1.0/(1.0 + Rho3);
     GSSSS->A2W2 = 1.0/(2.0*(1.0 + Rho3));
     GSSSS->A3W3 = 1.0/((1.0 + Rho1)*(1.0 + Rho2)*(1.0 + Rho3));
     GSSSS->A4W1 = 1.0/(1.0 + Rho3);
     GSSSS->A5W2 = (3.0 + Rho1 + Rho2 - Rho1*Rho2)/(2.0*(1.0 + Rho1)*(1.0 + Rho2)*(1.0 + Rho3));
     GSSSS->A6W1 = (2.0 + Rho1 + Rho2 + Rho3- Rho1*Rho2*Rho3)/((1.0 + Rho1)*(1.0 + Rho2)*(1.0 + Rho3));

     GSSSS->l1 = 1.0;
     GSSSS->l2 = 0.5;
     GSSSS->l3 = 1.0/((1.0 + Rho1)*(1.0 + Rho2));
     GSSSS->l4 = 1.0;
     GSSSS->l5 = (3.0 + Rho1 + Rho2 - Rho1*Rho2)/(2.0*(1.0 + Rho1)*(1.0 + Rho2));

     printf("A1W1 %lE A2W2 %lE A3W3 %lE A4W1 %lE A5W2 %lE A6W1 %lE W1 %lE l1 %lE l2 %lE l3 %lE l4 %lEl5 %lE\n", GSSSS->A1W1, GSSSS->A2W2, GSSSS->A3W3, GSSSS->A4W1, GSSSS->A5W2, GSSSS->A6W1, GSSSS->W1,
	    GSSSS->l1, GSSSS->l2, GSSSS->l3, GSSSS->l4, GSSSS->l5 );

}

void GSSSS_EffectiveForce_AForm( const MatrixVector_t *const Mass, const MatrixVector_t *const Damp,
				 const MatrixVector_t *const Stiff, const MatrixVector_t *const DispT,
				 const MatrixVector_t *const VelT, const MatrixVector_t *const AccT,
				 MatrixVector_t *const Tempvec, const TIntegration_GSSSS_t *const GSSSS,
				 const hysl_float_t DeltaT, MatrixVector_t *const Eff_ForceT )
{

     int incx = 1, incy = 1;  /* Stride in the vectors */
     char uplo = 'L';         /* The lower part (upper part in C) will be used and the upper part (lower part
			       * in C) will strictly not be referenced */
     hysl_float_t Alpha, Beta;  /* Constants for the BLAS routines */
     

    /* BLAS: tempvec = Disp */
     hysl_copy( &Tempvec->Rows, DispT->Array, &incx, Tempvec->Array, &incy );
     printf( "DispT %lE Tempvec %lE\n", DispT->Array[46], Tempvec->Array[46]);
     /* BLAS: tempvec = Disp + A1W1*Vel*DeltaT = tempvec + A1W1*Vel*DeltaT */
     Alpha = GSSSS->A1W1*DeltaT;
     hysl_axpy( &Tempvec->Rows, &Alpha, VelT->Array, &incx, Tempvec->Array, &incy );
     printf( "VelT %lE Tempvec %lE\n", VelT->Array[46], Tempvec->Array[46]);
     /* BLAS: tempvec = Disp + A1W1*Vel*DeltaT + (A2W2 - A3W3)*Acc*DeltaT^2 = tempvec + (A2W2
      * -A3W3)*a10*DeltaT^2 */
     Alpha = (GSSSS->A2W2 - GSSSS->A3W3)*DeltaT*DeltaT;
     hysl_axpy( &Tempvec->Rows, &Alpha, AccT->Array, &incx, Tempvec->Array, &incy );
     printf( "AccT %lE Tempvec %lE\n", AccT->Array[46], Tempvec->Array[46]);
     /* BLAS: Eff_Force = -Stiff*(Disp + A1W1*Vel*DeltaT + (A2W2 - A3W3)*Acc*DeltaT^2) = -Stiff*tempvec */
     Alpha = -1.0; Beta = 0.0;
     hysl_symv( &uplo, &Tempvec->Rows, &Alpha, Stiff->Array, &Tempvec->Rows, Tempvec->Array, &incx, &Beta,
		Eff_ForceT->Array, &incy );

     /* BLAS: tempvec = Vel */
     hysl_copy( &Tempvec->Rows, VelT->Array, &incx, Tempvec->Array, &incy );
     printf( "DispT %lE Tempvec %lE\n", VelT->Array[46], Tempvec->Array[46]);
     /* BLAS: tempvec = Vel + (A4W1 -A5W2)*DeltaT*Acc = tempvec + (A4W1 - A5W2)*DeltaT*Acc */
     Alpha = (GSSSS->A4W1 - GSSSS->A5W2)*DeltaT;
     hysl_axpy( &Tempvec->Rows, &Alpha, AccT->Array, &incx, Tempvec->Array, &incy );
     printf( "DispT %lE Tempvec %lE\n", AccT->Array[46], Tempvec->Array[46]);
     /* BLAS: Eff_Force = -Stiff*(Disp + A1W1*Vel*DeltaT + (A2W2 - A3W3)*Acc*DeltaT^2) -C*(Vel + (A4W1
      * -A5W2)*DeltaT*Acc)= Eff_Force - Damp*tempvec */
     Alpha = -1.0; Beta = 1.0;
     hysl_symv( &uplo, &Tempvec->Rows, &Alpha, Damp->Array, &Tempvec->Rows, Tempvec->Array, &incx, &Beta,
		Eff_ForceT->Array, &incy );

     /* BLAS: Eff_Force = Stiff*(Disp + A1W1*Vel*DeltaT + (A2W2 - A3W3)*Acc*DeltaT^2) -C*(Vel + (A4W1
      * -A5W2)*DeltaT*Acc) -(1-A6W1)*M*Acc= Eff_Force - (1-A6W1)*M*Acc */
     Alpha = -(1.0 - GSSSS->A6W1)*DeltaT;
     hysl_symv( &uplo, &Tempvec->Rows, &Alpha, Mass->Array, &Tempvec->Rows, AccT->Array, &incx, &Beta,
		Eff_ForceT->Array, &incy );
}

void GSSSS_EffectiveForce_AForm_PS( const MatrixVector_t *const Mass, const MatrixVector_t *const Damp,
				    const MatrixVector_t *const Stiff, const MatrixVector_t *const DispT,
				    const MatrixVector_t *const VelT, const MatrixVector_t *const AccT,
				    MatrixVector_t *const Tempvec, const TIntegration_GSSSS_t *const GSSSS,
				    const hysl_float_t DeltaT, MatrixVector_t *const Eff_ForceT )
{

     int incx = 1, incy = 1;  /* Stride in the vectors */
     char uplo = 'L';         /* The lower part (upper part in C) will be used and the upper part (lower part
			       * in C) will strictly not be referenced */
     hysl_float_t Alpha, Beta;  /* Constants for the BLAS routines */
     

    /* BLAS: tempvec = Disp */
     hysl_copy( &Tempvec->Rows, DispT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = Disp + A1W1*Vel*DeltaT = tempvec + A1W1*Vel*DeltaT */
     Alpha = GSSSS->A1W1*DeltaT;
     hysl_axpy( &Tempvec->Rows, &Alpha, VelT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = Disp + A1W1*Vel*DeltaT + (A2W2 - A3W3)*Acc*DeltaT^2 = tempvec + (A2W2
      * -A3W3)*a10*DeltaT^2 */
     Alpha = (GSSSS->A2W2 - GSSSS->A3W3)*DeltaT*DeltaT;
     hysl_axpy( &Tempvec->Rows, &Alpha, AccT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: Eff_Force = -Stiff*(Disp + A1W1*Vel*DeltaT + (A2W2 - A3W3)*Acc*DeltaT^2) = -Stiff*tempvec */
     Alpha = -1.0; Beta = 0.0;
     hysl_spmv( &uplo, &Tempvec->Rows, &Alpha, Stiff->Array, Tempvec->Array, &incx, &Beta,
		Eff_ForceT->Array, &incy );

     /* BLAS: tempvec = Vel */
     hysl_copy( &Tempvec->Rows, VelT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = Vel + (A4W1 -A5W2)*DeltaT*Acc = tempvec + (A4W1 - A5W2)*DeltaT*Acc */
     Alpha = (GSSSS->A4W1 - GSSSS->A5W2)*DeltaT;
     hysl_axpy( &Tempvec->Rows, &Alpha, AccT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: Eff_Force = -Stiff*(Disp + A1W1*Vel*DeltaT + (A2W2 - A3W3)*Acc*DeltaT^2) -C*(Vel + (A4W1
      * -A5W2)*DeltaT*Acc)= Eff_Force - Damp*tempvec */
     Alpha = -1.0; Beta = 1.0;
     hysl_spmv( &uplo, &Tempvec->Rows, &Alpha, Damp->Array, Tempvec->Array, &incx, &Beta,
		Eff_ForceT->Array, &incy );

     /* BLAS: Eff_Force = Stiff*(Disp + A1W1*Vel*DeltaT + (A2W2 - A3W3)*Acc*DeltaT^2) -C*(Vel + (A4W1
      * -A5W2)*DeltaT*Acc) -(1-A6W1)*M*Acc= Eff_Force - (1-A6W1)*M*Acc */
     Alpha = -(1.0 - GSSSS->A6W1)*DeltaT;
     hysl_spmv( &uplo, &Tempvec->Rows, &Alpha, Mass->Array, AccT->Array, &incx, &Beta,
		Eff_ForceT->Array, &incy );
}

void GSSSS_ComputeDisplacement_AForm( const MatrixVector_t *const DispT, const MatrixVector_t *const VelT,
				      const MatrixVector_t *const AccT, const MatrixVector_t *const AccTdT,
				      const TIntegration_GSSSS_t *const GSSSS, const hysl_float_t DeltaT,
				      MatrixVector_t *const DispTdT )
{
     int incx = 1, incy = 1;  /* Stride in the vectors */
     hysl_float_t Alpha;        /* Constant for the BLAS routines */

     /* BLAS: DispTdT = DispT */
     hysl_copy( &DispTdT->Rows, DispT->Array, &incx, DispTdT->Array, &incy ); 
     /* BLAS: DispTdT = DispT + l1*VelT*DeltaT = DispTdT + l1*VelT*DeltaT */
     Alpha = GSSSS->l1*DeltaT;
     hysl_axpy( &DispTdT->Rows, &Alpha, VelT->Array, &incx, DispTdT->Array, &incy );
     /* BLAS: DispTdT = DispT + l1*VelT*DeltaT + (l2 - l3)*AccT*DeltaT^2 = DispTdT + (l2 - l3)*AccT*DeltaT^2 */
     Alpha = (GSSSS->l2 - GSSSS->l3)*DeltaT*DeltaT;
     hysl_axpy( &DispTdT->Rows, &Alpha, AccT->Array, &incx, DispTdT->Array, &incy );
     /* BLAS: DispTdT = DispT + l1*VelT*DeltaT + (l2 - l3)*AccT*DeltaT^2 + l3*AccTdT*DeltaT^2 = DispTdT +
      * l3*AccTdT*DeltaT^2 */
     Alpha = GSSSS->l3*DeltaT*DeltaT;
     hysl_axpy( &DispTdT->Rows, &Alpha, AccTdT->Array, &incx, DispTdT->Array, &incy );
}

void GSSSS_ComputeVelocity_AForm( const MatrixVector_t *const VelT, const MatrixVector_t *const AccT,
				  const MatrixVector_t *const AccTdT, const TIntegration_GSSSS_t *const GSSSS,
				  const hysl_float_t DeltaT, MatrixVector_t *const VelTdT )
{
     int incx = 1, incy= 1;  /* Stride in the vectors */
     hysl_float_t Alpha;       /* Constant for the BLAS routines */

     /* BLAS: VelTdT = VelT */
     hysl_copy( &VelTdT->Rows, VelT->Array, &incx, VelTdT->Array, &incy );
     /* BLAS: VelTdT = VelT + (l4 - l5)*AccT*DeltaT = VelTdT + (l4 - l5)*AccT*DeltaT */
     Alpha = (GSSSS->l4 - GSSSS->l5)*DeltaT;
     hysl_axpy( &VelTdT->Rows, &Alpha, AccT->Array, &incx, VelTdT->Array, &incy );
     /* BLAS: VelTdT = VelT + (l4 - l5)*AccT*DeltaT + l5*AccTdT*DeltaT = VelTdT + l5*AccTdT*DeltaT */
     Alpha = GSSSS->l5*DeltaT;
     hysl_axpy( &VelTdT->Rows, &Alpha, AccTdT->Array, &incx, VelTdT->Array, &incy );
}

void GSSSS_Compute_NewState( const MatrixVector_t *const IGain, const MatrixVector_t *const Eff_ForceT,
			     const MatrixVector_t *const LoadTdT, const MatrixVector_t *const Err_ForceT,
			     const MatrixVector_t *const LoadT, const TIntegration_GSSSS_t *const GSSSS,
			     MatrixVector_t *const Tempvec, MatrixVector_t *const VecTdT_0 )
{
     int incx = 1, incy = 1;               /* Stride in the vectors */
     hysl_float_t Alpha = 0.0, Beta = 0.0;   /* Constants for the BLAS routines */
     char uplo = 'L';                      /* 
					    * The lower part (upper part in C) will be used and the upper part
					    * (lower part in C) will strictly not be referenced */

     /* BLAS: tempvec = Eff_ForceT */
     hysl_copy( &Tempvec->Rows, Eff_ForceT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = Eff_ForceT + W1*LoadTdT = tempvec + W1*LoadTdT */
     Alpha = GSSSS->W1;
     hysl_axpy( &Tempvec->Rows, &Alpha, LoadTdT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = Eff_ForceT + W1*LoadTdT + (1-W1)*LoadT = tempvec + (1 - W1)*LoadT */
     Alpha = (1.0 - GSSSS->W1);
     hysl_axpy( &Tempvec->Rows, &Alpha, LoadT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = Eff_ForceT + W1*LoadTdT + (1-W1)*LoadT - Err_ForceT = tempvec - Err_ForceT */
     Alpha = -1.0;
     hysl_axpy( &Tempvec->Rows, &Alpha, Err_ForceT->Array, &incx, Tempvec->Array, &incy );

     /* BLAS: VecTdT_0 = IGain*(Eff_ForceT + W1*LoadTdT + (1-W1)*LoadT - Err_ForceT) = IGain*tempvec */
     Alpha = 1.0;
     hysl_symv( &uplo, &Tempvec->Rows, &Alpha, IGain->Array, &Tempvec->Rows, Tempvec->Array, &incx, &Beta,
		VecTdT_0->Array, &incy );
}

void GSSSS_Compute_NewState_PS( const MatrixVector_t *const IGain, const MatrixVector_t *const Eff_ForceT,
				const MatrixVector_t *const LoadTdT, const MatrixVector_t *const Err_ForceT,
				const MatrixVector_t *const LoadT, const TIntegration_GSSSS_t *const GSSSS,
				MatrixVector_t *const Tempvec, MatrixVector_t *const VecTdT_0 )
{
     int incx = 1, incy = 1;               /* Stride in the vectors */
     hysl_float_t Alpha = 0.0, Beta = 0.0;   /* Constants for the BLAS routines */
     char uplo = 'L';                      /* 
					    * The lower part (upper part in C) will be used and the upper part
					    * (lower part in C) will strictly not be referenced */

     /* BLAS: tempvec = Eff_ForceT */
     hysl_copy( &Tempvec->Rows, Eff_ForceT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = Eff_ForceT + W1*LoadTdT = tempvec + W1*LoadTdT */
     Alpha = GSSSS->W1;
     hysl_axpy( &Tempvec->Rows, &Alpha, LoadTdT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = Eff_ForceT + W1*LoadTdT + (1-W1)*LoadT = tempvec + (1 - W1)*LoadT */
     Alpha = (1.0 - GSSSS->W1);
     hysl_axpy( &Tempvec->Rows, &Alpha, LoadT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = Eff_ForceT + W1*LoadTdT + (1-W1)*LoadT - Err_ForceT = tempvec - Err_ForceT */
     Alpha = -1.0;
     hysl_axpy( &Tempvec->Rows, &Alpha, Err_ForceT->Array, &incx, Tempvec->Array, &incy );

     /* BLAS: VecTdT_0 = IGain*(Eff_ForceT + W1*LoadTdT + (1-W1)*LoadT - Err_ForceT) = IGain*tempvec */
     Alpha = 1.0;
     hysl_spmv( &uplo, &Tempvec->Rows, &Alpha, IGain->Array, Tempvec->Array, &incx, &Beta, VecTdT_0->Array,
		&incy );
}

void GSSSS_ErrorForce_PID( const MatrixVector_t *const Mass, const MatrixVector_t *const Damp, const MatrixVector_t *Stiff,
			   const MatrixVector_t *const AccTdT, const MatrixVector_t *const AccT, const MatrixVector_t *const VelT,
			   const MatrixVector_t *const DispT, const MatrixVector_t *const fc, const MatrixVector_t *const LoadTdT,
			   const MatrixVector_t *const LoadT, const TIntegration_GSSSS_t *const GSSSS, const hysl_float_t DeltaT,
			   const MatrixVector_t *const Tempvec, const PID_t *const PID, MatrixVector_t *const fe )
{

     int incx = 1, incy = 1;  /* Stride in the vectors for BLAS routines */
     hysl_float_t Alpha, Beta;  /* Constants to use in the Sparse BLAS routines */
     char uplo = 'L';

     /* BLAS: tempvec = acct */
     hysl_copy( &Tempvec->Rows, AccT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = acct + A6W1*acctdt = tempvec + A6W1*acctdt */
     Alpha = GSSSS->A6W1;
     hysl_axpy( &Tempvec->Rows, &Alpha, AccTdT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = acct *A6W1*(acctdt - acct) = tempvec - A6W1*acct */
     Alpha = -GSSSS->A6W1;
     hysl_axpy( &Tempvec->Rows, &Alpha, AccT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: fe = -M*(acct *A6W1*(acctdt - acct)) = M*Tempvec */
     Alpha = -1.0; Beta = 0.0;
     hysl_symv( &uplo, &fe->Rows, &Alpha, Mass->Array, &fe->Rows, Tempvec->Array, &incx, &Beta, fe->Array, &incy );

     /* BLAS: tempvec = velt */
     hysl_copy( &Tempvec->Rows, AccT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = velt + A4W1*acct*DeltaT = tempvec + A4W1*DeltaT*acct */
     Alpha = GSSSS->A4W1*DeltaT;
     hysl_axpy( &Tempvec->Rows, &Alpha, AccT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = velt + A4W1*acct*DeltaT + A5W2*DeltaT*acctdt = tempvec + A5W2*DeltaT*acctdt */
     Alpha = GSSSS->A5W2*DeltaT;
     hysl_axpy( &Tempvec->Rows, &Alpha, AccTdT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = velt + A4W1*acct*DeltaT + A5W2*DeltaT*(acctdt - acct) = tempvec - A5W2*DeltaT*acct*/
     Alpha = -GSSSS->A5W2*DeltaT;
     hysl_axpy( &Tempvec->Rows, &Alpha, AccT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: fe = -M*(acct *A6W1*(acctdt - acct)) - C*(velt + A4W1*acct*DeltaT + A5W2*DeltaT*(acctdt - acct) =
      * fe + C*Tempvec */
     Alpha = -1.0; Beta = 1.0;
     hysl_symv( &uplo, &fe->Rows, &Alpha, Damp->Array, &fe->Rows, Tempvec->Array, &incx, &Beta, fe->Array, &incy );

     /* BLAS: tempvec = dispt */
     hysl_copy( &Tempvec->Rows, DispT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = dispt + A1W1*velt*DeltaT = tempvec + A1W1*DeltaT*velt */
     Alpha = GSSSS->A1W1*DeltaT;
     hysl_axpy( &Tempvec->Rows, &Alpha, VelT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = dispt + A1W1*acct*DeltaT + A2W2*DeltaT*DeltaT*acct = tempvec + A2W2*DeltaT*DeltaT*acct */
     Alpha = GSSSS->A2W2*DeltaT*DeltaT;
     hysl_axpy( &Tempvec->Rows, &Alpha, AccT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = dispt + A1W1*acct*DeltaT + A2W2*DeltaT*DeltaT*acct + A3W2*DeltaT*DeltaT*acctdt= tempvec + A3W3*DeltaT*DeltaT*acctdt */
     Alpha = GSSSS->A3W3*DeltaT*DeltaT;
     hysl_axpy( &Tempvec->Rows, &Alpha, AccTdT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = dispt + A1W1*acct*DeltaT + A2W2*DeltaT*DeltaT*acct + A3W2*DeltaT*DeltaT*(acctdt -
      * acct) = tempvec - A3W3*DeltaT*DeltaT*acct */
     Alpha = -GSSSS->A3W3*DeltaT*DeltaT;
     hysl_axpy( &Tempvec->Rows, &Alpha, AccT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: fe = -M*(acct *A6W1*(acctdt - acct)) - C*(velt + A4W1*acct*DeltaT + A5W2*DeltaT*(acctdt - acct) -
      * K*(dispt + A1W1*DeltaT*velt + A2W2*DeltaT*DeltaT*acct + A3W3*DeltaT*DeltaT*(acctdt - acct) =
      * fe + K*Tempvec */
     Alpha = -1.0; Beta = 1.0;
     hysl_symv( &uplo, &fe->Rows, &Alpha, Stiff->Array, &fe->Rows, Tempvec->Array, &incx, &Beta, fe->Array, &incy );

     /* BLAS: fe = M*(acct *A6W1*(acctdt - acct)) + C*(velt + A4W1*acct*DeltaT + A5W2*DeltaT*(acctdt - acct) +
      * K*(dispt + A1W1*DeltaT*velt + A2W2*DeltaT*DeltaT*acct + A3W3*DeltaT*DeltaT*(acctdt - acct) +
      * W1*LoadTdT = fe + W1*LoadTdT */
     Alpha = GSSSS->W1;
     hysl_axpy( &fe->Rows, &Alpha, LoadTdT->Array, &incx, fe->Array, &incy );

     /* BLAS: fe = M*(acct *A6W1*(acctdt - acct)) + C*(velt + A4W1*acct*DeltaT + A5W2*DeltaT*(acctdt - acct) +
      * K*(dispt + A1W1*DeltaT*velt + A2W2*DeltaT*DeltaT*acct + A3W3*DeltaT*DeltaT*(acctdt - acct) +
      * W1*LoadTdT + (1-W1)*LoadT = fe + (1 - W1)*LoadT */
     Alpha = (1.0 - GSSSS->W1);
     hysl_axpy( &fe->Rows, &Alpha, LoadT->Array, &incx, fe->Array, &incy );

     /* BLAS: fe = M*(acct *A6W1*(acctdt - acct)) + C*(velt + A4W1*acct*DeltaT + A5W2*DeltaT*(acctdt - acct) +
      * K*(dispt + A1W1*DeltaT*velt + A2W2*DeltaT*DeltaT*acct + A3W3*DeltaT*DeltaT*(acctdt - acct) +
      * W1*LoadTdT + (1-W1)*LoadT + fc = fe + fc */
     Alpha = 1.0;
     hysl_axpy( &fe->Rows, &Alpha, fc->Array, &incx, fe->Array, &incy );

     /* Apply the PID */
     /* BLAS: fe = P*fe */
     Alpha = -PID->P;
     hysl_scal( &fe->Rows, &Alpha, fe->Array, &incx );
}


void GSSSS_ErrorForce_PID_PS( const MatrixVector_t *const Mass, const MatrixVector_t *const Damp, const MatrixVector_t *Stiff,
			      const MatrixVector_t *const AccTdT, const MatrixVector_t *const AccT, const MatrixVector_t *const VelT,
			      const MatrixVector_t *const DispT, const MatrixVector_t *const fc, const MatrixVector_t *const LoadTdT,
			      const MatrixVector_t *const LoadT, const TIntegration_GSSSS_t *const GSSSS, const hysl_float_t DeltaT,
			      const MatrixVector_t *const Tempvec, const PID_t *const PID, MatrixVector_t *const fe )
{

     int incx = 1, incy = 1;  /* Stride in the vectors for BLAS routines */
     hysl_float_t Alpha, Beta;  /* Constants to use in the Sparse BLAS routines */
     char uplo = 'L';

     /* BLAS: tempvec = acct */
     hysl_copy( &Tempvec->Rows, AccT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = acct + A6W1*acctdt = tempvec + A6W1*acctdt */
     Alpha = GSSSS->A6W1;
     hysl_axpy( &Tempvec->Rows, &Alpha, AccTdT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = acct *A6W1*(acctdt - acct) = tempvec - A6W1*acct */
     Alpha = -GSSSS->A6W1;
     hysl_axpy( &Tempvec->Rows, &Alpha, AccT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: fe = -M*(acct *A6W1*(acctdt - acct)) = M*Tempvec */
     Alpha = -1.0; Beta = 0.0;
     hysl_spmv( &uplo, &fe->Rows, &Alpha, Mass->Array, Tempvec->Array, &incx, &Beta, fe->Array, &incy );

     /* BLAS: tempvec = velt */
     hysl_copy( &Tempvec->Rows, AccT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = velt + A4W1*acct*DeltaT = tempvec + A4W1*DeltaT*acct */
     Alpha = GSSSS->A4W1*DeltaT;
     hysl_axpy( &Tempvec->Rows, &Alpha, AccT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = velt + A4W1*acct*DeltaT + A5W2*DeltaT*acctdt = tempvec + A5W2*DeltaT*acctdt */
     Alpha = GSSSS->A5W2*DeltaT;
     hysl_axpy( &Tempvec->Rows, &Alpha, AccTdT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = velt + A4W1*acct*DeltaT + A5W2*DeltaT*(acctdt - acct) = tempvec - A5W2*DeltaT*acct*/
     Alpha = -GSSSS->A5W2*DeltaT;
     hysl_axpy( &Tempvec->Rows, &Alpha, AccT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: fe = -M*(acct *A6W1*(acctdt - acct)) - C*(velt + A4W1*acct*DeltaT + A5W2*DeltaT*(acctdt - acct) =
      * fe + C*Tempvec */
     Alpha = -1.0; Beta = 1.0;
     hysl_spmv( &uplo, &fe->Rows, &Alpha, Damp->Array, Tempvec->Array, &incx, &Beta, fe->Array, &incy );

     /* BLAS: tempvec = dispt */
     hysl_copy( &Tempvec->Rows, DispT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = dispt + A1W1*velt*DeltaT = tempvec + A1W1*DeltaT*velt */
     Alpha = GSSSS->A1W1*DeltaT;
     hysl_axpy( &Tempvec->Rows, &Alpha, VelT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = dispt + A1W1*acct*DeltaT + A2W2*DeltaT*DeltaT*acct = tempvec + A2W2*DeltaT*DeltaT*acct */
     Alpha = GSSSS->A2W2*DeltaT*DeltaT;
     hysl_axpy( &Tempvec->Rows, &Alpha, AccT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = dispt + A1W1*acct*DeltaT + A2W2*DeltaT*DeltaT*acct + A3W2*DeltaT*DeltaT*acctdt= tempvec + A3W3*DeltaT*DeltaT*acctdt */
     Alpha = GSSSS->A3W3*DeltaT*DeltaT;
     hysl_axpy( &Tempvec->Rows, &Alpha, AccTdT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = dispt + A1W1*acct*DeltaT + A2W2*DeltaT*DeltaT*acct + A3W2*DeltaT*DeltaT*(acctdt -
      * acct) = tempvec - A3W3*DeltaT*DeltaT*acct */
     Alpha = -GSSSS->A3W3*DeltaT*DeltaT;
     hysl_axpy( &Tempvec->Rows, &Alpha, AccT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: fe = -M*(acct *A6W1*(acctdt - acct)) - C*(velt + A4W1*acct*DeltaT + A5W2*DeltaT*(acctdt - acct) -
      * K*(dispt + A1W1*DeltaT*velt + A2W2*DeltaT*DeltaT*acct + A3W3*DeltaT*DeltaT*(acctdt - acct) =
      * fe + K*Tempvec */
     Alpha = -1.0; Beta = 1.0;
     hysl_spmv( &uplo, &fe->Rows, &Alpha, Stiff->Array, Tempvec->Array, &incx, &Beta, fe->Array, &incy );

     /* BLAS: fe = M*(acct *A6W1*(acctdt - acct)) + C*(velt + A4W1*acct*DeltaT + A5W2*DeltaT*(acctdt - acct) +
      * K*(dispt + A1W1*DeltaT*velt + A2W2*DeltaT*DeltaT*acct + A3W3*DeltaT*DeltaT*(acctdt - acct) +
      * W1*LoadTdT = fe + W1*LoadTdT */
     Alpha = GSSSS->W1;
     hysl_axpy( &fe->Rows, &Alpha, LoadTdT->Array, &incx, fe->Array, &incy );

     /* BLAS: fe = M*(acct *A6W1*(acctdt - acct)) + C*(velt + A4W1*acct*DeltaT + A5W2*DeltaT*(acctdt - acct) +
      * K*(dispt + A1W1*DeltaT*velt + A2W2*DeltaT*DeltaT*acct + A3W3*DeltaT*DeltaT*(acctdt - acct) +
      * W1*LoadTdT + (1-W1)*LoadT = fe + (1 - W1)*LoadT */
     Alpha = (1.0 - GSSSS->W1);
     hysl_axpy( &fe->Rows, &Alpha, LoadT->Array, &incx, fe->Array, &incy );

     /* BLAS: fe = M*(acct *A6W1*(acctdt - acct)) + C*(velt + A4W1*acct*DeltaT + A5W2*DeltaT*(acctdt - acct) +
      * K*(dispt + A1W1*DeltaT*velt + A2W2*DeltaT*DeltaT*acct + A3W3*DeltaT*DeltaT*(acctdt - acct) +
      * W1*LoadTdT + (1-W1)*LoadT + fc = fe + fc */
     Alpha = 1.0;
     hysl_axpy( &fe->Rows, &Alpha, fc->Array, &incx, fe->Array, &incy );

     /* Apply the PID */
     /* BLAS: fe = P*fe */
     Alpha = -PID->P;
     hysl_scal( &fe->Rows, &Alpha, fe->Array, &incx );
}
