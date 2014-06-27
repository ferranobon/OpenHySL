#include "GSSSS.h"
#include "MatrixVector.h"

#include "Definitions.h"

#if _MKL_
#include <mkl_blas.h>
#else
#include "Netlib.h"
#endif

void GSSSS_EffectiveForce_AForm( const MatrixVector_t *const Mass, const MatrixVector_t *const Damp,
				 const MatrixVector_t *const Stiff, const MatrixVector_t *const DispT,
				 const MatrixVector_t *const VelT, const MatrixVector_t *const AccT,
				 MatrixVector_t *const Tempvec, const TIntegration_GSSSS_t *const GSSSS,
				 const HYSL_FLOAT DeltaT, MatrixVector_t *const Eff_ForceT )
{

     int incx = 1, incy = 1;  /* Stride in the vectors */
     char uplo = 'L';         /* The lower part (upper part in C) will be used and the upper part (lower part
			       * in C) will strictly not be referenced */
     HYSL_FLOAT Alpha, Beta;  /* Constants for the BLAS routines */
     

    /* BLAS: tempvec = Disp */
     hysl_copy( &Tempvec->Rows, DispT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = Disp + A1W1*Vel*DeltaT = tempvec + A1W1*Vel*DeltaT */
     Alpha = GSSSS->A1W1*DeltaT;
     hysl_axpy( &Tempvec->Rows, &Alpha, VelT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = Disp + A1W1*Vel*DeltaT + (A2W2 - A3W3)*Acc*DeltaT^2 = tempvec + (A2W2 -A3W3)*a10*DeltaT^2 */
     Alpha = (GSSSS->A2W2 - GSSSS->A3W3)*DeltaT*DeltaT;
     hysl_axpy( &Tempvec->Rows, &Alpha, AccT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: Eff_Force = -Stiff*(Disp + A1W1*Vel*DeltaT + (A2W2 - A3W3)*Acc*DeltaT^2) = -Stiff*tempvec */
     Alpha = -1.0; Beta = 0.0;
     hysl_symv( &uplo, &Tempvec->Rows, &Alpha, Stiff->Array, &Tempvec->Rows, Tempvec->Array, &incx, &Beta,
		Eff_ForceT->Array, &incy );

     /* BLAS: tempvec = Vel */
     hysl_copy( &Tempvec->Rows, VelT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = Vel + (A4W1 -A5W2)*DeltaT*Acc = tempvec + (A4W1 - A5W2)*DeltaT*Acc */
     Alpha = (GSSSS->A4W1 - GSSSS->A5W2)*DeltaT;
     hysl_axpy( &Tempvec->Rows, &Alpha, AccT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: Eff_Force = -Stiff*(Disp + A1W1*Vel*DeltaT + (A2W2 - A3W3)*Acc*DeltaT^2) -C*(Vel + (A4W1 -A5W2)*DeltaT*Acc)= Eff_Force - Damp*tempvec */
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
				    const HYSL_FLOAT DeltaT, MatrixVector_t *const Eff_ForceT )
{

     int incx = 1, incy = 1;  /* Stride in the vectors */
     char uplo = 'L';         /* The lower part (upper part in C) will be used and the upper part (lower part
			       * in C) will strictly not be referenced */
     HYSL_FLOAT Alpha, Beta;  /* Constants for the BLAS routines */
     

    /* BLAS: tempvec = Disp */
     hysl_copy( &Tempvec->Rows, DispT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = Disp + A1W1*Vel*DeltaT = tempvec + A1W1*Vel*DeltaT */
     Alpha = GSSSS->A1W1*DeltaT;
     hysl_axpy( &Tempvec->Rows, &Alpha, VelT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = Disp + A1W1*Vel*DeltaT + (A2W2 - A3W3)*Acc*DeltaT^2 = tempvec + (A2W2 -A3W3)*a10*DeltaT^2 */
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
     /* BLAS: Eff_Force = -Stiff*(Disp + A1W1*Vel*DeltaT + (A2W2 - A3W3)*Acc*DeltaT^2) -C*(Vel + (A4W1 -A5W2)*DeltaT*Acc)= Eff_Force - Damp*tempvec */
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
				      const TIntegration_GSSSS_t *const GSSSS, const HYSL_FLOAT DeltaT,
				      MatrixVector_t *const DispTdT )
{
     int incx = 1, incy = 1;  /* Stride in the vectors */
     HYSL_FLOAT Alpha;        /* Constant for the BLAS routines */

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
				  const HYSL_FLOAT DeltaT, MatrixVector_t *const VelTdT )
{
     int incx = 1, incy= 1;  /* Stride in the vectors */
     HYSL_FLOAT Alpha;       /* Constant for the BLAS routines */

     /* BLAS: VelTdT = VelT */
     hysl_copy( &VelTdT->Rows, VelT->Array, &incx, VelTdT->Array, &incy );
     /* BLAS: VelTdT = VelT + (l4 - l5)*AccT*DeltaT = VelTdT + (l4 - l5)*AccT*DeltaT */
     Alpha = (GSSSS->l4 - GSSSS->l5)*DeltaT;
     hysl_axpy( &VelTdT->Rows, &Alpha, AccT->Array, &incx, VelTdT->Array, &incy );
     /* BLAS: VelTdT = VelT + (l4 - l5)*AccT*DeltaT + l5*AccTdT*DeltaT = VelTdT + l5*AccTdT*DeltaT */
     Alpha = GSSSS->l5*DeltaT;
     hysl_axpy( &VelTdT->Rows, &Alpha, AccTdT->Array, &incx, VelTdT->Array, &incy );
}
