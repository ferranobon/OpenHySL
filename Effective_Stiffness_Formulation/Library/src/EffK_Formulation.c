#include "EffK_Formulation.h"

void EffK_EffectiveForce( const MatrixVector_t *const Mass, const MatrixVector_t *const Damp, const MatrixVector_t *const DispT,
			  const MatrixVector_t *const VelT, const MatrixVector_t *const AccT, MatrixVector_t *const Tempvec,
			  const double a0, const double a1, const double a2, const double a3, const double a4,
			  const double a5, MatrixVector_t *const Eff_ForceT )
{

     static int incx = 1, incy = 1;  /* Stride in the vectors */
     static char uplo = 'L';         /* The lower part (upper part in C) will be used and the upper part (lower part in C)
				      *	will strictly not be referenced */
     static double Alpha, Beta;      /* Constants for the BLAS routines */

     /* BLAS: tempvec = Disp */
     dcopy_( &Tempvec->Rows, DispT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = a0*Disp = a0*tempvec */
     Alpha = a0;
     dscal_( &Tempvec->Rows, &Alpha, Tempvec->Array, &incx );
     /* BLAS: tempvec = a0*Disp + a2*Vel = tempvec + a2*Vel */
     Alpha = a2;
     daxpy_( &Tempvec->Rows, &Alpha, VelT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = a0*Disp + a2*Vel + a3*Acc = tempvec + a3*Acc */
     Alpha = a3;
     daxpy_( &Tempvec->Rows, &Alpha, AccT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: Eff_Force = Mass*(a0*Disp + a2*Vel + a3*Acc) = Mass*tempvec */
     Alpha = 1.0; Beta = 0.0;
     dsymv_( &uplo, &Tempvec->Rows, &Alpha, Mass->Array, &Tempvec->Rows, Tempvec->Array, &incx, &Beta,
	     Eff_ForceT->Array, &incy );

     /* BLAS: tempvec = Disp */
     dcopy_( &Tempvec->Rows, DispT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = a1*Disp = a0*tempvec */
     Alpha = a1;
     dscal_( &Tempvec->Rows, &Alpha, Tempvec->Array, &incx );
     /* BLAS: tempvec = a1*Disp + a4*Vel = tempvec + a4*Vel */
     Alpha = a4;
     daxpy_( &Tempvec->Rows, &Alpha, VelT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = a1*Disp + a4*Vel + a5*Acc = tempvec + a5*Acc */
     Alpha = a5;
     daxpy_( &Tempvec->Rows, &Alpha, AccT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: Eff_Force = Mass*(a0*Disp + a2*Vel + a3*Acc) + Damp*(a1*Disp + a4*Vel + a5*Acc) = Eff_Force + Damp*tempvec */
     Alpha = 1.0; Beta = 1.0;
     dsymv_( &uplo, &Tempvec->Rows, &Alpha, Damp->Array, &Tempvec->Rows, Tempvec->Array, &incx, &Beta,
	     Eff_ForceT->Array, &incy );
     
}

void EffK_ComputeAcceleration( const MatrixVector_t *const DispTdT, const MatrixVector_t *const DispT,
			       const MatrixVector_t *const VelT, const MatrixVector_t *const AccT,
			       const double a0, const double a2, const double a3, MatrixVector_t *const AccTdT )
{
     static int incx = 1, incy = 1;  /* Stride in the vectors */
     static double Alpha;            /* Constant for the BLAS routines */

     /* BLAS: AccTdT = DispTdT */
     dcopy_( &AccTdT->Rows, DispTdT->Array, &incx, AccTdT->Array, &incy ); 
     /* BLAS: AccTdT = DispTdT - DispT = AccTdT - DispT */
     Alpha = -1.0;
     daxpy_( &AccTdT->Rows, &Alpha, DispT->Array, &incx, AccTdT->Array, &incy );
     /* BLAS: AccTdT = a0*(DispTdT - DispT) = a0*AccTdT */
     Alpha = a0;
     dscal_( &AccTdT->Rows, &Alpha, AccTdT->Array, &incx );
     /* BLAS: AccTdT = a0*(DispTdT - DispT) - a2*VelT = AccTdT - a2*VelT */
     Alpha = -a2;
     daxpy_( &AccTdT->Rows, &Alpha, VelT->Array, &incx, AccTdT->Array, &incy );
     /* BLAS: AccTdT = a0*(DispTdT - DispT) - a2*VelT - a3*AccT = AccTdT - a3*AccT */
     Alpha = -a3;
     daxpy_( &AccTdT->Rows, &Alpha, AccT->Array, &incx, AccTdT->Array, &incy );
}

void EffK_ComputeVelocity( const MatrixVector_t *const VelT, const MatrixVector_t *const AccT, const MatrixVector_t *const AccTdT,
			   const double a6, const double a7, MatrixVector_t *const VelTdT )
{
     static int incx = 1, incy= 1;  /* Stride in the vectors */
     static double Alpha;           /* Constant for the BLAS routines */

     /* BLAS: VelTdT = VelT */
     dcopy_( &VelTdT->Rows, VelT->Array, &incx, VelTdT->Array, &incy );
     /* BLAS: VelTdT = VelT + a6*AccT = VelTdT + a6*AccT */
     Alpha = a6;
     daxpy_( &VelTdT->Rows, &Alpha, AccT->Array, &incx, VelTdT->Array, &incy );
     /* BLAS: VelTdT = VelT + a6*AccT + a7*AccTdT = VelTdT + a7*AccTdT */
     Alpha = a7;
     daxpy_( &VelTdT->Rows, &Alpha, AccTdT->Array, &incx, VelTdT->Array, &incy );
}
