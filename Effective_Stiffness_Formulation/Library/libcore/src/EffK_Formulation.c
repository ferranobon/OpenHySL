#include "EffK_Formulation.h"
#include "MatrixVector.h"

#include "Definitions.h"

#if _MKL_
#include <mkl_blas.h>
#else
#include "Netlib.h"
#endif

void EffK_EffectiveForce( const MatrixVector_t *const Mass, const MatrixVector_t *const Damp,
			  const MatrixVector_t *const DispT, const MatrixVector_t *const VelT,
			  const MatrixVector_t *const AccT, MatrixVector_t *const Tempvec,
			  const HYSL_FLOAT a0, const HYSL_FLOAT a1, const HYSL_FLOAT a2, const HYSL_FLOAT a3, const HYSL_FLOAT a4,
			  const HYSL_FLOAT a5, MatrixVector_t *const Eff_ForceT )
{

     int incx = 1, incy = 1;  /* Stride in the vectors */
     char uplo = 'L';         /* The lower part (upper part in C) will be used and the upper part (lower part
			       * in C) will strictly not be referenced */
     HYSL_FLOAT Alpha, Beta;      /* Constants for the BLAS routines */

     /* BLAS: tempvec = Disp */
     hysl_copy( &Tempvec->Rows, DispT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = a0*Disp = a0*tempvec */
     Alpha = a0;
     hysl_scal( &Tempvec->Rows, &Alpha, Tempvec->Array, &incx );
     /* BLAS: tempvec = a0*Disp + a2*Vel = tempvec + a2*Vel */
     Alpha = a2;
     hysl_axpy( &Tempvec->Rows, &Alpha, VelT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = a0*Disp + a2*Vel + a3*Acc = tempvec + a3*Acc */
     Alpha = a3;
     hysl_axpy( &Tempvec->Rows, &Alpha, AccT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: Eff_Force = Mass*(a0*Disp + a2*Vel + a3*Acc) = Mass*tempvec */
     Alpha = 1.0; Beta = 0.0;
     hysl_symv( &uplo, &Tempvec->Rows, &Alpha, Mass->Array, &Tempvec->Rows, Tempvec->Array, &incx, &Beta,
	    Eff_ForceT->Array, &incy );

     /* BLAS: tempvec = Disp */
     hysl_copy( &Tempvec->Rows, DispT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = a1*Disp = a0*tempvec */
     Alpha = a1;
     hysl_scal( &Tempvec->Rows, &Alpha, Tempvec->Array, &incx );
     /* BLAS: tempvec = a1*Disp + a4*Vel = tempvec + a4*Vel */
     Alpha = a4;
     hysl_axpy( &Tempvec->Rows, &Alpha, VelT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = a1*Disp + a4*Vel + a5*Acc = tempvec + a5*Acc */
     Alpha = a5;
     hysl_axpy( &Tempvec->Rows, &Alpha, AccT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: Eff_Force = Mass*(a0*Disp + a2*Vel + a3*Acc) + Damp*(a1*Disp + a4*Vel + a5*Acc) = Eff_Force +
      * Damp*tempvec */
     Alpha = 1.0; Beta = 1.0;
     hysl_symv( &uplo, &Tempvec->Rows, &Alpha, Damp->Array, &Tempvec->Rows, Tempvec->Array, &incx, &Beta,
	    Eff_ForceT->Array, &incy );     
}

void EffK_EffectiveForce_PS( const MatrixVector_t *const Mass, const MatrixVector_t *const Damp,
			     const MatrixVector_t *const DispT, const MatrixVector_t *const VelT,
			     const MatrixVector_t *const AccT, MatrixVector_t *const Tempvec,
			     const HYSL_FLOAT a0, const HYSL_FLOAT a1, const HYSL_FLOAT a2, const HYSL_FLOAT a3,
			     const HYSL_FLOAT a4, const HYSL_FLOAT a5, MatrixVector_t *const Eff_ForceT )
{

     int incx = 1, incy = 1;  /* Stride in the vectors */
     char uplo = 'L';         /* The lower part (upper part in C) will be used and the upper part (lower part
			       * in C) will strictly not be referenced */
     HYSL_FLOAT Alpha, Beta;      /* Constants for the BLAS routines */

     /* BLAS: tempvec = Disp */
     hysl_copy( &Tempvec->Rows, DispT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = a0*Disp = a0*tempvec */
     Alpha = a0;
     hysl_scal( &Tempvec->Rows, &Alpha, Tempvec->Array, &incx );
     /* BLAS: tempvec = a0*Disp + a2*Vel = tempvec + a2*Vel */
     Alpha = a2;
     hysl_axpy( &Tempvec->Rows, &Alpha, VelT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = a0*Disp + a2*Vel + a3*Acc = tempvec + a3*Acc */
     Alpha = a3;
     hysl_axpy( &Tempvec->Rows, &Alpha, AccT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: Eff_Force = Mass*(a0*Disp + a2*Vel + a3*Acc) = Mass*tempvec */
     Alpha = 1.0; Beta = 0.0;
     hysl_spmv( &uplo, &Tempvec->Rows, &Alpha, Mass->Array, Tempvec->Array, &incx, &Beta,
	     Eff_ForceT->Array, &incy );

     /* BLAS: tempvec = Disp */
     hysl_copy( &Tempvec->Rows, DispT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = a1*Disp = a0*tempvec */
     Alpha = a1;
     hysl_scal( &Tempvec->Rows, &Alpha, Tempvec->Array, &incx );
     /* BLAS: tempvec = a1*Disp + a4*Vel = tempvec + a4*Vel */
     Alpha = a4;
     hysl_axpy( &Tempvec->Rows, &Alpha, VelT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = a1*Disp + a4*Vel + a5*Acc = tempvec + a5*Acc */
     Alpha = a5;
     hysl_axpy( &Tempvec->Rows, &Alpha, AccT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: Eff_Force = Mass*(a0*Disp + a2*Vel + a3*Acc) + Damp*(a1*Disp + a4*Vel + a5*Acc) = Eff_Force +
      * Damp*tempvec */
     Alpha = 1.0; Beta = 1.0;
     hysl_spmv( &uplo, &Tempvec->Rows, &Alpha, Damp->Array, Tempvec->Array, &incx, &Beta,
	    Eff_ForceT->Array, &incy );
     
}

void EffK_ComputeAcceleration( const MatrixVector_t *const DispTdT, const MatrixVector_t *const DispT,
			       const MatrixVector_t *const VelT, const MatrixVector_t *const AccT,
			       const HYSL_FLOAT a0, const HYSL_FLOAT a2, const HYSL_FLOAT a3, MatrixVector_t *const AccTdT )
{
     int incx = 1, incy = 1;  /* Stride in the vectors */
     HYSL_FLOAT Alpha;            /* Constant for the BLAS routines */

     /* BLAS: AccTdT = DispTdT */
     hysl_copy( &AccTdT->Rows, DispTdT->Array, &incx, AccTdT->Array, &incy ); 
     /* BLAS: AccTdT = DispTdT - DispT = AccTdT - DispT */
     Alpha = -1.0;
     hysl_axpy( &AccTdT->Rows, &Alpha, DispT->Array, &incx, AccTdT->Array, &incy );
     /* BLAS: AccTdT = a0*(DispTdT - DispT) = a0*AccTdT */
     Alpha = a0;
     hysl_scal( &AccTdT->Rows, &Alpha, AccTdT->Array, &incx );
     /* BLAS: AccTdT = a0*(DispTdT - DispT) - a2*VelT = AccTdT - a2*VelT */
     Alpha = -a2;
     hysl_axpy( &AccTdT->Rows, &Alpha, VelT->Array, &incx, AccTdT->Array, &incy );
     /* BLAS: AccTdT = a0*(DispTdT - DispT) - a2*VelT - a3*AccT = AccTdT - a3*AccT */
     Alpha = -a3;
     hysl_axpy( &AccTdT->Rows, &Alpha, AccT->Array, &incx, AccTdT->Array, &incy );
}

void EffK_ComputeVelocity( const MatrixVector_t *const VelT, const MatrixVector_t *const AccT,
			   const MatrixVector_t *const AccTdT, const HYSL_FLOAT a6, const HYSL_FLOAT a7,
			   MatrixVector_t *const VelTdT )
{
     int incx = 1, incy= 1;  /* Stride in the vectors */
     HYSL_FLOAT Alpha;           /* Constant for the BLAS routines */

     /* BLAS: VelTdT = VelT */
     hysl_copy( &VelTdT->Rows, VelT->Array, &incx, VelTdT->Array, &incy );
     /* BLAS: VelTdT = VelT + a6*AccT = VelTdT + a6*AccT */
     Alpha = a6;
     hysl_axpy( &VelTdT->Rows, &Alpha, AccT->Array, &incx, VelTdT->Array, &incy );
     /* BLAS: VelTdT = VelT + a6*AccT + a7*AccTdT = VelTdT + a7*AccTdT */
     Alpha = a7;
     hysl_axpy( &VelTdT->Rows, &Alpha, AccTdT->Array, &incx, VelTdT->Array, &incy );
}
