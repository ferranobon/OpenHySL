#include "EffK_Formulation.h"
#include "MatrixVector_MPI.h"

#if _MKL_
#include <mkl_pblas.h>
#else
#include "Netlib.h"
#endif

void EffK_EffectiveForce_MPI( const PMatrixVector_t *const Mass, const PMatrixVector_t *const Damp,
			      const PMatrixVector_t *const Disp, const PMatrixVector_t *const Vel,
			      const PMatrixVector_t *const Acc, const PMatrixVector_t *const Tempvec,
			      const double a0, const double a1, const double a2, const double a3,
			      const double a4, const double a5, PMatrixVector_t *const Eff_Force )
{
     int ione = 1;
     int incx = 1, incy = 1;  /* Stride in the vectors */
     char uplo = 'L';         /* The lower part (upper part in C) will be used and the upper part (lower part
			       * in C) will strictly not be referenced */
     double Alpha, Beta;      /* Constants for the BLAS routines */

     /* PBLAS: tempvec = Disp */
     pdcopy_( &Tempvec->GlobalSize.Row, Disp->Array, &ione, &ione, Disp->Desc, &ione, Tempvec->Array, &ione,
	      &ione, Tempvec->Desc, &ione );
     /* PBLAS: tempvec = a0*Disp = a0*tempvec */
     Alpha = a0;
     pdscal_( &Tempvec->GlobalSize.Row, &Alpha, Tempvec->Array, &ione, &ione, Tempvec->Desc, &incx );
     /* BLAS: tempvec = a0*Disp + a2*Vel = tempvec + a2*Vel */
     Alpha = a2;
     pdaxpy_( &Tempvec->GlobalSize.Row, &Alpha, Vel->Array, &ione, &ione, Vel->Desc, &incx, Tempvec->Array,
	      &ione, &ione, Tempvec->Desc, &incy );
     /* PBLAS: tempvec = a0*Disp + a2*Vel + a3*Acc = tempvec + a3*Acc */
     Alpha = a3;
     pdaxpy_( &Tempvec->GlobalSize.Row, &Alpha, Acc->Array, &ione, &ione, Acc->Desc, &incx, Tempvec->Array,
	      &ione, &ione, Tempvec->Desc, &incy );
     /* PBLAS: Eff_Force = Mass*(a0*Disp + a2*Vel + a3*Acc) = Mass*tempvec */
     Alpha = 1.0; Beta = 0.0;
     pdsymv_( &uplo, &Tempvec->GlobalSize.Row, &Alpha, Mass->Array, &ione, &ione, Mass->Desc, Tempvec->Array,
	      &ione, &ione, Tempvec->Desc, &incx, &Beta, Eff_Force->Array, &ione, &ione, Eff_Force->Desc,
	      &incy );

     /* PBLAS: tempvec = Disp */
     pdcopy_( &Tempvec->GlobalSize.Row, Disp->Array, &ione, &ione, Disp->Desc, &ione, Tempvec->Array, &ione,
	      &ione, Tempvec->Desc, &ione );
     /* PBLAS: tempvec = a1*Disp = a0*tempvec */
     Alpha = a1;
     pdscal_( &Tempvec->GlobalSize.Row, &Alpha, Tempvec->Array, &ione, &ione, Tempvec->Desc, &incx );
     /* PBLAS: tempvec = a1*Disp + a4*Vel = tempvec + a4*Vel */
     Alpha = a4;
     pdaxpy_( &Tempvec->GlobalSize.Row, &Alpha, Vel->Array, &ione, &ione, Vel->Desc, &incx, Tempvec->Array,
	      &ione, &ione, Tempvec->Desc, &incy );
     /* PBLAS: tempvec = a1*Disp + a4*Vel + a5*Acc = tempvec + a5*Acc */
     Alpha = a5;
     pdaxpy_( &Tempvec->GlobalSize.Row, &Alpha, Acc->Array, &ione, &ione, Acc->Desc, &incx, Tempvec->Array,
	      &ione, &ione, Tempvec->Desc, &incy );
     /* PBLAS: Eff_Force = Mass*(a0*Disp + a2*Vel + a3*Acc) + Damp*(a1*Disp + a4*Vel + a5*Acc) = Eff_Force +
      * Damp*tempvec */
     Alpha = 1.0; Beta = 1.0;
     pdsymv_( &uplo, &Tempvec->GlobalSize.Row, &Alpha, Damp->Array, &ione, &ione, Damp->Desc, Tempvec->Array,
	      &ione, &ione, Tempvec->Desc, &incx, &Beta, Eff_Force->Array, &ione, &ione, Eff_Force->Desc,
	      &incy );
}

void EffK_ComputeAcceleration_MPI( const PMatrixVector_t *const DispTdT, const PMatrixVector_t *const DispT,
				   const PMatrixVector_t *const VelT, const PMatrixVector_t *const AccT,
				   const double a0, const double a2, const double a3,
				   PMatrixVector_t *const AccTdT )
{

     int ione = 1;
     int incx = 1, incy = 1;  /* Stride in the vectors */
     double Alpha;            /* Constant for the BLAS routines */

     /* PBLAS: AccTdT = DispTdT */
     pdcopy_( &AccTdT->GlobalSize.Row, DispTdT->Array, &ione, &ione, DispTdT->Desc, &ione, AccTdT->Array,
	      &ione, &ione, AccTdT->Desc, &ione );
     /* PBLAS: AccTdT = DispTdT - DispT = AccTdT - DispT */
     Alpha = -1.0;
     pdaxpy_( &AccTdT->GlobalSize.Row, &Alpha, DispT->Array, &ione, &ione, DispT->Desc, &incx, AccTdT->Array,
	      &ione, &ione, AccTdT->Desc, &incy );
     /* PBLAS: AccTdT = a0*(DispTdT - DispT) = a0*AccTdT */
     Alpha = a0;
     pdscal_( &AccTdT->GlobalSize.Row, &Alpha, AccTdT->Array, &ione, &ione, AccTdT->Desc, &incx );
     /* PBLAS: AccTdT = a0*(DispTdT - DispT) - a2*VelT = AccTdT - a2*VelT */
     Alpha = -a2;
     pdaxpy_( &AccTdT->GlobalSize.Row, &Alpha, VelT->Array, &ione, &ione, VelT->Desc, &incx, AccTdT->Array,
	      &ione, &ione, AccTdT->Desc, &incy );
     /* PBLAS: AccTdT = a0*(DispTdT - DispT) - a2*VelT - a3*AccT = AccTdT - a3*AccT */
     Alpha = -a3;
     pdaxpy_( &AccTdT->GlobalSize.Row, &Alpha, AccT->Array, &ione, &ione, AccT->Desc, &incx, AccTdT->Array,
	      &ione, &ione, AccTdT->Desc, &incy );
}

void EffK_ComputeVelocity_MPI( const PMatrixVector_t *const VelT, const PMatrixVector_t *const AccT,
			       const PMatrixVector_t *const AccTdT, const double a6, const double a7,
			       PMatrixVector_t *const VelTdT )
{
     int ione = 1;
     int incx = 1, incy= 1;  /* Stride in the vectors */
     double Alpha;           /* Constant for the BLAS routines */

     /* PBLAS: VelTdT = VelT */
     pdcopy_( &VelTdT->GlobalSize.Row, VelT->Array, &ione, &ione, VelT->Desc, &incx, VelTdT->Array, &ione,
	      &ione, VelTdT->Desc, &incy );
     /* BLAS: VelTdT = VelT + a6*AccT = VelTdT + a6*AccT */
     Alpha = a6;
     pdaxpy_( &VelTdT->GlobalSize.Row, &Alpha, AccT->Array, &ione, &ione, AccT->Desc, &incx, VelTdT->Array,
	      &ione, &ione, VelTdT->Desc, &incy );
     /* BLAS: VelTdT = VelT + a6*AccT + a7*AccTdT = VelTdT + a7*AccTdT */
     Alpha = a7;
     pdaxpy_( &VelTdT->GlobalSize.Row, &Alpha, AccTdT->Array, &ione, &ione, AccTdT->Desc, &incx,
	      VelTdT->Array, &ione, &ione, VelTdT->Desc, &incy );
}