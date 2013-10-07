#include "EffM_Formulation.h"
#include "MatrixVector_MPI.h"

#if _MKL_
#include <mkl_pblas.h>
#include <Cblacs.h>
#else
#include "Netlib.h"
#endif

void EffM_EffectiveForce_MPI( PMatrixVector_t *const Stiff, PMatrixVector_t *const Damp,
			      PMatrixVector_t *const DispT, PMatrixVector_t *const VelT,
			      PMatrixVector_t *const AccT, PMatrixVector_t *const Tempvec,
			      const double a6, const double a9, const double a10,
			      PMatrixVector_t *const Eff_ForceT )
{
     int ione = 1;
     int incx = 1, incy = 1;  /* Stride in the vectors */
     char uplo = 'L';         /* The lower part (upper part in C) will be used and the upper part (lower part
			       * in C) will strictly not be referenced */
     double Alpha, Beta;      /* Constants for the BLAS routines */

     /* PBLAS: tempvec = Disp */
     pdcopy_( &Tempvec->GlobalSize.Row, DispT->Array, &ione, &ione, DispT->Desc, &ione, Tempvec->Array, &ione,
	      &ione, Tempvec->Desc, &ione );
     /* PBLAS: tempvec = tempvec = Disp + a9*Vel = tempvec + a9*Vel */
     Alpha = a9;
     pdaxpy_( &Tempvec->GlobalSize.Row, &Alpha, VelT->Array, &ione, &ione, VelT->Desc, &incx, Tempvec->Array,
	      &ione, &ione, Tempvec->Desc, &incy );
     /* PBLAS: tempvec = Disp + a9*Vel + a10*Acc = tempvec + a10*Acc */
     Alpha = a10;
     pdaxpy_( &Tempvec->GlobalSize.Row, &Alpha, AccT->Array, &ione, &ione, AccT->Desc, &incx, Tempvec->Array,
	      &ione, &ione, Tempvec->Desc, &incy );
     /* PBLAS: Eff_Force = -Stiff*(Disp + a2*Vel + a10*Acc) = -Stiff*tempvec */
     Alpha = -1.0; Beta = 0.0;
     pdsymv_( &uplo, &Tempvec->GlobalSize.Row, &Alpha, Stiff->Array, &ione, &ione, Stiff->Desc, Tempvec->Array,
	      &ione, &ione, Tempvec->Desc, &incx, &Beta, Eff_ForceT->Array, &ione, &ione, Eff_ForceT->Desc,
	      &incy );

     /* PBLAS: tempvec = Vel */
     pdcopy_( &Tempvec->GlobalSize.Row, VelT->Array, &ione, &ione, VelT->Desc, &ione, Tempvec->Array, &ione,
	      &ione, Tempvec->Desc, &ione );
     /* PBLAS: tempvec = Vel + a6*Acc = tempvec + a6*Acc */
     Alpha = a6;
     pdaxpy_( &Tempvec->GlobalSize.Row, &Alpha, AccT->Array, &ione, &ione, AccT->Desc, &incx, Tempvec->Array,
	      &ione, &ione, Tempvec->Desc, &incy );
     /* PBLAS: Eff_Force = -Stiff*(Disp + a9*Vel + a10*Acc) - Damp*(Vel + 6*Acc) = Eff_Force - Damp*tempvec */
     Alpha = 1.0; Beta = 1.0;
     pdsymv_( &uplo, &Tempvec->GlobalSize.Row, &Alpha, Damp->Array, &ione, &ione, Damp->Desc, Tempvec->Array,
	      &ione, &ione, Tempvec->Desc, &incx, &Beta, Eff_ForceT->Array, &ione, &ione, Eff_ForceT->Desc,
	      &incy );
}

void EffM_ComputeDisplacement_MPI( PMatrixVector_t *const DispT, PMatrixVector_t *const VelT,
				   PMatrixVector_t *const AccT, PMatrixVector_t *const AccTdT,
				   const double a8, const double a9, const double a10,
				   PMatrixVector_t *const DispTdT )
{

     int ione = 1;
     int incx = 1, incy = 1;  /* Stride in the vectors */
     double Alpha;            /* Constant for the BLAS routines */

     /* PBLAS: DispTdT = DispT */
     pdcopy_( &DispTdT->GlobalSize.Row, DispT->Array, &ione, &ione, DispT->Desc, &ione, DispTdT->Array,
	      &ione, &ione, DispTdT->Desc, &ione );
     /* PBLAS: DispTdT = DispT + a9*VelT = DispTdT + a9*VelT */
     Alpha = a9;
     pdaxpy_( &DispTdT->GlobalSize.Row, &Alpha, VelT->Array, &ione, &ione, VelT->Desc, &incx, DispTdT->Array,
	      &ione, &ione, DispTdT->Desc, &incy );
     /* PBLAS: DispTdT = DispT + a9*VelT + a10*AccT = DispTdT + a10*AccT */
     Alpha = a10;
     pdaxpy_( &DispTdT->GlobalSize.Row, &Alpha, AccT->Array, &ione, &ione, AccT->Desc, &incx, DispTdT->Array,
	      &ione, &ione, AccTdT->Desc, &incy );
     /* PBLAS: DispTdT = DispT + a9*VelT + a10*AccT + a8*AccTdT = DispTdT + a8*AccT */
     Alpha = a8;
     pdaxpy_( &DispTdT->GlobalSize.Row, &Alpha, AccTdT->Array, &ione, &ione, AccTdT->Desc, &incx, DispTdT->Array,
	      &ione, &ione, DispTdT->Desc, &incy );
}

void EffM_ComputeVelocity_MPI( PMatrixVector_t *const VelT, PMatrixVector_t *const AccT,
			       PMatrixVector_t *const AccTdT, const double a6, const double a7,
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