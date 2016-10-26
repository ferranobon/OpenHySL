#include "EffK_Formulation.h"
#include "MatrixVector_MPI.h"

#include "Definitions.h"

#if _MKL_
#include <mkl_pblas.h>
#include <Cblacs.h>
#else
#include "Netlib.h"
#endif

void EffK_EffectiveForce_MPI( PMatrixVector_t *const Mass, PMatrixVector_t *const Damp,
			      PMatrixVector_t *const DispT, PMatrixVector_t *const VelT,
			      PMatrixVector_t *const AccT, PMatrixVector_t *const Tempvec,
			      const HYSL_FLOAT a0, const HYSL_FLOAT a1, const HYSL_FLOAT a2, const HYSL_FLOAT a3,
			      const HYSL_FLOAT a4, const HYSL_FLOAT a5, PMatrixVector_t *const Eff_ForceT )
{
     int ione = 1;
     int incx = 1, incy = 1;  /* Stride in the vectors */
     char uplo = 'L';         /* The lower part (upper part in C) will be used and the upper part (lower part
			       * in C) will strictly not be referenced */
     HYSL_FLOAT Alpha, Beta;      /* Constants for the BLAS routines */

     /* PBLAS: tempvec = Disp */
     hysl_pcopy( &Tempvec->GlobalSize.Row, DispT->Array, &ione, &ione, DispT->Desc, &ione, Tempvec->Array, &ione,
	      &ione, Tempvec->Desc, &ione );
     /* PBLAS: tempvec = a0*Disp = a0*tempvec */
     Alpha = a0;
     hysl_pscal( &Tempvec->GlobalSize.Row, &Alpha, Tempvec->Array, &ione, &ione, Tempvec->Desc, &incx );
     /* BLAS: tempvec = a0*Disp + a2*Vel = tempvec + a2*Vel */
     Alpha = a2;
     hysl_paxpy( &Tempvec->GlobalSize.Row, &Alpha, VelT->Array, &ione, &ione, VelT->Desc, &incx, Tempvec->Array,
	      &ione, &ione, Tempvec->Desc, &incy );
     /* PBLAS: tempvec = a0*Disp + a2*Vel + a3*Acc = tempvec + a3*Acc */
     Alpha = a3;
     hysl_paxpy( &Tempvec->GlobalSize.Row, &Alpha, AccT->Array, &ione, &ione, AccT->Desc, &incx, Tempvec->Array,
	      &ione, &ione, Tempvec->Desc, &incy );
     /* PBLAS: Eff_Force = Mass*(a0*Disp + a2*Vel + a3*Acc) = Mass*tempvec */
     Alpha = 1.0; Beta = 0.0;
     hysl_psymv( &uplo, &Tempvec->GlobalSize.Row, &Alpha, Mass->Array, &ione, &ione, Mass->Desc, Tempvec->Array,
	      &ione, &ione, Tempvec->Desc, &incx, &Beta, Eff_ForceT->Array, &ione, &ione, Eff_ForceT->Desc,
	      &incy );

     /* PBLAS: tempvec = Disp */
     hysl_pcopy( &Tempvec->GlobalSize.Row, DispT->Array, &ione, &ione, DispT->Desc, &ione, Tempvec->Array, &ione,
	      &ione, Tempvec->Desc, &ione );
     /* PBLAS: tempvec = a1*Disp = a0*tempvec */
     Alpha = a1;
     hysl_pscal( &Tempvec->GlobalSize.Row, &Alpha, Tempvec->Array, &ione, &ione, Tempvec->Desc, &incx );
     /* PBLAS: tempvec = a1*Disp + a4*Vel = tempvec + a4*Vel */
     Alpha = a4;
     hysl_paxpy( &Tempvec->GlobalSize.Row, &Alpha, VelT->Array, &ione, &ione, VelT->Desc, &incx, Tempvec->Array,
	      &ione, &ione, Tempvec->Desc, &incy );
     /* PBLAS: tempvec = a1*Disp + a4*Vel + a5*Acc = tempvec + a5*Acc */
     Alpha = a5;
     hysl_paxpy( &Tempvec->GlobalSize.Row, &Alpha, AccT->Array, &ione, &ione, AccT->Desc, &incx, Tempvec->Array,
	      &ione, &ione, Tempvec->Desc, &incy );
     /* PBLAS: Eff_Force = Mass*(a0*Disp + a2*Vel + a3*Acc) + Damp*(a1*Disp + a4*Vel + a5*Acc) = Eff_Force +
      * Damp*tempvec */
     Alpha = 1.0; Beta = 1.0;
     hysl_psymv( &uplo, &Tempvec->GlobalSize.Row, &Alpha, Damp->Array, &ione, &ione, Damp->Desc, Tempvec->Array,
	      &ione, &ione, Tempvec->Desc, &incx, &Beta, Eff_ForceT->Array, &ione, &ione, Eff_ForceT->Desc,
	      &incy );
}

void EffK_ComputeAcceleration_MPI( PMatrixVector_t *const DispTdT, PMatrixVector_t *const DispT,
				   PMatrixVector_t *const VelT, PMatrixVector_t *const AccT,
				   const HYSL_FLOAT a0, const HYSL_FLOAT a2, const HYSL_FLOAT a3,
				   PMatrixVector_t *const AccTdT )
{

     int ione = 1;
     int incx = 1, incy = 1;  /* Stride in the vectors */
     HYSL_FLOAT Alpha;            /* Constant for the BLAS routines */

     /* PBLAS: AccTdT = DispTdT */
     hysl_pcopy( &AccTdT->GlobalSize.Row, DispTdT->Array, &ione, &ione, DispTdT->Desc, &ione, AccTdT->Array,
	      &ione, &ione, AccTdT->Desc, &ione );
     /* PBLAS: AccTdT = DispTdT - DispT = AccTdT - DispT */
     Alpha = -1.0;
     hysl_paxpy( &AccTdT->GlobalSize.Row, &Alpha, DispT->Array, &ione, &ione, DispT->Desc, &incx, AccTdT->Array,
	      &ione, &ione, AccTdT->Desc, &incy );
     /* PBLAS: AccTdT = a0*(DispTdT - DispT) = a0*AccTdT */
     Alpha = a0;
     hysl_pscal( &AccTdT->GlobalSize.Row, &Alpha, AccTdT->Array, &ione, &ione, AccTdT->Desc, &incx );
     /* PBLAS: AccTdT = a0*(DispTdT - DispT) - a2*VelT = AccTdT - a2*VelT */
     Alpha = -a2;
     hysl_paxpy( &AccTdT->GlobalSize.Row, &Alpha, VelT->Array, &ione, &ione, VelT->Desc, &incx, AccTdT->Array,
	      &ione, &ione, AccTdT->Desc, &incy );
     /* PBLAS: AccTdT = a0*(DispTdT - DispT) - a2*VelT - a3*AccT = AccTdT - a3*AccT */
     Alpha = -a3;
     hysl_paxpy( &AccTdT->GlobalSize.Row, &Alpha, AccT->Array, &ione, &ione, AccT->Desc, &incx, AccTdT->Array,
	      &ione, &ione, AccTdT->Desc, &incy );
}

void EffK_ComputeVelocity_MPI( PMatrixVector_t *const VelT, PMatrixVector_t *const AccT,
			       PMatrixVector_t *const AccTdT, const HYSL_FLOAT a6, const HYSL_FLOAT a7,
			       PMatrixVector_t *const VelTdT )
{
     int ione = 1;
     int incx = 1, incy= 1;  /* Stride in the vectors */
     HYSL_FLOAT Alpha;           /* Constant for the BLAS routines */

     /* PBLAS: VelTdT = VelT */
     hysl_pcopy( &VelTdT->GlobalSize.Row, VelT->Array, &ione, &ione, VelT->Desc, &incx, VelTdT->Array, &ione,
	      &ione, VelTdT->Desc, &incy );
     /* BLAS: VelTdT = VelT + a6*AccT = VelTdT + a6*AccT */
     Alpha = a6;
     hysl_paxpy( &VelTdT->GlobalSize.Row, &Alpha, AccT->Array, &ione, &ione, AccT->Desc, &incx, VelTdT->Array,
	      &ione, &ione, VelTdT->Desc, &incy );
     /* BLAS: VelTdT = VelT + a6*AccT + a7*AccTdT = VelTdT + a7*AccTdT */
     Alpha = a7;
     hysl_paxpy( &VelTdT->GlobalSize.Row, &Alpha, AccTdT->Array, &ione, &ione, AccTdT->Desc, &incx,
	      VelTdT->Array, &ione, &ione, VelTdT->Desc, &incy );
}
