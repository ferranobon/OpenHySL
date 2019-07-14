#include "New_State.h"
#include "MatrixVector_MPI.h"
#include "Definitions.h"

#if _MKL_
#include <mkl_pblas.h>
#else
#include "Netlib.h"
#endif

void Compute_NewState_MPI( PMatrixVector_t *const IGain, PMatrixVector_t *const Eff_ForceT,
			   PMatrixVector_t *const In_LoadT, PMatrixVector_t *const Err_ForceT,
			   PMatrixVector_t *const Tempvec, PMatrixVector_t *const VecTdT_0 )
{
     int ione;
     int incx, incy;                   /* Stride in the vectors */
     HYSL_FLOAT Alpha = 1.0, Beta = 0.0;   /* Constants for the BLAS routines */
     char uplo = 'L';                  /* The lower part (upper part in C) will be used and the upper part
					* (lower part in C) will strictly not be referenced */

     ione = 1;
     incx = 1; incy = 1;
     Alpha = 1.0; Beta = 0.0;
     uplo = 'L';

     /* PBLAS: tempvec = Eff_Force */
     hysl_pcopy( &Tempvec->GlobalSize.Row, Eff_ForceT->Array, &ione, &ione, Eff_ForceT->Desc, &incx,
	      Tempvec->Array, &ione, &ione, Tempvec->Desc, &incy );

     /* PBLAS: tempvec = Eff_Force + LoadTdT = tempvec + LoadTdT */
     hysl_paxpy( &Tempvec->GlobalSize.Row, &Alpha, In_LoadT->Array, &ione, &ione, In_LoadT->Desc, &incx,
	      Tempvec->Array, &ione, &ione, Tempvec->Desc, &incy );

     /* PBLAS: tempvec = Eff_Force + LoadTdT - Err_Force = tempvec - Err_Force. */
     Alpha = -1.0;
     hysl_paxpy( &Tempvec->GlobalSize.Row, &Alpha, Err_ForceT->Array, &ione, &ione, Err_ForceT->Desc, &incx,
	      Tempvec->Array,  &ione, &ione, Tempvec->Desc, &incy );

     /* PBLAS: VecTdT_0 = IGain*(Eff_Force + LoadTdT - Err_Force) = IGain*Tempvec */
     Alpha = 1.0;
     hysl_psymv( &uplo, &Tempvec->GlobalSize.Row, &Alpha, IGain->Array, &ione, &ione, IGain->Desc, 
	      Tempvec->Array, &ione, &ione, Tempvec->Desc, &incx,
	      &Beta, VecTdT_0->Array, &ione, &ione, VecTdT_0->Desc, &incy );
}
