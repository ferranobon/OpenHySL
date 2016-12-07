#include "Error_Compensation.h"
#include "MatrixVector_MPI.h"
#include "Definitions.h"

#if _MKL_
#include <mkl_pblas.h>
#else
#include "Netlib.h"
#endif

void ErrorForce_PID_MPI( PMatrixVector_t *const Mass, PMatrixVector_t *const Damp,
			 PMatrixVector_t *Stiff, PMatrixVector_t *const AccTdT,
			 PMatrixVector_t *const VelTdT, PMatrixVector_t *const DispTdT,
			 PMatrixVector_t *const fc, PMatrixVector_t *const LoadTdT,
			 const PID_t *const PID, PMatrixVector_t *const fe )
{

     int ione = 1;
     int incx = 1, incy = 1;
     HYSL_FLOAT Alpha, Beta;
     char uplo = 'L';

     /* PBLAS: fe = fc */
     hysl_pcopy( &fe->GlobalSize.Row, fc->Array, &ione, &ione, fc->Desc, &ione, fe->Array, &ione, &ione,
	      fe->Desc, &incy );
     /* PBLAS: fe = fc + LoadTdT */
     Alpha = 1.0;
     hysl_paxpy( &fe->GlobalSize.Row, &Alpha, LoadTdT->Array, &ione, &ione, LoadTdT->Desc, &incx, fe->Array,
	      &ione, &ione, fe->Desc, &incy );

     /* PBLAS: fe = fc + LoadTdT - Mass*AccTdT = fe - Mass*AccTdT */
     Alpha = -1.0; Beta = 1.0;
     hysl_psymv( &uplo, &fe->GlobalSize.Row, &Alpha, Mass->Array, &ione, &ione, Mass->Desc, AccTdT->Array,
	      &ione, &ione, AccTdT->Desc, &incx, &Beta, fe->Array, &ione, &ione, fe->Desc, &incy );
     /* PBLAS: fe = fc + LoadTdT - Mass*AccTdT - Damp*VelTdT = fe - Damp*VelTdT */
     hysl_psymv( &uplo, &fe->GlobalSize.Row, &Alpha, Damp->Array, &ione, &ione, Damp->Desc, VelTdT->Array,
	      &ione, &ione, VelTdT->Desc, &incx, &Beta, fe->Array, &ione, &ione, fe->Desc, &incy );
     /* PBLAS: fe = fc + LoadTdT - Mass*AccTdT - Damp*VelTdT - Stiff*DispTdT = fe - Stiff*DispTdT */
     hysl_psymv( &uplo, &fe->GlobalSize.Row, &Alpha, Stiff->Array, &ione, &ione, Stiff->Desc, DispTdT->Array,
	      &ione, &ione, DispTdT->Desc, &incx, &Beta, fe->Array, &ione, &ione, fe->Desc, &incy );

     /* Apply the PID */
     /* PBLAS: fe = P*fe */
     Alpha = -PID->P;
     hysl_pscal( &fe->GlobalSize.Row, &Alpha, fe->Array, &ione, &ione, fe->Desc, &incx );

}