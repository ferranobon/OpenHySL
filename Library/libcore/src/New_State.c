#include "MatrixVector.h"
#include "New_State.h"
#include "Definitions.h"

#if _MKL_
#include <mkl_blas.h>
#else
#include "Netlib.h"
#endif

void Compute_NewState( const MatrixVector_t *const IGain, const MatrixVector_t *const Eff_ForceT,
		       const MatrixVector_t *const In_LoadT, const MatrixVector_t *const Err_ForceT,
		       MatrixVector_t *const Tempvec, MatrixVector_t *const VecTdT_0 )
{
     int incx = 1, incy = 1;             /* Stride in the vectors */
     HYSL_FLOAT Alpha = 1.0, Beta = 0.0; /* Constants for the BLAS routines */
     char uplo = 'L';                    /* The lower part (upper part in C) will be used and the upper part
					  * (lower part in C) will strictly not be referenced */

     /* BLAS: tempvec = Eff_Force */
     hysl_copy( &Tempvec->Rows, Eff_ForceT->Array, &incx, Tempvec->Array, &incy );
     
     /* BLAS: tempvec = Eff_Force + LoadTdT = tempvec + LoadTdT */
     hysl_axpy( &Tempvec->Rows, &Alpha, In_LoadT->Array, &incx, Tempvec->Array, &incy );
     
     /* BLAS: tempvec = Eff_Force + LoadTdT - Err_Force = tempvec - Err_Force. */
     Alpha = -1.0;
     hysl_axpy( &Tempvec->Rows, &Alpha, Err_ForceT->Array, &incx, Tempvec->Array, &incy );

     /* BLAS: VecTdT_0 = Keinv*(Eff_Force + LoadTdT - Err_Force) = Keinv*Tempvec */
     Alpha = 1.0;
     hysl_symv( &uplo, &Tempvec->Rows, &Alpha, IGain->Array, &Tempvec->Rows, Tempvec->Array, &incx, &Beta,
	    VecTdT_0->Array, &incy );
}

void Compute_NewState_PS( const MatrixVector_t *const IGain, const MatrixVector_t *const Eff_ForceT,
			  const MatrixVector_t *const In_LoadT, const MatrixVector_t *const Err_ForceT,
			  MatrixVector_t *const Tempvec, MatrixVector_t *const VecTdT_0 )
{
     int incx = 1, incy = 1;             /* Stride in the vectors */
     HYSL_FLOAT Alpha = 1.0, Beta = 0.0; /* Constants for the BLAS routines */
     char uplo = 'L';                    /* The lower part (upper part in C) will be used and the upper part
					  * (lower part in C) will strictly not be referenced */

     /* BLAS: tempvec = Eff_Force */
     hysl_copy( &Tempvec->Rows, Eff_ForceT->Array, &incx, Tempvec->Array, &incy );
     
     /* BLAS: tempvec = Eff_Force + LoadTdT = tempvec + LoadTdT */
     hysl_axpy( &Tempvec->Rows, &Alpha, In_LoadT->Array, &incx, Tempvec->Array, &incy );
     
     /* BLAS: tempvec = Eff_Force + LoadTdT - Err_Force = tempvec - Err_Force. */
     Alpha = -1.0;
     hysl_axpy( &Tempvec->Rows, &Alpha, Err_ForceT->Array, &incx, Tempvec->Array, &incy );

     /* BLAS: Disp0 = Keinv*(Eff_Force + LoadTdT + Err_Force) = Keinv*Tempvec */
     Alpha = 1.0;
     hysl_spmv( &uplo, &Tempvec->Rows, &Alpha, IGain->Array, Tempvec->Array, &incx, &Beta,
	    VecTdT_0->Array, &incy );
}

void Compute_NewState_HHT( const MatrixVector_t *const IGain, const MatrixVector_t *const Eff_ForceT,
			   const MatrixVector_t *const In_LoadT, const MatrixVector_t *const In_LoadTdT,
			   const MatrixVector_t *const Err_ForceT, MatrixVector_t *const Tempvec,
			   const HYSL_FLOAT Alpha_H, MatrixVector_t *const VecTdT_0 )
{
     int incx = 1, incy = 1;             /* Stride in the vectors */
     HYSL_FLOAT Alpha = 1.0, Beta = 0.0; /* Constants for the BLAS routines */
     char uplo = 'L';                    /* The lower part (upper part in C) will be used and the upper part
					  * (lower part in C) will strictly not be referenced */
     
     /* BLAS: tempvec = Eff_Force */
     hysl_copy( &Tempvec->Rows, Eff_ForceT->Array, &incx, Tempvec->Array, &incy );

     /* BLAS: tempvec = Eff_Force + (1 + alphaH)*LoadTdT = tempvec + (1 + alphaH)*LoadTdT */
     Alpha = 1.0 + Alpha_H;
     hysl_axpy( &Tempvec->Rows, &Alpha, In_LoadTdT->Array, &incx, Tempvec->Array, &incy );

     /* BLAS: tempvec = Eff_Force + (1 + alphaH)*LoadTdT - alphaH*LoadT = tempvec - alphaH*LoadT */
     Alpha = -Alpha_H;
     hysl_axpy( &Tempvec->Rows, &Alpha, In_LoadTdT->Array, &incx, Tempvec->Array, &incy );
     
     /* BLAS: tempvec = Eff_Force + (1 + alphaH)*LoadTdT - alphaH*LoadT - Err_Force = tempvec - Err_Force. */
     Alpha = -1.0;
     hysl_axpy( &Tempvec->Rows, &Alpha, Err_ForceT->Array, &incx, Tempvec->Array, &incy );

     /* BLAS: VecTdT_0 = Keinv*(Eff_Force + (1 + alphaH)*LoadTdT - alphaH*LoadT - Err_Force) = Keinv*Tempvec */
     Alpha = 1.0;
     hysl_symv( &uplo, &Tempvec->Rows, &Alpha, IGain->Array, &Tempvec->Rows, Tempvec->Array, &incx, &Beta,
		VecTdT_0->Array, &incy );
}

void Compute_NewState_Zienkiewicz( const MatrixVector_t *const Meff, const MatrixVector_t *const MatA,
				   const MatrixVector_t *const MatB, const MatrixVector_t *const DispT,
				   const MatrixVector_t *const DispT0, const MatrixVector_t *const In_LoadT,
				   const MatrixVector_t *const Err_ForceT, const MatrixVector_t *const ForceT,
				   const MatrixVector_t *const ForceT0, const HYSL_FLOAT a8, const HYSL_FLOAT a17,
				   const HYSL_FLOAT a18, MatrixVector_t *const Tempvec, MatrixVector_t *const VecTdT_0 )
{
     int incx = 1, incy = 1;              /* Stride in the vectors */
     HYSL_FLOAT Alpha = 1.0, Beta = 1.0;  /* Constants for the BLAS routines */
     char uplo = 'L';                     /* The lower part (upper part in C) will be used and the upper part
					   * (lower part in C) will strictly not be referenced */

     /* BLAS: Tempvec = ForceT */
     hysl_copy( &Tempvec->Rows, ForceT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: Tempvec = (1/2 - 2*beta + gamma)*ForceT = a17*ForceT */
     Alpha = a17;
     hysl_scal( &Tempvec->Rows, &Alpha, Tempvec->Array, &incx );
     /* BLAS: Tempvec = (1/2 - 2*beta + gamma)*Deltat^2*ForceT + (1/2 + beta - gamma)*ForceT0 = Tempvec + a18*ForceT0 */
     Alpha = a18;
     hysl_axpy( &Tempvec->Rows, &Alpha, ForceT0->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: Tempvec = 1/2*beta*Deltat^2*In_LoadT +(1/2 - 2*beta + gamma)*ForceT + (1/2 + beta - gamma)*ForceT0 = Tempvec + a8*In_LoadT */
     Alpha = a8;
     hysl_axpy( &Tempvec->Rows, &Alpha, In_LoadT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: Tempvec = MatA*DispT + 1/2*beta*Deltat^2*In_LoadT +(1/2 - 2*beta + gamma)*ForceT + (1/2 + beta - gamma)*ForceT0 = Tempvec + MatA*DispT */
     Alpha = 1.0;
     hysl_symv( &uplo, &Tempvec->Rows, &Alpha, MatA->Array, &Tempvec->Rows, DispT->Array, &incx, &Beta,
		Tempvec->Array, &incy );
     /* BLAS: Tempvec = MatA*DispT - MatB*DispT0 + 1/2*beta*Deltat^2*In_LoadT +(1/2 - 2*beta + gamma)*ForceT
      * + (1/2 + beta - gamma)*ForceT0 = Tempvec - MatB*DispT0 */
     /* Note: MatB has already the minus sign applied */
     hysl_symv( &uplo, &Tempvec->Rows, &Alpha, MatB->Array, &Tempvec->Rows, DispT0->Array, &incx, &Beta,
		Tempvec->Array, &incy );
     /* BLAS: Tempvec = MatA*DispT - MatB*DispT0 + 1/2*beta*Deltat^2*In_LoadT +(1/2 - 2*beta + gamma)*ForceT
      * + (1/2 + beta - gamma)*ForceT0 = Tempvec - MatB*DispT0 + 1/2*beta*Deltat^2*Err_ForceT = Tempvec +
      * a8*Err_ForceT */
     Alpha = a8;
     hysl_axpy( &Tempvec->Rows, &Alpha, Err_ForceT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: VecTdT_0 = Meff*Tempvec */
     Alpha = 1.0;
     Beta = 0.0;
     hysl_symv( &uplo, &VecTdT_0->Rows, &Alpha, Meff->Array, &Tempvec->Rows, Tempvec->Array, &incx, &Beta,
		VecTdT_0->Array, &incy );
}

void Compute_NewState_Zienkiewicz_PS( const MatrixVector_t *const Meff, const MatrixVector_t *const MatA,
				      const MatrixVector_t *const MatB, const MatrixVector_t *const DispT,
				      const MatrixVector_t *const DispT0, const MatrixVector_t *const In_LoadT,
				      const MatrixVector_t *const Err_ForceT, const MatrixVector_t *const ForceT,
				      const MatrixVector_t *const ForceT0, const HYSL_FLOAT a8, const HYSL_FLOAT a17,
				      const HYSL_FLOAT a18, MatrixVector_t *const Tempvec, MatrixVector_t *const VecTdT_0 )
{
     int incx = 1, incy = 1;              /* Stride in the vectors */
     HYSL_FLOAT Alpha = 1.0, Beta = 1.0;  /* Constants for the BLAS routines */
     char uplo = 'L';                     /* The lower part (upper part in C) will be used and the upper part
					   * (lower part in C) will strictly not be referenced */

     /* Blas: Tempvec = ForceT */
     hysl_copy( &Tempvec->Rows, ForceT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: Tempvec = (1/2 - 2*beta + gamma)*ForceT = a17*ForceT */
     Alpha = a17;
     hysl_scal( &Tempvec->Rows, &Alpha, Tempvec->Array, &incx );
     /* BLAS: Tempvec = (1/2 - 2*beta + gamma)*Deltat^2*ForceT + (1/2 + beta - gamma)*ForceT0 = Tempvec + a18*ForceT0 */
     Alpha = a18;
     hysl_axpy( &Tempvec->Rows, &Alpha, ForceT0->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: Tempvec = 1/2*beta*Deltat^2*In_LoadT +(1/2 - 2*beta + gamma)*ForceT + (1/2 + beta - gamma)*ForceT0 = Tempvec + a8*In_LoadT */
     Alpha = a8;
     hysl_axpy( &Tempvec->Rows, &Alpha, In_LoadT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: Tempvec = MatA*DispT + 1/2*beta*Deltat^2*In_LoadT +(1/2 - 2*beta + gamma)*ForceT + (1/2 + beta - gamma)*ForceT0 = Tempvec + MatA*DispT */
     Alpha = 1.0;
     hysl_spmv( &uplo, &Tempvec->Rows, &Alpha, MatA->Array, DispT->Array, &incx, &Beta, Tempvec->Array,
		&incy );
     /* BLAS: Tempvec = MatA*DispT + MatB*DispT0 + 1/2*beta*Deltat^2*In_LoadT +(1/2 - 2*beta + gamma)*ForceT
      * + (1/2 + beta - gamma)*ForceT0 = Tempvec - MatB*DispT0 */
     /* Note: MatB has already the minus sign applied */
     hysl_spmv( &uplo, &Tempvec->Rows, &Alpha, MatB->Array, DispT0->Array, &incx, &Beta, Tempvec->Array,
		&incy );
     /* BLAS: Tempvec = MatA*DispT - MatB*DispT0 + 1/2*beta*Deltat^2*In_LoadT +(1/2 - 2*beta + gamma)*ForceT
      * + (1/2 + beta - gamma)*ForceT0 = Tempvec - MatB*DispT0 + 1/2*beta*Deltat^2*Err_ForceT = Tempvec +
      * a8*Err_ForceT */
     Alpha = a8;
     hysl_axpy( &Tempvec->Rows, &Alpha, Err_ForceT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: VecTdT_0 = Meff*Tempvec */
     Alpha = 1.0;
     Beta = 0.0;
     hysl_spmv( &uplo, &VecTdT_0->Rows, &Alpha, Meff->Array, Tempvec->Array, &incx, &Beta, VecTdT_0->Array,
		&incy );
}
