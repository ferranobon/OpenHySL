#include "GSSSS.h"
#include "MatrixVector_Sp.h"
#include "Error_Compensation.h"

#include "Definitions.h"

#include <mkl_blas.h>
#include <mkl_spblas.h>


void GSSSS_EffectiveForce_AForm_Sp( const MatrixVector_Sp_t *const Mass, const MatrixVector_Sp_t *const Damp,
				    const MatrixVector_Sp_t *const Stiff, const MatrixVector_t *const DispT,
				    const MatrixVector_t *const VelT, const MatrixVector_t *const AccT,
				    MatrixVector_t *const Tempvec, const TIntegration_GSSSS_t *const GSSSS,
				    const HYSL_FLOAT DeltaT, MatrixVector_t *const Eff_ForceT )
{

     int incx = 1, incy = 1;  /* Stride in the vectors */
     char trans = 'N';
     char matdescra[6] = {'S',  /* The matrix is symmetric */
			  'U',  /* The upper part is referenced */
			  'N',  /* Non-unit values in the diagonal */
			  'F'}; /* One based index */
     HYSL_FLOAT Alpha, Beta;    /* Constants for the BLAS routines */
     

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
     hysl_mkl_csrmv( &trans, &Tempvec->Rows, &Tempvec->Rows, &Alpha, matdescra, Stiff->Values, Stiff->Columns,
		     Stiff->RowIndex, &Stiff->RowIndex[1], Tempvec->Array, &Beta, Eff_ForceT->Array );

     /* BLAS: tempvec = Vel */
     hysl_copy( &Tempvec->Rows, VelT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: tempvec = Vel + (A4W1 -A5W2)*DeltaT*Acc = tempvec + (A4W1 - A5W2)*DeltaT*Acc */
     Alpha = (GSSSS->A4W1 - GSSSS->A5W2)*DeltaT;
     hysl_axpy( &Tempvec->Rows, &Alpha, AccT->Array, &incx, Tempvec->Array, &incy );
     /* BLAS: Eff_Force = -Stiff*(Disp + A1W1*Vel*DeltaT + (A2W2 - A3W3)*Acc*DeltaT^2) -C*(Vel + (A4W1
      * -A5W2)*DeltaT*Acc)= Eff_Force - Damp*tempvec */
     Alpha = -1.0; Beta = 1.0;
     hysl_mkl_csrmv( &trans, &Tempvec->Rows, &Tempvec->Rows, &Alpha, matdescra, Damp->Values, Damp->Columns,
		     Damp->RowIndex, &Damp->RowIndex[1], Tempvec->Array, &Beta, Eff_ForceT->Array );

     /* BLAS: Eff_Force = Stiff*(Disp + A1W1*Vel*DeltaT + (A2W2 - A3W3)*Acc*DeltaT^2) -C*(Vel + (A4W1
      * -A5W2)*DeltaT*Acc) -(1-A6W1)*M*Acc= Eff_Force - (1-A6W1)*M*Acc */
     Alpha = -(1.0 - GSSSS->A6W1)*DeltaT;
     hysl_mkl_csrmv( &trans, &Tempvec->Rows, &Tempvec->Rows, &Alpha, matdescra, Mass->Values, Mass->Columns,
		     Mass->RowIndex, &Mass->RowIndex[1], AccT->Array, &Beta, Eff_ForceT->Array );
}


void GSSSS_ErrorForce_PID_Sp( const MatrixVector_Sp_t *const Mass, const MatrixVector_Sp_t *const Damp,
			      const MatrixVector_Sp_t *Stiff, const MatrixVector_t *const AccTdT,
			      const MatrixVector_t *const AccT, const MatrixVector_t *const VelT,
			      const MatrixVector_t *const DispT, const MatrixVector_t *const fc,
			      const MatrixVector_t *const LoadTdT, const MatrixVector_t *const LoadT,
			      const TIntegration_GSSSS_t *const GSSSS, const HYSL_FLOAT DeltaT,
			      const MatrixVector_t *const Tempvec, const PID_t *const PID,
			      MatrixVector_t *const fe )
{

     int incx = 1, incy = 1;  /* Stride in the vectors */
     char trans = 'N';
     char matdescra[6] = {'S',  /* The matrix is symmetric */
			  'U',  /* The upper part is referenced */
			  'N',  /* Non-unit values in the diagonal */
			  'F'}; /* One based index */
     HYSL_FLOAT Alpha, Beta;    /* Constants for the BLAS routines */

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
     hysl_mkl_csrmv( &trans, &Tempvec->Rows, &Tempvec->Rows, &Alpha, matdescra, Mass->Values, Mass->Columns,
		     Mass->RowIndex, &Mass->RowIndex[1], Tempvec->Array, &Beta, fe->Array );

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
     hysl_mkl_csrmv( &trans, &Tempvec->Rows, &Tempvec->Rows, &Alpha, matdescra, Damp->Values, Damp->Columns,
		     Damp->RowIndex, &Damp->RowIndex[1], Tempvec->Array, &Beta, fe->Array );

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
     hysl_mkl_csrmv( &trans, &Tempvec->Rows, &Tempvec->Rows, &Alpha, matdescra, Stiff->Values, Stiff->Columns,
		     Stiff->RowIndex, &Stiff->RowIndex[1], Tempvec->Array, &Beta, fe->Array );

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
