#include "Common_Formulation.h"
#include "MatrixVector.h"

#if _MKL_
#include <mkl_blas.h>
#else
#include "Netlib.h"
#endif

void InputLoad_AbsValues( const MatrixVector_t *const Stiff, const MatrixVector_t *const Damp,
			  const MatrixVector_t *const GDisp, const MatrixVector_t *const GVel,
			  MatrixVector_t *const InLoad )
{

     int incx, incy;       /* Stride in the vectors for BLAS library */
     double Alpha, Beta;   /* Constants to use in the BLAS library */
     char uplo;            /* Character to use in the BLAS library */

     incx = 1; incy = 1;
     Alpha = 1.0; Beta = 0.0;
     uplo = 'L';                  /* The lower part (FORTRAN) of the symmetric matrix is
				   * referenced (upper part in C) */


     /* BLAS: li = K*ug */
     dsymv( &uplo, &InLoad->Rows, &Alpha, Stiff->Array, &InLoad->Rows, GDisp->Array, &incx, &Beta,
	     InLoad->Array, &incy );

     /* BLAS: li = K*ug + C*vg = li + C*vg */
     Beta = 1.0;
     dsymv( &uplo, &InLoad->Rows, &Alpha, Damp->Array, &InLoad->Rows, GVel->Array, &incx, &Beta,
	     InLoad->Array, &incy );
}

void InputLoad_AbsValues_PS( const MatrixVector_t *const Stiff, const MatrixVector_t *const Damp,
			     const MatrixVector_t *const GDisp, const MatrixVector_t *const GVel,
			     MatrixVector_t *const InLoad )
{

     int incx, incy;       /* Stride in the vectors for BLAS library */
     double Alpha, Beta;   /* Constants to use in the BLAS library */
     char uplo;            /* Character to use in the BLAS library */

     incx = 1; incy = 1;
     Alpha = 1.0; Beta = 0.0;
     uplo = 'L';                  /* The lower part (FORTRAN) of the symmetric matrix is
				   * referenced (upper part in C) */


     /* BLAS: li = K*ug */
     dspmv( &uplo, &InLoad->Rows, &Alpha, Stiff->Array, GDisp->Array, &incx, &Beta,
	    InLoad->Array, &incy );

     /* BLAS: li = K*ug + C*vg = li + C*vg */
     Beta = 1.0;
     dspmv( &uplo, &InLoad->Rows, &Alpha, Damp->Array, GVel->Array, &incx, &Beta,
	    InLoad->Array, &incy );
}

void InputLoad_RelValues( const MatrixVector_t *const Mass, const MatrixVector_t *const GAcc,
			  MatrixVector_t *const InLoad )
{

     int incx = 1, incy = 1;           /* Stride in the vectors for BLAS library */
     double Alpha = -1.0, Beta = 0.0;  /* Constants to use in the BLAS library */
     char uplo = 'L';                  /* Character defining that the lower part (FORTRAN)
					* of the symmetric matrix is referenced (upper
					* part in C) */

     /* BLAS: li = -M*ag */
     dsymv( &uplo, &InLoad->Rows, &Alpha, Mass->Array, &InLoad->Rows, GAcc->Array, &incx, &Beta,
	     InLoad->Array, &incy );
}

void InputLoad_RelValues_PS( const MatrixVector_t *const Mass, const MatrixVector_t *const GAcc,
			  MatrixVector_t *const InLoad )
{

     int incx = 1, incy = 1;           /* Stride in the vectors for BLAS library */
     double Alpha = -1.0, Beta = 0.0;  /* Constants to use in the BLAS library */
     char uplo = 'L';                  /* Character defining that the lower part (FORTRAN)
					* of the symmetric matrix is referenced (upper
					* part in C) */

     /* BLAS: li = -M*ag */
     dspmv( &uplo, &InLoad->Rows, &Alpha, Mass->Array, GAcc->Array, &incx, &Beta,
	    InLoad->Array, &incy );
}

void InputLoad_Generate_LoadVectorForm( int *DOF, MatrixVector_t *const LoadVectorForm )
{
     int i, j;

     i = 0;
     while( i < LoadVectorForm->Rows ){	  
	  for ( j = 1; j < DOF[0]; j++ ){
	       LoadVectorForm->Array[i] = (double) DOF[j];
	       i = i + 1;
	  }
     }
}

void InputLoad_Apply_LoadVectorForm( const MatrixVector_t *const LoadForm, const double Value, MatrixVector_t *const LoadVector )
{
     int incx = 1;
     int incy = 1;
     double Scalar;

     Scalar = Value;
     dcopy( &LoadVector->Rows, LoadForm->Array, &incx, LoadVector->Array, &incy );
     dscal( &LoadVector->Rows, &Scalar, LoadVector->Array, &incx );
}
