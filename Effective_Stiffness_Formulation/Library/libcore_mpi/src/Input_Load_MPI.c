#include "Input_Load.h"
#include "MatrixVector_MPI.h"

#if _MKL_
#include <mkl_pblas.h>
#include "Cblacs.h"
#include "Scalapack_Aux.h"
#else
#include "Netlib.h"
#endif

void InputLoad_AbsValues_MPI( PMatrixVector_t *const Stiff, PMatrixVector_t *const Damp,
			      PMatrixVector_t *const GDisp, PMatrixVector_t *const GVel,
			      PMatrixVector_t *const InLoad )
{

     int incx, incy;     /* Stride in the vectors for PBLAS library */
     int ione;           /* Integer variable of value 1 for PBLAS library */
     double Alpha, Beta; /* Constants to use in the PBLAS library */
     char uplo;          /* Character to use in the PBLAS library */

     incx = 1; incy = 1;
     Alpha = 1.0; Beta = 0.0;
     uplo = 'L';
     ione = 1;

     /* PBLAS: li = K*ug */
     pdsymv_( &uplo, &InLoad->GlobalSize.Row, &Alpha, Stiff->Array, &ione, &ione, Stiff->Desc, GDisp->Array,
	      &ione, &ione, GDisp->Desc, &incx, &Beta, InLoad->Array, &ione, &ione, InLoad->Desc, &incy );

     /* BLAS: li = K*ug + C*vg = li + C*vg */
     Beta = 1.0;
     pdsymv_( &uplo, &InLoad->GlobalSize.Row, &Alpha, Damp->Array, &ione, &ione, Damp->Desc, GVel->Array,
	      &ione, &ione, GVel->Desc, &incx, &Beta, InLoad->Array, &ione, &ione, InLoad->Desc, &incy );
}

void InputLoad_RelValues_MPI( PMatrixVector_t *const Mass, PMatrixVector_t *const GAcc,
			      PMatrixVector_t *const InLoad )
{
     int incx, incy;     /* Stride in the vectors for PBLAS library */
     int ione;           /* Integer variable of value 1 for PBLAS library */
     double Alpha, Beta; /* Constants to use in the PBLAS library */
     char uplo;          /* Character to use in the PBLAS library */

     incx = 1; incy = 1;
     Alpha = -1.0; Beta = 0.0;
     uplo = 'L';

     pdsymv_( &uplo, &InLoad->GlobalSize.Row, &Alpha, Mass->Array, &ione, &ione, Mass->Desc, GAcc->Array,
	      &ione, &ione, GAcc->Desc, &incx, &Beta, InLoad->Array, &ione, &ione, InLoad->Desc, &incy );
}

void InputLoad_Generate_LoadVectorForm_MPI( int *DOF, PMatrixVector_t *const LoadVectorForm )
{
     int i, j, ione = 1;
     int myrow, mycol, nprow, npcol;
     int RowProcess, ColProcess;
     int LRowIndex, LColIndex;

     /* Get grid information */
     Cblacs_gridinfo( LoadVectorForm->Desc[1], &nprow, &npcol, &myrow, &mycol );

     i = 0;
     while( i < LoadVectorForm->GlobalSize.Row ){	  
	  for ( j = 1; j < DOF[0]; j++ ){

	       /* Get the local coordinates */
	       infog2l_( &i, &ione, LoadVectorForm->Desc, &nprow, &npcol, &myrow, &mycol, &LRowIndex,
			 &LColIndex, &RowProcess, &ColProcess );
	       /* If its the right process, assign the value */
	       if ( (Cblacs_pnum( LoadVectorForm->Desc[1], RowProcess, ColProcess ) == 0) &&
		    (myrow == RowProcess) && (mycol == ColProcess) ){ 
		    LoadVectorForm->Array[i] = (double) DOF[j];
	       }
	       i = i + 1;
	  }
     }
}

void InputLoad_Apply_LoadVectorForm_MPI( PMatrixVector_t *const LoadForm, const double Value,
					 PMatrixVector_t *const LoadVector )
{
     int incx = 1;
     int incy = 1;
     int ione = 1;
     double Scalar;

     Scalar = Value;
     
     pdcopy_( &LoadVector->GlobalSize.Row, LoadForm->Array, &ione, &ione, LoadForm->Desc, &ione,
	      LoadVector->Array, &ione, &ione, LoadVector->Desc, &incy );
     pdscal_( &LoadVector->GlobalSize.Row, &Scalar, LoadVector->Array, &ione, &ione, LoadVector->Desc, &incx );
}
