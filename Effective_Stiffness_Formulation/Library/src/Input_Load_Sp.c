#include "Common_Formulation_Sp.h"
#include "MatrixVector.h"
#include "MatrixVector_Sp.h"

#include <mkl_spblas.h>

void InputLoad_AbsValues_Sp( const MatrixVector_Sp_t *const Stiff, const MatrixVector_Sp_t *const Damp,
			     const MatrixVector_t *const GDisp, const MatrixVector_t *const GVel,
			     MatrixVector_t *const InLoad )
{

     static double Alpha, Beta;        /* Constants to use in the Sparse BLAS routines */
     static char trans = 'N';          /* No transpose operation */
     static char matdescra[6] = {'S',  /* The matrix is symmetric */
				 'U',  /* The upper part is referenced */
				 'N',  /* Non-unit values in the diagonal */
				 'F'}; /* One based index */

     Alpha = 1.0; Beta = 0.0;

     /* Sparse BLAS: li = K*ug */
     mkl_dcsrmv( &trans, &InLoad->Rows, &InLoad->Rows, &Alpha, matdescra, Stiff->Values, Stiff->Columns,
		 Stiff->RowIndex, &Stiff->RowIndex[1], GDisp->Array, &Beta, InLoad->Array );

     /* Sparse BLAS: li = K*ug + C*vg = li + C*vg */
     Beta = 1.0;
     mkl_dcsrmv( &trans, &InLoad->Rows, &InLoad->Rows, &Alpha, matdescra, Damp->Values, Damp->Columns,
		 Damp->RowIndex, &Damp->RowIndex[1], GVel->Array, &Beta, InLoad->Array );
}

void InputLoad_RelValues_Sp( const MatrixVector_Sp_t *const Mass, const MatrixVector_t *const GAcc, 
			     MatrixVector_t *const InLoad )
{
     
     static double Alpha, Beta;        /* Constants to use in the Sparse BLAS routines */
     static char trans = 'N';          /* No transpose operation */
     static char matdescra[6] = {'S',  /* The matrix is symmetric */
				 'U',  /* The upper part is referenced */
				 'N',  /* Non-unit values in the diagonal */
				 'F'}; /* One based index */

     Alpha = -1.0; Beta = 0.0;

     /* Sparse BLAS: li = -M*ag */
     mkl_dcsrmv( &trans, &InLoad->Rows, &InLoad->Rows, &Alpha, matdescra, Mass->Values, Mass->Columns,
		 Mass->RowIndex, &Mass->RowIndex[1], GAcc->Array, &Beta, InLoad->Array );
}
