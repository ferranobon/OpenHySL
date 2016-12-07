#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "MatrixVector.h"

#include "Substructure.h"
#include "Substructure_Auxiliary.h"

#include "Auxiliary_Math.h"   /* For Max( ) */
#include "Definitions.h"

#if _MKL_
#include <mkl_blas.h>
#else
#include "Netlib.h"
#endif

void Substructure_JoinNonCouplingPart( MatrixVector_t *const VecTdT_m, const MatrixVector_t *const Gain_m,
				       const MatrixVector_t *const fcprevsub,
				       const CouplingNode_t *const CNodes, MatrixVector_t *const VecTdT )			  
{
     int icoup;                 /* Counter for the coupling nodes */
     int incx, incy;            /* Stride in the vectors */
     HYSL_FLOAT Alpha, Beta;    /* Constants for the BLAS routines */
     char trans;                /* Use or not the transpose */
     int Rows, Cols;            /* Number of Rows and columns */
     int lda;                   /* Leading dimension */
     int Length, PosX, PosXm;   /* Length and position counters */
     
     incx = 1; incy = 1;
     trans = 'N';
     Alpha = 1.0; Beta = 1.0;
     Rows = Gain_m->Rows;
     Cols = Gain_m->Cols;
     lda = Max( 1, Gain_m->Rows);

     /* Update the VecTdT_m displacments to include the effects of the coupling force */
     /* BLAS: VecTdT_m = Gain_m*fcprevsub */
     hysl_gemv( &trans, &Rows, &Cols, &Alpha, Gain_m->Array, &lda,
	     fcprevsub->Array, &incx, &Beta, VecTdT_m->Array, &incy );

     /* Copy the updated values into the complete displacement vector */
     PosX = 0; PosXm = 0;
     for ( icoup = 0; icoup < CNodes->Order; icoup++ ){
	  Length = CNodes->Array[icoup] - PosX -1;
	  hysl_copy( &Length, &VecTdT_m->Array[PosXm], &incx, &VecTdT->Array[PosX], &incy );
	  PosX = CNodes->Array[icoup];
	  PosXm = PosXm + Length;
     }

     /* Add the elements between the final coupling node and the final element of the complete displacement
      * vector */
     Length = VecTdT->Rows - CNodes->Array[CNodes->Order -1];
     hysl_copy( &Length, &VecTdT_m->Array[PosXm], &incx, &VecTdT->Array[PosX], &incy );	
}

void Substructure_MatrixXc( const MatrixVector_t *const Mat, const CouplingNode_t *const CNodes,
			    MatrixVector_t *const MatCouple )
{

     int icoup;    /* Counter for the coupling nodes */
     int jcoup;    /* Another counter */
     
     for( icoup = 0; icoup < CNodes->Order; icoup++ ){

	  for (jcoup = icoup; jcoup < CNodes->Order; jcoup++){
	       
	       MatCouple->Array[icoup*CNodes->Order + jcoup] = Mat->Array[(CNodes->Array[icoup]-1)*Mat->Cols + CNodes->Array[jcoup] - 1];
	       /* Now add the elements belonging to the same row as the current coupling node but also
		* belonging to the same column as the rest of the coupling nodes */
	       if( icoup != jcoup ){
		    MatCouple->Array[jcoup*CNodes->Order + icoup] = MatCouple->Array[icoup*CNodes->Order + jcoup];
	       }
	  }
     }
}

void Substructure_MatrixXc_PS( const MatrixVector_t *const Mat, const CouplingNode_t *const CNodes,
			       MatrixVector_t *const MatCouple )
{

     int icoup;       /* Counter for the coupling nodes */
     int jcoup;       /* Another counter */
     unsigned int Index_PS; /* Index of the matrix in packed storage */    


     for( icoup = 0; icoup < CNodes->Order; icoup++ ){

	  for (jcoup = icoup; jcoup < CNodes->Order; jcoup++){
	       Index_PS = MatrixVector_ReturnIndex_UPS( (unsigned int) CNodes->Array[icoup],
							(unsigned int) CNodes->Array[jcoup], Mat->Rows );

	       MatCouple->Array[icoup*CNodes->Order + jcoup] = Mat->Array[Index_PS];
	       /* Now add the elements belonging to the same row as the current coupling node but also
		* belonging to the same column as the rest of the coupling nodes */
	       if( icoup != jcoup ){
		    MatCouple->Array[jcoup*CNodes->Order + icoup] = Mat->Array[Index_PS];
	       }
	  }
     }
}

void Substructure_MatrixXcm( const MatrixVector_t *const Mat, const CouplingNode_t *const CNodes,
			     MatrixVector_t *const MatXcm )
{

     int Length;
     int icoup;      /* Counter for the coupling nodes */
     int jcoup;
     int incx, incy;
     int PosXcm, Acumulated_Length;

     incx = Mat->Rows;
     incy = 1;          /* The values in the Xm matrix are stored in columns. Therefore the stride has to be 1
			 * because of the way FORTRAN handles arrays (major column ordering) */

     /* Since the matrix is symmetric and only a part of it is stored, this routine has to be splitted into
      * two parts. The first will copy the elements above the coupling node, while the second will focus on
      * the other part */

     /* Copy until the first coupling node */
     PosXcm = 0;
     Length = CNodes->Array[0] - 1;
     for ( jcoup = 0; jcoup < CNodes->Order; jcoup++ ){
	  PosXcm = jcoup*MatXcm->Rows;
	  hysl_copy( &Length, &Mat->Array[CNodes->Array[jcoup] - 1], &incx, &MatXcm->Array[PosXcm], &incy );
     }

     /* Copy until the last coupling node */
     Acumulated_Length = Length;
     for( icoup = 1; icoup < CNodes->Order; icoup++ ){
	     
	  Length = CNodes->Array[icoup] - CNodes->Array[icoup-1] - 1;
	  for ( jcoup = icoup; jcoup < CNodes->Order; jcoup++ ){
	 
	       PosXcm = jcoup*MatXcm->Rows + Acumulated_Length;
	       
	       hysl_copy( &Length, &Mat->Array[(CNodes->Array[jcoup] - 1) + (CNodes->Array[icoup-1])*Mat->Rows],
		      &incx, &MatXcm->Array[PosXcm], &incy );
	  }
	  Acumulated_Length = Acumulated_Length + Length;
     }

     /* Do the same but starting from the opposite side */
     incx = 1;
     Length = Mat->Rows - CNodes->Array[CNodes->Order -1];
     for ( jcoup = CNodes->Order - 1; jcoup >= 0; jcoup = jcoup - 1 ){
	  PosXcm = MatXcm->Rows*(jcoup+1) - Length;
	  hysl_copy( &Length, &Mat->Array[(CNodes->Array[jcoup] - 1)*Mat->Rows + (CNodes->Array[CNodes->Order-1])],
		 &incx, &MatXcm->Array[PosXcm], &incy );
     }
     Acumulated_Length = Length;
     incx = 1;
     for( icoup = CNodes->Order -2; icoup >= 0; icoup = icoup -1 ){
	  Length = CNodes->Array[icoup + 1] - CNodes->Array[icoup] - 1;
	  for( jcoup = icoup; jcoup >= 0; jcoup = jcoup - 1){
	       PosXcm = (jcoup+1)*MatXcm->Rows - Acumulated_Length - Length;
	       hysl_copy( &Length, &Mat->Array[(CNodes->Array[jcoup] - 1)*Mat->Rows + (CNodes->Array[icoup])],
		      &incx, &MatXcm->Array[PosXcm], &incy );
	  }
	  Acumulated_Length = Acumulated_Length + Length;
     }
}

void Substructure_MatrixXcm_PS( const MatrixVector_t *const Mat, const CouplingNode_t *const CNodes,
				MatrixVector_t *const MatXcm )
{
     int icoup, jcoup, kcoup;  /* Counters for the coupling nodes */
     int IndexXcm;             /* Index of the Xcm matrix */
     unsigned int RowIndex, ColIndex;   /* Row and column indixes for the full matrix in packed storage */
     unsigned int Index;       /* Index of the matrix in packed storage: Mat */

     /* Since the matrix is symmetric and only a part of it is stored, this routine has to be splitted into
      * two parts. The first will copy the elements above the coupling node, while the second will focus on
      * the other part */


     IndexXcm = 0;
     for( icoup = 0; icoup < CNodes->Order; icoup++ ){
	  ColIndex = (unsigned int) CNodes->Array[icoup];
	  RowIndex = 1;
	  jcoup = 0; kcoup = 0;
	  while( RowIndex <= (unsigned int) Mat->Rows ){
	       if ( RowIndex == (unsigned int) CNodes->Array[kcoup] ){
		    kcoup = kcoup + 1;
	       } else {
		    IndexXcm = icoup*MatXcm->Rows + jcoup;
		    Index = MatrixVector_ReturnIndex_UPS( RowIndex, ColIndex, Mat->Rows );
		    MatXcm->Array[IndexXcm] = Mat->Array[Index];
		    jcoup = jcoup + 1;
	       }
	       RowIndex = RowIndex + 1;
	  }
     }
}

void Substructure_VectorXm( const MatrixVector_t *const VectorX, const CouplingNode_t *const CNodes,
			    MatrixVector_t *const VectorXm )
{

	int incx, incy;
	int Length;
	int PosX, PosXm;
	int icoup;    /* Counter for the coupling nodes */

	incx = 1; incy = 1;
	PosX = 0; PosXm = 0;
	for( icoup = 0; icoup < CNodes->Order; icoup++ ){

	     /* Copy the part of the vector between twwo positions */
	     Length = CNodes->Array[icoup] - PosX - 1;
	     hysl_copy( &Length, &VectorX->Array[PosX], &incx, &VectorXm->Array[PosXm], &incy );
	     /* Update the values of the position in the vectors */
	     PosX = CNodes->Array[icoup];
	     PosXm = PosXm + Length;
	}

	/* Copy the elements from the last position until the end of the vector */
	Length = VectorX->Rows - CNodes->Array[CNodes->Order-1];
	hysl_copy( &Length, &VectorX->Array[PosX], &incx, &VectorXm->Array[PosXm], &incy );
}

void Substructure_VectorXc( const MatrixVector_t *const VecX, const CouplingNode_t *const CNodes,
			    MatrixVector_t *const VecXc )
{
     int icoup;    /* Counter for the coupling nodes */

#pragma omp parallel for
     for( icoup = 0; icoup < CNodes->Order; icoup++ ){
	  VecXc->Array[icoup] = VecX->Array[CNodes->Array[icoup] - 1];
     }
}

int Substructure_FindPosition( const int Pos_in_VecX, const CouplingNode_t *const CNodes )
{
     int pos;
     bool found;

     pos = 0;
     found = false;
     
     while (!found && (pos < CNodes->Order)) {
	  if (CNodes->Array[pos] == Pos_in_VecX){
	       found = true;
	  } else {
	       pos = pos + 1;
	  }
     }

     if (!found){
	  printf( "Position not found" );
	  exit( EXIT_FAILURE );
     }

     return pos;	  
}