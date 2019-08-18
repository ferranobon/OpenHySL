#include "MatrixVector.h"
#include "MatrixVector_MPI.h"

#include "Substructure.h"
#include "Substructure_Auxiliary.h"

#include "mpi.h"

#if _MKL_
#include <mkl_pblas.h>
#include "mkl_blacs.h"
#include "mkl_blacs.h"
#include "Cblacs.h"
#include "Scalapack_Aux.h"
#else
#include "Netlib.h"
#endif

void Substructure_JoinNonCouplingPart_MPI( PMatrixVector_t *const VecTdT_m,
					   PMatrixVector_t *const Gain_m,
					   PMatrixVector_t *const fcprevsub,
					   const CouplingNode_t *const CNodes, PMatrixVector_t *const VecTdT )
{

     int incx, incy, ione;
     int icoup, Length;
     hysl_float_t Alpha, Beta;
     char trans = 'N';
     int Rows, Cols;
     int PosX_Row, PosX_Col, PosXm_Row, PosXm_Col;

     incx = 1; incy = 1;
     ione = 1;
     Rows = Gain_m->GlobalSize.Row;
     Cols = Gain_m->GlobalSize.Col;
     Alpha = 1.0; Beta = 1.0;
     trans = 'N';

     /* Update the VecTdT_m displacments to include the effects of the coupling force */
     /* PBLAS: VecTdT_m = Gain_m*fcprevsub */
     hysl_pgemv( &trans, &Rows, &Cols, &Alpha, Gain_m->Array, &ione, &ione, Gain_m->Desc, fcprevsub->Array,
	      &ione, &ione, fcprevsub->Desc, &incx, &Beta, VecTdT_m->Array, &ione, &ione, VecTdT_m->Desc,
	      &incy );

     PosX_Col = 1;
     PosXm_Col = 1;
     
     PosX_Row = 1;
     PosXm_Row = 1;
     for( icoup = 0; icoup < CNodes->Order; icoup++ ){
	  Length = CNodes->Array[icoup] - PosX_Row;
	  
	  /* Copy the part of the vector between two positions */
	  hysl_pcopy( &Length, VecTdT_m->Array, &PosXm_Row, &PosXm_Col, VecTdT_m->Desc, &incx, VecTdT->Array,
		   &PosX_Row, &PosX_Col, VecTdT->Desc, &incy );
	  
	  /* Update the values of the position in the vectors */
	  PosX_Row = CNodes->Array[icoup] + 1; /* 1 based index */
	  PosXm_Row = PosXm_Row + Length;
	}

     /* Copy the elements from the last position until the end of the vector */
     Length = VecTdT->GlobalSize.Row - CNodes->Array[CNodes->Order-1];
     hysl_pcopy( &Length, VecTdT_m->Array, &PosXm_Row, &PosXm_Col, VecTdT_m->Desc, &incx, VecTdT->Array,
	      &PosX_Row, &PosX_Col, VecTdT->Desc, &incy );

}

void Substructure_MatrixXc_MPI( const MPI_Comm Comm, const CouplingNode_t *const CNodes,
				PMatrixVector_t *const Mat, MatrixVector_t *const MatCouple )
{

     int icoup, jcoup;
     int nprow, npcol, myrow, mycol;
     int rank;
     int GRowIndex, GColIndex;
     int LRowIndex, LColIndex;
     int RowProcess, ColProcess;

     MPI_Status status;

     MPI_Comm_rank( Comm, &rank );

     /* Get grid info */
     Cblacs_gridinfo( Mat->Desc[1], &nprow, &npcol, &myrow, &mycol );
     
     for( icoup = 0; icoup < CNodes->Order; icoup++ ){
	  GColIndex = CNodes->Array[icoup];
	  for (jcoup = icoup; jcoup < CNodes->Order; jcoup++){
	       GRowIndex = CNodes->Array[jcoup];

	       /* Given the global index of an element (GRowIndex, GColIndex) returns the local index of the
		  element (LRowIndex, LColIndex) and the coordinates of the process (Row Process,
		  ColProcess) */
	       infog2l_( &GRowIndex, &GColIndex, Mat->Desc, &nprow, &npcol, &myrow, &mycol, &LRowIndex,
			 &LColIndex, &RowProcess, &ColProcess );

	       if ( Cblacs_pnum( Mat->Desc[1], RowProcess, ColProcess ) == 0 && myrow == RowProcess &&
		    mycol == ColProcess ){ 
		    MatCouple->Array[icoup*CNodes->Order + jcoup] = Mat->Array[(LRowIndex - 1) + Mat->LocalSize.Row*(LColIndex - 1)];

		    /* Now add the elements belonging to the same row as the current coupling node but also
		     * belonging to the same column as the rest of the coupling nodes */
		    MatCouple->Array[jcoup*CNodes->Order + icoup] = MatCouple->Array[icoup*CNodes->Order + jcoup];
	       } else {
		    if ( myrow == RowProcess && mycol == ColProcess ){			 
			 MPI_Send( &Mat->Array[(LRowIndex - 1)+Mat->LocalSize.Row*(LColIndex - 1)], 1,
				   MPI_HYSL_FLOAT, 0, MATRIX_XC, Comm );
		    }

		    if ( rank == 0 ){
			 /* Store the diagonal elements */
			 MPI_Recv( &MatCouple->Array[icoup*CNodes->Order + jcoup], 1, MPI_HYSL_FLOAT,
				   Cblacs_pnum( Mat->Desc[1], RowProcess, ColProcess ), MATRIX_XC, Comm,
				   &status );
			 /* Now add the elements belonging to the same row as the current coupling node but
			  * also belonging to the same column as the rest of the coupling nodes */
			 MatCouple->Array[jcoup*CNodes->Order + icoup] = MatCouple->Array[icoup*CNodes->Order + jcoup];
		    }
	       }
	  }

     }     
}

void Substructure_MatrixXcm_MPI( const MPI_Comm Comm, PMatrixVector_t *const Mat,
				 const CouplingNode_t *const CNodes, PMatrixVector_t *const MatXcm )
{


     int Length, Accumulated_Length;
     int icoup;      /* Counter for the coupling nodes */
     int jcoup;
     int incx, incy;
     int PosXcm_Rows, PosXcm_Cols;
     int PosX_Rows, PosX_Cols;
     int rank;

     MPI_Comm_rank( Comm, &rank );

     incx = Mat->GlobalSize.Row;
     incy = 1;          /* The values in the Xm matrix are stored in columns. Therefore the stride has to be 1
			 * because of the way FORTRAN handles arrays (major column ordering) */

     /* Since the matrix is symmetric and only a part of it is stored, this routine has to be splitted into
      * two parts. The first will copy the elements above the coupling node, while the second will focus on
      * the other part */

     /* Copy until the first coupling node */
     Length = CNodes->Array[0] - 1;
     Accumulated_Length = 0;
     for( icoup = 0; icoup < CNodes->Order; icoup++ ){
	  PosX_Rows = CNodes->Array[icoup];     /* 1 based index */
	  PosX_Cols = 1; /* 1 based index */
	  PosXcm_Rows = 1;
	  PosXcm_Cols = icoup + 1;
	  hysl_pcopy( &Length, Mat->Array, &PosX_Rows, &PosX_Cols, Mat->Desc, &incx, MatXcm->Array, &PosXcm_Rows,
		   &PosXcm_Cols, MatXcm->Desc, &incy );
     }
     Accumulated_Length = Accumulated_Length + Length;
     for( icoup = 1; icoup < CNodes->Order; icoup++ ){
	  Length = CNodes->Array[icoup] - CNodes->Array[icoup - 1] - 1;
	  if ( Length > 0 ){
	       for( jcoup = icoup; jcoup < CNodes->Order; jcoup++ ){
		    
		    PosX_Rows = CNodes->Array[jcoup];     /* 1 based index */
		    PosX_Cols = CNodes->Array[icoup-1] + 1; /* 1 based index */
		    PosXcm_Rows = Accumulated_Length + 1;
		    PosXcm_Cols = jcoup + 1;		    
		    hysl_pcopy( &Length, Mat->Array, &PosX_Rows, &PosX_Cols, Mat->Desc, &incx, MatXcm->Array,
			     &PosXcm_Rows, &PosXcm_Cols, MatXcm->Desc, &incy );
	       }
	       Accumulated_Length = Accumulated_Length + Length;
	  }
	  
     }

     incx = 1;
     Length = Mat->GlobalSize.Row - CNodes->Array[CNodes->Order - 1];
     if( Length > 0 ){
	  for ( icoup = CNodes->Order - 1; icoup >= 0; icoup = icoup -1 ){
	       PosX_Rows = CNodes->Array[CNodes->Order - 1] + 1;
	       PosX_Cols = CNodes->Array[icoup];
	       PosXcm_Rows = MatXcm->GlobalSize.Row - Length + 1;
	       PosXcm_Cols = icoup + 1;
	       hysl_pcopy( &Length, Mat->Array, &PosX_Rows, &PosX_Cols, Mat->Desc, &incx, MatXcm->Array,
			&PosXcm_Rows, &PosXcm_Cols, MatXcm->Desc, &incy );
	  }
     }

     Accumulated_Length = Length;
     for( icoup = CNodes->Order -2; icoup >= 0; icoup = icoup -1 ){
	  Length = CNodes->Array[icoup + 1] - CNodes->Array[icoup] - 1;
	  if( Length > 0 ){
	       for( jcoup = icoup; jcoup >= 0; jcoup = jcoup - 1){
		    PosX_Rows = CNodes->Array[icoup] + 1;
		    PosX_Cols = CNodes->Array[jcoup];
		    PosXcm_Rows = MatXcm->GlobalSize.Row - Accumulated_Length - Length + 1;
		    PosXcm_Cols = jcoup + 1;	           
		    hysl_pcopy( &Length, Mat->Array, &PosX_Rows, &PosX_Cols, Mat->Desc, &incx, MatXcm->Array,
			     &PosXcm_Rows, &PosXcm_Cols, MatXcm->Desc, &incy );
	       }
	       Accumulated_Length = Accumulated_Length + Length;
	  }
     }
}

void Substructure_VectorXm_MPI( PMatrixVector_t *const VectorX, const CouplingNode_t *const CNodes,
				PMatrixVector_t *const VectorXm )
{

	int incx, incy;
	int Length;
	int PosX_Row, PosX_Col, PosXm_Row, PosXm_Col;
	int icoup;    /* Counter for the coupling nodes */

	incx = 1; incy = 1;
	PosX_Col = 1;
	PosXm_Col = 1;

	PosX_Row = 1;
	PosXm_Row = 1;
	for( icoup = 0; icoup < CNodes->Order; icoup++ ){
	     Length = CNodes->Array[icoup] - PosX_Row;

	     /* Copy the part of the vector between two positions */
	     hysl_pcopy( &Length, VectorX->Array, &PosX_Row, &PosX_Col, VectorX->Desc, &incx, VectorXm->Array,
		      &PosXm_Row, &PosXm_Col, VectorXm->Desc, &incy );

	     /* Update the values of the position in the vectors */
	     PosX_Row = CNodes->Array[icoup] + 1;
	     PosXm_Row = PosXm_Row + Length;
	}

	/* Copy the elements from the last position until the end of the vector */
	Length = VectorX->GlobalSize.Row - CNodes->Array[CNodes->Order-1];
	hysl_pcopy( &Length, VectorX->Array, &PosX_Row, &PosX_Col, VectorX->Desc, &incx, VectorXm->Array,
		 &PosXm_Row, &PosXm_Col, VectorXm->Desc, &incy );
}


void Substructure_VectorXc_MPI( const MPI_Comm Comm, PMatrixVector_t *const VecX,
				const CouplingNode_t *const CNodes, MatrixVector_t *const VecXc )
{
     int icoup;
     int nprow, npcol, myrow, mycol;
     int rank;
     int GRowIndex, GColIndex;
     int LRowIndex, LColIndex;
     int RowProcess, ColProcess;

     MPI_Status status;

     MPI_Comm_rank( Comm, &rank );

     /* Get grid info */
     Cblacs_gridinfo( VecX->Desc[1], &nprow, &npcol, &myrow, &mycol );
     
     GColIndex = 1;

#pragma omp parallel for
     for( icoup = 0; icoup < CNodes->Order; icoup++ ){

	  GRowIndex = CNodes->Array[icoup];   /* 1 based */
	  /* Given the global index of an element (GRowIndex, GColIndex) returns the local index of the
	     element (LRowIndex, LColIndex) and the coordinates of the process (Row Process, ColProcess) */
	  infog2l_( &GRowIndex, &GColIndex, VecX->Desc, &nprow, &npcol, &myrow, &mycol, &LRowIndex, &LColIndex,
		    &RowProcess, &ColProcess );

	  if ( Cblacs_pnum( VecX->Desc[1], RowProcess, ColProcess ) == 0 && myrow == RowProcess
	       && mycol == ColProcess ){ 
	       VecXc->Array[icoup] = VecX->Array[(LRowIndex - 1)*VecX->LocalSize.Col + (LColIndex - 1)];    
	  } else {
	       if ( myrow == RowProcess && mycol == ColProcess ){
			 
		    MPI_Send( &VecX->Array[(LRowIndex - 1)*VecX->LocalSize.Col + (LColIndex - 1)], 1,
			      MPI_HYSL_FLOAT, 0, VECTOR_XC, Comm );
	       }
	       
	       if ( rank == 0 ){
		    /* Store the diagonal elements */
		    MPI_Recv( &VecXc->Array[icoup], 1, MPI_HYSL_FLOAT,
			      Cblacs_pnum( VecX->Desc[1], RowProcess, ColProcess ), VECTOR_XC, Comm, &status );
	       }
	  }	  
     }     
}
