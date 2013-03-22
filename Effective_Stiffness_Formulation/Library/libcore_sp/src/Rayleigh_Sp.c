#include <stdio.h>            /* For printf(), fprintf() */
#include <stdlib.h>           /* For exit() */

#include "MatrixVector.h"     /* MatrixVector definition */
#include "MatrixVector_Sp.h"
#include "Print_Messages.h"   /* For Print_Header() */
#include "Rayleigh.h"         /* Rayleigh damping routines */

#include <mkl_blas.h>
#include <mkl_spblas.h>

void Rayleigh_Damping_Sp( const MatrixVector_Sp_t *const Mass, const MatrixVector_Sp_t *const Stif, MatrixVector_Sp_t *const Damp,
			  const Rayleigh_t *const Rayleigh )
{
     MatrixVector_Sp_t Temp;  /* Temporal matrix */
     int  i;                  /* A counter */
     int Length;
     int incx, incy;          /* Stride in the operations */
     double alpha, beta;      /* Constants */
     char trans;
     int job, sort, info;     /* MKL variables */

     alpha = Rayleigh->Alpha;
     beta = Rayleigh->Beta;

     MatrixVector_Create_Sp( Mass->Rows, Mass->Cols, Mass->Num_Nonzero, &Temp );

     incx = 1; incy = 1;
     Length = Temp.Num_Nonzero;
     dcopy( &Length, Mass->Values, &incx, Temp.Values, &incy );

     /* Copy the column array */
#pragma omp parallel for
     for (i = 0; i < Length; i++ ){
	  Temp.Columns[i] = Mass->Columns[i];
     }

     /* Scal the Values array */
     dscal( &Length, &alpha, Temp.Values, &incx );

     /* Copy the RowIndex array */
     Length = Temp.Rows + 1;
#pragma omp parallel for
     for (i = 0; i < Length; i++ ){
	  Temp.RowIndex[i] = Mass->RowIndex[i];
     }

     trans = 'N';  /* The operation C = Temp + beta*B is performed */
     job = 0;      /* The routine computes the addition */
     sort = 0;     /* The routine does not perform any reordering */
     mkl_dcsradd( &trans, &job, &sort, &Temp.Rows, &Temp.Cols, Temp.Values, Temp.Columns, Temp.RowIndex,
		  &beta, Stif->Values, Stif->Columns, Stif->RowIndex,
		  Damp->Values, Damp->Columns, Damp->RowIndex, &Damp->Num_Nonzero, &info );

     /* Delete the previously allocated sparse matrix */
     MatrixVector_Destroy_Sp( &Temp );

     if ( info > 0){
	  Print_Header( ERROR );
	  fprintf( stderr, "Number of elements exceeded while calculating the Damping matrix.\n" );
	  exit( EXIT_FAILURE );
     } else if ( info < 0 ){
	  Print_Header( ERROR);
	  fprintf( stderr, "I do not understand.\n" );
	  exit( EXIT_FAILURE );
     } else {
	  Print_Header( SUCCESS );
	  printf( "Damping matrix successfully calculated.\n" );
     }
}
