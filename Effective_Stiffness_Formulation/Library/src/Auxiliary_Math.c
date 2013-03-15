#include <stdio.h>              /* For printf(), fprintf() */
#include <stdlib.h>             /* For exit() */

#include "Auxiliary_Math.h"
#include "Print_Messages.h"

int Max ( const int a, const int b )
{
     if ( a >= b ){
	  return a;
     } else return b;
}

int Min ( const int a, const int b )
{
     if ( a <= b ){
	  return a;
     } else return b;
}

MatrixVector_t Generate_IdentityMatrix( int Rows, int Cols )
{
     MatrixVector_t Identity;
     unsigned short int i;

     if( Rows != Cols ){
	  Print_Header( ERROR );
	  fprintf( stderr, "Generate_IdentityMatrix: The number of rows and columns must be the same.\n" );
	  exit( EXIT_FAILURE );
     }

     MatrixVector_Create( Rows, Cols, &Identity );
     
     for( i = 0; i < Rows; i++ ){
	  Identity.Array[i + Rows*i] = 1.0;
     }

     return Identity;
}

unsigned int MatrixVector_ReturnIndex_UPS( const unsigned int RowIndex, const unsigned int ColIndex, const int n )
{
     unsigned int Index;

     if( RowIndex >= ColIndex ){
	  Index = RowIndex + (2*(unsigned int) n - ColIndex)*(ColIndex - 1)/2 - 1;
     } else {
	  Index = ColIndex + (2*(unsigned int) n - RowIndex)*(RowIndex - 1)/2 - 1;
     }
     return Index;
}

unsigned int MatrixVector_ReturnIndex_LPS( const unsigned int RowIndex, const unsigned int ColIndex )
{
     unsigned int Index;

     if( ColIndex >= RowIndex ){
	  Index = RowIndex + ColIndex*(ColIndex - 1)/2 - 1;
     } else {
	  Index = ColIndex + RowIndex*(RowIndex - 1)/2 - 1;
     }
     return Index;
}
