#include <stdio.h>
#include <stdlib.h>

int main( int argc, char **argv )
{

     float fone = 1.0;
     int i, Max = 22200;
     FILE *TheFile;

     TheFile = fopen( "22200LV.txt", "w" );
     for( i = 0; i < Max; i++ ){
	  fprintf( TheFile, "%f\n", fone );
     }

     fclose( TheFile );

     return 0;
}
