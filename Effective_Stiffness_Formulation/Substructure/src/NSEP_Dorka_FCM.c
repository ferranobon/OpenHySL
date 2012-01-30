#include <stdio.h>       /* For printf() and fprintf() */
#include <stdlib.h>      /* For atoi() and exit( ) */
#include <sys/socket.h>  /* For socket(), bind() and connect() */
#include <arpa/inet.h>   /* For sockaddr_in and inet_ntoa() */
#include <string.h>      /* For memset() */
#include <unistd.h>      /* For close() */
#include <sys/time.h>    /* For gettimeoffday()  */

#include "NSEP_Communication_Sync.h"    /* Prototypes of several functions involved in the CGM */
#include "NSEP_Definitions.h"  /* Definition of various NSEP constants */
#include "Send_Receive_Data.h" /* Send and receive data routines */
#include "ErrorHandling.h"     /* Error handling routines */
#include "RoutinesADwin.h"     /* Routines to communicate with ADwin */


int main ( int argc, char **argv )
{

     int i, j;                                /* Counters */
     int Is_Not_Finished;

     int OrderCouple = 1;
     const int Num_Sub = 4;
     const float DeltaT_Sub = 0.01/(float)Num_Sub;

     float *Gc;
     float *u0c, *uc;

     float *fcprev, *fc;
     float Send[3*OrderCouple], Recv[OrderCouple];

     float *ADWIN_DATA;

     Send[0] = 0.0;
     Send[1] = 0.0;
     Send[2] = 0.0;

     Recv[0] = 0.0;

     Gc = calloc( OrderCouple*OrderCouple, sizeof( float ) );
 
     u0c = calloc( OrderCouple, sizeof( float ) );
     uc = calloc( OrderCouple, sizeof( float ) );

     fcprev = calloc( OrderCouple, sizeof( float ) );
     fc = calloc( OrderCouple, sizeof( float ) );

     ADWIN_DATA = calloc( 4*4096*22, sizeof(float));

     /* Using NSEP */
     /* Open the Socket */
     printf( "Establishing connection with PNSE server.\n" );
     Communicate_With_PNSE( 0, 0.0, Send, Recv, OrderCouple );

     /* Receive the matrix Gc from the CGM facility */
     for ( i = 0; i < OrderCouple; i++ ){
	  Communicate_With_PNSE( 3, 0.0, Send, Recv, OrderCouple*OrderCouple );
	  for ( j = 0; j < OrderCouple; j++ ){
	       Gc[i*OrderCouple + j] = Recv[j];
	  }
	  /* This is in done so that the PNSE don't overtake the first step */
	  Communicate_With_PNSE( 4, 0.0, Send, Recv, 3*OrderCouple );
     }
     
     /* Set Gc to ADwin */
     ADWIN_SetGc( Gc, OrderCouple*OrderCouple );

     for( i = 0; i< OrderCouple*OrderCouple; i++ ){
	  printf("Gc = %e\t", Gc[i]);
     }
     printf("\n");

     Is_Not_Finished = 1;
     while ( Is_Not_Finished ){

	  /* Receive the displacement from the CGM facility */
	  Communicate_With_PNSE( 3, 0.0, Send, Recv, OrderCouple );
	  for ( j = 0; j < OrderCouple; j++ ){
	       u0c[j] =  Recv[j];
	  }

	  /* Perform the substepping process */
	  ADWIN_Substep( u0c, uc, fcprev, fc, OrderCouple, Num_Sub, DeltaT_Sub );

	  /* Compose the data to send */
	  for (i = 0; i < OrderCouple; i++) {
	       Send[i] = uc[i];
	       Send[i+OrderCouple] = fcprev[i];
	       Send[i+2*OrderCouple] = fc[i];
	  }
	  /* Send the NSEP_CSIG package. WhatToDo = 4 */
	  Communicate_With_PNSE( 4, 0.0, Send, Recv, 3*OrderCouple );
          /* Receive the state of the experiment. WhatToDo = 5 */
	  Communicate_With_PNSE( 5, 0.0, Send, Recv, OrderCouple );
	  if ( Recv[0] == -9999.0 ){
	       Is_Not_Finished = 0;
	  }
     }

     /* Say to PNSE server that the process has finished. WhatToDo = 6 */
     Communicate_With_PNSE( 6, 0.0, Send, Recv, 0 );

     GetDataADwin( 4096, Num_Sub, ADWIN_DATA );

     /* Free the dinamically allocated memory */
     free( ADWIN_DATA );

     free( Gc );

     free( u0c );
     free( uc );

     free( fcprev );
     free( fc );

     return 0;
}
