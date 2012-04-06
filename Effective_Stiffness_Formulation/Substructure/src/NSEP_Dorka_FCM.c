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
#include "Substructure.h"


int main ( int argc, char **argv )
{

     int i, j;                                /* Counters */
     int Is_Not_Finished;

     ConstSub Cnst;

     float *Gc;
     float *u0c, *uc;

     float *fcprev, *fc;
     float *Send, *Recv;

     float *ADWIN_DATA;

     /* Dynamically allocate memory */
     Gc = calloc( Cnst.Order_Couple*Cnst.Order_Couple, sizeof( float ) );
 
     u0c = calloc( Cnst.Order_Couple, sizeof( float ) );
     uc = calloc( Cnst.Order_Couple, sizeof( float ) );

     fcprev = calloc( Cnst.Order_Couple, sizeof( float ) );
     fc = calloc( Cnst.Order_Couple, sizeof( float ) );

     Send = calloc( 3*Cnst.Order_Couple, sizeof( float ) );
     Recv = calloc( Cnst.Order_Couple, sizeof( float ) );

#if SIMULATE_SUB_
     /* Do nothing */
#else
     ADWIN_DATA = calloc( Cnst.Num_Sub*Cnst.Num_Steps*NUM_CHANNELS, sizeof( float ) );
#endif

     /* Using NSEP */
     /* Open the Socket */
     printf( "Establishing connection with PNSE server.\n" );
     Communicate_With_PNSE( 0, 0.0, Send, Recv, Cnst.Order_Couple );

     /* Receive the matrix Gc from the CGM facility */
     for ( i = 0; i < Cnst.Order_Couple; i++ ){
	  Communicate_With_PNSE( 3, 0.0, Send, Recv, Cnst.Order_Couple*Cnst.Order_Couple );
	  for ( j = 0; j < Cnst.Order_Couple; j++ ){
	       Gc[i*Cnst.Order_Couple + j] = Recv[j];
	  }
	  /* This is in done so that the PNSE don't overtake the first step */
	  Communicate_With_PNSE( 4, 0.0, Send, Recv, 3*Cnst.Order_Couple );
     }
#if SIMULATE_SUB_  /* Run this without ADwin */
     /* Do nothing */
     printf("Simulating the substructure\n");
#else
     /* Set Gc to ADwin */
     ADWIN_SetGc( Gc, Cnst.Order_Couple*Cnst.Order_Couple );
#endif

     for( i = 0; i< Cnst.Order_Couple*Cnst.Order_Couple; i++ ){
	  printf("Gc = %e\t", Gc[i]);
     }
     printf("\n");

     Is_Not_Finished = 1;
     while ( Is_Not_Finished ){

	  /* Receive the displacement from the CGM facility */
	  Communicate_With_PNSE( 3, 0.0, Send, Recv, Cnst.Order_Couple );
	  for ( j = 0; j < Cnst.Order_Couple; j++ ){
	       u0c[j] =  Recv[j];
	  }

#if SIMULATE_SUB_  /* Run this without ADwin */
	  Simulate_Substructure( u0c, uc, fcprev, fc, Cnst.Order_Couple, Cnst.Num_Sub, Cnst.DeltaT_Sub );
#else              /* Run using ADwin */
	  /* Perform the substepping process */
	  ADWIN_Substep( u0c, uc, fcprev, fc, Cnst.Order_Couple, Cnst.Num_Sub, Cnst.DeltaT_Sub );
#endif	  

	  /* Compose the data to send */
	  for (i = 0; i < Cnst.Order_Couple; i++) {
	       Send[i] = uc[i];
	       Send[i+Cnst.Order_Couple] = fcprev[i];
	       Send[i+2*Cnst.Order_Couple] = fc[i];
	  }
	  /* Send the NSEP_CSIG package. WhatToDo = 4 */
	  Communicate_With_PNSE( 4, 0.0, Send, Recv, 3*Cnst.Order_Couple );
          /* Receive the state of the experiment. WhatToDo = 5 */
	  Communicate_With_PNSE( 5, 0.0, Send, Recv, Cnst.Order_Couple );
	  if ( Recv[0] == -9999.0 ){
	       Is_Not_Finished = 0;
	  }
     }

     /* Say to PNSE server that the process has finished. WhatToDo = 6 */
     Communicate_With_PNSE( 6, 0.0, Send, Recv, 0 );

#if SIMULATE_SUB_  /* Run this without ADwin */
     printf("The simulatiovn has finished\n");
#else
     /* Get the Data from ADwin */
     printf("Getting the data from ADwin...");
     GetDataADwin( Cnst.Num_Steps, Cnst.Num_Sub, ADWIN_DATA );
     printf(" DONE!\n");

     free( ADWIN_DATA );
#endif

     /* Free the dinamically allocated memory */     

     free( Gc );

     free( u0c );
     free( uc );

     free( fcprev );
     free( fc );

     free( Send );
     free( Recv );

     return 0;
}
