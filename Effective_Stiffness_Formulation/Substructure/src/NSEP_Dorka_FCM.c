#include <stdio.h>       /* For printf() and fprintf() */
#include <stdlib.h>      /* For atoi() and exit( ) */
#include <sys/socket.h>  /* For socket(), bind() and connect() */
#include <arpa/inet.h>   /* For sockaddr_in and inet_ntoa() */
#include <string.h>      /* For memset() */
#include <unistd.h>      /* For close() */
#include <sys/time.h>    /* For gettimeoffday()  */
#include <getopt.h>      /* For getopt_long() */

#include "NSEP_Communication_Sync.h"    /* Prototypes of several functions involved in the CGM */
#include "NSEP_Definitions.h"  /* Definition of various NSEP constants */
#include "Send_Receive_Data.h" /* Send and receive data routines */
#include "ErrorHandling.h"     /* Error handling routines */

#include "Substructure.h"

#if ADWIN_
#include "RoutinesADwin.h"
#endif

void Print_Help( const char *Program_Name );

int main ( int argc, char **argv )
{

     int i, j;            /* Counters */

     /* Variables for the sub-structure testing/simulation */
     int Is_Not_Finished;

     ConstSub Cnst;

     float *Gc;
     float *u0c, *uc;

     float *fcprev, *fc;
     float *Send, *Recv;

#if ADWIN_
     /* Array where the data from ADwin will be stored */
     float *ADWIN_DATA;
#endif

     TMD_Sim Num_TMD;

     /* Variables to deal with arguments */
     int Mode, Selected_Option;

     struct option long_options[] = {
	  {"help", no_argument, 0, 'h'},
	  {"mode", required_argument, 0, 'm'},
	  {0, 0, 0, 0}
     };

     /* Set the default value for Port and Mode before the user input. */
     Mode = 1;     /* Simulate the sub-structure using an exact solution. */

    /* This is only used if there are no arguments */
     if ( argc == 1 ){	  
	  printf("Defaulting on mode 1.\n");
     }

     /* Assign each argument to the correct variable */
     while( (Selected_Option = getopt_long( argc, argv, "m:h", long_options, NULL )) != -1 ){
	  switch( Selected_Option ){
	  case 'm':
	       Mode = atoi( optarg );

	       if ( Mode < 0 || Mode > 2 ){
		    fprintf( stderr, "Mode %d is not a valid mode value.\n", Mode );
		    Print_Help( argv[0] );
		    return EXIT_FAILURE;
	       }
	       break;
	  case 'h':
	       Print_Help( argv[0] );
	       return EXIT_FAILURE;
	       break;
	  case '?':
	       /* Long options already prints an error message telling that there is an unrecognised option */
	       Print_Help( argv[0] );
	       return EXIT_FAILURE;
	  case ':':
	       /* Long options already prints an error message telling that the option requieres an argument */
	       Print_Help( argv[0] );
	       return EXIT_FAILURE;
	  }
     }

     /* Dynamically allocate memory */
     Gc = calloc( Cnst.Order_Couple*Cnst.Order_Couple, sizeof( float ) );
 
     u0c = calloc( Cnst.Order_Couple, sizeof( float ) );
     uc = calloc( Cnst.Order_Couple, sizeof( float ) );

     fcprev = calloc( Cnst.Order_Couple, sizeof( float ) );
     fc = calloc( Cnst.Order_Couple, sizeof( float ) );

     Send = calloc( 3*Cnst.Order_Couple, sizeof( float ) );
     Recv = calloc( Cnst.Order_Couple, sizeof( float ) );

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

     if ( Mode == USE_ADWIN ){
#if ADWIN_
	  /* Run with ADwin */
	  ADWIN_SetGc( Gc, Cnst.Order_Couple*Cnst.Order_Couple );
	  printf( "Using ADwin to perform the sub-stepping process.\n" );
	  ADWIN_DATA = calloc( Cnst.Num_Sub*Cnst.Num_Steps*NUM_CHANNELS, sizeof( float ) );
#else
	  fprintf(stderr, "The program was not compiled with ADwin support.\n");
	  exit( EXIT_FAILURE );
#endif
     } else if ( Mode == USE_EXACT ){
	  /* Simulate the substructure numerically */
	  printf( "Simulating the sub-structure using an exact integration method.\n");
	  ExactSolution_Init( 285, 352.18177, 68000, Cnst.DeltaT_Sub, &Num_TMD );
     } else {
	  printf( "Simulating the sub-structure using measured values as an input.\n");
	  /* Do nothing for the moment */
     }

     Is_Not_Finished = 1;
     while ( Is_Not_Finished ){

	  /* Receive the displacement from the CGM facility */
	  Communicate_With_PNSE( 3, 0.0, Send, Recv, Cnst.Order_Couple );
	  for ( j = 0; j < Cnst.Order_Couple; j++ ){
	       u0c[j] =  Recv[j];
	  }

	  /* Perform the substepping process */
	  if ( Mode == USE_ADWIN ){
#if ADWIN_
	       /* Run using ADwin */
	       ADWIN_Substep( u0c, uc, fcprev, fc, Cnst.Order_Couple, Cnst.Num_Sub, Cnst.DeltaT_Sub );
#endif
	  } else if ( Mode == USE_EXACT ){
	       /* Run without ADwin and simulating the substructure using an exact
		* solution.
		*/
	       Simulate_Substructure( &Num_TMD, Gc, u0c, uc, fcprev, fc, Cnst.Order_Couple, Cnst.Num_Sub, Cnst.DeltaT_Sub );
	  } else {
	       /* Run without ADwin and simulating the substructure using measured
		* values of the coupling force.
		*/
	       Simulate_Substructure_Measured_Values( "fc.txt", Gc, u0c, uc, fcprev, fc, Cnst.Order_Couple, Cnst.Num_Sub, Cnst.DeltaT_Sub );
	  }

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


     if ( Mode == USE_ADWIN ){
#if ADWIN_
	  /* Get the Data from ADwin */
	  printf("Getting the data from ADwin...");
	  GetDataADwin( Cnst.Num_Steps, Cnst.Num_Sub, ADWIN_DATA );
	  printf(" DONE!\n");
     
	  free( ADWIN_DATA );
#endif
     } else {
	  printf("The simulatiovn has finished\n");
     }

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


void Print_Help( const char *Program_Name )
{

     fprintf( stderr, "Usage: %s [-h] -m <Mode>", Program_Name );
     fprintf( stderr,
	      "  -h  --help    This help text.\n"
	      "  -m  --mode    The mode used by the program. Default value 1.\n"
	      "                  0 - ADwin will be used to perform the sub-stepping process.\n"
	      "                  1 - The Substructure will be simulated using an exact solution.\n"
	      "                  2 - The Substructure will be simulated using measured values.\n" );
}
