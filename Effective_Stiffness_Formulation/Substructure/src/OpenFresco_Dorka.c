#include <stdio.h>
#include <stdlib.h>       /* For exit( ) and atoi( ) */
#include <string.h>       /* For memset( ) */

#include <getopt.h>      /* For getopt_long() */
#include "Substructure.h" /* For NUM_CHANNELS, ConstSub and Init_Constants_Substructure( ) */
#include "OPSocket.h"

#if ADWIN_
#include "RoutinesADwin.h"
#endif


void Print_Help( const char *Program_Name );

int main ( int argc, char **argv )
{
     int i;
     /* Variables concerning data communication */
     int Server_Socket;      /* Socket for the server */
     unsigned int Port;      /* Port */

     int DataTypeSize;       /* Size of the data type to be transfered */
     char *Data;             /* To Send and receive data */
     int Length;             /* Length of the data to be transfered */
     int ierr;               /* Error. OpenFresco routines */
     int Is_Not_Finished;    /* To check if the process is finished or not */

     /* Variables required by OpenFresco */
     int iData[11];

     ConstSub Cnst;

     float *Gc;
     float *u0c, *uc;
     float *fcprev, *fc;
     double *Send, *Recv;

     TMD_Sim Num_TMD;
     UHYDE_Sim Num_UHYDE;

#if ADWIN_
     /* Array where the data from ADwin will be stored */
     float *ADWIN_DATA;
#endif

     /* Variables to deal with arguments */
     int Mode, Selected_Option;

     struct option long_options[] = {
	  {"help", no_argument, 0, 'h'},
	  {"mode", required_argument, 0, 'm'},
	  {"server-port", required_argument, 0, 'p'},
	  {0, 0, 0, 0}
     };

     /* Set the default value for Port and Mode before the user input. */
     Port = 3333;  /* Default port */
     Mode = 1;     /* Simulate the sub-structure using an exact solution. */

    /* This is only used if there are no arguments */
     if ( argc == 1 ){	  
	  printf("Defaulting on mode 1.\n");
     }

     /* Assign each argument to the correct variable */
     while( (Selected_Option = getopt_long( argc, argv, "m:p:h", long_options, NULL )) != -1 ){
	  switch( Selected_Option ){
	  case 'm':
	       Mode = atoi( optarg );

	       if ( Mode < 0 || Mode > 3 ){
		    fprintf( stderr, "Mode %d is not a valid mode value.\n", Mode );
		    Print_Help( argv[0] );
		    return EXIT_FAILURE;
	       }
	       break;
	  case 'p':
	       Port = (unsigned int) atoi( optarg );
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


     /* Create a TCP/IP socket for the server using OpenFresco routine */
     printf( "Waiting for the ECGeneric experimental facility to connect.\n" );
     setupconnectionserver( &Port, &Server_Socket );
     if ( Server_Socket >= 0 ){
	  printf("Connection successfully established\n" );
     } else {
	  exit( EXIT_FAILURE );
     }

     /* Receive an ID from OpenFresco with sizeCtrl, sizeDaq and dataSize */
     Data = (char *) iData;
     DataTypeSize = sizeof( int );
     Length = 11;
     recvdata( &Server_Socket, &DataTypeSize, Data, &Length, &ierr );
     printf("Received ID vector: \n" );

     for ( i = 0; i < 11; i++ ){
	  printf("%d\t",iData[i] );
     }
     printf("\n");

     /* Initialise the constants of the substructure */
     Init_Constants_Substructure( &Cnst );

     Gc = calloc( Cnst.Order_Couple*Cnst.Order_Couple, sizeof( float ) );

     u0c = calloc( Cnst.Order_Couple, sizeof( float ) );
     uc = calloc( Cnst.Order_Couple, sizeof( float ) );

     fcprev = calloc( Cnst.Order_Couple, sizeof( float ) );
     fc = calloc( Cnst.Order_Couple, sizeof( float ) );

     /* The size of the data to be exchanged is given by the last element of iData */
     Length = iData[10];
     printf("%d\n", Length);
     Send = (double*) calloc( Length, sizeof(double) );
     Recv = (double*) calloc( Length, sizeof(double) );   
     

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
     } else if ( Mode == USE_UHYDE ){
	  printf( "Simulating the friction device UHYDE-fbr.\n" );
	  Simulate_UHYDE_1D_Init( 500.0, 1.0, 9.0, &Num_UHYDE );
     } else {
	  printf( "Simulating the sub-structure using measured values as an input.\n");
	  /* Do nothing for the moment */
     }

     /* The process is not finished. */
     Is_Not_Finished = 1;  
     DataTypeSize = sizeof( double );
     while( Is_Not_Finished ){
     
	  Data = (char *) Recv;
	  recvdata( &Server_Socket, &DataTypeSize, Data, &Length, &ierr );

	  if ( Recv[0] == 2.0 ){
	       /* TODO */
	  } else if ( Recv[0] == 3.0 ){
	       /* Check if what has been received are the control values */
	       for ( i = 0; i < Cnst.Order_Couple; i++ ){
		    printf("Received control values values\n");
		    u0c[i] = (float) Recv[1+i];
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
		    Simulate_Substructure( &Num_TMD, Mode, Gc, u0c, uc, fcprev, fc, Cnst.Order_Couple, Cnst.Num_Sub, Cnst.DeltaT_Sub );
	       } else if ( Mode == USE_UHYDE ){
		    /*
		     * Simulate the UHYDE-fbr without ADwin
		     */
		    Simulate_Substructure( &Num_UHYDE, Mode, Gc, u0c, uc, fcprev, fc, Cnst.Order_Couple, Cnst.Num_Sub, Cnst.DeltaT_Sub );
	       } else {
		    /* Run without ADwin and simulating the substructure using measured
		     * values of the coupling force.
		     */
		    Simulate_Substructure_Measured_Values( "fc.txt", Gc, u0c, uc, fcprev, fc, Cnst.Order_Couple, Cnst.Num_Sub, Cnst.DeltaT_Sub );
	       }

	  } else if ( Recv[0] == 6.0 ){
	       printf("Received query to send DAQ values\n");
	       /* Compose the data */
	       for ( i = 0; i < Cnst.Order_Couple; i++ ){
		    Send[i] = (double) uc[i];
		    Send[i + Cnst.Order_Couple] = (double) fcprev[i];
		    Send[i + 2*Cnst.Order_Couple] = (double) fc[i];
	       }
	       Data = (char *) Send;
	       senddata( &Server_Socket, &DataTypeSize, Data, &Length, &ierr );
	  } else if ( Recv[0] == 99.0 ){
	       Is_Not_Finished = 0;
	       /* End the connection with OpenFresco */
	       closeconnection( &Server_Socket, &ierr );
	  }
     }


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

     fprintf( stderr, "Usage: %s [-h] -m <Mode> -p <Port>", Program_Name );
     fprintf( stderr,
	      "  -h  --help    This help text.\n"
	      "  -m  --mode    The mode used by the program. Default value 1.\n"
	      "                  0 - ADwin will be used to perform the sub-stepping process.\n"
	      "                  1 - The Substructure will be simulated using an exact solution.\n"
	      "                  2 - The Substructure will be simulated using measured values.\n"
	      "  -p  --port    Port used for communication with the client program. Default value 3333.\n" );
}
