#include <stdio.h>
#include <stdlib.h>      /* For atoi( ) and exit() */
#include <string.h>      /* For memset( ) */
#include <getopt.h>      /* For getopt_long() */

#include <netdb.h>

#include "ErrorHandling.h"
#include "Substructure.h"
#include "Send_Receive_Data.h"
#include "PIndependent.h"

#if ADWIN_
#include "RoutinesADwin.h"
#endif

int main( int argc, char **argv )
{
     unsigned int i,j;

     int Server_Socket;                /* Socket for the server */
     int Client_Socket;                /* Socket for the client */

     int Socket_Type;
     const char *Port;
     struct sockaddr_storage Client_Addr;

     int Is_Not_Finished;
     unsigned int Length;

     ConstSub Cnst;

     float *Gc;
     float *u0c0, *u0c, *uc;
     float *fcprev, *fc;
     float *Send;

     TMD_Sim *Num_TMD;
     UHYDE_Sim *Num_UHYDE;


     /* Array where the data from ADwin will be stored */
     float *ADWIN_DATA = NULL;

     /* Variables to deal with arguments */
     int Mode, Selected_Option;

     struct option long_options[] = {
	  {"help", no_argument, 0, 'h'},
	  {"mode", required_argument, 0, 'm'},
	  {"server-port", required_argument, 0, 'p'},
	  {"socket-type", required_argument, 0, 's'},
	  {0, 0, 0, 0}
     };

     /* Initialise the constants of the substructure */
     Init_Constants_Substructure( &Cnst, "ConfFile.conf" );

     /* Set the default value for Port, socket type and Mode before the user input. */
     Port = "3333";  /* Default port */
     Mode = 1;       /* Simulate the sub-structure using an exact solution. */
     Socket_Type = PROTOCOL_TCP;

    /* This is only used if there are no arguments */
     if ( argc == 1 ){
	  printf("Defaulting to TCP communication protocol.\n");	  
	  printf("Defaulting on mode 1.\n");
     }

     /* Assign each argument to the correct variable */
     while( (Selected_Option = getopt_long( argc, argv, "m:p:s:h", long_options, NULL )) != -1 ){
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
	       Port = optarg;
	       break;
	  case 's':
	       Socket_Type = atoi( optarg );
	       
	       if ( Socket_Type < 1 || Socket_Type > 4 ){
		    fprintf( stderr, "Socket type %d is not a valid.\n", Socket_Type );
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

     /* Accept the connection from the client in case of TCP sockets */
     if ( Socket_Type == PROTOCOL_TCP || Socket_Type == PROTOCOL_UDP ){
	  /* Create a TCP/IP socket for the server */
	  Server_Socket = Setup_Server_Socket( Port, Socket_Type );
	  if( Server_Socket < -1 ){
	       PrintErrorAndExit( "Setup_Server_Socket() failed" );	  
	  } else {
	       if( Socket_Type == PROTOCOL_TCP ){
		    printf( "TCP Server created.\nListening on port %s for incoming client connections...\n", Port );
		    Client_Socket = Accept_TCP_Client_Connection( Server_Socket );
	       } else {
		    printf( "UDP Server created.\nListening on port %s for incoming client connections...\n", Port );
	       }
	  }
     } else if ( Socket_Type == PROTOCOL_NSEP ){
	  /* Using NSEP */
	  /* Open the Socket */
	  printf( "Establishing connection with PNSE server.\n" );
	  Communicate_With_PNSE( 0, 0.0, Send, Recv, Cnst.Order_Couple );
     } else if ( Socket_Type == PROTOCOL_OF ){
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
     }

     /* Dynamically allocate memory */
     Gc = (float *) calloc( (size_t) Cnst.Order_Couple*Cnst.Order_Couple, sizeof( float ) );
     u0c0 = (float *) calloc( (size_t) Cnst.Order_Couple, sizeof( float ) );
     u0c = (float *) calloc( (size_t) Cnst.Order_Couple, sizeof( float ) );
     uc = (float *) calloc( (size_t) Cnst.Order_Couple, sizeof( float ) );

     fcprev = (float *) calloc( (size_t) Cnst.Order_Couple, sizeof( float ) );
     fc = (float *) calloc( (size_t) Cnst.Order_Couple, sizeof( float ) );

     Send = (float *) calloc( (size_t) 3*Cnst.Order_Couple, sizeof( float ) );

     /* Receive matrix Gc */
     if( Socket_Type == PROTOCOL_TCP ){
	  Length = Cnst.Order_Couple*Cnst.Order_Couple;
	  Receive_Data( Gc, Length, Client_Socket );
	  for( i = 0; i < Cnst.Order_Couple; i++ ){
	       for( j = 0; j < Cnst.Order_Couple; j++ ){
		    printf("%e\t", Gc[i*Cnst.Order_Couple + j] );
	       }
	       printf("\n");
	  }

     } else if ( Socket_Type == PROTOCOL_UDP ) {
	  Length = Cnst.Order_Couple*Cnst.Order_Couple;
	  Client_AddrLen = sizeof(Client_Addr);
	  if ( recvfrom( Server_Socket, Gc, Length*sizeof(float), 0, (struct sockaddr *) &Client_Addr, &Client_AddrLen) < 0 ){
	       PrintErrorAndExit("Recvfrom failed" );
	  }
     } else if ( Socket_Type == PROTOCOL_NSEP ){
	      /* Receive the matrix Gc from the CGM facility */
	  for ( i = 0; i < Cnst.Order_Couple; i++ ){
	       Communicate_With_PNSE( 3, 0.0, Send, Recv, Cnst.Order_Couple*Cnst.Order_Couple );
	       for ( j = 0; j < Cnst.Order_Couple; j++ ){
		    Gc[i*Cnst.Order_Couple + j] = Recv[j];
	       }
	       /* This is in done so that the PNSE don't overtake the first step */
	       Communicate_With_PNSE( 4, 0.0, Send, Recv, 3*Cnst.Order_Couple );
	  }
     }

     if ( Mode == USE_ADWIN ){
#if ADWIN_
	  /* Run with ADwin */
	  ADWIN_SetGc( Gc, Cnst.Order_Couple*Cnst.Order_Couple );
	  printf( "Using ADwin to perform the sub-stepping process.\n" );
	  ADWIN_DATA = (float *) calloc( (size_t) Cnst.Num_Sub*Cnst.Num_Steps*NUM_CHANNELS, sizeof( float ) );
#else
	  fprintf(stderr, "The program was not compiled with ADwin support.\n");
	  exit( EXIT_FAILURE );
#endif
     } else if ( Mode == USE_EXACT ){
	  /* Simulate the substructure numerically */
	  printf( "Simulating the sub-structure using an exact integration method.\n");
	  Num_TMD = (TMD_Sim *) malloc( sizeof(TMD_Sim)*(size_t) Cnst.Order_Couple );
	  for( i = 0; i < Cnst.Order_Couple; i++ ){
	       ExactSolution_Init( 285.0f, 352.18177f, 68000.0f, Cnst.DeltaT_Sub, &Num_TMD[i] );
	  }
     } else if ( Mode == USE_UHYDE ){
	  printf( "Simulating the friction device UHYDE-fbr.\n" );
	  /* Allocate memory for the number of sub-structures */
	  Num_UHYDE = (UHYDE_Sim *) malloc( sizeof(UHYDE_Sim)*(size_t) Cnst.Order_Couple );
	  for( i = 0; i < Cnst.Order_Couple; i++ ){
	       Simulate_UHYDE_1D_Init( 0.0002f, 0.9f, 500.0f, &Num_UHYDE[i] );
	  }
     } else {
	  printf( "Simulating the sub-structure using measured values as an input.\n");
	  /* Do nothing for the moment */
     }

     Is_Not_Finished = 1;
     while( Is_Not_Finished ){

	  /* Receive the displacement */
	  Length = Cnst.Order_Couple;
	  if ( Socket_Type == PROTOCOL_TCP ){
	       Receive_Data( u0c, Length, Client_Socket );

	  } else if ( Socket_Type ==  PROTOCOL_UDP ){
	       Client_AddrLen = sizeof(Client_Addr);
	       if ( recvfrom( Server_Socket, u0c, Length*sizeof(float), 0, (struct sockaddr *) &Client_Addr, &Client_AddrLen) < 0 ){
		    PrintErrorAndExit("recvfrom() failed" );
	       }
	  } else if ( Socket_Type == PROTOCOL_NSEP ){
	       /* Receive the state of the experiment. WhatToDo = 5 */
	       Communicate_With_PNSE( 5, 0.0, Send, Recv, Cnst.Order_Couple );
	       u0c[0] = Recv[0];
	       if( u0c != -9999.0f ){
		    Communicate_With_PNSE( 3, 0.0, Send, Recv, Cnst.Order_Couple );
		    for ( j = 0; j < Cnst.Order_Couple; j++ ){
			 u0c[j] =  Recv[j];
		    }
	       }
	  }
  
	  if ( u0c[0] == -9999.0f ){
	       Is_Not_Finished = 0;
	  } else if ( ){
	       /* Perform the substepping process */
	       if ( Mode == USE_ADWIN ){
#if ADWIN_
		    /* Run using ADwin */
		    ADWIN_Substep( u0c, uc, fcprev, fc, Cnst.Order_Couple );
#endif
	       } else if ( Mode == USE_EXACT ){
		    /* Run without ADwin and simulating the substructure using an exact
		     * solution.
		     */
		    Simulate_Substructure( Num_TMD, Mode, Gc, u0c0, u0c, uc, fcprev, fc, Cnst.Order_Couple, Cnst.Num_Sub, Cnst.DeltaT_Sub );
	       } else if ( Mode == USE_UHYDE ){
		    /*
		     * Simulate the UHYDE-fbr without ADwin
		     */
		    Simulate_Substructure( Num_UHYDE, Mode, Gc, u0c0, u0c, uc, fcprev, fc, Cnst.Order_Couple, Cnst.Num_Sub, Cnst.DeltaT_Sub );
	       } else {
		    /* Run without ADwin and simulating the substructure using measured
		     * values of the coupling force.
		     */
		    Simulate_Substructure_Measured_Values( "fc.txt", Gc, u0c, uc, fcprev, fc, Cnst.Order_Couple, Cnst.Num_Sub );
	       }

	       /* Compose the data to send */
	       for (i = 0; i < Cnst.Order_Couple; i++) {
		    Send[i] = uc[i];
		    Send[i+Cnst.Order_Couple] = fcprev[i];
		    Send[i+2*Cnst.Order_Couple] = fc[i];
	       }

	       /* Send the response */
	       Length = 3*Cnst.Order_Couple;
	       if( Socket_Type == PROTOCOL_TCP ){
		    Send_Data( Send, Length, Client_Socket );
	       } else if ( PROTOCOL_UDP ){
		    if( sendto( Server_Socket, Send, Length*sizeof(float), 0, (struct sockaddr *) &Client_Addr, sizeof(Client_Addr) ) != (int) sizeof(float)*Length ){  /* Sizeof() returns unsigned int */
			 PrintErrorAndExit("sendto() failed or send a different ammount of data");
		    }
	       } else if ( PROTOCOL_NSEP ){
		    /* Send the NSEP_CSIG package. WhatToDo = 4 */
		    Communicate_With_PNSE( 4, 0.0, Send, Recv, 3*Cnst.Order_Couple );
	       }
	  }
     }

     /* Close the connection with the Client */
     if ( Socket_Type == PROTOCOL_TCP ){
	  Close_Socket( &Client_Socket );
     } else if ( Socket_Type == PROTOCOL_UDP ){
	  Close_Socket( &Server_Socket );
     } else if ( Socket_Type == PROTOCOL_NSEP ){
	  Communicate_With_PNSE( 6, 0.0, Send, Recv, 0 );
     }

     if ( Mode == USE_ADWIN ){
#if ADWIN_
	  /* Get the Data from ADwin */
	  printf("Getting the data from ADwin...");
	  GetDataADwin( Cnst.Num_Steps, Cnst.Num_Sub, ADWIN_DATA );
	  printf(" DONE!\n");
     
	  free( ADWIN_DATA );
#endif
     } else if ( Mode == USE_EXACT ){
	  printf("The simulatiovn has finished\n");
	  free( Num_TMD );
     } else if ( Mode == USE_UHYDE ){
	  printf("The simulatiovn has finished\n");
	  free( Num_UHYDE );
     } else {
	  printf("The simulatiovn has finished\n");
     }
     
     /* Free the dinamically allocated memory */
     free( Gc );
     
     free( u0c0 );
     free( u0c );
     free( uc );

     free( fcprev );
     free( fc );

     free( Send );

     return 0;
}



void Print_Help( const char *Program_Name )
{

     fprintf( stderr, "Usage: %s [-h] -m <Mode> -p <Port>", Program_Name );
     fprintf( stderr,
	      "  -h  --help           This help text.\n"
	      "  -m  --mode           The mode used by the program. Default value 1.\n"
	      "                            0 - ADwin will be used to perform the sub-stepping process.\n"
	      "                            1 - The Substructure will be simulated using an exact solution (default).\n"
	      "                            2 - The Substructure will be simulated using measured values.\n"
	      "  -p  --server-port    Port used for communication with the client program. Default value 3333.\n"
	      "  -s  --socket-type    Type of socket to use. Default value 1.\n"
	      "                            1 - Use TCP socket.\n"
	      "                            2 - Use UDP socket.\n" );
}