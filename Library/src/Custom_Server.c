#include <stdio.h>
#include <stdlib.h>      /* For atoi( ) and exit() */
#include <stdbool.h>
#include <getopt.h>      /* For getopt_long() */
#include <netdb.h>

#include "Algorithm_Aux.h"

#include "Definitions.h"
#include "Print_Messages.h"
#include "Substructure.h"
#include "Substructure_Exact.h"
#include "Substructure_Remote.h"
#include "Substructure_UHYDEfbr.h"
#include "Substructure_SimMeasured.h"
#include "Substructure_Experimental.h"
#include "Substructure_CouplingNodes.h"

#include "Custom_Server.h"

#if _ADWIN_
#include "ADwin_Routines.h"
#endif

int main( int argc, char **argv )
{
     int i, pos
	  ;
     const char* FileConf;
     AlgConst_t InitCnt;
     CouplingNode_t CNodes;

     /* Variables concerning data communication */
     int Server_Socket;                /* Socket for the server */
     int Client_Socket;                /* Socket for the client */

     int Socket_Type;
     const char *Port;
     struct sockaddr_storage Client_Addr;
     socklen_t Client_AddrLen;

     /* Variables for the sub-structure testing/simulation */
     int Is_Not_Finished;
     int Length;

     hysl_float_t *IGain;
     bool Called_Sub = false, Called_ADwin = false;
     hysl_float_t *Recv = NULL, *Recv_ADwin = NULL, *VecTdT0_c_ADwin = NULL;
     hysl_float_t *Send = NULL;
     hysl_float_t *VecTdT0_c;
     hysl_float_t *fcprev, *fc;

     hysl_float_t GAcc = 0.0;

     /* Array where the data from ADwin will be stored */
     hysl_float_t *ADWIN_DATA = NULL;

     /* Variables to deal with arguments */
     int Selected_Option;

     struct option long_options[] = {
	  {"help", no_argument, 0, 'h'},
	  {"config-file", required_argument, 0, 'c'},
	  {"server-port", required_argument, 0, 'p'},
	  {"socket-type", required_argument, 0, 's'},
	  {0, 0, 0, 0}
     };

     Socket_Type = REMOTE_TCP;
     Port = "3333";
     FileConf = "ConfFile_Remote.conf";

     /* This is only used if there are no arguments */
     if( argc == 1 ){
	  Print_Header( INFO );
	  printf( "Assuming the configuration file to be: %s.\n", FileConf );
	  Print_Header( INFO );
	  printf( "Defaulting to TCP/IP protocol on port 3333\n" );
     }

     /* Assign each argument to the correct variable */
     while( (Selected_Option = getopt_long( argc, argv, "c:p:s:h", long_options, NULL )) != -1 ){
	  switch( Selected_Option ){
	  case 'c':
	       FileConf = optarg;
	       break;
	  case 'p':
	       Port = optarg;
	       break;
	  case 's':
	       Socket_Type = atoi( optarg );
	       
	       if ( Socket_Type != REMOTE_TCP && Socket_Type != REMOTE_UDP ){
		    Print_Header( ERROR );
		    fprintf( stderr, "Socket type %d is not a valid.\n", Socket_Type );
		    CustomServer_PrintHelp( argv[0] );
		    return EXIT_FAILURE;
	       }
	       break;		    
	  case 'h':
	       CustomServer_PrintHelp( argv[0] );
	       return EXIT_FAILURE;
	       break;
	  case '?':
	       /* Long options already prints an error message telling that there is an unrecognised option */
	       CustomServer_PrintHelp( argv[0] );
	       return EXIT_FAILURE;
	  case ':':
	       /* Long options already prints an error message telling that the option requieres an argument */
	       CustomServer_PrintHelp( argv[0] );
	       return EXIT_FAILURE;
	  }
     }

     /* Create a TCP/IP socket for the server */
     Server_Socket = Substructure_Remote_SetupServer( Port, Socket_Type );
     if( Server_Socket < -1 ){
	  Print_Header( ERROR );
	  fprintf( stderr, "Setup_Server_Socket() failed" );
	  return EXIT_FAILURE;
     } else {
	  if( Socket_Type == REMOTE_TCP ){
	  Print_Header( SUCCESS );
	  printf( "TCP Server created.\n" );
	  Print_Header( INFO );
	  printf( "Listening on port %s for incoming client connections...\n", Port );
	  } else {
	       Print_Header( SUCCESS );
	       printf( "UDP Server created.\nListening on port %s for incoming client connections...\n", Port );
	  }
     }

     /* Accept the connection from the client in case of TCP sockets */
     if ( Socket_Type == REMOTE_TCP ){
	  Client_Socket = Substructure_Remote_AcceptTCPClientConnection( Server_Socket );
     }

     /* Initialise the constants of the substructure */
     CustomServer_Init( FileConf, &InitCnt );

     /* Read the coupling nodes from a file */
     Substructure_ReadCouplingNodes( &InitCnt, &CNodes );


     /* Dynamically allocate memory */
     IGain = (hysl_float_t *) calloc( (size_t) InitCnt.OrderSub*(size_t) InitCnt.OrderSub, sizeof( hysl_float_t ) );

     VecTdT0_c = (hysl_float_t *) calloc( (size_t) InitCnt.OrderSub, sizeof( hysl_float_t ) );

     fcprev = (hysl_float_t *) calloc( (size_t) InitCnt.OrderSub, sizeof( hysl_float_t ) );
     fc = (hysl_float_t *) calloc( (size_t) InitCnt.OrderSub, sizeof( hysl_float_t ) );

     Send = (hysl_float_t *) calloc( (size_t) 3*(size_t) InitCnt.OrderSub, sizeof( hysl_float_t ) );
     Recv = (hysl_float_t *) calloc( (size_t) 1+(size_t) InitCnt.OrderSub, sizeof( hysl_float_t ) );

     /* Receive matrix IGain */
     Length = InitCnt.OrderSub*InitCnt.OrderSub;

     if( Socket_Type == REMOTE_TCP ){
	  Substructure_Remote_Receive( Client_Socket, (unsigned int) Length, sizeof(hysl_float_t), (char *const) IGain );
     } else {
	  Client_AddrLen = sizeof(Client_Addr);
	  if ( recvfrom( Server_Socket, IGain, (size_t) Length*sizeof(hysl_float_t), 0, (struct sockaddr *) &Client_Addr, &Client_AddrLen) < 0 ){
	       Print_Header( ERROR );
	       fprintf( stderr, "Recvfrom failed" );
	       return EXIT_FAILURE;
	  }
     }

     for( i = 0; i < CNodes.Order; i++ ){
	  if( CNodes.Sub[i].Type == EXP_ADWIN ){
	       Substructure_SendGainMatrix( IGain, (unsigned int) CNodes.Order, &CNodes.Sub[i] );
	       /* Allocate the memory for reading the ADwin data */
	       if( ADWIN_DATA == NULL ){
		    ADWIN_DATA = (hysl_float_t *) calloc( (size_t) InitCnt.OrderSub*InitCnt.NSubstep*NUM_CHANNELS, sizeof( hysl_float_t ) );
	       }
	  }
     }

     Is_Not_Finished = 1;
     while( Is_Not_Finished ){

	  Called_Sub = false;

	  /* Receive the displacement */
	  Length = InitCnt.OrderSub;
	  if ( Socket_Type == REMOTE_TCP ){
	       Substructure_Remote_Receive( Client_Socket, (unsigned int) Length + 1, sizeof(hysl_float_t), (char *const) Recv );
	  } else {
	       Client_AddrLen = sizeof(Client_Addr);
	       if ( recvfrom( Server_Socket, Recv, (size_t) (Length+1)*sizeof(hysl_float_t), 0, (struct sockaddr *) &Client_Addr, &Client_AddrLen) < 0 ){
		    Print_Header( ERROR );
		    fprintf( stderr, "recvfrom() failed" );
		    return EXIT_FAILURE;
	       }
	  }
  
	  for( i = 0; i < CNodes.Order; i++ ){
	       VecTdT0_c[i] = Recv[i];
	  }
	  GAcc = Recv[CNodes.Order];

	  if ( VecTdT0_c[0] == -9999.0 ){
	       Is_Not_Finished = 0;
	  } else {
	       for( i = 0; i < CNodes.Order; i++ ){
		    switch ( CNodes.Sub[i].Type ){
		    case SIM_EXACT_MDOF:
			 /* This is the same case as SIM_MEASURED. All the simulated substructures are treated together
			  * in the same routine.*/
		    case SIM_EXACT_SDOF:
			 /* This is the same case as SIM_MEASURED. All the simulated substructures are treated together
			  * in the same routine.*/
		    case SIM_EXACT_ESP:
			 /* This is the same case as SIM_MEASURED. All the simulated substructures are treated together
			  * in the same routine.*/
		    case SIM_NEWMARK:
			 /* This is the same case as SIM_MEASURED. All the simulated substructures are treated together
			  * in the same routine.*/
		    case SIM_BOUCWEN:
			 /* This is the same case as SIM_MEASURED. All the simulated substructures are treated together
			  * in the same routine.*/
		    case SIM_UHYDE:
			 /* This is the same case as SIM_MEASURED. All the simulated substructures are treated together
			  * in the same routine.*/
		    case SIM_STONEDRUMS:
			 /* This is the same case as SIM_MEASURED. All the simulated substructures are treated together
			  * in the same routine.*/
		    case SIM_MEASURED:
			 /* Call the Simulate_Substructures() function only once. All the simulated substructures are
			  * handled together in this routine */
			 if( !Called_Sub ){
			      Substructure_Simulate( IGain, VecTdT0_c, GAcc, InitCnt.NSubstep, InitCnt.DeltaT_Sub, &CNodes, &Recv[0],
						     &Recv[CNodes.Order], &Recv[2*CNodes.Order] );
			      Called_Sub = true;
			 }
			 break;
#if _ADWIN_
		    case EXP_ADWIN:
			 /* Tell ADwin to perform the substepping process */
			 if( !Called_ADwin ){
			      ADwin_Substep_Pre( VecTdT0_c_ADwin, (unsigned int) CNodes.OrderADwin );
			      Called_ADwin = true;
			 }
			 break;
#endif
		    }
	       }

#if _ADWIN_
	       if( CNodes.OrderADwin > 0 ){
		    ADwin_Substep_Post( (unsigned int) CNodes.OrderADwin, 75.0, &Recv_ADwin[0], &Recv_ADwin[CNodes.OrderADwin], &Recv_ADwin[2*CNodes.OrderADwin] );
	       }
#endif

	       pos = 0;
	       for ( i = 0; i < CNodes.Order; i++ ){
		    if (CNodes.Sub[i].Type == EXP_ADWIN){
			 Send[i] = Recv_ADwin[pos];
			 Send[i + InitCnt.OrderSub] = Recv_ADwin[CNodes.OrderADwin + pos];
			 Send[i + 2*InitCnt.OrderSub] = Recv_ADwin[2*CNodes.OrderADwin + pos];
			 pos = pos + 1;
		    } else {
			 Send[i] = Recv[i];
			 Send[i + InitCnt.OrderSub]= Recv[CNodes.Order + i];
			 Send[i + 2*InitCnt.OrderSub] = Recv[2*CNodes.Order + i];
		    }
	       }

	       /* Send the response */
	       Length = 3*InitCnt.OrderSub;
	       if( Socket_Type == REMOTE_TCP ){
		    Substructure_Remote_Send( Client_Socket, (unsigned int) Length, sizeof(hysl_float_t), (char *const) Send );
	       } else{
		    if( sendto( Server_Socket, Send, (size_t) Length*sizeof(hysl_float_t), 0, (struct sockaddr *) &Client_Addr, sizeof(Client_Addr) ) != (int) sizeof(hysl_float_t)*Length ){  /* Sizeof() returns unsigned int */
			 Print_Header( ERROR );
			 fprintf( stderr, "sendto() failed or send a different ammount of data");
			 return EXIT_FAILURE;
		    }
	       }
	  }
     }

     Substructure_Remote_CloseSocket( &Client_Socket );

     /* Free initiation values */
     CustomServer_Destroy( &InitCnt );

     /* Free the coupling nodes memory and close sockets if appropiate */
     Substructure_DeleteCouplingNodes( &CNodes );
     
     /* Free the dinamically allocated memory */
     free( IGain );
     
     free( VecTdT0_c );

     free( fcprev );
     free( fc );

     free( Send );
     free( Recv );

     free( ADWIN_DATA );

     return 0;
}



void CustomServer_PrintHelp( const char *Program_Name )
{

     fprintf( stderr, "Usage: %s [-h] -c <ConfFile> -p <Port> -s <SoketType>\n", Program_Name );
     fprintf( stderr,
	      "  -h  --help           This help text.\n"
	      "  -c  --config-file    The name of the configuration file. Default value: ConfFile_Remote.conf\n"
	      "  -p  --server-port    Port used for communication with the client program. Default value 3333.\n"
	      "  -s  --socket-type    Type of socket to use. Default value 1.\n"
	      "                            1 - Use TCP socket.\n"
	      "                            2 - Use UDP socket.\n" );
}


void CustomServer_Init( const char *FileName, AlgConst_t *const InitConst )
{

     ConfFile_t *Config;
     
     Config = ConfFile_Create( 70 );

     ConfFile_ReadFile( FileName, Config );

#if _FLOAT_
     InitConst->Delta_t = ConfFile_GetFloat( Config, "General:Delta" );
#else
     InitConst->Delta_t = ConfFile_GetDouble( Config, "General:Delta" );
#endif
     if ( InitConst->Delta_t <= 0.0 ){
	  Print_Header( ERROR );
	  fprintf( stderr, "Invalid time step.\n" );
	  exit( EXIT_FAILURE );
     }

     /* Newmark integration constants */
#if _FLOAT_
     InitConst->TIntConst.Gamma = ConfFile_GetFloat( Config, "Newmark:Gamma" );
     InitConst->TIntConst.Beta = ConfFile_GetFloat( Config, "Newmark:Beta" );
#else
     InitConst->TIntConst.Gamma = ConfFile_GetDouble( Config, "Newmark:Gamma" );
     InitConst->TIntConst.Beta = ConfFile_GetDouble( Config, "Newmark:Beta" );
#endif

     /* Constants for Ending Step */
#if _FLOAT_
     InitConst->a0 = 1.0f/(InitConst->TIntConst.Beta*InitConst->Delta_t*InitConst->Delta_t);
     InitConst->a2 = 1.0f/(InitConst->TIntConst.Beta*InitConst->Delta_t);
     InitConst->a3 = 1.0f/(2.0f*InitConst->TIntConst.Beta) - 1.0f;
     InitConst->a4 = InitConst->TIntConst.Gamma/InitConst->TIntConst.Beta - 1.0f;
     InitConst->a5 = (InitConst->Delta_t/2.0f)*(InitConst->TIntConst.Gamma/InitConst->TIntConst.Beta - 2.0f);
     InitConst->a6 = (1.0f - InitConst->TIntConst.Gamma)*InitConst->Delta_t;
#else
     InitConst->a0 = 1.0/(InitConst->TIntConst.Beta*InitConst->Delta_t*InitConst->Delta_t);
     InitConst->a2 = 1.0/(InitConst->TIntConst.Beta*InitConst->Delta_t);
     InitConst->a3 = 1.0/(2.0*InitConst->TIntConst.Beta) - 1.0;
     InitConst->a4 = InitConst->TIntConst.Gamma/InitConst->TIntConst.Beta - 1.0;
     InitConst->a5 = (InitConst->Delta_t/2.0)*(InitConst->TIntConst.Gamma/InitConst->TIntConst.Beta - 2.0);
     InitConst->a6 = (1.0 - InitConst->TIntConst.Gamma)*InitConst->Delta_t;
#endif

     InitConst->a1 = InitConst->TIntConst.Gamma/(InitConst->TIntConst.Beta*InitConst->Delta_t);
     InitConst->a7 = InitConst->TIntConst.Gamma*InitConst->Delta_t;

     /* File Names */
     InitConst->FileCNodes = strdup( ConfFile_GetString( Config, "FileNames:Coupling_Nodes" ) );

     /* Number of substructures */
     InitConst->OrderSub = ConfFile_GetInt( Config, "Substructure:Order" );
     if ( InitConst->OrderSub < 0 ){
	  Print_Header( ERROR );
	  fprintf( stderr, "Invalid option for the number of sub-structures of the matrices.\n" );
	  exit( EXIT_FAILURE );
     }
     
     /* Number of substructures */
     InitConst->NSubstep = (unsigned int) ConfFile_GetInt( Config, "Substructure:Num_Substeps" );

     InitConst->DeltaT_Sub = InitConst->Delta_t/(hysl_float_t) InitConst->NSubstep;

     ConfFile_Destroy( Config );

     Print_Header( SUCCESS );
     printf( "Initialisation succcessfully completed.\n" );
}

void CustomServer_Destroy( AlgConst_t *const InitConst )
{

     free( InitConst->FileCNodes );

}
