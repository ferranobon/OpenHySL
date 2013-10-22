#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h> /* For gettimeofday() */
#include <assert.h>   /* For assert() */

#include <libadwin.h>        /* ADwin routines Boot(), Set_DeviceNo(), ... */
#include <libadwin/errno.h>  /* ADwin error handling */

#include "ADwin_Routines.h"
#include "Print_Messages.h"  /* For Print_Header() */

#include "Definitions.h"

void ADwin_Boot( const int32_t Device_Number, const char *Boot_Path  )
{
     /* Set the device into 0x150 (hexadecimal) or 336 */
     Set_DeviceNo( Device_Number );

     /* Boot path should be something like "/opt/adwin/share/btl/ADwin11.btl" */
     Boot( Boot_Path );
 
     /* Check if there is any response from ADwin and if the correct version has
      * been loaded
      */
     ADwin_TestVersion( );
     Print_Header( SUCCESS );
     printf( "ADwin: Boot successful.\n" );  
}

void ADwin_CheckProcessStatus( const int ProcessNumber )
{
     int ADwinStatus;

     Set_DeviceNo( (int32_t) 336 );
     ADwinStatus = Process_Status ( ProcessNumber );

     if ( ADwinStatus == 1){
	  Print_Header( INFO );
	  printf( "Process %i is running.\n", ProcessNumber );
     } else if ( ADwinStatus == 0 ){
	  Print_Header( ERROR );
	  fprintf( stderr, "Process %i is not running.\n", ProcessNumber );
	  exit( EXIT_FAILURE );
     } else if ( ADwinStatus == -1 ){
	  Print_Header( INFO );
	  printf( "Process %i is being stopped and is waiting for the last event.\n", ProcessNumber );
     } else assert(0);
}

void ADwin_ManageProcess( const char* PName, const int PNum, const int dowhat )
{
     Set_DeviceNo( (int32_t) 336 );

     if ( dowhat == 1 ){
	  /* Load the process on ADwin*/
	  Load_Process( PName );
     } else if ( dowhat == 2 ){
	  /* Start the process on ADwin */
	  Start_Process( (int32_t) PNum );
	  /* Check the status of the process */
	  ADwin_CheckProcessStatus( PNum );
     } else if ( dowhat == 3 ){
	  /* Stop the process on ADwin */
	  Stop_Process( (int32_t) PNum );
	  /* Check the status of the process */
	  ADwin_CheckProcessStatus( PNum );
     } else assert( 0 );
}

void ADwin_SendArray( const unsigned int Index, const HYSL_FLOAT *const Array, const unsigned int Length )
{
   
     Set_DeviceNo( (int32_t) 336 );

#if _FLOAT_
     SetData_Float( (int32_t) Index, Array, 1, (int32_t) Length );
#else
     SetData_Double( (int32_t) Index, Array, 1, (int32_t) Length );
#endif

}

void ADwin_Substep( const HYSL_FLOAT *const VecTdT_0c, const unsigned int OrderC, const HYSL_FLOAT Time_To_Wait, HYSL_FLOAT *const VecTdT_c,
		    HYSL_FLOAT *const fcprev_c, HYSL_FLOAT *const fc_c )
{

     unsigned int i;
     unsigned int Length_Receive, Length_Send;
     
     HYSL_FLOAT *Send_ADwin = NULL, *ReceiveADwin = NULL;
     int ADWinReady;

     struct timeval t1;
     struct timeval t2;

     double ElapsedTime;

     Length_Receive = 3*OrderC + 1;
     Length_Send = OrderC + 1;

     Send_ADwin = (HYSL_FLOAT *) calloc( (size_t) Length_Send, sizeof(HYSL_FLOAT) );
     ReceiveADwin = (HYSL_FLOAT *) calloc( (size_t) Length_Receive, sizeof(HYSL_FLOAT) );

     for ( i = 0; i < Length_Receive; i++ ){
	  ReceiveADwin[i] = 0.0;
     }
     ADWinReady = 0;

     /* Prepare the data to be sent. The first value indicates that the data has been uploaded and
      * ADwin can start the substepping process */
     Send_ADwin[0] = 1.0;
     for ( i = 0; i < OrderC; i++ ){
	  Send_ADwin[i + 1] = VecTdT_0c[i];
     }

     /* Set the displacement. In ADwin the array storing the displacement is DATA_2 */
#if _FLOAT_
     SetData_Float( 2, Send_ADwin, 1, (int32_t) Length_Send );
#else
     SetData_Double( 2, Send_ADwin, 1, (int32_t) Length_Send );
#endif

     /* Do nothing until a certain time has passed to avoid overloading adwin system */
     gettimeofday( &t1, NULL );
     ElapsedTime = 0.0;
          
     while ( ElapsedTime < Time_To_Wait ){ /*(NSub - 1.0)*Deltat_Sub*1000.0 - 0.5){*/
	  gettimeofday(&t2, NULL );
	  ElapsedTime = (double) (t2.tv_sec - t1.tv_sec)*1000.0;
	  ElapsedTime += (double) (t2.tv_usec -t1.tv_usec)/1000.0;
     }

     /* Get the displacement when substep is over */
     while( ADWinReady == 0 ){
#if _FLOAT_
	  GetData_Float( 3, ReceiveADwin, 1, (int32_t) Length_Receive );
#else
	  GetData_Double( 3, ReceiveADwin, 1, (int32_t) Length_Receive );
#endif
	  if ( ReceiveADwin[0] == -1.0){

	       /* Get the data from ADWIN (DATA_3) */
	       ADWinReady = 1;
	       for( i = 0; i < OrderC; i++ ){
		    VecTdT_c[i] = ReceiveADwin[i+1];
		    fcprev_c[i] = ReceiveADwin[OrderC + i+1];
		    fc_c[i] = ReceiveADwin[2*OrderC + i+1];
	       }
	  }
     }

     /* Free the dynamically allocated memory */
     free( Send_ADwin );
     free( ReceiveADwin );
}

void ADwin_TestVersion( void )
{
     int Status;
     Set_DeviceNo( (int32_t) 336 );

     Status = Test_Version( );

     if ( Status == 0 ){
	  Print_Header( SUCCESS );
	  printf( "ADwin: Everything is OK.\n" );
     } else if ( Status == 1 ){
	  Print_Header( ERROR );
	  fprintf( stderr, "ADwin: Wrong driver version, processor continues working.\n" );
	  exit( EXIT_FAILURE );
     } else if ( Status == 2 ){
	  Print_Header( ERROR );
	  fprintf( stderr, "ADwin: Wrong driver version, processor stops.\n" );
	  exit( EXIT_FAILURE );
     } else if ( Status == 3 ){
	  Print_Header( ERROR);
	  fprintf( stderr, "ADwin: No response from the ADwin system.\n" );
	  exit( EXIT_FAILURE );
     } else {
	  Print_Header( ERROR );
	  fprintf( stderr, "ADwin: Unrecognised return value.\n" );
	  exit( EXIT_FAILURE );
     }
}

void ADwin_SaveData_ASCII( const char *FileName, const unsigned int Num_Steps, const unsigned int Num_Sub,
			   const unsigned short int Num_Channels, const char **Chan_Names, const int DataIndex )
{
     unsigned int i, j, Length;
     HYSL_FLOAT *Data;
     FILE *OutFile;

     OutFile = fopen( FileName, "w" );

     if( OutFile == NULL ){
	  Print_Header( WARNING );
	  fprintf( stderr, "ADwin_SaveData_TXT: Could not open %s.\n", FileName );
     }

     Length = Num_Sub*Num_Steps*Num_Channels;
     Data = (HYSL_FLOAT *) calloc( (size_t) Length, sizeof( HYSL_FLOAT ) );
     if( Data == NULL ){
	  Print_Header( WARNING );
	  fprintf( stderr, "ADwin_SaveData_HDF5: Out of memory. Manual extraction of the data required.\n" );
     }

     /* Get the data from ADwin */
#if _FLOAT_
     GetData_Float( (int32_t) DataIndex, Data, 1, (int32_t) Length );
#else
     GetData_Double( (int32_t) DataIndex, Data, 1, (int32_t) Length );
#endif

     for( i = 0; i < Num_Channels; i++ ){
	  fprintf( OutFile, "%s\t", Chan_Names[i] );
     }
     fprintf( OutFile, "\n" );

     for( i = 0; i < Num_Sub*Num_Steps; i++ ){
	  for( j = 0; j < Num_Channels; j++ ){
	       fprintf( OutFile, "%.8lE\t", Data[i*Num_Channels + j] );
	  }
	  fprintf( OutFile, "\n" );
     }
     
     /* Free the data and close the file */
     free( Data );
     fclose( OutFile );
}
