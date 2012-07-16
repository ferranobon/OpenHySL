#include <stdio.h>
#include <stdlib.h>
#include <assert.h>    /* For assert() */
#include <sys/time.h>  /* For gettimeoffday() */

/* ADwin header files */
#include <libadwin.h>
#include <libadwin/errno.h>

#include "RoutinesADwin.h"
#include "ErrorHandling.h"


int BootADWIN( const int Device_Number, const char *Boot_Path  )
{

     /* Set the device into 0x150 (hexadecimal) or 336 */
     Set_DeviceNo( (int32_t) Device_Number );
 
     /* Boot path should be something like "/opt/adwin/share/btl/ADwin11.btl" */
     Boot( Boot_Path );
 
     /* Check if there is any response from ADwin and if the correct version has
      * been loaded
      */
     if ( ADWIN_TestVersion( ) ){
	  return 1;
     } else return 0;

}

int ADWIN_TestVersion( void )
{
     int Status;
     Set_DeviceNo( (int32_t) 336 );

     Status = Test_Version( );

     if ( Status == 0 ){
	  printf( "ADWIN: Everything is OK.\n" );
	  return 1;
     } else if ( Status == 1 ){
	  PrintErrorAndExit( "ADWIN: Wrong driver version, processor continues working." );
	  return 0;
     } else if ( Status == 2 ){
	  PrintErrorAndExit( "ADWIN: Wrong driver version, processor stops." );
	  return 0;
     } else if ( Status == 3 ){
	  PrintErrorAndExit( "ADWIN: No response from the ADwin system." );
	  return 0;
     } else {
	  PrintErrorAndExit( "ADWIN: Unrecognised return value." );
	  return 0;
     }

}

void ADWIN_ManageProcess( const char* PName, const int PNum, const int dowhat )
{
     Set_DeviceNo( (int32_t) 336 );

     if ( dowhat == 1 ){
	  /* Load the process on ADwin*/
	  Load_Process( PName );
     } else if ( dowhat == 2 ){
	  /* Start the process on ADwin */
	  Start_Process( (int32_t) PNum );
	  /* Check the status of the process */
	  ADWIN_CheckProcessStatus( 2 );
     } else if ( dowhat == 3 ){
	  /* Stop the process on ADwin */
	  Stop_Process( (int32_t) PNum );
	  /* Check the status of the process */
	  ADWIN_CheckProcessStatus( 2 );
     } else assert( 0 );

}

void ADWIN_CheckProcessStatus( int ProcessNumber )
{
     int ADwinStatus;

     Set_DeviceNo( (int32_t) 336 );
     ADwinStatus = Process_Status ( ProcessNumber );

     if ( ADwinStatus == 1){
	  printf( "Process is running.\n" );
     } else if ( ADwinStatus == 0 ){
	  printf("Process %i is not running\n", ProcessNumber );
	  PrintErrorAndExit( "Process is not running." );
     } else if ( ADwinStatus == -1 ){
	  printf( "Process is being stopped and is waiting for the last event.\n" );
     }
}

void ADWIN_SetGc( const float *const Gc, const int length )
{
     /* In the code of ADWIN, the variable G is stored in DATA_1 */
     Set_DeviceNo( (int32_t) 336 );

     SetData_Float( 1, (float *)Gc, 1, length );

}

void ADWIN_Substep( const float *const u0c, float *const uc, float *const fcprev, float *const fc, const int OrderC )
{

     static int i;
     
     float *Send_ADwin, *ReceiveADwin;
     static int ADWinReady;

     static struct timeval t1;
     static struct timeval t2;
     static struct timeval t3;
     static double ElapsedTime;

     Send_ADwin = (float *) calloc( OrderC + 1, sizeof(float) );
     ReceiveADwin = (float *) calloc( 3*OrderC + 1, sizeof(float) );

     for ( i = 0; i < 3*OrderC + 1; i++ ){
	  ReceiveADwin[i] = 0.0;
     }
     ADWinReady = 0;

     /* Prepare the data to be sent. The first value indicates that the data has been uploaded and
      * ADwin can start the substepping process */
     Send_ADwin[0] = 1.0;
     for ( i = 0; i < OrderC; i++ ){
	  Send_ADwin[i + 1] = u0c[i];
     }

     /* Set the displacement. In ADwin the array storing the displacement is DATA_2 */
     SetData_Float( 2, Send_ADwin, 1, OrderC + 1 );

     /* Do nothing until a certain time has passed to avoid overloading adwin system */
     gettimeofday( &t1, NULL );
     ElapsedTime = 0.0;
          

     while ( ElapsedTime < 8.0 ){ //(NSub - 1.0)*Deltat_Sub*1000.0 - 0.5){
	  gettimeofday(&t2, NULL );
	  ElapsedTime = (t2.tv_sec - t1.tv_sec)*1000.0;
	  ElapsedTime += (t2.tv_usec -t1.tv_usec)/1000.0;
     }
     
     /* Get the displacement when substep is over */
     while( ADWinReady == 0 ){
	  GetData_Float( 3, ReceiveADwin, 1, 3*OrderC + 1);
	  if ( ReceiveADwin[0] == -1.0){

	       gettimeofday(&t3, NULL );
	       ElapsedTime = (t3.tv_sec - t1.tv_sec)*1000.0;
	       ElapsedTime += (t3.tv_usec -t1.tv_usec)/1000.0;
	       /* Get the data from ADWIN (DATA_3) */

	       ADWinReady = 1;
	       for( i = 0; i < OrderC; i++ ){
		    uc[i] = ReceiveADwin[i+1];
		    fcprev[i] = ReceiveADwin[OrderC + i+1];
		    fc[i] = ReceiveADwin[2*OrderC + i+1];
	       }
	  }
     }

     /* Free the dynamically allocated memory */
     free( Send_ADwin );
     free( ReceiveADwin );
}

void GetDataADwin( const int Num_Steps, const int Num_Sub, float *const Data )
{
     FILE *OutFile;
     int i,j;

     GetData_Float( 97, Data, 1, Num_Sub*Num_Steps*NUM_CHANNELS );
     OutFile = fopen( "Data.txt", "w");
  
     /* Print the header file */
     fprintf( OutFile, "substep\tdctr1\tdctr2\tdmeas1\tdmeas2\taccctrl1\taccctrl2\tameas1\tameas2\taccTMDy\taccTMDx\tAccFrame\tFcY1\tFcY2\tfc\tdTMDx\tdTMDy\tPress_pctrl\tPress_pmeas\tuc\tTime do substep\ttime substep\ttimesynPC\ttime_from_sub1_to_newPCdata\n" );

     for (i = 0; i < Num_Sub*Num_Steps; i++ ){
	  for ( j = 0; j < NUM_CHANNELS; j++ ){
	       fprintf(OutFile, "%.10le\t", Data[i*NUM_CHANNELS+j] );
	  }
	  fprintf(OutFile,"\n");
     }

     fclose(OutFile);
}
