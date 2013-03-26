#include <string>
/* ADwin header files */
#include <libadwin.h>
#include <libadwin/errno.h>

#include "ADwin.hpp"

int32_t ADwin_Class::Boot_ADwin( const int Device_Number, const char *Boot_Path  )
{
     /* Set the device into 0x150 */
     Set_DeviceNo( (int32_t) Device_Number );
 
     Boot( Boot_Path );
 
     /* Check if there has been any error during the boot process. */
     return (Get_Last_Error( ));
}     

void ADwin_Class::Load_Process_ADwin( const char *Process_Path, int *Process_Is_Loaded, int32_t *const Error )
{
     *Error = 0;  /* Assuming everything is OK */

     /* If the process is not yet loaded, do so */
     if ( !(*Process_Is_Loaded) ){
	  Load_Process( Process_Path );

	  /* Check if there is some error while loading the process.
	   * If everything is ok, set the variable to Process_Is_Loaded to 1
	   */
	  *Error = Get_Last_Error( );
	  if ( *Error == 0 ){
	       *Process_Is_Loaded = 1;
	  }
	  /* If the Process is already loaded, the Error is -1 and do nothing.*/
     } else *Error = -1;
}

void ADwin_Class::Start_Process_ADwin( const int32_t Process_Num, int *Process_Is_Started, int32_t *const Error )
{
     int32_t PStatus;

     *Error = 0; /* Assuming everything is OK */

     /* If the process is not yet started, do so */
     if ( !(*Process_Is_Started) ){
	  Start_Process( Process_Num );
    
	  /* Check if there is some error while starting the process. */
	  PStatus = Process_Status( Process_Num );
	  if ( PStatus == 0 ){
	       /* Figure out why the process is not running */
	       *Error = Get_Last_Error( );
	  } else if ( PStatus == 1 ){
	       /* The process is running */
	       *Process_Is_Started = 1;
	  } else if ( PStatus == -1 ){
	       *Process_Is_Started = -1;
	  }
	  /* If the Process is already running, set Error to -1 and do nothing. */
     } else *Error = -1;
}
	       
void ADwin_Class::Stop_Process_ADwin( const int32_t Process_Num, int *Process_Is_Running, int32_t *const Error ){

     int32_t PStatus;

     *Error = 0;  /* Assuming everything is OK */

     /* If the process is running, stop it */
     if ( Process_Is_Running ){
	  Stop_Process( Process_Num );

	  /* Check if the process has really stoped */
	  PStatus = Process_Status( Process_Num );
	  /* If the status is -1, that is the process is running but it is waiting until the
	   * last Event in order to finish, keep reading until it status changes */
	  while ( PStatus == -1 ){
	       PStatus = Process_Status( Process_Num );
	  }

	  if ( PStatus == 1 ){
	       /* The process is still running. Check what is wrong. */
	       *Error = Get_Last_Error( );
	  } else Process_Is_Running = 0;
     } else {
	  /* The process has already been stoped do nothing. */
	  *Error = -1;
     }
}

void ADwin_Class::Set_Par_ADwin( const int32_t Index, const int32_t Value, int32_t *const Error )
{
     Set_Par( Index, Value );
     *Error = Get_Last_Error( );  /* Check if there has been any error during this process */
}

void ADwin_Class::Get_Par_ADwin( const int32_t Index, int32_t *const Value, int32_t *const Error )
{
     *Value = Get_Par( Index );
     *Error = Get_Last_Error( );  /* Check if there has been any error during this process */
}

void ADwin_Class::Set_FPar_ADwin( const int32_t Index, const float Value, int32_t *const Error )
{
     Set_FPar( Index, Value );
     *Error = Get_Last_Error( );  /* Check if there has been any error during this process */
}

void ADwin_Class::Get_FPar_ADwin( const int32_t Index, float *const Value, int32_t *const Error )
{
     *Value = Get_FPar( Index );
     *Error = Get_Last_Error( );  /* Check if there has been any error during this process */
}

void ADwin_Class::Get_DataLength_ADwin( const int DataNo, int32_t *const Length, int32_t *const Error )
{
     *Length = Data_Length( DataNo );
     *Error = Get_Last_Error( );  /* Check if there has been any error during this process */
}

void ADwin_Class::Get_DataFloat_ADwin( const int DataNo, float *const Data, const int32_t Start_Index, const int32_t Count, int32_t *const Error )
{

     GetData_Float( DataNo, Data, Start_Index, Count );
     *Error = Get_Last_Error( );  /* Check if there has been any error during this process */
}

char* ADwin_Class::Get_Error_Text( const int32_t Error )
{
     char* Error_Text;

     Error_Text = Get_Last_Error_Text( Error );

     return Error_Text;
}

