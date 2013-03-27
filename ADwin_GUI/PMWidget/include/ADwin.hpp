#ifndef ADWIN_HPP
#define ADWIN_HPP

#include <string.h>

class ADwin_Class
{
public:

     /* TODO, is this correct??? */
     static const int32_t NEW_PRESSURE_COMMAND = 5;
     static const int32_t CURRENT_PRESSURE = 6;
     static const int32_t EMERGENCY_STOP = 10;
     static const int32_t ADWIN_READY = 11;
     static const int32_t TESTING = 12;
     static const int32_t CONTROLBOX_READY = 14;
     static const int32_t LIMIT_CHECK = 18;

     /* Other variables */
     static const int32_t NO_PRESSURE = 0;
     static const int32_t LOW_PRESSURE = 1;
     static const int32_t HIGH_PRESSURE = 2;
     static const int32_t STOP_PRESSURE = 3;

     /* PAR - variables */
     static const int32_t MESSAGE_FLAG = 15;
     static const int OVER_DMIN = 1;
     static const int OVER_DMAX = 2;
     static const int OVER_VMIN = 3;
     static const int OVER_VMAX = 4;
     static const int OVER_AMIN = 5;
     static const int OVER_AMAX = 6;
     static const int OVER_FMAX = 7;


     /* FPAR - Limit variables */
     static const int32_t DMIN1 = 51;
     static const int32_t DMIN2 = 52;
     static const int32_t DMAX1 = 53;
     static const int32_t DMAX2 = 54;

     /* FPAR - Control variables */
     static const int32_t DCNTRL1 = 3;
     static const int32_t DCNTRL2 = 4;

     /* FPAR - PID variables */
     static const int32_t PID_P = 12;
     static const int32_t PID_I = 13;
     static const int32_t PID_D = 14;

     /* Measurement variables */
     static const int32_t DMEAS1 = 72;
     static const int32_t DMEAS2 = 73;

     /* Limit default values */
     static const float Default_DMin1 = -0.149;
     static const float Default_DMax1 = 0.149;
     static const float Default_DMin2 = -0.149;
     static const float Default_DMax2 = 0.149;
     /* Limit control values */
     static const float DCntrlMin1 = -0.15;
     static const float DCntrlMax1 = 0.15;
     static const float DCntrlMin2 = -0.15;
     static const float DCntrlMax2 = 0.15;

     /* PID default values */
     static const float Default_P = 3000.0;
     static const float Default_I = 0.0;
     static const float Default_D = 0.0;
     /* Limit PID values */
     static const float P_Min = 1000.0;
     static const float P_Max = 20000.0;
     static const float I_Min = 0.0;
     static const float I_Max = 500.0;
     static const float D_Min = 0.0;
     static const float D_Max = 100.0;



     int32_t Boot_ADwin( const int Device_Number, const char *Boot_Path );
     void Load_Process_ADwin( const char *Process_Path, int *Process_Is_Loaded, int32_t *const Error );
     void Start_Process_ADwin( const int32_t Process_Num, int *Process_Is_Started, int32_t *const Error );
     void Stop_Process_ADwin( const int32_t Process_Num, int *Process_Is_Running, int32_t *const Error );
     void Set_Par_ADwin( const int32_t Index, const int32_t Value, int32_t *const Error );
     void Get_Par_ADwin( const int32_t Index, int32_t *const Value,  int32_t *const Error );
     void Set_FPar_ADwin( const int32_t Index, const float Value, int32_t *const Error );
     void Get_FPar_ADwin( const int32_t Index, float *const Value,  int32_t *const Error );
     void Get_DataLength_ADwin( const int DataNo, int32_t *const Length, int32_t *const Error );
     void Get_DataFloat_ADwin( const int DataNo, float *const Data, const int32_t Start_Index, const int32_t Count, int32_t *const Error );
     char* Get_Error_Text( const int32_t Error );
     
};

#endif
