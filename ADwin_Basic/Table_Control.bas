'<ADbasic Header, Headerversion 001.001>
' Process_Number                 = 1
' Initial_Processdelay           = 3000
' Eventsource                    = Timer
' Control_long_Delays_for_Stop   = No
' Priority                       = High
' Version                        = 1
' ADbasic_Version                = 5.0.8
' Optimize                       = Yes
' Optimize_Level                 = 1
' Info_Last_Save                 = BASILISC  Basilisc\ferran
'<Header End>
#include ADwinPRO_ALL.inc ' Include-file for Pro-system t11

' ----------------------------------------------------------------------------------
' ---------------------------- Definition of constants -----------------------------
' ----------------------------------------------------------------------------------
#define Num_Channel            2                   ' Number of channels. It is equal to the number of directions that
Rem                                                  the shake table can handle.
#define Offset_DA_AD           32768               ' Zero offset of the DA and AD converter
#define Max_Error_Int          0.1                 ' Limit value of the integration error
#define Clock_To_Time_Conv     ((10.0/3.0)*1.0e-9) ' Conversion from clock to time units [seconds]
#define Default_dtcontrol      0.00025             ' Default value of the time increment for the control loop. Possible values are
Rem                                                         - 0.00010 in the case of real time tests.
Rem                                                         - 0.0025 in the case of distributed tests with a step time of 2s
#define Default_dtdata         0.5                 ' Time increment for displacement/new data. Ask Van Thuan
#define Default_Freq_Filter    35.0                ' Default frequency filter
#define DtcompVA               0.005               ' Time interval for computing velocity and acceleration from measured values
#define Length_His_Disp        1000                ' This value should be larger than DtCompVA/dtcontrol
#define Length_His_Vel         1000                ' This value should be larger than DtCompVA/dtcontrol

' Different types of error flags. Used in combination with ERROR_FLAG variable.
#define NO_ERROR               0                   ' No error
#define ERROR_UNDER_DMIN       1                   ' The measured value is below the minimum displacement value
#define ERROR_OVER_DMAX        2                   ' The measured value is over the maximum displacement value
#define ERROR_UNDER_VMIN       3                   ' The measured value is below the minimum velocity value
#define ERROR_OVER_VMAX        4                   ' The measured value is over the maximum velocity value
#define ERROR_UNDER_AMIN       5                   ' The measured value is below the minimum acceleration value
#define ERROR_OVER_AMAX        6                   ' The measured value is over the maximum acceleration value
#define ERROR_OVER_FMAX        7                   ' The measured value is over the maximum force value

' ----------------------------------------------------------------------------------
' ------------------------- Definition of global variables -------------------------
' ----------------------------------------------------------------------------------
#define NEW_DATA_AVAILABLE     PAR_4               ' Indicates whether there is new data available (value of 1) or not (value of 0).
#define NEW_PRESSURE_COMMAND   PAR_5               ' Used to change the status of the hydraulic pressure in Check_PressureControl().
Rem                                                        0 = Do nothing;
Rem                                                        1 = Require Low pressure;
Rem                                                        2 = Require High Pressure;
Rem                                                        3 = Stop pressure.
#define CURRENT_PRESSURE       PAR_6               ' Indicates the status of the pressure. 0 = No Pressure; 1 = Low Pressure; 2 = High Pressure
Rem                                                  It follows this combination:
Rem                                                  No Pressure: Y1 = 0 and Y2 = 0
Rem                                                  Low Pressure: Y1 = 0 and Y2 = 1
Rem                                                  High Pressure: Y1 = 1 and Y2 = 0
Rem                                                  Transition time (High to low and viceversa): Y1 = 1 and Y2 = 1
#define Y1                     PAR_7               ' Variable used to control the hydraulic pressure in the hydraulic block together with Y2.
#define Y2                     PAR_8               ' Variable used to control the hydraulic pressure in the hydraulic block together with Y1.
#define NEW_STATUS_DOUT        PAR_9               ' Variable used to identify whether a digital output signal has changed or not.
#define EMERGENCY_STOP         PAR_10              ' Variable to stop the process suddenly
#define ADWIN_READY            PAR_11              ' ADwin is ready to begin the test
#define TESTING                PAR_12              ' The substepping process has started and it is currently running a test. This is an inter-process communication variable.
#define CONTROLBOX_READY       PAR_14              ' The hydraulic pressure can be set
#define ERROR_FLAG             PAR_15              ' Variable used to store the error during the check of measured values against limit values.
#define CHECK_LIMITS_VAR       PAR_18              ' Variable to enable (value equal to 1) or disable (value equal to 0) the limits in order to
Rem                                                  operate in the safety zone.

#define dtcontrol              FPAR_1              ' time increment of the control loop
#define dtdata                 FPAR_2              ' Data to be controlled. Ask Van Thuan

#define dctrl1                 FPAR_3              ' Variable containing the control displacement in direction 1
#define dctrl2                 FPAR_4              ' Variable containing the control displacement in direction 2

#define PID_P                  FPAR_12             ' Proportional part of the PID controller for the hydraulic cylinder
#define PID_I                  FPAR_13             ' Integration constant of the PID controller for the hydraulic cylinder
#define PID_D                  FPAR_14             ' Derivative constant of the PID controller for the hydraulic cylinder

#define EVENT_LOOP_TIME_P1     FPAR_26             ' Time consumed during an Event in process 1.

#define vmin1                  FPAR_35             ' Minimum velocity in direction 1
#define vmax1                  FPAR_36             ' Maximum velocity in direction 1
#define vmin2                  FPAR_37             ' Minimum velocity in direction 2
#define vmax2                  FPAR_38             ' Maximum velocity in direction 2
#define amin1                  FPAR_39             ' Minimum acceleration in direction 1
#define amax1                  FPAR_40             ' Maximum acceleration in direction 1
#define amin2                  FPAR_41             ' Minimum acceleration in direction 2
#define amax2                  FPAR_42             ' Maximum acceleration in direction 2
#define dmin1                  FPAR_51             ' Minimum displacement in direction 1
#define dmax1                  FPAR_52             ' Maximum displacement in direction 1
#define dmin2                  FPAR_53             ' Minimum displacement in direction 2
#define dmax2                  FPAR_54             ' Maximum displacement in direction 2

#define dmeas1                 FPAR_72             ' Measured displacement in direction 1. Used for inter-process communication with the sub-stepping routine
#define dmeas2                 FPAR_73             ' Measured displacement in direction 2. Used for inter-process communication with the sub-stepping routine
#define ameas1                 FPAR_74             ' Measured acceleration in direction 1. Used for inter-process communication with the sub-stepping routine
#define ameas2                 FPAR_75             ' Measured acceleration in direction 2. Used for inter-process communication with the sub-stepping routine
#define vmeas1                 FPAR_76             ' Measured velocity in direction 1. Used for inter-process communication with the sub-stepping routine
#define vmeas2                 FPAR_77             ' Measured velocity in direction 2. Used for inter-process communication with the sub-stepping routine

dim isubcontrolstep as long                   ' Counter
dim fsubcontrolstep as long                   ' Counter

dim Disp0_CylT[Num_Channel] as float          ' Backup the different
dim Disp_CylT[Num_Channel] as float           ' Array to store the current displacement of the actuators (after filtering).
dim Acc_CylT[Num_Channel] as float            ' Array to store the current acceleration of the actuators. This is a calculated value and it is used for safety reasons mainly.
dim Vel_CylT[Num_Channel] as float            ' Array to store the current velocity of the actuators. This is a calculated value and it is used for safety reasons mainly.
dim Error_Prop_Old[Num_Channel] as float
dim Error_Prop[Num_Channel] as float
dim Error_Int[Num_Channel] as float
dim Error_Dif[Num_Channel] as float
dim Cntrl_Current[Num_Channel] as float

dim His_Disp[Length_His_Disp] as float
dim His_Vel[Length_His_Vel] as float
dim Num_Stored_Data as long

dim Disp_Conversion[Num_Channel] as float     ' Conversion factor for displacement
dim Current_Conversion[Num_Channel] as float  ' Conversion factor for current

dim dmmin[Num_Channel] as float               ' Array to store the minimum displacement values in order to operate in the safe zone
dim dmmax[Num_Channel] as float               ' Array to store the maximum displacement values in order to operate in the safe zone
dim vmmin[Num_Channel] as float               ' Array to store the minimum velocity values in order to operate in the safe zone
dim vmmax[Num_Channel] as float               ' Array to store the maximum velocity values in order to operate in the safe zone
dim ammin[Num_Channel] as float               ' Array to store the minimum acceleration values in order to operate in the safe zone
dim ammax[Num_Channel] as float               ' Array to store the maximum acceleration values in order to operate in the safe zone

dim Target_Disp[Num_Channel] as float         ' Array to store the target displacements of the cylinders
dim Cntrl_Disp_Old[Num_Channel] as float      ' Array to store the old control displacement
dim Cntrl_Disp[Num_Channel] as float          ' Array to store the control displacement. Value modified by the sub-stepping process

dim AD_Input[Num_Channel] as short            ' Array to store the AD input module
dim Channel_Input[Num_Channel] as short       ' Array to store the channels for input
dim DA_Output[Num_Channel] as short           ' Array to store the DA output module
dim Channel_Output[Num_Channel] as short      ' Array to store the channels for output

dim ramp as float

' PID Control variables
dim KPdis[Num_Channel] as float
dim KIdis[Num_Channel] as float
dim KDdis[Num_Channel] as float
dim Tidis[Num_Channel] as float
dim Tddis[Num_Channel] as float

' Variables to control the process time, working time etc
dim Time_Event_Init as long     ' Time at the start of an EVENT
dim Pressure_Transit_Time, Pressure_Transit_Time_0, Pressure_Transit_Time_1 as long

INIT:
  Init_Variables()
  Setled(1,ad,1)
  Setled(2,ad,1)
  Setled(3,ad,1)
  Setled(1,da,1)
  Setled(2,da,0)
EVENT:
  ' Get the time when the event starts
  Time_Event_Init = read_timer()
  
  isubcontrolstep = isubcontrolstep + 1
  
  ' Check if ADwin is Ready
  Check_ADwinReady()
  
  If(EMERGENCY_STOP = 1) Then
    EmergencyStop()
  Else
    ' Check the pressure control
    Check_PressureControl() ' Partial
    
    'Update Limits each event to match those set in the C++ interface
    Update_Limits()
    
    ' Measure the response of the system
    Measure_Response() ' Partial
    
    ' Check the limits in order to operate in the safe zone if the variable Check_Limits_Var is set to true.
    if(CHECK_LIMITS_VAR = 1) then
      Check_Limits()
    endif
    
    ' Update the new measured values
    Update_New_Values()
    
    Get_New_Input_Data()
    
    PID_Control_Cylinder()
      
    EVENT_LOOP_TIME_P1 = read_timer() - Time_Event_Init
    
  endif
FINISH:
  Setled(1,ad,0)
  Setled(2,ad,0)
  Setled(3,ad,0)
  Setled(1,da,0)
  Setled(2,da,0)
  
SUB Init_Variables()
  
  dim k as Long
  dtcontrol = Default_dtcontrol                ' The dtcontrol value is set to be the default value
  dtdata = Default_dtdata                      ' The dtdata value is set to be the default value, which is a large value since at the
  Rem                                            beginning of the test everything should move slowly. This value is modified to the
  Rem                                            substep time value once the Test is started in the sub-stepping process.
  PROCESSDELAY = 300000*Default_dtcontrol*1000 ' 1 ms = 300000 Clocks
  
  Num_Stored_Data = (dtcompVA / dtcontrol) + 1
  if(Num_Stored_Data > Length_His_Disp) then
    Num_Stored_Data = Length_His_Disp
  endif
  if(Num_Stored_Data < 2) then
    Num_Stored_Data = 2
  endif
  
  isubcontrolstep = 0                           ' Initialise the counter
  
  ' Set the default maximum/minimum acceptable values of displacement, velocity and acceleration.
  ' Units are all in International system.
  dmin1 = -0.14                                 ' [m]
  dmax1 = 0.14                                  ' [m]
  vmin1 = -0.5                                  ' [m/s]  the correct velocity of the cylinder is around 0.2 m/s.
  vmax1 = 0.5                                   ' [m/s]
  amin1 = -50.0                                 ' [m/s^2] 
  amax1 = 50.0                                  ' [m/s^2]
  
  dmin2 = dmin1
  dmax2 = dmax1
  vmin2 = vmin1
  vmax2 = vmax1
  amin2 = amin1
  amax2 = amax1
  
  
  Pressure_Transit_Time = 0
  DIGPROG1(1,11111111b)                         '(byte 1 out, byte 2 in)
  DIGPROG2(1,0b)        '(byte 3 in, byte 4 in)
  DIGOUT_F(1,0,0)
  
  ' Initialise various flags of the program.
  ' The values of the switches (Y1 and Y2) to control the pressure in the hydraulic block are set to 0. No pressure.
  Y1 = 0                                       ' Value of switch 1 set to 0.
  Y2 = 0                                       ' Value of switch 2 set to 0.
  DIGOUT_F(1,1,Y1)  
  DIGOUT_F(1,2,Y2)
  CURRENT_PRESSURE = 0                         ' No hydraulic pressure.
  NEW_PRESSURE_COMMAND = 0                     ' Initial value. Nothing is done.
  
  ADWIN_READY = 0                              ' ADwin is not yet ready. It should be readied through the control interface
  CONTROLBOX_READY = 0                         ' The control box is not yet ready.
  ERROR_FLAG = 0                               ' There are no errors during initialisation... and hopefully none during the test
 
  NEW_DATA_AVAILABLE = 1                       ' There are new data available at the beginning of the process
  
  TESTING = 0                                  ' No tests are performed yet. This value will be changed by the sub-stepping process
  
  ' Initialisation of control variables
  dctrl1 = 0.0
  dctrl2 = 0.0
  
  For k=1 To Num_Channel
    Disp0_CylT[k] = 0.0
    Disp_CylT[k] = 0.0
    Vel_CylT[k] = 0.0
    Acc_CylT[k] = 0.0
    Target_Disp[k] = 0.0
    Cntrl_Disp_Old[k] = 0.0
    Cntrl_Disp[k] = 0.0
    Error_Prop_Old[k] = 0.0
    Error_Prop[k] = 0.0
    Error_Int[k] = 0.0
    Error_Dif[k] = 0.0
    Cntrl_Current[k] = 0.0

    AD_Input[k] = 1                            ' AD input modules
    Channel_Input[k] = k                       ' Channel input
    DA_Output[k] = 1                           ' DA output modules
    Channel_Output[k] = k                      ' Channel Output

    ' Sensitivities    
    ' For mere information Acceleration convertion would be 390.69/32768 +/- 390.69m/s2 for 10V = 32768, 250mv/g or 9.8065m/s2 
    Disp_Conversion[k] = 0.2/32768.0           ' +/- 0.2m = 10V = 32768
    Current_Conversion[k] = 30.0/32768.0       ' +/- 30mA  = 10V = 32768
  Next k 
  
  For k = 1 To Length_His_Disp
    His_Disp[k] = 0.0
    His_Vel[k] = 0.0
  Next k
ENDSUB
SUB Check_ADwinReady()
  If(ADWIN_READY = 1) Then
    DIGOUT_F(1,0,1)                            ' Set channel 0 of module 1 to high level. Channel 0 is used to allow the control box to manually be set to ready (green LED)
    EMERGENCY_STOP = 0                         ' Set emergency stop to 0
    ADWIN_READY = 0                            ' ADWIN_READY set to 0 so that this function is not called anymore 
    CONTROLBOX_READY = 1                       ' The Control box is ready.
  EndIf
ENDSUB

' Emergency stop routine
SUB EmergencyStop()
  ' Set the hydraulic pressure to 0. System is stopped suddenly although the hydraulic cylinder
  ' will still move a little due to inercy.
 
  ' The values of the switches to control the pressure in the hydraulic block are set to 0. No pressure.
  Y1 = 0                                       ' Set the value of switch 1 in the hydraulic block to 0.
  Y2 = 0                                       ' Set the value of switch 2 in the hydraulic block to 0.
  
  DIGOUT_F(1,0,0)                              ' Set the channel to the control box to low level
  DIGOUT_F(1,1,Y1)                             ' Set the output of the first hydraulic block switch to Y1 (low level in this case)
  DIGOUT_F(1,2,Y2)                             ' Set the output of the second hydraulic block switch to Y2 (low level in this case)
  
  CURRENT_PRESSURE = 0                         ' As a result the pressure is set to 0.
ENDSUB

' Routine to control the pressure in the hydraulic block.
SUB Check_PressureControl()
  SelectCase NEW_PRESSURE_COMMAND
    Case 0                                     ' This is the do nothing case
      NEW_PRESSURE_COMMAND = 0
    Case 1                                     ' In this case low pressure is required.
      SelectCase CURRENT_PRESSURE              ' Evaluate what is the current pressure in the hydraulic block and act accordingly.
        Case 0                                 ' The pressure should change from 'No pressure' to 'Low pressure'
          Y2 = 1                               ' The second switch of the hydraulic block is set to 1.
          Y1 = 0                               ' The first switch should remain with a value of 0.
          CURRENT_PRESSURE = 1                 ' As a result the pressure is now set to 'Low pressure'
          NEW_PRESSURE_COMMAND = 0             ' Do nothing in the next events unless manually specified in the C++ interface.
          NEW_STATUS_DOUT = 1                  ' One of the digital output signals has changed.
        Case 1                                 ' Keep Low pressure
          NEW_PRESSURE_COMMAND = 0             ' This should not usually happen. Nevertheless, do nothing in this case except setting the NEW_PRESSURE_COMMAND variable to 0.
        Case 2                                 ' Change from high pressure to low pressure
          ' Since the pressure cannot be changed immediately due to the characteristics of the hydraulic block, timing
          ' functions have to be used to smooth the transition.
          If (Pressure_Transit_Time = 0) Then
            Pressure_Transit_Time_0 = read_timer()
          EndIf
          
          Y2 = 1                             ' Tell the hydraulic block that we want to change the pressure from high to low. The other
          Rem                                  switch Y1 remains unmodified until a certain time (2 seconds) has passed
          NEW_STATUS_DOUT = 1                ' One of the digital output signals has changed
          
          Pressure_Transit_Time_1 = read_timer()
          ' Make sure that the long variable does not reset to 0 due to reaching its upper limit
          If(Pressure_Transit_Time_1 > Pressure_Transit_Time_0) Then
            Pressure_Transit_Time = Pressure_Transit_Time + ((Pressure_Transit_Time_1 - Pressure_Transit_Time_0)*Clock_To_Time_Conv)
          Else
            Pressure_Transit_Time = Pressure_Transit_Time + ((Pressure_Transit_Time_1 - Pressure_Transit_Time_0)*Clock_To_Time_Conv) + 14.30
          EndIf
          Pressure_Transit_Time_0 = Pressure_Transit_Time_1 ' Backup the new value for further use in the next event
          
          If (Pressure_Transit_time < 2.0) Then
            Y1 = Y1                            ' Do nothing
          Else                                 ' Once the correct amount of time has been reached, the pressure can effectively
            Rem                                               be switched to 'Low pressure'
            Y1 = 0                             ' Tell the hydraulic block that we effectively are in low pressure by switching Y1 to 0
            CURRENT_PRESSURE = 1               ' The new pressure is set to 'Low pressure'
            NEW_STATUS_DOUT = 1                ' One of the digital output signals has changed
            NEW_PRESSURE_COMMAND = 0           ' There is not a NEW_PRESSURE_COMMAND in the next event unless manually specified
            Pressure_Transit_Time = 0          ' Reset the time counter
          Endif
      EndSelect
    Case 2                                     ' In this case high pressure is required
      SelectCase CURRENT_PRESSURE
        Case 0                                 ' The pressure cannot be changed from no pressure to high pressure in a single step.
          NEW_PRESSURE_COMMAND = 0
        Case 1                                 ' Change the pressure from 'Low pressure' to 'High Pressure'
          Rem                                    Since the pressure cannot be changed immediately due to the characteristics of the hydraulic block, timing
          Rem                                    functions have to be used to smooth the transition.
          If (Pressure_Transit_Time = 0) Then
            Pressure_Transit_Time_0 = read_timer()
          Endif
          
          Y1 = 1                             ' Tell the hydraulic block that we want to change the pressure from high to low. The other
          Rem                                  switch Y2 remains unmodified until a certain time (2 seconds) has passed
          NEW_STATUS_DOUT = 1                ' One of the digital output signals has changed
                    
          Pressure_Transit_Time_1 = read_timer()
          ' Make sure that the long variable does not reset to 0 due to reaching its upper limit
          If(Pressure_Transit_Time_1 > Pressure_Transit_Time_0) Then
            Pressure_Transit_Time = Pressure_Transit_Time + ((Pressure_Transit_Time_1 - Pressure_Transit_Time_0)*Clock_To_Time_Conv)
          Else
            Pressure_Transit_Time = Pressure_Transit_Time + ((Pressure_Transit_Time_1 - Pressure_Transit_Time_0)*Clock_To_Time_Conv) + 14.30
          EndIf
          Pressure_Transit_Time_0 = Pressure_Transit_Time_1 ' Backup the new value for further use in the next event
          
          If (Pressure_Transit_Time < 2.0) Then
            Y2 = Y2  ' Do nothing
          Else   ' Once the correct amount of time has been reached, the pressure can effectively
            Rem                                     be switched to 'High pressure'
            Y2 = 0                                ' Tell the hydraulic block that we effectively are in high pressure by switching Y2 to 0
            CURRENT_PRESSURE = 2                  ' The new pressure is set to 'High pressure'
            NEW_STATUS_DOUT = 1                   ' One of the digital output signals has changed
            NEW_PRESSURE_COMMAND = 0              ' There is not a NEW_PRESSURE_COMMAND in the next event unless manually specified
            Pressure_Transit_Time = 0             ' Reset the time counter
          endif
        Case 2                                    ' Do nothing, we are already in high pressure.
          NEW_PRESSURE_COMMAND = 0
      EndSelect
    Case 3 ' Stop the pressure.
      Y1 = 0                                      ' The first switch of the hydraulic block is set to 0.
      Y2 = 0                                      ' The second switch of the hydraulic block is set to 0.
      CURRENT_PRESSURE = 0                        ' The actual pressure variable is set to 0 since there is none.
      NEW_STATUS_DOUT = 1                         ' There is one action to be taken. Change the pressure from 'High' or 'Low' to none.
      NEW_PRESSURE_COMMAND = 0                    ' There is no new pressure command during this step.
  EndSelect
  
  ' If there has been a change in the pressure control, tell the digital signal to do so.
  If (NEW_STATUS_DOUT = 1) Then
    DIGOUT_F(1,1,Y1)                           ' Set the output of the first hydraulic block switch to Y1
    DIGOUT_F(1,2,Y2)                           ' Set the output of the first hydraulic block switch to Y1
    NEW_STATUS_DOUT = 0                        ' Reset the value to 0.
  Endif
  
ENDSUB
' Update the limits if they are changed in the C++ interface
SUB Update_Limits()
  
  ' Update displacement limits
  dmmax[1] = dmax1
  dmmax[2] = dmax2
  dmmin[1] = dmin1
  dmmin[2] = dmin2
  
  ' Update velocity limits
  vmmax[1] = vmax1
  vmmax[2] = vmax2
  vmmin[1] = vmin1
  vmmin[2] = vmin2
  
  ' Update acceleration limits
  ammax[1] = amax1
  ammax[2] = amax2
  ammin[1] = amin1
  ammin[2] = amin2

ENDSUB
Rem The original routine had the feature of geometry control. It has not been included for the moment to simplify the ADbasic code.
SUB Measure_Response()
  dim alpha as float
  dim k as Long
  
  alpha = dtcontrol/(1.0/(2*3.14159*Default_Freq_Filter) + dtcontrol)
  
  For k=1 To Num_Channel
    ' Backup the values
    Disp0_CylT[k] = Disp_CylT[k]
    ' Measure the displacement in the actuators
    Disp_CylT[k] = (ADCF(AD_Input[k],Channel_Input[k]) - Offset_DA_AD) * Disp_Conversion[k]
    
    ' If the displacement of the cylinder exceeds the maximum displacement of the table, restrict it to 20 cm
    If (Disp_CylT[k] > 0.2)  Then
      Disp_CylT[k] = 0.2 
    EndIf
    If (Disp_CylT[k] < -0.2) Then
      Disp_CylT[k] = -0.2
    EndIf
    
    ' Apply the frequency filter
    Disp_CylT[k] = Disp0_CylT[k] + alpha*(Disp_CylT[k] - Disp0_CylT[k])
  Next k
  
  ' Calculate Velocity and acceleration
  Calculate_Vel_Acc()
  
    
endsub
SUB Calculate_Vel_Acc()
  
  dim k,l as long
  
  For k=1 To Num_Channel
    For l = 1 To (Num_Stored_Data-1)
      His_Disp[((k-1)*Num_Stored_Data) + l] = His_Disp[((k-1)*Num_Stored_Data) + l + 1]
    Next l
    His_Disp[k*Num_Stored_Data] = Disp_CylT[k]
  Next k
  For k=1 To Num_Channel
    Vel_CylT[k]  = (His_Disp[((k-1)*Num_Stored_Data) + Num_Stored_Data] - His_Disp[((k-1)*Num_Stored_Data) + 1]) / (dtcontrol*(Num_Stored_Data-1.0))
  Next k
  
  For k=1 To Num_Channel
    For l = 1 To (Num_Stored_Data-1)
      His_Vel[((k-1)*Num_Stored_Data) + l] = His_Vel[((k-1)*Num_Stored_Data) + l + 1]
    Next l
    His_Vel[k*Num_Stored_Data] = Vel_CylT[k]
  Next k
  
  for k=1 to Num_Channel
    Acc_CylT[k]  = (His_Vel[((k-1)*Num_Stored_Data) + Num_Stored_Data] - His_Vel[((k-1)*Num_Stored_Data) + 1]) / (dtcontrol*(Num_Stored_Data-1.0))
  next k
  
Endsub
SUB Check_Limits()
  dim k as long
  
  ' There is no need to repeat this routine if some error in the previous events is found
  If(ERROR_FLAG = NO_ERROR) Then
    ' Check for maximum displacements, velocities and accelerations of the cylinder
    for k=1 to Num_Channel
      If((Disp_CylT[k] > dmmax[k]) Or (Disp_CylT[k] < dmmin[k])) Then
        EMERGENCY_STOP = 1
        If(Disp_CylT[k] > dmmax[k]) Then
          ERROR_FLAG = ERROR_OVER_DMAX
        Else
          ERROR_FLAG = ERROR_UNDER_DMIN 
        EndIf
      EndIf
      
      If((Vel_CylT[k] > vmmax[k]) Or (Vel_CylT[k] < vmmin[k])) Then
        EMERGENCY_STOP = 1
        If(Vel_CylT[k] > vmmax[k]) Then
          ERROR_FLAG = ERROR_OVER_VMAX
        Else
          ERROR_FLAG = ERROR_UNDER_VMIN 
        EndIf
      EndIf
      
      If((Acc_CylT[k] > ammax[k]) Or (Acc_CylT[k] < ammin[k])) Then
        EMERGENCY_STOP = 1
        If(Acc_CylT[k] > ammax[k]) Then
          ERROR_FLAG = ERROR_OVER_AMAX
        Else
          ERROR_FLAG = ERROR_UNDER_AMIN 
        EndIf
      EndIf
    Next k
    
  EndIf    
ENDSUB
SUB Update_New_Values()  
  dmeas1 = Disp_CylT[1]
  dmeas2 = Disp_CylT[2]
  
  ameas1 = Acc_CylT[1]
  ameas2 = Acc_CylT[2]
  
  vmeas1 = Vel_CylT[1]
  vmeas2 = Vel_CylT[2]
EndSub
SUB Get_New_Input_Data()
  dim k as long
  
  If (isubcontrolstep > fsubcontrolstep) Then
    NEW_DATA_AVAILABLE = 1
  Endif
  
  If (NEW_DATA_AVAILABLE = 1) Then
    ' Backup data
    For k = 1 To Num_Channel
      Cntrl_Disp_Old[k] = Cntrl_Disp[k]
    Next k
    
    ' Assign new values
    Cntrl_Disp[1] = dctrl1                     ' Assign the displacement values coming from the sub-stepping process
    Cntrl_Disp[2] = dctrl2

    isubcontrolstep = 1                        ' Start a new control step
    NEW_DATA_AVAILABLE = 0                     ' Reset the value. This value is changed in the sub-stepping process.   
  Endif
  
  ' Assign the new displacements to the controller
  fsubcontrolstep = dtdata/dtcontrol
  ramp = isubcontrolstep/fsubcontrolstep
  for k = 1 to Num_Channel
    Target_Disp[k] = Cntrl_Disp_Old[k] + (Cntrl_Disp[k] - Cntrl_Disp_Old[k])*ramp
  next k
  
endsub
SUB PID_Control_Cylinder()
  dim k as long
  
  For k = 1 To Num_Channel
    KPdis[k] = PID_P
    KIdis[k] = PID_I
    KDdis[k] = PID_D
    Tidis[k] = 1
    Tddis[k] = 1
  Next k
  
  
  For k = 1 To Num_Channel
    ' Backup data
    Error_Prop_Old[k] = Error_Prop[k]
    
    ' Calculate the new proportional error
    Error_Prop[k] = Target_Disp[k] - Disp_CylT[k]
    
    ' Calculate the integration error and check that it is between the valid range
    Error_Int[k] = Error_Int[k] + Error_Prop[k]
    
    If (Error_Int[k] < -Max_Error_Int) Then
      Error_Int[k] = -Max_Error_Int
    EndIf
    
    If (Error_Int[k] > Max_Error_Int) Then
      Error_Int[k] = Max_Error_Int
    EndIf
    
    ' Calculate differentiation error
    Error_Dif[k] = (Error_Prop[k] - Error_Prop_Old[k])/dtcontrol
    
    ' Send the response    
    Cntrl_Current[k] = (KPdis[k] * Error_Prop[k] + (KIdis[k]/Tidis[k])*Error_Int[k] + KDdis[k]*Tddis[k]*Error_Dif[k])
    If (Cntrl_Current[k] < -30.0) Then
      Cntrl_Current[k] = -30.0
    EndIf
    
    If (Cntrl_Current[k] > 30) Then
      Cntrl_Current[k] = 30
    EndIf
    
    DAC(DA_Output[k], Channel_Output[k], Cntrl_Current[k]/Current_Conversion[k] + Offset_DA_AD)
  Next k
  
ENDSUB
