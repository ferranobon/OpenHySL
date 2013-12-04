'<ADbasic Header, Headerversion 001.001>
' Process_Number                 = 2
' Initial_Processdelay           = 3000
' Eventsource                    = Timer
' Control_long_Delays_for_Stop   = No
' Priority                       = Low
' Priority_Low_Level             = 1
' Version                        = 1
' ADbasic_Version                = 5.0.8
' Optimize                       = Yes
' Optimize_Level                 = 1
' Info_Last_Save                 = BASILISC  Basilisc\ferran
'<Header End>
#include ADwinPRO_ALL.inc         ' Include-file for Pro-system t11


#define timeoffset 14316.5576170    ' Time offset for correcting time calculation 
Rem                                   when the sign of clock counter changes from max to min

'======================================================================================================
'====================================                           =======================================
'==================================== GENERAL CONTROL VARIABLES =======================================
'====================================                           =======================================
'======================================================================================================
#define Num_Step               4096     ' Number of steps in the sub-structure test
#define dtstep                 0.01     ' Time increment of a step
#define Num_Substep            4        ' Number of sub-steps
#define Num_Event_Substep      250      ' Number of events per sub-step. Initially 500
#define Select_Air_Pressure    0        ' Pressure in the air chamber of the friction device
Rem                                             0                 No pressure in the air chamber
Rem                                             1                 0.5 bar air pressure
Rem                                             2                 1 bar air pressure
Rem                                             Any other value   No pressure in the air chamber    
#define UseAccelerations_As_Fc 0        ' Use accelerations for measuring the coupling force
Rem                                             0 Use load cells
Rem                                             1 Use accelerometers
#define Freq_Cut_Fc            65.0     ' Frequency cut for the coupling forces in Hz*2Pi
#define Freq_Cut_Acc           15.0     ' Frequency cut for the coupling forces in Hz*2Pi

' Order of matrices and vectors
#define Order                  1        ' Number of DOF of the sub-structure
#define Order_Input            2        ' Order of the input array = Order + 1
#define Order_Output           4        ' Order of the output array = 3*Order + 1
#define Order_Gain             1        ' Order of the gain matrix = Order*Order

' ADlog related
#define MaxSubStep             16384    ' number of of substeps in whole test (Num_Step*Num_Substep)
#define lenghtDataPC           393216   ' length of data for data record that will be read by the PC at
Rem                                     the end of the process.

'======================================================================================================
'====================================                           =======================================
'====================================    DELAY COMPENSATION     =======================================
'====================================                           =======================================
'======================================================================================================
#define WithActuator           1
#define MaxThetaOder           99
#define Maxndelay              30
#define dcompCompensation      1
#define dcomphsmax             1.0
'#define dcomphsmax            0.0
#define dcompnu                5
#define dcomplamdaValue        0.99
#define dcompndelay            12
'#define dcompUsingUrUsingforwardstep   0
#define dcomppumless           1
#define dCompStartPoint        125
#define dCompFullPoint         300

'#define dcompDecreasingStep   2000
'#define dcompStopStep         2500

#define dcompDecreasingStep    3700
#define dcompStopStep          3900

'#define dcomp_devide_Po       50.0
#define dcomp_devide_Po        35.0
#define dcompFilter            0        ' 0 or 1
'#define dcomp_fc              30.0
#define dcomp_fc               15.0

#define dcompMaxSubstepDelay   200    
#define dcompMaxHis            200      ' at least = nu*dcompndelay + 1
#define dcompfcut              20.0

'======================================================================================================
'====================================      LIMITATIONS FOR      =======================================
'====================================      EMERGENCY STOP       =======================================
'====================================      AND ERROR FLAGS      =======================================
'======================================================================================================
' Limit check
#define Avoid_Limit_Check      1        ' 0 = normal, 1 = avoiding limit checking

' Limitations for emergency stop
#define dtabmax                0.09     ' Maximum displacement of the shake table in m
#define accmax                 50.0     ' Maximum acceleration in m/s^2
#define dTMDmax                0.095    ' Maximum displacement of the TMD in m
#define Fmax                   9600.0   ' Maximum force in the load cells in N (load cells of the DFG project)
'#define Fmax                  6000.0   ' Maximum force in the load cells in N (New load cells)

' Different types of error flags. Used in combination with ERROR_FLAG variable.
#define NO_ERROR               0        ' The no error case
#define ERROR_UNDER_DMIN       1        ' The measured value is below the minimum displacement value
#define ERROR_OVER_DMAX        2        ' The measured value is over the maximum displacement value
#define ERROR_OVER_VMIN        3        ' (Not used) The measured velocity is below the minimum value
#define ERROR_OVER_VMAX        4        ' (Not used) The measured velocity is above the maximum value
#define ERROR_OVER_AMIN        5        ' (Not used) The measured acceleration is below the minimum value
#define ERROR_OVER_AMAX        6        ' (Not used) The measured acceleration is above the maximum value
#define ERROR_OVER_FMAX        7        ' The measured force is over the maximum allowed force

' Zero of the different accelerometers. Check at the begining of each test through adlog, so that it is set
' to zero.
#define acc1_zero              -0.0
#define acc2_zero              -0.0
#define acc3_zero              -0.0

' The filter has to be used when reading the pressure from the device since otherwise it shows some
' shocks. This was experienced by Van Thuan during some of the E-FAST tests.
#define compensateFilter       1.0225 

'======================================================================================================
'====================================                           =======================================
'====================================     DECLARATION OF PAR    =======================================
'====================================     AND FPAR VARIABLES    =======================================
'====================================                           =======================================
'======================================================================================================

#define NEW_DATA_AVAILABLE     PAR_4    ' Indicates whether there is new data available (value of 1) or not (value of 0).
#define EMERGENCY_STOP         PAR_10   ' Variable to stop the process suddenly
#define TESTING                PAR_12   ' The substepping process has started (value of 1) and it is currently running a test.
Rem                                       This is an inter-process communication variable.
#define ERROR_FLAG             PAR_15   ' Variable used to store the error during the check of measured values against limit values.
#define Current_Step           PAR_24   ' Counter for the number of steps
#define NumRecordedData        PAR_29   ' for mornitoring number of tested substeps Ask Van Thuan
#define Num_Data_Received      PAR_80
' Control variables
#define dtdata                 FPAR_2   ' Data to be controlled. Ask Van Thuan

#define dctrl1                 FPAR_3   ' Variable containing the control displacement in direction 1
#define dctrl2                 FPAR_4   ' Variable containing the control displacement in direction 2
#define accctrl1               FPAR_5
#define accctrl2               FPAR_6
#define velctrl1               FPAR_7
#define velctrl2               FPAR_8

' Synchronisation flags
#define SYNC_PC                FPAR_31  ' This flag is used to synchronise the uploaded data by the PC. It can have two values:
Rem                                       SYNC_PC = 0. Do nothing or the substep is performed
Rem                                             SYNC_PC = 1. New data has been uploaded from the PC.
#define SYNC_ADWIN             FPAR_32  ' This flag is used to synchronise ADwin and the PC. It can have three possible values:
Rem                                             SYNC_ADWIN = 0.0:  Do nothing
Rem                                             SYNC_ADWIN = 1.0:  Perform the sub-stepping process since new data has been uploaded
Rem                                             SYNC_ADWIN = -1.0: The sub-step process is finished and the PC can collect the new
Rem                                                                data (see Output_Data).

#define Press_pctrl            FPAR_48  ' Control pressure
#define Press_pmeas            FPAR_47  ' Measured pressure

'#define F12x                  FPAR_68  ' Measured Force in the load cell (new load cells)
'#define F12y                  FPAR_69  ' Measured Force in the load cell (new load cells)
'#define F34x                  FPAR_70  ' Measured Force in the load cell (new load cells)
'#define F34y                  FPAR_71  ' Measured Force in the load cell (new load cells)
#define dmeas1                 FPAR_72  ' Measured displacement in direction 1. Used for inter-process communication with the sub-stepping routine
#define dmeas2                 FPAR_73  ' Measured displacement in direction 2. Used for inter-process communication with the sub-stepping routine
#define ameas1                 FPAR_74  ' Measured acceleration in direction 1. Used for inter-process communication with the sub-stepping routine
#define ameas2                 FPAR_75  ' Measured acceleration in direction 2. Used for inter-process communication with the sub-stepping routine
#define vmeas1                 FPAR_76  ' Measured velocity in direction 1. Used for inter-process communication with the sub-stepping routine
#define vmeas2                 FPAR_77  ' Measured velocity in direction 2. Used for inter-process communication with the sub-stepping routine

#define Time                   FPAR_78

#define Gain_Matrix            DATA_1   ' The gain matrix.
#define Input_Data             DATA_2   ' Array to store the u0 vector from the PC. When Input_Data[1] = 1, the PC has uploaded
Rem                                       a new vector starting at Input_Data[2]. The Input_Data[1] is switched back to 0.0 when
Rem                                       ADwin reads it and start the sub-stepping.
#define Output_Data            DATA_3   ' Array to store the result of the sub-step process. Output_Data[1] can have different values (see SYNC_ADWIN).
Rem                                             Output_Data[1] = -1.0      The PC will pick the available values and continue with the stepping process.
Rem                                                                        Any other valuewill make the PC to scan this value until it is -1.0 in order
Rem                                                                        to proceed with the stepping process.
Rem                                             Output_Data[2]             Stores the uc vector resulting from the sub-stepping process
Rem                                             Output_Data[2+Order]       Stores the coupling force of the substep 'Num_Substeps -1'
Rem                                             Output_Data[2*(1 + Order)] Stores the coupling force of the last sub-step

'======================================================================================================
'====================================                           =======================================
'====================================      ADLOG DEFINITION     =======================================
'====================================      AND STORING DATA     =======================================
'====================================                           =======================================
'======================================================================================================

#define DataPC                 DATA_97
#define dataDAQ                DATA_199 ' data accquisition: substep, data data data ... substep, data, data, data, ...

#define NumberLONGsPerStep     20       ' number of longs per each step of record set for ADlog
#DEFINE ADlogData1             DATA_180
#DEFINE ADlogData2             DATA_181
#DEFINE WritePointer           ADlogData2[1]
#DEFINE LoopCounter            ADlogData2[2]
#DEFINE BufferSize             ADlogData2[3]
#DEFINE ValuesCount            ADlogData2[4]
#DEFINE Flags                  ADlogData2[5]
#DEFINE BUFSIZE                4000000 ' whole-numbered multiple of values count

DIM ADlogData1[BUFSIZE] AS LONG AT DRAM_EXTERN
DIM ADlogData2[200] AS LONG AT DM_LOCAL

dim TimeInit, TimeEnd as float

dim Gain_Matrix[Order_Gain] as float at DM_LOCAL    ' Allocate the necesssary memory for the Gain matrix. The values are
Rem                                                   suplied by the PC.
dim Input_Data[Order_Input] as float at DM_LOCAL    ' Allocate the necessary space
dim Output_Data[Order_Output] as float at DM_LOCAL  ' Allocate the necessary space

dim DataPC[lenghtDataPC] as float

dim alphaRecursiveFilter as float

dim i as long                           ' A simple counter

dim dtevent as float                    ' Time increment of each event in ADwin
dim dtsub as float                      ' Time increment of the sub-step

dim Current_Substep as long             ' Counter for the sub-steps
dim Current_Event   as long             ' Counter for the number of events in a sub-step

dim ramp  as float                      ' Variable to store the operation Current_Substep/Num_Substep.
dim ramp0 as float                      ' Variable to store the operation 1.0 - Current_Substep/Num_Substep

dim DO_RAMP_FUNCTION as short           ' Variable to identify whether the ramp function should be performed (value of 1) or not
Rem                                       (value of 0)
dim DID_RAMP_FUNCTION as short
dim DO_PREVIOUS_RAMP as short           ' Variable to identify whether the ramp function from the previous step should be performed
Rem                                       (value of 1) or not (value of 0). This only shoud be true when there has been no update from
Rem                                       the PC at the desired time.

dim vi[order] as float                  ' Stores the velocity within a sub-step. Used in the delay compensation.
dim u0[order] as float                  ' To store the data received from the PC at the previous step
dim u01[order] as float                 ' To store the data received from the PC at the beginning of the sub-step (current step)
dim uc[order] as float                  ' Vector to be sent to the PC at the end of each sub-stepping process.
dim ucprev[order] as float              ' Backup the displacement vector to compute velocity
dim fcprev[order] as float              ' Vector containing the coupling force measured at Num_Substep - 1
dim fc[order] as float                  ' Vector containing the coupling force measured at the last sub-step

dim Disp_TMDx, Disp_TMDy as float       ' Displacement of the TMD in the x and y directions respectively
dim FcoupY1, FcoupY2, FcoupX as float   ' Coupling force in the load cells (DFG load cells).
dim FcoupX1, FcoupX1Old as float        ' Auxiliary variables for filtering the coupling force in the X direction.

dim accExt1,accExt2,accExt3 as float    ' Variables to store the measurement of external accelerometers
dim alphaRecursiveFilter as float       ' Recursive filter for accelerations and coupling forces

dim Press_pctrl_InitV  as float    ' Initial value for the air pressure in the friction device

dim beginning_time,ending_time,time_between_substep as float
dim clock,clock_synPC,clock_synPC_old,clock_do_substep,clock_do_substep_old,clock_beginning_step as long
dim time_synPC,time_substep,time_step,time_do_substep,time_from_sub1_to_newPCdata as float


'delay compensation
dim dcompr,dcompdU1,dcompdU2,dcompdU1Old,dcompdu,dcompY,dcompr,dcomperr as float   
dim dcompphi2,dcompphi,dcompphiphi,dcomplamda,dcompP,dcompPold,dcompPinvertStart as float
dim dcompnTheta as long
dim dcompPHI[MaxThetaOder],dcompPHI2[MaxThetaOder],dcompPHIOLD[MaxThetaOder],dcompTHETA[MaxThetaOder],dcompTHETAOLD[MaxThetaOder],dcompK[MaxThetaOder] as float
dim dcompur1value,dcompur2value as float
dim dcompdelaydata[dcompMaxSubstepDelay] as float
dim dcompdctrl[dcompMaxSubstepDelay] as float
dim dcompdis[dcompMaxSubstepDelay] as float
dim dcompvel[dcompMaxSubstepDelay] as float
dim dcompdUold[dcompMaxSubstepDelay] as float
dim dcompHisdU[dcompMaxHis] as float
dim dcompHisU[dcompMaxHis] as float
dim dcompHisV[dcompMaxHis] as float
dim dcompHisE[dcompMaxHis] as float
dim idcomp,jdcomp,dcompnewpos as long
'dim v0vec[num],d0vec[num],a0vec[num] as float
dim dcomphs,dcompuijold as float
dim dcomperrfiltered as float

INIT:
  
  Init_Variables()
  ADlogInit()
EVENT:
  
  dtdata = dtsub                      ' Set the time increment for new data/measurements to the sub-step time. This variable
  Rem                                   is also used in the high priority process.
  
  If (TESTING = 0)  Then           ' This is only accessed before the start of the test.
    If (Input_Data[1] = 1.0) Then  ' If the PC has uploaded some data, start the test and perform the first sub-stepping process
      Current_Substep = 1          ' Initialise the sub-step counter
      Current_Step = 1             ' Initialise the step counter counter
      Current_Event = 1            ' Initialise the event counter counter
                  
      TESTING = 1                  ' Start the testing process since new data has been received from the computer.
      
      TimeInit = Read_Timer()
    Else
      ADlogPre()                    ' Acquire sensor data for further calibration.
      if((Avoid_Limit_Check=0)) then
        Check_Limitation()
      endif
      ' Do nothing else until the PC uploads the first data.
    Endif
  Else                             ' The test has started.
    If (Input_Data[1] = 1.0) Then  ' The new data is ready (SYNC_ADWIN = 1.0) and the sub-stepping process will run
      Rem                            when Current_Event = 1
      SYNC_ADWIN = 1.0             ' ADwin is synchronised since this is the first step
      Output_Data[1] = SYNC_ADWIN  ' The data is not yet ready (ADwin - PC synchronisation)
      SYNC_PC = 0.0                ' The values have been received and now the sub-step process will be performed
      Input_Data[1] = SYNC_PC      ' Inform the PC about this change in the process

      Num_Data_Received = Num_Data_Received + 1
      For i = 1 to Order
        u01[i] = Input_Data[i+1]   ' Copy the new data coming from the PC into u01, the new target displacement
      Next i
    EndIf
    
    If (SYNC_ADWIN = 1.0) Then
      If (Current_Event = 1) Then  ' Only do the ramp function if the number of events is 1 for the current sub-step.
        DO_RAMP_FUNCTION = 1       ' The sub-stepping should only take place at the begining of the scanning loop in order
        Rem                          to make sure that the process runs at the required time. Otherwise the test results
        Rem                          would not be coherent since the substeps will not take place in the specified time
        Rem                          increments.
      EndIf
    EndIf
    
    If ((SYNC_ADWIN = -1.0) And (Current_Substep = 1)) Then
      If (Current_Event = 1) Then  ' If the PC has not uploaded any data at the begining of the step, then the previous ramp
        Rem                          function will be used for the first step. If after this the data is not available, there
        Rem                          will not be any more displacement updates until a new value arrives.
        DO_RAMP_FUNCTION = 1
        DO_PREVIOUS_RAMP = 1
      EndIf
    EndIf 
    
    If (DO_RAMP_FUNCTION = 1) Then
      
      If(dcompCompensation=1) Then
        Process_dcomp_Before()
      EndIf
      
      ' Update the ramp function values
      ramp = Current_Substep / Num_Substep
      ramp0= 1.0 - ramp
      
      ucprev[1] = uc[1]               ' Backup the displacement in order to calculate the velocity.
      If (DO_PREVIOUS_RAMP = 0) Then       
        'mul_const_vector(ramp0,u0,order,vec1)
        'mul_const_vector(ramp,u01,order,vec2)
        'sum_vector_vector(vec1,order,vec2,vec3)
        'mul_matrix_vector(Gain,order,order,fc,vec1)
        'sum_vector_vector(vec3,order,vec1,uc)
        uc[1] = ramp0 * u0[1] + ramp * u01[1] + Gain_Matrix[1]*fc[1]
      Else
        ' Continue with the previous ramp function for one more substep
        ramp = (Num_Substep + 1)/Num_Substep
        ramp0 = 1.0 - ramp
        uc[1] = ramp0 * u0[1] + ramp * u01[1] + Gain_Matrix[1]*fc[1]
        DO_PREVIOUS_RAMP = 0
      EndIf
           
      vi[1] = (uc[1]-ucprev[1])/dtsub  ' Velocity for delay compensation 
       
      Write_Displacement_Signal()
      
      DO_RAMP_FUNCTION = 0             ' Do not perform any more displacement updates until the next cycle of events to assure that the
      Rem                                specified amount of time has passed between to consecutive sub-steps
      DID_RAMP_FUNCTION = 1
    EndIf
 
    If( DID_RAMP_FUNCTION = 1 ) Then
      ' Measure the coupling forces at the end of substep
      If (Current_Event = Num_Event_Substep) Then
        ' Use the load cells to measure the force
        ForceMeasurement()
        if (UseAccelerations_As_Fc = 0) then        
          Coupling_Force_LoadCells()
          fc[1] =  -FcoupX1
        else
          Coupling_Force_Accelerations()
          fc[1] = -FcoupX1
        endif
      
        ' Backup the coupling force vector at the Num_Substep -1 substep.
        if ((Current_Substep = (Num_Substep -1)) or (Num_Substep = 1)) then
          scopy(fc, Order, fcprev)
        endif
      
        ' Check the values for safety reasons
        Check_Limitation()
      
        If(dcompCompensation=1) Then
          Process_DComp_After()
        endif      
      
        ADlog()
      
        If (Current_Substep = Num_Substep) Then
          ' Copy the data to Output_Data so that the PC can read it
          for i = 1 to order
            Output_Data[1+i] = uc[i]
            Output_Data[1+ (1*order)+i] = fcprev[i]
            Output_Data[1 +(2*order)+i] = fc[i]
          next i
        
          SYNC_ADWIN = -1.0
          Output_Data[1] = SYNC_ADWIN
        
          PressureControl()
        
          Current_Substep = 1                    ' Reset the sub-step counter
          Current_Step = Current_Step + 1
          scopy(u01,Order,u0)                    ' Copy the vector received from the PC during this step to u0 in order to compute
          Rem                                      uc in the next step         
        Else
          Current_Substep = Current_Substep + 1  ' Increase the sub-step counter
        Endif
      
        Current_Event = 1                        ' Reset the event counter
        DID_RAMP_FUNCTION = 0
      Else
        Current_Event = Current_Event + 1        ' Increase the event counter
      EndIf
    EndIf
  EndIf
  
  If (Current_Step > Num_Step) Then
    TimeEnd = Read_Timer()
    Time = (TimeEnd - TimeInit)*((10.0/3.0)*1.0e-6)
    End   ' Finish the Process
  EndIf
  
  
FINISH:
  
  dtdata = 5.0    ' Reset the new data/measurement time increment to a larger value in seconds
  Rem               to move the hydraulic cylinders slowlier.

  DAC(1,3,((0.0 / 2.0)*32768) + 32768)  ' Set the pressure in the friction device to zero
  
  Rem -----------------------------------------------------------------------------------------------------------------------
  Rem ------------------------------------------------- SUB-ROUTINES SECTION ------------------------------------------------
  Rem -----------------------------------------------------------------------------------------------------------------------
  
SUB Init_Variables()
  dtsub = dtstep/Num_Substep         ' Defining the sub-step time increment
  dtevent = dtsub/Num_Event_Substep  ' Defining the time increment of each event
  PROCESSDELAY = 300000000*dtevent   ' Define the process delay
  
  TESTING = 0                        ' Not yet performing a substructure test
  Num_Data_Received = 0
  Current_Substep = 1                ' Set the sub-step counter to its starting value
  Current_Event = 1                  ' Set the event counter to its starting value
  Current_Step = 1
  
  ERROR_FLAG = NO_ERROR              ' No error
  TESTING = 0
  
  DO_RAMP_FUNCTION = 0               ' The ramp function should not be performed until the first update arrives
  DO_PREVIOUS_RAMP = 0               ' There must not be any actuator movement at the begining.
  
  ' Initialise the vectors to 0.0
  Set2Value(Gain_Matrix, Order_Gain, 0.0)
  Set2Value(Input_Data, Order_Input, 0.0)
  Set2Value(Output_Data, Order_Output, 0.0)
  Set2Value(u0, Order, 0.0)
  Set2Value(u01, Order, 0.0)
  Set2Value(uc, Order, 0.0)
  Set2Value(fcprev, Order, 0.0)
  Set2Value(fc, Order, 0.0)
  
  ' Measuring variables are set to 0
  Disp_TMDx = 0.0
  Disp_TMDy = 0.0
  
  FcoupX = 0.0
  FcoupX1Old = 0.0
  FcoupY1 = 0.0
  FcoupY2 = 0.0

  fc[1] = 0.0
  fcprev[1] = 0.0
  ucprev[1] = 0.0
  
  accExt1 = 0.0
  accExt2 = 0.0
  accExt3 = 0.0
  
  idcomp= 1
  jdcomp= 1
  dcomperrfiltered = 0
  
  dcompInit()

  NumRecordedData = 0

  SelectCase Select_Air_Pressure
    Case 0
      Press_pctrl_InitV = 0.0    ' No pressure in the friction device
    Case 1
      Press_pctrl_InitV = 0.05   ' Pressure of 0.5bar in the friction device
    Case 2
      Press_pctrl_InitV = 0.25   ' Pressure of 0.25bar in the friction device
    CaseElse
      Press_pctrl_InitV = 0.0    ' For any other value, the pressure should be 0 for safety reasons
  EndSelect
  
EndSub
SUB ADlogInit()
  'Init ADLog
  ' write pointer
  WritePointer = 1 'ADLogData2[1]
  ' buffer overleap counter
  LoopCounter = 0 'ADLogData2[2]
  ' buffer size
  BufferSize = BUFSIZE 'ADLogData2[3]
  ' values count (number of LONGs per step)
  ValuesCount = NumberLONGsPerStep 'ADLogData2[4]
  ' Bit 0 : use LoopCounter
  Flags = 1 'ADLogData2[5]
ENDSUB
SUB Write_Displacement_Signal()
  
  dctrl1 = 0               ' Direction 1 is never used in these tests 
  dctrl2 = uc[1]' + dcompdU * dcomphs  ' Set the control displacement to the value of uc unless it exceeds the table displacement
  
  if(dctrl2 > dtabmax) then
    dctrl2 = dtabmax
  endif
  if(dctrl2 < -dtabmax) then
    dctrl2 = -dtabmax
  endif
  
  NEW_DATA_AVAILABLE = 1        ' Tell the other process that new data is available
EndSub
SUB ForceMeasurement()
  'MeasureCouplingForce
  FcoupY1 = -((ADCF(1,3) - 32768.0) / 32768.0)* 10000 'N
  FcoupY2 = -((ADCF(1,4) - 32768.0) / 32768.0)* 10000 'N
  FcoupX  =  ((ADCF(1,5) - 32768.0) / 32768.0)* 10000 'N  
ENDSUB
SUB Coupling_Force_LoadCells()
  alphaRecursiveFilter  = (dtstep/Num_Substep) / ((1.0/Freq_Cut_Fc) + (dtstep/Num_Substep))
  FcoupX1 =  FcoupX1Old + alphaRecursiveFilter*(FcoupX - FcoupX1Old)
  FcoupX1Old = FcoupX1
ENDSUB
SUB Coupling_Force_Accelerations()
  dim temp, temp1 as float
  dim AccTMD as float
  dim AccFrame as float
  
  alphaRecursiveFilter  = (dtstep/Num_Substep) / ((1.0/(Freq_Cut_Acc*2*3.14)) + (dtstep/Num_Substep))
  
  ' TMD. x-direction
  temp = - (((ADCF(3, 1)-32768)/32768.0)*10.0 / 0.936) * 9.8065 - acc2_zero       ' PCB accelerometer - TMD
  temp1 =  accExt2 + alphaRecursiveFilter*(temp - accExt2)
  accExt2 = - (((ADCF(3, 1)-32768)/32768.0)*10.0 / 0.936) * 9.8065 - acc2_zero       ' PCB accelerometer - TMD temp1
  
  AccTMD = accExt2
  
  ' x direction of the frame
  temp = - (((ADCF(2, 8)-32768)/32768.0)*10.0 / 0.507 ) * 9.8065 - acc3_zero     ' Kistler 8640A
  temp1 =  accExt3 + alphaRecursiveFilter*(temp - accExt3)
  accExt3 = - (((ADCF(2, 8)-32768)/32768.0)*10.0 / 0.507 ) * 9.8065 - acc3_zero     ' Kistler 8640A temp1
  
  AccFrame = accExt3
  
  FcoupX1 = AccFrame*600.0 + AccTMD*298
  
ENDSUB
SUB Check_Limitation()
  
  if(ERROR_FLAG = 0) then
    ' Forces are checked in the high priority process
  
    ' Displacements
    if((dmeas1 > dtabmax) or (dmeas1 < -dtabmax)) then
      EMERGENCY_STOP = 1
      if(dmeas1> dtabmax) then 
        ERROR_FLAG = ERROR_OVER_DMAX
      else
        ERROR_FLAG = ERROR_UNDER_DMIN
      endif
    endif
    
    if((dmeas2 > dtabmax) or (dmeas2 < -dtabmax)) then
      EMERGENCY_STOP = 1
      if(dmeas2 > dtabmax) then 
        ERROR_FLAG = ERROR_OVER_DMAX
      else
        ERROR_FLAG = ERROR_UNDER_DMIN
      endif
    endif
    
    if((Disp_TMDx > dTMDmax) or ( Disp_TMDx < -dTMDmax)) then
      EMERGENCY_STOP = 1
      if(Disp_TMDx > dtabmax) then 
        ERROR_FLAG = ERROR_OVER_DMAX
      else
        ERROR_FLAG = ERROR_UNDER_DMIN
      endif
    endif
    
    if((Disp_TMDy > dTMDmax) or (Disp_TMDy < -dTMDmax)) then
      EMERGENCY_STOP = 1
      if(Disp_TMDy > dtabmax) then 
        ERROR_FLAG = ERROR_OVER_DMAX
      else
        ERROR_FLAG = ERROR_UNDER_DMIN
      endif
    endif
    
    if(absf(FcoupY1) > Fmax) then 
      EMERGENCY_STOP = 1
      ERROR_FLAG = ERROR_OVER_FMAX
    endif
    
    if(absf(FcoupY2) > Fmax) then 
      EMERGENCY_STOP = 1
      ERROR_FLAG = ERROR_OVER_FMAX
    endif
    
    if(absf(FcoupX) > Fmax) then 
      EMERGENCY_STOP = 1
      ERROR_FLAG = ERROR_OVER_FMAX
    endif  
  endif
EndSub
SUB  ADlogPre()
  
  dim data as LONG
  dim temp1 as float
  
  data = 0                                     'trigger
  ADlogData1[WritePointer] = data
  
  Forcemeasurement()

  GetAndStoreData_Pre()
  
  WritePointer = WritePointer + NumberLONGsPerStep
  
  ' check writepointer
  IF (WritePointer > BUFSIZE) THEN
    WritePointer = 1
    INC LoopCounter
  ENDIF

  DAC(1,3,((Press_pctrl_InitV / 2.0)*32768) + 32768)  ' Initialise the pressure in the device
EndSub
SUB  ADlog()
  dim data as LONG
  dim temp1 as float
  
  data = 2147483647                                 'trigger
  ADlogData1[WritePointer] = data
  
  GetAndStoreData()
  'GetAndStoreData_withDAQ()
  
  WritePointer = WritePointer + NumberLONGsPerStep
  
  ' check writepointer
  IF (WritePointer > BUFSIZE) THEN
    WritePointer = 1
    INC LoopCounter
  ENDIF
  
  NumRecordedData = NumRecordedData+1
endsub
SUB GetAndStoreData_Pre()
  dim data,dtab,temp as float

  dim i,j,length as long
  dim numchannelDAQ as long ' number of channels for DAQ
  
  
  ' TMD. y-direction
  temp =  - (((ADCF(2, 7)-32768)/32768.0)*10.0 / 0.104) * 9.8065 - acc1_zero      ' Endevco 61-100 Acc tabe                 
  accExt1 = temp

  ' TMD. x-direction
  temp = - (((ADCF(3, 1)-32768)/32768.0)*10.0 / 0.936) * 9.8065 - acc2_zero       ' PCB accelerometer - TMD
  accExt2 = temp
  
  ' x direction of the frame
  temp = - (((ADCF(2, 8)-32768)/32768.0)*10.0 / 0.507 ) * 9.8065 - acc3_zero     ' Kistler 8640A
  accExt3 = temp    
  
  Disp_TMDx = - ((ADCF(3, 2)-32768)/32768.0) * 0.1 
  Disp_TMDy = - ((ADCF(2, 5)-32768)/32768.0) * 0.1 
  
  data =  (dctrl1 /0.2)* 2147483647                'dctrl1
  ADlogData1[WritePointer+1] = data
  
  data =  (dctrl2 /0.2)* 2147483647                'dctrl2
  ADlogData1[WritePointer+2] = data
  
  data =  (dmeas1 /0.2)* 2147483647                'dmeas1
  ADlogData1[WritePointer+3] = data
  
  data =  (dmeas2 /0.2)* 2147483647                'dmeas2
  ADlogData1[WritePointer+4] = data
  
  data =  (accctrl1 / (9.8065 * 10.0)) * 2147483647     'accctrl1
  ADlogData1[WritePointer+5] = data
  
  data =  (accctrl2 / (9.8065 * 10.0)) * 2147483647     'accctrl2
  ADlogData1[WritePointer+6] = data
  
  data =  (ameas1 / (9.8065 * 10.0)) * 2147483647     'ameas1 by computing dis
  ADlogData1[WritePointer+7] = data
  
  data =  (ameas2 / (9.8065 * 10.0)) * 2147483647     'ameas2 by computing dis
  ADlogData1[WritePointer+8] = data
  
  data =  (accExt1 / (9.8065 * 10.0)) * 2147483647     'acc1
  ADlogData1[WritePointer+9] = data
  
  data =  (accExt2 / (9.8065 * 10.0)) * 2147483647     'acc2
  ADlogData1[WritePointer+10] = data
  
  data =  (accExt3 / (9.8065 * 10.0)) * 2147483647     'acc3
  ADlogData1[WritePointer+11] = data
  
  data = (FcoupY1 / 10000.0) * 2147483647               'FcY1 (N)
  ADlogData1[WritePointer+12] = data
  
  data = (FcoupY2 / 10000.0) * 2147483647               'FcY2 (N)
  ADlogData1[WritePointer+13] = data
  
  data = (FcoupX1  / 10000.0) * 2147483647                'FcX
  ADlogData1[WritePointer+14] = data
  
  data = (Disp_TMDx / 0.1) * 2147483647               'D tmdx (mm)
  ADlogData1[WritePointer+15] = data
  
  data = (Disp_TMDy / 0.1) * 2147483647               'D tmdy (mm)
  ADlogData1[WritePointer+16] = data
  
  data = Press_pctrl/2.0 * 2147483647  'PressureCtrl
  ADlogData1[WritePointer+17] = data
  
  Press_pmeas = (((ADCF(3,7) - 32768) * compensateFilter) / 32768.0) * 2.0
  data = (Press_pmeas/2.0 ) * 2147483647 
  ADlogData1[WritePointer+18] = data
  
  data =  (uc[1] / 0.1) * 2147483647     ' disi1
  ADlogData1[WritePointer+19] = data
  
ENDSUB
SUB GetAndStoreData()
  dim data,dtab,temp as float
  dim dpisint as long
  dim i,j,length as long
  dim numchannelDAQ as long ' number of channels for DAQ
  
  
  'alphaRecursiveFilter  = (dtstep/Num_Substep) / ((1.0/(Freq_Cut_Acc*2*3.14)) + (dtstep/Num_Substep))
  
  if (UseAccelerations_As_Fc = 0) then
    ' TMD. y-direction
    'temp =     
    'temp1 =  accExt1 + alphaRecursiveFilter*(temp - accExt1)  
    accExt1 = - (((ADCF(2, 7)-32768)/32768.0)*10.0 / 0.104) * 9.8065 - acc1_zero      ' Endevco 61-100 Acc tabe temp1

    ' TMD. x-direction
    'temp = - (((ADCF(3, 1)-32768)/32768.0)*10.0 / 0.936) * 9.8065 - acc2_zero       ' PCB accelerometer - TMD
    'temp1 =  accExt2 + alphaRecursiveFilter*(temp - accExt2)
    accExt2 = - (((ADCF(3, 1)-32768)/32768.0)*10.0 / 0.936) * 9.8065 - acc2_zero       ' PCB accelerometer - TMD temp1
  
    ' x direction of the frame
    'temp = - (((ADCF(2, 8)-32768)/32768.0)*10.0 / 0.507 ) * 9.8065 - acc3_zero     ' Kistler 8640A
    'temp1 =  accExt3 + alphaRecursiveFilter*(temp - accExt3)
    accExt3 = - (((ADCF(2, 8)-32768)/32768.0)*10.0 / 0.507 ) * 9.8065 - acc3_zero     ' Kistler 8640A temp1
  else
    ' TMD. y-direction
    'temp =  - (((ADCF(2, 7)-32768)/32768.0)*10.0 / 0.104) * 9.8065 - acc1_zero      ' Endevco 61-100 Acc tabe   
    ' temp1 =  accExt1 + alphaRecursiveFilter*(temp - accExt1)  
    accExt1 = - (((ADCF(2, 7)-32768)/32768.0)*10.0 / 0.104) * 9.8065 - acc1_zero      ' Endevco 61-100 Acc tabe    temp1
  endif
      
  Disp_TMDx = - ((ADCF(3, 2)-32768)/32768.0) * 0.1 
  Disp_TMDy = - ((ADCF(2, 5)-32768)/32768.0) * 0.1
  

  length = 24
  i = (Current_Step-1)*Num_Substep + Current_Substep-1
  j = (((Current_Step-1)*Num_Substep + Current_Substep)  - 1)*length
  
  DataPC[j+1] = i

  'data =  (dctrl1 /0.2)* 2147483647                'dctrl1
  data = dcompdU
  ADlogData1[WritePointer+1] = data
  DataPC[j+2] = data

  data =  (dctrl2 /0.2)* 2147483647                'dctrl2
  ADlogData1[WritePointer+2] = data
  DataPC[j+3] = dctrl2

  data =  (dmeas1 /0.2)* 2147483647                'dmeas1
  ADlogData1[WritePointer+3] = data
  DataPC[j+4] = dmeas1

  data =  (dmeas2 /0.2)* 2147483647                'dmeas2
  ADlogData1[WritePointer+4] = data
  DataPC[j+5] = dmeas2
  
  data =  (accctrl1 / (9.8065 * 10.0)) * 2147483647     'accctrl1
  ADlogData1[WritePointer+5] = data
  DataPC[j + 6] = accctrl1
  
  data =  (accctrl2 / (9.8065 * 10.0)) * 2147483647     'accctrl2
  ADlogData1[WritePointer+6] = data
  DataPC[j+ 7] = accctrl2
  
  data =  (ameas1 / (9.8065 * 10.0)) * 2147483647     'ameas1 by computing dis
  ADlogData1[WritePointer+7] = data
  DataPC[j+8] = ameas1
  
  data =  (ameas2 / (9.8065 * 10.0)) * 2147483647     'ameas2 by computing dis
  ADlogData1[WritePointer+8] = data
  DataPC[j+9] = ameas2
  
  data =  (accExt1 / (9.8065 * 10.0)) * 2147483647     'acc1
  ADlogData1[WritePointer+9] = data
  DataPC[j+10] = accExt1
  
  data =  (accExt2 / (9.8065 * 10.0)) * 2147483647     'acc2
  ADlogData1[WritePointer+10] = data
  DataPC[j+11] = accExt2
  
  data =  (accExt3 / (9.8065 * 10.0)) * 2147483647     'acc3
  ADlogData1[WritePointer+11] = data
  DataPC[j + 12] = accExt3
  
  data = (FcoupY1 / 10000.0) * 2147483647               'FcY1 (N)
  ADlogData1[WritePointer+12] = data
  DataPC[j+13] = FcoupY1
  
  data = (FcoupY2 / 10000.0) * 2147483647               'FcY2 (N)
  ADlogData1[WritePointer+13] = data
  DataPC[j+14] = FcoupY2
  
  data = (FcoupX1 / 10000.0) * 2147483647                'FcX
  ADlogData1[WritePointer+14] = data
  DataPC[j+15] = FcoupX
  
  data = (Disp_TMDx / 0.1) * 2147483647               'D tmdx (mm)
  ADlogData1[WritePointer+15] = data
  DataPC[j+16] = Disp_TMDx
  
  data = (Disp_TMDy / 0.1) * 2147483647               'D tmdy (mm)
  ADlogData1[WritePointer+16] = data
  DataPC[j+17] = Disp_TMDy
  
  data = Press_pctrl/2.0 * 2147483647  'PressureCtrl
  ADlogData1[WritePointer+17] = data
  DataPC[j+18] = Press_pctrl
  
  Press_pmeas = (((ADCF(3,7) - 32768) * compensateFilter) / 32768.0) * 2.0
  data = (Press_pmeas/2.0 ) * 2147483647 
  ADlogData1[WritePointer+18] = data
  DataPC[j+19] = Press_pmeas
  
  data =  (uc[1] / 0.1) * 2147483647     ' disi1
  ADlogData1[WritePointer+19] = data
  DataPC[j+20] = uc[1]
  
  DataPC[j+21] = time_do_substep
  DataPC[j+22] = time_substep
  DataPC[j+23] = time_synPC
  DataPC[j+24] = time_from_sub1_to_newPCdata
  
endsub
SUB PressureControl()
  SelectCase Select_Air_Pressure
    Case 0
      Press_pctrl = Press_pctrl_InitV
    Case 1
      Press_pctrl = Press_pctrl_InitV
    Case 2
      Press_pctrl = Press_pctrl_InitV
    CaseElse
      Press_pctrl = 0
  EndSelect
  DAC(1,3,((Press_pctrl / 2.0)*32768) + 32768)
ENDSUB
SUB dCompInit()
  
  dcompdU1old= 0.0
  dcompdU = 0.0
  dcompdU1   = 0.0
  dcompdU2   = 0.0
  dcomperr  = 0.0
  dcompPinvertStart = 0.0
  dcomplamda= dcomplamdaValue
  dcompnTheta=3*dcompnu
  for i = 1 to MaxThetaOder
    dcompPHIOLD[i]= 0.0
    dcompTHETA[i] = 0.0
    dcompTHETAOLD[i]=0.0
    dcompK[i]=0.0
  next i
        
  for i = 1 to dcompMaxSubstepDelay
    dcompdelaydata[i] = 0.0
    dcompdctrl[i]=0.0
    dcompdis[i]=0.0
    dcompvel[i]=0.0
    dcompdUold[i]=0.0
  next i
  
  for i = 1 to dcompMaxSubstepDelay
    dcompHisdU[i]=0.0
    dcompHisU[i]=0.0
    dcompHisV[i]=0.0
    dcompHisE[i]=0.0
  next i
  dcompnewpos = (dcompndelay * dcompnu) + 1

  dcompuijold = 0.0
  dcomphs = 0
ENDSUB

SUB dCompSetPHI()
  'Put new vectors to phi[]
  for i = 1 to dCompnu
    dCompPHI[i]           = -dCompHisdU[dcompnewpos-(i*dcompndelay)]
    'dCompPHI[i+dCompnu]   = dCompHisU[dcompnewpos-((i-1)*dcompndelay)]
    dCompPHI[i+1*dCompnu] = (dCompHisV[dcompnewpos-(i*dcompndelay)] )
    dCompPHI[i+2*dCompnu] = dCompHisE[dcompnewpos-(i*dcompndelay)]
  next i
ENDSUB
'=================================
SUB  dCompEstimatedU()
  dCompdU = 0.0
  if(dcomppumless = 0) then
    for i = 1 to dCompnTheta
      dCompdU = dCompdU + (dCompPHI[i] * dCompTHETA[i])
    next i
  else
    dcompr = jdcomp
    dcompr = dcompr / dcompndelay 
    for i = 1 to dCompnTheta
      dCompdU = dCompdU + (dCompPHI[i] * (dCompTHETAold[i]*(1-dcompr)+dCompTHETA[i]*dcompr))
    next i
  endif
ENDSUB
'=================================
SUB dCompUpdateModel()
  
  'backup old values
  for i = 1 to dCompnTheta
    dCompTHETAOLD[i] = dCompTHETA[i]
  next i

  'new theta
  dCompphiphi = 0.0
  for i = 1 to dcompnTheta
    dcompphiphi = dcompphiphi + (dcompPHI[i]*dcompPHI[i])
  next i
  dcompP = (dcompPold - (dcompPold*dcompPold*dcompphiphi) /(dComplamda + dCompPold*dCompphiphi)) / dComplamda
  for i = 1 to dCompnTheta
    dCompK[i] = dCompP*dCompPHI[i]
    dCompTHETA[i] = dCompTHETAOLD[i] + dCompK[i]*dComperr
  next i

  'backup old values
  for i = 1 to dCompnTheta
    dCompPHIOLD[i]   = dCompPHI[i]
  next i
  
ENDSUB
'====================================================================================================================
SUB dCompCheckInitationP()
  if(idcomp<(dCompStartPoint-1)) then
    'dcompPinvertStart = dcompPinvertStart + dcomperr*dcomperr + velvecj[1]*velvecj[1]
    dcompPinvertStart = dcompPinvertStart + dcomperr*dcomperr + vi[1]*vi[1]
    'dcompPinvertStart = dcompPinvertStart + dcomperr*dcomperr + uc[1]*uc[1]
  else
    if(idcomp=(dCompStartPoint-1))then
      dcompPold = (1.0/dcompPinvertStart) / dcomp_devide_Po
      'dcompPold = (1.0/dcompPinvertStart)
      ' larger value dcompP will more sensitive and fast convergense
      ' the devide value from 1.0 to 20.0 has been tested, ok 
      ' devide 20, result no fluctiation 
      dcompP = dcompPold
    endif
  endif
ENDSUB
SUB Process_Dcomp_Before()
  'calculate displacement compensation dcopmdU
  if(dcompCompensation = 1) then
    'Measure the initation value Po for time delay compensation
    if((jdcomp = 1) and  (idcomp<dCompStartPoint)) then
      dCompCheckInitationP()
    endif
  
    if(Current_Step<dcompStartPoint) then 
      dcomphs = 0.0
    else
      if(Current_Step<dcompFullPoint) then 
        dcomphs = (Current_Step-dcompStartPoint)*dcomphsmax / (dcompFullPoint-dcompStartPoint)
      else
        if(Current_Step<dcompDecreasingStep) then
          dcomphs = dcomphsmax
        else
          if(Current_Step<dcompStopStep) then
            dcomphs = dcomphsmax * (1.0 - (Current_Step-dcompDecreasingStep) / (dcompStopStep - dcompDecreasingStep))
          else 
            dcomphs = 0.0
          endif
        endif
      endif
    endif
  endif
  if((dcompCompensation=1) and (WithActuator=1)) then
    dCompSetPHI()
    dCompEstimatedU()
  else
    dcompdU = 0.0
  endif
ENDSUB

SUB Process_Dcomp_After() 
  
  dim temp1,temp2,fc,kfc as float
  
  'Backup History of data for dcomp
  if(dcompFilter = 0) then 
    dcomperr  = (uc[1] - dmeas2) - (1.0-dcomphs)*dcompdU
  else
    fc = dcomp_fc
    kfc = dtdata / ((1.0/(2*3.14*fc))+ dtdata)
    temp1 = dcomperrfiltered
    temp2 = (uc[1] - dmeas2) - (1.0-dcomphs)*dcompdU
    dcomperrfiltered = temp1 + kfc*(temp2 - temp1)
    dcomperr = dcomperrfiltered
  endif
  
  
  
  for i = 1 to dcompnewpos-1
    dcompHisdU[i]=  dcompHisdU[i+1]
    dcompHisU[i] =  dcompHisU[i+1]
    dcompHisV[i] =  dcompHisV[i+1]
    dcompHisE[i] =  dcompHisE[i+1]
  next i
  'dcompHisdU[(dcompndelay*dcompnu)+1]=  dcompdUfiltered
  dcompHisdU[(dcompndelay*dcompnu)+1]= dcompdU
  
  'Put new displacement and new error to the history ***
  dcompHisU[(dcompndelay*dcompnu)+1] =  uc[1]
  dcompHisV[(dcompndelay*dcompnu)+1] =  (uc[1] - dcompuijold)/dtdata
  dcompHisE[(dcompndelay*dcompnu)+1] =  dcomperr
  
  dcompuijold = uc[1]
  
  'Update model at each delay time period
  if(WithActuator=1) then
    if(jdcomp = dcompndelay) then
      if(((dcompCompensation = 1) and (Current_Step>dCompStartPoint)))then
        dCompUpdatemodel()
      endif
      jdcomp = 0
      idcomp = idcomp + 1
    endif
    jdcomp=jdcomp+1
  endif
  
ENDSUB
' Set all the elements of the vector to a specified value
SUB Set2Value(vec, Num_Elements, Value)
  dim k as long
  
  for k = 1 to Num_Elements
    vec[k] = Value
  next k
EndSub
' Copies the values of a vector to another vector
SUB scopy(vec_in, Num_Elements, vec_out)
  dim k as long
  
  for k = 1 to Num_Elements
    vec_out[k] = vec_in[k]
  next k
EndSub
