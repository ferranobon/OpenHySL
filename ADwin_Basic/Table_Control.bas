'<ADbasic Header, Headerversion 001.001>
' Process_Number                 = 2
' Initial_Processdelay           = 3000
' Eventsource                    = Timer
' Control_long_Delays_for_Stop   = No
' Priority                       = Low
' Priority_Low_Level             = 1
' Version                        = 1
' ADbasic_Version                = 5.0.6
' Optimize                       = Yes
' Optimize_Level                 = 1
' Info_Last_Save                 = BASILISC  Basilisc\Ferran
'<Header End>
#include ADwinPRO_ALL.inc         ' Include-file for Pro-system t11

' todo check the adlog related routines (variables etc)

#define timeoffset 14316.5576170  'time offset for correcting time calculation 
'when the sign of clock counter changes from max to min

#define Num_Step          4096    ' Number of steps in the sub-structure test
#define dtstep            0.5     ' Time increment of a step
#define Num_Substep       4       ' Number of sub-steps

' ADlog related
#define MaxSubStep    32768   ' number of of substeps in whole test
#define lenghtDataAll 3276801 'lenght of data for data record that will be read by VC software 
#define lenghtDataDAQ 3276801 'lengh of data array for accquisition

' Order of matrices and vectors
#define Order             1       ' Number of DOF of the sub-structure
#define Order_Input       2       ' Order of the input array = Order + 1
#define Order_Output      4       ' Order of the output array = 3*Order + 1
#define Order_Gain        1       ' Order of the gain matrix = Order*Order

#define Num_Event_Substep 500     ' Number of events per sub-step. Initially 500

' Limitations for emergency stop
#define dtabmax 							0.09
#define accmax 								50.0
#define dTMDmax               0.095
#define Fmax                  6000.0

' Different types of error flags. Used in combination with ERROR_FLAG variable.
#define ERROR_UNDER_DMIN          1       ' The measured value is below the minimum displacement value
#define ERROR_OVER_DMAX           2       ' The measured value is over the maximum displacement value

' Zero of the different accelerometers. Check at the begining of each test through adlog, so that it is set
' to zero.
#define acc1_zero          -2.0
#define acc2_zero          -0.219
#define acc3_zero          -0.42
' The filter has to be used when reading the pressure from the device since otherwise it shows some
' shocks. This was experienced by Van Thuan during some of the E-FAST tests.
#define compensateFilter   1.0225 

' Limit check
#define Avoid_Limit_Check       1     ' 0 = normal, 1 = avoiding limit checking

#define NEW_DATA_AVAILABLE      PAR_4   ' Indicates whether there is new data available (value of 1) or not (value of 0).
#define EMERGENCY_STOP          PAR_10  ' Variable to stop the process suddenly
#define TESTING                 PAR_12  ' The substepping process has started (value of 1) and it is currently running a test.
Rem                                       This is an inter-process communication variable.
#define ERROR_FLAG              PAR_15  ' Variable used to store the error during the check of measured values against limit values.
#define Current_Step            PAR_24  ' Counter for the number of steps
#define NumRecordedData					PAR_29  ' for mornitoring number of tested substeps Ask Van Thuan

' Control variables
#define dtdata                  FPAR_2  ' Data to be controlled. Ask Van Thuan

#define dctrl1                  FPAR_3  ' Variable containing the control displacement in direction 1
#define dctrl2                  FPAR_4  ' Variable containing the control displacement in direction 2
#define accctrl1 	              FPAR_5
#define accctrl2 	              FPAR_6
#define velctrl1 	              FPAR_7
#define velctrl2 	              FPAR_8

' Synchronisation flags
#define SYNC_PC                 FPAR_31 ' This flag is used to synchronise the uploaded data by the PC. It can have two values:
Rem                                       - SYNC_PC = 0. Do nothing or the substep is performed
Rem                                       - SYNC_PC = 1. New data has been uploaded from the PC.
#define SYNC_ADWIN              FPAR_32 ' This flag is used to synchronise ADwin and the PC. It can have three possible values:
Rem                                       - SYNC_ADWIN = 0.0. Do nothing
Rem                                       - SYNC_ADWIN = 1.0. Perform the sub-stepping process since new data has been uploaded
Rem                                       - SYNC_ADWIN = -1.0. The sub-step process is finished and the PC can collect the new
Rem                                         data (see Output_Data).

#define Press_pctrl							FPAR_48 ' Control pressure
#define Press_pmeas							FPAR_47 ' Measured pressure

#define F12x                    FPAR_68 ' Measured Force in the load cell. Used for inter-process communication with the sub-stepping routine
#define F12y                    FPAR_69 ' Measured Force in the load cell. Used for inter-process communication with the sub-stepping routine
#define F34x                    FPAR_70 ' Measured Force in the load cell. Used for inter-process communication with the sub-stepping routine
#define F34y                    FPAR_71 ' Measured Force in the load cell. Used for inter-process communication with the sub-stepping routine
#define dmeas1                  FPAR_72 ' Measured displacement in direction 1. Used for inter-process communication with the sub-stepping routine
#define dmeas2                  FPAR_73 ' Measured displacement in direction 2. Used for inter-process communication with the sub-stepping routine
#define ameas1                  FPAR_74 ' Measured acceleration in direction 1. Used for inter-process communication with the sub-stepping routine
#define ameas2                  FPAR_75 ' Measured acceleration in direction 2. Used for inter-process communication with the sub-stepping routine
#define vmeas1                  FPAR_76 ' Measured velocity in direction 1. Used for inter-process communication with the sub-stepping routine
#define vmeas2                  FPAR_77 ' Measured velocity in direction 2. Used for inter-process communication with the sub-stepping routine

#define Gain_Matrix             DATA_1  ' The gain matrix.
#define Input_Data              DATA_2  ' Array to store the u0 vector from the PC. When Input_Data[1] = 1, the PC has uploaded
Rem                                       a new vector starting at Input_Data[2]. The Input_Data[1] is switched back to 0.0 when
Rem                                       ADwin reads it and start the sub-stepping.
#define Output_Data             DATA_3  ' Array to store the result of the sub-step process.
Rem                                       - Output_Data[1] can have different values (see SYNC_ADWIN). When Output_Data[1] = -1.0
Rem                                         the PC will pick the available values and continue with the stepping process. Any
Rem                                         other valuewill make the PC to scan this value until it is -1.0 in order to proceed
Rem                                         with the stepping process.
Rem                                       - Output_Data[2] stores the uc vector resulting from the sub-stepping process
Rem                                       - Output_Data[2+Order] stores the coupling force of the substep 'Num_Substeps -1'
Rem                                       - Output_Data[2*(1 + Order)] stores the coupling force of the last sub-step

#define UVADATA1 	DATA_92
#define DataAll 	DATA_97
#define dataDAQ   DATA_199   ' data accquisition: substep, data data data ... substep, data, data, data, ...

'define of ADlog
#define NumberLONGsPerStep              51    ' number of longs per each step of record set for ADlog
#DEFINE ADlogData1 DATA_180
#DEFINE ADlogData2 DATA_181
#DEFINE WritePointer ADlogData2[1]
#DEFINE LoopCounter	ADlogData2[2]
#DEFINE BufferSize ADlogData2[3]
#DEFINE ValuesCount ADlogData2[4]
#DEFINE Flags ADlogData2[5]
#DEFINE BUFSIZE 4000000 ' whole-numbered multiple of values count
DIM ADlogData1[BUFSIZE] AS LONG AT DRAM_EXTERN
DIM ADlogData2[200] AS LONG AT DM_LOCAL
dim byte1,byte2 as long
'end of define ADLog

dim Gain_Matrix[Order_Gain] as float at DM_LOCAL    ' Allocate the necesssary memory for the Gain matrix. The values are
Rem                                                   suplied by the PC.
dim Input_Data[Order_Input] as float at DM_LOCAL    ' Allocate the necessary space
dim Output_Data[Order_Output] as float at DM_LOCAL  ' Allocate the necessary space

dim UVADATA1[MaxSubStep] as float
dim DataAll[lenghtDataAll] as float
dim DataDAQ[lenghtDataDAQ] as float

dim i as long                           ' A simple counter

dim dtevent as float                    ' Time increment of each event in ADwin
dim dtsub as float                      ' Time increment of the sub-step

dim Current_Substep as long             ' Counter for the sub-steps
dim Current_Event   as long             ' Counter for the number of events in a sub-step

dim ramp  as float                      ' Variable to store the operation Current_Substep/Num_Substep.
dim ramp0 as float                      ' Variable to store the operation 1.0 - Current_Substep/Num_Substep

dim u0[order] as float                  ' To store the data received from the PC at the previous step
dim u01[order] as float                 ' To store the data received from the PC at the beginning of the sub-step (current step)
dim uc[order] as float                  ' Vector to be sent to the PC at the end of each sub-stepping process.
dim fcprev[order] as float              ' Vector containing the coupling force measured at Num_Substep - 1
dim fc[order] as float                  ' Vector containing the coupling force measured at the last sub-step

dim dTMD as float                       ' Displacement of the TMD

dim IS_FIRST_STEP as short
dim DO_SUB_STEP   as short              ' Variable to identify whether the sub-stepping process should be performed (value of 1) or not
Rem                                       (value of 0)

dim accExt1,accExt2,accExt3 as float    ' Variables to store the measurement of external accelerometers

dim beginning_time,ending_time,time_between_substep as float
dim clock,clock_synPC,clock_synPC_old,clock_do_substep,clock_do_substep_old,clock_beginning_step as long
dim time_synPC,time_substep,time_step,time_do_substep,time_from_sub1_to_newPCdata as float

dim T,TcompSub,TcompSubMax as float
INIT:
  
  Init_Variables()
  ADlogInnit()
EVENT:
  
  dtdata = dtsub                      ' Set the time increment for new data/measurements to the sub-step time. This variable
  Rem                                   is also used in the high priority process.
  
  if(Input_Data[1] = 1.0) then        ' The PC has uploaded some data and it is ready to be used
    
    'Time checking
    clock_synPC_old = clock_synPC
    clock_synPC = Read_Timer()
    time_synPC = (clock_synPC - clock_synPC_old) * (10.0/3.0)/1000000.0
    if(time_synPC<0.0) then
      time_synPC = time_synPC + timeoffset
    endif
    if(Current_Step>1) then
      time_from_sub1_to_newPCdata = (clock_synPC - clock_beginning_step) * (10.0/3.0)/1000000.0
    else
      time_from_sub1_to_newPCdata = 0.0
    endif
    if(time_from_sub1_to_newPCdata<0.0) then
      time_from_sub1_to_newPCdata = time_from_sub1_to_newPCdata + timeoffset
    endif
    'End time checking
    
    if(IS_FIRST_STEP = 1) then        ' In the first step, ADwin should start inmediately and without waiting.
      Current_Substep = 1             ' Initialise the counter of sub-steps
      Current_Event = 1               ' Initialise the event counter
      Current_Step = 1
      
      SYNC_ADWIN = 1.0                ' ADwin is synchronised since this is the first step
      Output_Data[1] = SYNC_ADWIN     ' The data is not yet ready (ADwin - PC synchronisation)
      SYNC_PC = 0.0                   ' The values have been received and now the sub-step process will be performed
      Input_Data[1] = SYNC_PC         ' Inform the PC about this change in the process      
      
      IS_FIRST_STEP = 0               ' No more first step case.
      TESTING = 1                     ' The testing process is now oficially running
    else                              ' ADwin notes that the new data is ready since 
      SYNC_ADWIN = 1.0
      Output_Data[1] = SYNC_ADWIN
      SYNC_PC = 0.0
      Input_Data[1] = SYNC_PC
    endif
  endif
  
  
  if(SYNC_ADWIN = 1.0) then              ' Perform sub-stepping
    if(Current_Event = 1) then        ' The sub-stepping should only take place at the begining of the scanning loop in order
      Rem                               to make sure that the process runs at the required time. Otherwise the test results
      Rem                               would not be coherent since the sub-steps will not take place in the specified time
      Rem                               increments
      DO_SUB_STEP = 1                 ' Tell ADwin that the sub-step process should be started/continued in the next events
    endif
    
    if(DO_SUB_STEP = 1) then
      
      if((Current_Step = 1) and (Current_Substep=1)) then
        clock_do_substep = Read_Timer()
        time_do_substep = 0.0
      else
        clock_do_substep_old = clock_do_substep
        clock_do_substep = Read_Timer()
        time_do_substep = (clock_do_substep - clock_do_substep_old) * (10.0/3.0)/1000000.0
        if(time_do_substep<0.0) then
          time_do_substep = time_do_substep + timeoffset
        endif
      endif
      
      if(Current_Substep=1) then
        clock_beginning_step = clock_do_substep
      endif
      
      ramp = Current_Substep/Num_Substep         ' Update the ramp function values
      ramp0 = 1.0 - ramp
      
      if(Current_Substep = 1) then
        for i = 1 to Order
          u01[i] = Input_Data[i+1]               ' Copy the new data coming from the PC at the begining of the sub-stepping process into u01
        next i
      endif
          
      uc[1] = ramp0*u0[1] + ramp*u01[1]
        
      Write_Displacement_Signal()
        
      ' The forces are measured in the high priority process
      fc[1] = F12y + F34y
      fc[1] = 0.0  
      ' Backup the coupling force vector at the Num_Substep -1 substep.
      if ((Current_Substep = (Num_Substep -1)) or (Num_Substep = 1)) then
        scopy(fc, Order, fcprev)
      endif
        
      ' Check the values for safety reasons
      Check_Limitation() 
        
      if(Current_Substep < Num_Substep) then
        'Time checking
        clock = Read_Timer()
        time_substep = (clock - clock_synPC)*(10.0/3.0)/1000000.0
        if(time_substep<0.0) then
          time_substep = time_substep + timeoffset
        endif
        'ADlog()
      endif      
      
      if(Current_Substep = Num_Substep) then
        ' Copy the data to Output_Data so that the PC can read it
        for i = 1 to Order
          Output_Data[1+i] = uc[i]
          Output_Data[1+Order+i] = fcprev[i]
          Output_Data[1+2*Order + i] = fc[i]
        next i
        SYNC_ADWIN = -1.0
        Output_Data[1] = SYNC_ADWIN       ' Say the PC that the Data is ready. Output_Data[1] = -1.0
        
        'Time checking
        clock = Read_Timer()
        time_substep = (clock - clock_synPC)*(10.0/3.0)/1000000.0
        if(time_substep<0.0) then
          time_substep = time_substep + timeoffset
        endif
          
        'ADlog()
        Current_Step = Current_Step + 1
        Current_Substep = 0               ' Reset the sub-step counter
        scopy(u01,Order,u0)               ' Copy the vector received from the PC during this step to u0 in order to compute
        Rem                                 uc in the next step   
      endif
      
      Current_Substep = Current_Substep + 1 ' Increase the substep counter
      
      if(TcompSubmax < TcompSub) then
        TcompSubmax = TcompSub
      endif
      
      DO_SUB_STEP = 0                       ' Do not perform any sub-step until the next cycle of events to assure that the
      Rem                                     specified amount of time has passed between to consecutive sub-steps
    endif    
  endif
  
  if (Current_Event = Num_Event_Substep) then
    Current_Event = 0                       ' Reset the counter in case it has reached its maximum value
  endif
  
  Current_Event = Current_Event + 1         ' Increase the Event counter 
  
  if((SYNC_ADWIN = 0.0) and (Current_Event = 1)) then
    'ADlogPre()
    if((Avoid_Limit_Check=0)) then
      Check_Limitation()
    endif
    clock_synPC = Read_Timer()
  endif
  
FINISH:
  
  dtdata = 5.0    ' Reset the new data/measurement time increment to a larger value in seconds
  Rem               to move the hydraulic cylinders slowlier.
  
  Rem -----------------------------------------------------------------------------------------------------------------------
  Rem ------------------------------------------------- SUB-ROUTINES SECTION ------------------------------------------------
  Rem -----------------------------------------------------------------------------------------------------------------------
  
SUB Init_Variables()
  dtsub = dtstep/Num_Substep        ' Defining the sub-step time increment
  dtevent = dtsub/Num_Event_Substep ' Defining the time increment of each event
  PROCESSDELAY = 300000000*dtevent  ' Define the process delay
  
  Current_Substep = 1               ' Set the sub-step counter to its starting value
  Current_Event = 1                 ' Set the event counter to its starting value
  Current_Step = 1
  
  ERROR_FLAG = 0              ' No error
  TESTING = 0                 ' We are not yet running a sub-structure test
  IS_FIRST_STEP = 1           ' The first step will take place once the PC is synchronised
  
  ' Initialise the vectors to 0.0
  Set2Value(Gain_Matrix, Order_Gain, 0.0)
  Set2Value(Input_Data, Order_Input, 0.0)
  Set2Value(Output_Data, Order_Output, 0.0)
  Set2Value(u0, Order, 0.0)
  Set2Value(u01, Order, 0.0)
  Set2Value(uc, Order, 0.0)
  Set2Value(fcprev, Order, 0.0)
  Set2Value(fc, Order, 0.0)
  
  T = 0.0
  TcompSub = 0.0
  TcompSubMax = 0.0
  NumRecordedData = 0
  
EndSub
SUB ADlogInnit()
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
  
  dctrl1 = 0          ' Direction 1 is never used in these tests 
  dctrl2 = uc[1]      ' Set the control displacement to the value of uc unless it exceeds the table displacement
  
  if(dctrl2 > dtabmax) then
    dctrl2 = dtabmax
  endif
  if(dctrl2 < -dtabmax) then
    dctrl2 = -dtabmax
  endif
  
  NEW_DATA_AVAILABLE = 1        ' Tell the other process that new data is available
EndSub
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
    
    if((dTMD > dTMDmax) or (dTMD< -dTMDmax)) then
      EMERGENCY_STOP = 1
      if(dTMD> dtabmax) then 
        ERROR_FLAG = ERROR_OVER_DMAX
      else
        ERROR_FLAG = ERROR_UNDER_DMIN
      endif
    endif
  endif
EndSub
SUB	ADlogPre()
  
  dim data as LONG
  dim temp1 as float
  
  data = 0                                     'trigger
  ADlogData1[WritePointer] = data
  
  ' Forces are measured in the high priority process
  
  GetAndStoreData_Pre()
  
  WritePointer = WritePointer + NumberLONGsPerStep
  
  ' check writepointer
  IF (WritePointer > BUFSIZE) THEN
    WritePointer = 1
    INC LoopCounter
  ENDIF
EndSub
SUB	ADlog()
  dim data as LONG
  dim temp1 as float
  
  data = 2147483647                                 'trigger
  ADlogData1[WritePointer] = data
  
  GetAndStoreData_withDAQ()
  
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
  dim dpisint as long
  dim i,j,lenght as long
  dim numchannelDAQ as long ' number of channels for DAQ
  
  temp =  - (((ADCF(2, 7)-32768)/32768.0)*10.0 / 0.104) * 9.8065 - acc1_zero      ' Endevco 61-100 Acc tabe                 
  accExt1 = temp

  temp = - (((ADCF(3, 1)-32768)/32768.0)*10.0 / 0.936) * 9.8065 - acc2_zero       ' PCB accelerometer - TMD
  accExt2 = temp
  
  temp = - (((ADCF(2, 8)-32768)/32768.0)*10.0 / 0.507 ) * 9.8065 - acc3_zero       ' Kistler Acc Frame                     
  accExt3 = temp
  
  dpisint = (ADCF(1,5) -  32768)
  dtab = ((ADCF(2, 5)-32768)/32768.0) * 0.1  
  dTMD = - ((ADCF(3, 2)-32768)/32768.0) * 0.1 
  
        
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

  temp = ((ADCF(1, 3)-32768)/32768.0) * 250000.0        
  data = (temp / 250000.0) * 2147483647               'F cylinder 1
  ADlogData1[WritePointer+11] = data
  
  temp = ((ADCF(1, 4)-32768)/32768.0) * 250000.0        
  data = (temp / 250000.0) * 2147483647               'F cylinder 2
  ADlogData1[WritePointer+12] = data
  
  
  data = (dtab / 0.1) * 2147483647               'D table (mm)
  ADlogData1[WritePointer+13] = data
  
 
  data = (dTMD / 0.1) * 2147483647               'D tmd (mm)
  ADlogData1[WritePointer+14] = data
  
  data = (F12x / 20000.0) * 2147483647               'F12x (N)
  ADlogData1[WritePointer+15] = data
  
  data = (F34x / 20000.0) * 2147483647               'F34x (N)
  ADlogData1[WritePointer+16] = data
  
  data = (F12y / 20000.0) * 2147483647               'F12y (N)
  ADlogData1[WritePointer+17] = data
  
  data = (F34y / 20000.0) * 2147483647               'F34y (N)
  ADlogData1[WritePointer+18] = data
  
  data = Press_pctrl/2.0 * 2147483647  'PressureCtrl
  ADlogData1[WritePointer+19] = data
  
  Press_pmeas = (((ADCF(3,8) - 32768) * compensateFilter) / 32768.0) * 2.0
  data = (Press_pmeas/2.0 ) * 2147483647 
  ADlogData1[WritePointer+20] = data
  
  data =  (accExt3 / (9.8065 * 10.0)) * 2147483647     'acc3  (servo Acc)
  ADlogData1[WritePointer+21] = data
  
  data = (dpisint / 32768.0) * 2147483647  ' dcylmeas
  ADlogData1[WritePointer+22] = data
  
  data =  ( UVADATA1[Current_Step] / (0.1)) * 2147483647     'dbase 
  ADlogData1[WritePointer+23] = data
  
  data =  (uc[1] / 0.1) * 2147483647     ' disi1
  ADlogData1[WritePointer+24] = data
  
  data =  (0.0 / 1.0) * 2147483647     ' veli1
  ADlogData1[WritePointer+25] = data
  
  data =  (0.0 /(9.8065 * 10.0)) * 2147483647     ' acci1
  ADlogData1[WritePointer+26] = data
  
  data =  (fc[1] /20000) * 2147483647     ' fc
  ADlogData1[WritePointer+27] = data
  
  'data =  (fu[1] /20000) * 2147483647     ' fui1
  data = 0
  ADlogData1[WritePointer+28] = data
  
  'data =  (Ff /20000) * 2147483647     ' Ff
  data = 0
  ADlogData1[WritePointer+29] = data
  
  data =  (0 * 0 / 0.1) * 2147483647     ' diss
  ADlogData1[WritePointer+30] = data
  
  data =  (0 * 0 / 1.0) * 2147483647     ' vels
  ADlogData1[WritePointer+31] = data
  
  data =  (0 * 0 /(9.8065 * 10.0)) * 2147483647     ' accs
  ADlogData1[WritePointer+32] = data
  
  data =  (0.0 /1.0) * 2147483647     ' 
  ADlogData1[WritePointer+33] = data
  
  data =  (0.0 /1.0) * 2147483647     ' 
  ADlogData1[WritePointer+34] = data
  
  'data =  (FucompdF /20000.0) * 2147483647     ' 
  data =  0
  ADlogData1[WritePointer+35] = data
  
  'data =  (FucompTHETA[1] /1000.0) * 2147483647     ' 
  data = 0
  ADlogData1[WritePointer+36] = data
  
  'data =  (FucompTHETA[2] /1000.0) * 2147483647     ' 
  data = 0
  ADlogData1[WritePointer+37] = data
  
  'data =  (FucompTHETA[3] /1000.0) * 2147483647     ' 
  data = 0
  ADlogData1[WritePointer+38] = data
  
  'data =  (FucompTHETA[4] /1000.0) * 2147483647     '
  data = 0 
  ADlogData1[WritePointer+39] = data
  
  'data =  (FucompTHETA[5] /1000.0) * 2147483647     ' 
  data = 0
  ADlogData1[WritePointer+40] = data
  
  'data =  (FucompTHETA[6] /1000.0) * 2147483647     ' 
  data = 0
  ADlogData1[WritePointer+41] = data
  
  'data =  (FucompTHETA[7] /1000.0) * 2147483647     ' 
  data = 0
  ADlogData1[WritePointer+42] = data
  
  'data =  (dcompdU /1.0) * 2147483647     ' 
  data = 0
  ADlogData1[WritePointer+43] = data
  
  'data =  (dcompTHETA[1] /1000.0) * 2147483647     ' 
  data = 0
  ADlogData1[WritePointer+44] = data
  
  'data =  (dcompTHETA[2] /1000.0) * 2147483647     ' 
  data = 0
  ADlogData1[WritePointer+45] = data
  
  ' data =  (dcompTHETA[3] /1000.0) * 2147483647     ' 
  data = 0
  ADlogData1[WritePointer+46] = data
  
  'data =  (dcompTHETA[4] /1000.0) * 2147483647     ' 
  data = 0
  ADlogData1[WritePointer+47] = data
  
  'data =  (dcompTHETA[5] /1000.0) * 2147483647     ' 
  data = 0
  ADlogData1[WritePointer+48] = data
  
  'data =  (dcompTHETA[6] /1000.0) * 2147483647     ' 
  data = 0
  ADlogData1[WritePointer+49] = data
  
  'data =  (dcompTHETA[7] /1000.0) * 2147483647     ' 
  'data =  (q/0.1) * 2147483647     ' 
  data = 0
  ADlogData1[WritePointer+50] = data
ENDSUB
SUB GetAndStoreData_withDAQ()
  dim data,dtab,temp as float
  dim dpisint as long
  dim i,j,lenght as long
  dim numchannelDAQ as long ' number of channels for DAQ
  
  temp =  - (((ADCF(2, 7)-32768)/32768.0)*10.0 / 0.104) * 9.8065 - acc1_zero      ' Endevco 61-100 Acc tabe                 
  accExt1 = temp

  temp = - (((ADCF(3, 1)-32768)/32768.0)*10.0 / 0.936) * 9.8065 - acc2_zero       ' PCB accelerometer - TMD
  accExt2 = temp
  
  temp = - (((ADCF(2, 8)-32768)/32768.0)*10.0 / 0.507 ) * 9.8065 - acc3_zero       ' Kistler Acc Frame                     
  accExt3 = temp
  
  dpisint = (ADCF(1,5) -  32768)
  dtab = ((ADCF(2, 5)-32768)/32768.0) * 0.1  
  dTMD = - ((ADCF(3, 2)-32768)/32768.0) * 0.1 
  
        
  data =  (dctrl1 /0.2)* 2147483647                'dctrl1
  ADlogData1[WritePointer+1] = data
   
  
  data =  (dctrl2 /0.2)* 2147483647                'dctrl2
  ADlogData1[WritePointer+2] = data
  
  data =  (dmeas1 /0.2)* 2147483647                'dmeas1
  ADlogData1[WritePointer+3] = data
  
  data =  (dmeas2 /0.2)* 2147483647                'dmeas2
  ADlogData1[WritePointer+4] = data
  
  data =  (accctrl1 / (9.8065 * 10.0)) * 2147483647     'accctrl1 really used ?
  ADlogData1[WritePointer+5] = data
  
  data =  (accctrl2 / (9.8065 * 10.0)) * 2147483647     'accctrl2 really used? Ask van thuan
  ADlogData1[WritePointer+6] = data
  
  data =  (ameas1 / (9.8065 * 10.0)) * 2147483647     'ameas1 by computing dis
  ADlogData1[WritePointer+7] = data
  
  data =  (ameas2 / (9.8065 * 10.0)) * 2147483647     'ameas2 by computing dis
  ADlogData1[WritePointer+8] = data
  
  data =  (accExt1 / (9.8065 * 10.0)) * 2147483647     'acc1
  ADlogData1[WritePointer+9] = data
  
  data =  (accExt2 / (9.8065 * 10.0)) * 2147483647     'acc2
  ADlogData1[WritePointer+10] = data

  temp = ((ADCF(1, 3)-32768)/32768.0) * 250000.0        
  data = (temp / 250000.0) * 2147483647               'F cylinder 1
  ADlogData1[WritePointer+11] = data
  
  temp = ((ADCF(1, 4)-32768)/32768.0) * 250000.0        
  data = (temp / 250000.0) * 2147483647               'F cylinder 2
  ADlogData1[WritePointer+12] = data
  
  
  data = (dtab / 0.1) * 2147483647               'D table (mm)
  ADlogData1[WritePointer+13] = data
  
 
  data = (dTMD / 0.1) * 2147483647               'D tmd (mm)
  ADlogData1[WritePointer+14] = data
  
  data = (F12x / 20000.0) * 2147483647               'F12x (N)
  ADlogData1[WritePointer+15] = data
  
  data = (F34x / 20000.0) * 2147483647               'F34x (N)
  ADlogData1[WritePointer+16] = data
  
  data = (F12y / 20000.0) * 2147483647               'F12y (N)
  ADlogData1[WritePointer+17] = data
  
  data = (F34y / 20000.0) * 2147483647               'F34y (N)
  ADlogData1[WritePointer+18] = data
  
  data = Press_pctrl/2.0 * 2147483647  'PressureCtrl
  ADlogData1[WritePointer+19] = data
  
  Press_pmeas = (((ADCF(3,8) - 32768) * compensateFilter) / 32768.0) * 2.0
  data = (Press_pmeas/2.0 ) * 2147483647 
  ADlogData1[WritePointer+20] = data
  
  data =  (accExt3 / (9.8065 * 10.0)) * 2147483647     'acc3  (servo Acc)
  ADlogData1[WritePointer+21] = data
  
  data = (dpisint / 32768.0) * 2147483647  ' dcylmeas
  ADlogData1[WritePointer+22] = data
  
  data =  ( UVADATA1[Current_Step] / (0.1)) * 2147483647     'dbase 
  ADlogData1[WritePointer+23] = data
  
  data =  (uc[1] / 0.1) * 2147483647     ' disi1
  ADlogData1[WritePointer+24] = data
  
  data =  (0.0 / 1.0) * 2147483647     ' veli1
  ADlogData1[WritePointer+25] = data
  
  data =  (0.0 /(9.8065 * 10.0)) * 2147483647     ' acci1
  ADlogData1[WritePointer+26] = data
  
  data =  (fc[1] /20000) * 2147483647     ' fc
  ADlogData1[WritePointer+27] = data
  
  'data =  (fu[1] /20000) * 2147483647     ' fui1
  data = 0
  ADlogData1[WritePointer+28] = data
  
  'data =  (Ff /20000) * 2147483647     ' Ff
  data = 0
  ADlogData1[WritePointer+29] = data
  
  data =  (0 * 0 / 0.1) * 2147483647     ' diss
  ADlogData1[WritePointer+30] = data
  
  data =  (0 * 0 / 1.0) * 2147483647     ' vels
  ADlogData1[WritePointer+31] = data
  
  data =  (0 * 0 /(9.8065 * 10.0)) * 2147483647     ' accs
  ADlogData1[WritePointer+32] = data
  
  data =  (0.0 /1.0) * 2147483647     ' 
  ADlogData1[WritePointer+33] = data
  
  data =  (0.0 /1.0) * 2147483647     ' 
  ADlogData1[WritePointer+34] = data
  
  'data =  (FucompdF /20000.0) * 2147483647     ' 
  data =  0
  ADlogData1[WritePointer+35] = data
  
  'data =  (FucompTHETA[1] /1000.0) * 2147483647     ' 
  data = 0
  ADlogData1[WritePointer+36] = data
  
  'data =  (FucompTHETA[2] /1000.0) * 2147483647     ' 
  data = 0
  ADlogData1[WritePointer+37] = data
  
  'data =  (FucompTHETA[3] /1000.0) * 2147483647     ' 
  data = 0
  ADlogData1[WritePointer+38] = data
  
  'data =  (FucompTHETA[4] /1000.0) * 2147483647     '
  data = 0 
  ADlogData1[WritePointer+39] = data
  
  'data =  (FucompTHETA[5] /1000.0) * 2147483647     ' 
  data = 0
  ADlogData1[WritePointer+40] = data
  
  'data =  (FucompTHETA[6] /1000.0) * 2147483647     ' 
  data = 0
  ADlogData1[WritePointer+41] = data
  
  'data =  (FucompTHETA[7] /1000.0) * 2147483647     ' 
  data = 0
  ADlogData1[WritePointer+42] = data
  
  'data =  (dcompdU /1.0) * 2147483647     ' 
  data = 0
  ADlogData1[WritePointer+43] = data
  
  'data =  (dcompTHETA[1] /1000.0) * 2147483647     ' 
  data = 0
  ADlogData1[WritePointer+44] = data
  
  'data =  (dcompTHETA[2] /1000.0) * 2147483647     ' 
  data = 0
  ADlogData1[WritePointer+45] = data
  
  ' data =  (dcompTHETA[3] /1000.0) * 2147483647     ' 
  data = 0
  ADlogData1[WritePointer+46] = data
  
  'data =  (dcompTHETA[4] /1000.0) * 2147483647     ' 
  data = 0
  ADlogData1[WritePointer+47] = data
  
  'data =  (dcompTHETA[5] /1000.0) * 2147483647     ' 
  data = 0
  ADlogData1[WritePointer+48] = data
  
  'data =  (dcompTHETA[6] /1000.0) * 2147483647     ' 
  data = 0
  ADlogData1[WritePointer+49] = data
  
  'data =  (dcompTHETA[7] /1000.0) * 2147483647     ' 
  'data =  (q/0.1) * 2147483647     ' 
  data = 0
  ADlogData1[WritePointer+50] = data
    
  'these indexes for recording dataDAQ[]
  numchannelDAQ = 16
  j = (((Current_Step-1)*Num_Substep + Current_Substep)  - 1) * numchannelDAQ
  DataDAQ[j+1] = (Current_Step-1)*Num_Substep + Current_Substep
  DataDAQ[j+2] = uc[1]
  DataDAQ[j+3] = dctrl2
  DataDAQ[j+4] = dmeas2
  DataDAQ[j+5] = (dpisint/32768.0)*0.1
  DataDAQ[j+6] = accExt1 ' atab
  DataDAQ[j+7] = accExt2 ' atmd
  DataDAQ[j+8] = accExt3 ' aframe
  DataDAQ[j+9] = dtab ' dtab
  DataDAQ[j+10] = dtmd ' atmd
  DataDAQ[j+11] = F12x ' Fc
  DataDAQ[j+12] = F34x ' Fc
  DataDAQ[j+13] = F12y ' Fc
  DataDAQ[j+14] = F34y ' Fc
  DataDAQ[j+15] = Press_pctrl '
  DataDAQ[j+16] = Press_pmeas '
    
  lenght = Num_Step*Num_Substep
  i = (Current_Step-1)*Num_Substep + Current_Substep  
  j = 0
  DataAll[j*lenght + i] = i
  j = j+1
  DataAll[j*lenght + i] = time_do_substep
  j = j+1
  DataAll[j*lenght + i] = time_substep
  j = j+1
  DataAll[j*lenght + i] = time_synPC
  j = j+1
  DataAll[j*lenght + i] =  time_from_sub1_to_newPCdata
  
  j = j+1
  DataAll[j*lenght + i] =  u01[1]
  j = j+1
  DataAll[j*lenght + i] =  uc[1]
  j = j+1
  DataAll[j*lenght + i] =  fc[1]
  j = j+1
  DataAll[j*lenght + i] =  dctrl2
  j = j+1
  DataAll[j*lenght + i] =  dmeas2
  
  j = j+1
  DataAll[j*lenght + i] =  (dpisint/32768.0)*0.1
  j = j+1
  DataAll[j*lenght + i] =  accExt1 ' atab
  j = j+1
  DataAll[j*lenght + i] =  accExt2 ' atmd
  j = j+1
  DataAll[j*lenght + i] =  accExt3 ' aframe
  j = j+1
  DataAll[j*lenght + i] =  dtab ' dtab
  
  j = j+1
  DataAll[j*lenght + i] =  dtmd ' atmd
  j = j+1
  DataAll[j*lenght + i] =  F12x ' Fc
  j = j+1
  DataAll[j*lenght + i] =  F34x ' Fc
  j = j+1
  DataAll[j*lenght + i] =  F12y ' Fc
  j = j+1
  DataAll[j*lenght + i] =  F34y ' Fc
  
  j = j+1
  DataAll[j*lenght + i] =  Press_pctrl '
  j = j+1
  DataAll[j*lenght + i] =  Press_pmeas '
  
  'Total number of channels for data recoreding: 22
  
endsub
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
