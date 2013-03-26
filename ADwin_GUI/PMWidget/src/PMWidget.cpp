#include <QFileDialog>
#include <QString>
#include <QTimer>
#include <QDateTime>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cassert>

#include "PMWidget.hpp"
#include "Support_Routines.hpp"

PMWidget::PMWidget( QWidget *parent )
     : QWidget( parent )
{
     setupUi( this );

     ADwin_is_Booted = 0;

     CntrlP_is_Loaded = 0;
     CntrlP_is_Started = 0;

     TestP_is_Loaded = 0;
     TestP_is_Started = 0;

     Configure_Cntrl_SpinBoxes( );
     Configure_PID_SpinBoxes( );

     /* Boot file management */
     connect( BFiletoolButton, SIGNAL( clicked() ), this, SLOT( On_BFiletoolButton_Clicked() ) );
     /* Boot ADwin*/
     connect( BootButton, SIGNAL(clicked() ), this, SLOT(On_BootButton_Clicked() ));

     /* Displacement control process file */
     connect( PFiletoolButton, SIGNAL( clicked() ), this, SLOT( On_PFiletoolButton_Clicked() ) );
     /* Load Displacement control process into ADwin */
     connect( LoadCtrlPButton, SIGNAL( clicked() ), this, SLOT( On_LoadCtrlPButton_Clicked() ) );
     /* Start Displacement control process in ADwin */
     connect( StartCntrlPButton, SIGNAL( clicked() ), this, SLOT( On_StartCtrlPButton_Clicked() ) );
     /* Stop displacement control process */
     connect( StopCntrlPButton, SIGNAL( clicked() ), this, SLOT( On_StopCntrlPButton_Clicked() ) );

     /* Set ADwin to be ready */
     connect( ADwinReadyButton, SIGNAL( clicked() ), this, SLOT( On_ADwinReadyButton_Clicked() ) );
     /* Emergency button */
     connect( EmergencyButton, SIGNAL( clicked() ), this, SLOT( On_EmergencyButton_Clicked() ) );

     /* Low pressure button */
     connect( LowPressureButton, SIGNAL( clicked() ), this, SLOT( On_LowPressureButton_Clicked() ) );
     /* High pressure button */
     connect( HighPressureButton, SIGNAL( clicked() ), this, SLOT( On_HighPressureButton_Clicked() ) );
     /* Stop pressure button */
     connect( StopPressureButton, SIGNAL( clicked() ), this, SLOT( On_StopPressureButton_Clicked() ) );

     /* Substructure process file */
     connect( TFiletoolButton, SIGNAL( clicked() ), this, SLOT( On_TFiletoolButton_Clicked() ) );
     /* Load Test process into ADwin */
     connect( LoadTestPButton, SIGNAL( clicked() ), this, SLOT( On_LoadTestPButton_Clicked() ) );
     /* Start the test process in ADwin */
     connect( StartTestButton, SIGNAL( clicked() ), this, SLOT( On_StartTestButton_Clicked() ) );
     /* Stop the test process in ADwin */
     connect( StopTestButton, SIGNAL( clicked() ), this, SLOT( On_StopTestButton_Clicked() ) );
     /* Read the data */
     connect( ReadDataButton, SIGNAL( clicked() ), this, SLOT( On_ReadDataButton_Clicked() ) );

     /* Input File */
     connect( IFiletoolButton, SIGNAL( clicked() ), this, SLOT( On_IFiletoolButton_Clicked() ) );
     /* Output Data File */
     connect( OFiletoolButton, SIGNAL( clicked() ), this, SLOT( On_OFiletoolButton_Clicked() ) );

     /* Connect Spin boxes in Displacement control */
     connect( DCtrl1_DSpinBox, SIGNAL( valueChanged( double ) ), this, SLOT( On_DCtrl1_DSpinBox_valueChanged( double ) ) );
     connect( DCtrl2_DSpinBox, SIGNAL( valueChanged( double ) ), this, SLOT( On_DCtrl2_DSpinBox_valueChanged( double ) ) );

     /* Connect the spin boxes for PID control */
     connect( P_DSpinBox, SIGNAL( valueChanged( double ) ), this, SLOT( On_P_DSpinBox_valueChanged( double ) ) );
     connect( I_DSpinBox, SIGNAL( valueChanged( double ) ), this, SLOT( On_I_DSpinBox_valueChanged( double ) ) );
     connect( D_DSpinBox, SIGNAL( valueChanged( double ) ), this, SLOT( On_D_DSpinBox_valueChanged( double ) ) );

     /* Connect Limits Check Box */
     connect( EnabLimCheckBox, SIGNAL( stateChanged( int ) ), this, SLOT( On_SetLimCheckBox_Checked( int ) ) );

     /* Reset to default values */
     connect( DefValPButton, SIGNAL( clicked() ), this, SLOT( On_DefValPButton_Clicked() ) );

}

void PMWidget::Configure_Cntrl_SpinBoxes( )
{
     DCtrl1_DSpinBox->setDecimals( 5 );
     DCtrl1_DSpinBox->setSingleStep( 0.01 );
     DCtrl1_DSpinBox->setRange( DCntrlMin1, DCntrlMax1 );

     DCtrl2_DSpinBox->setDecimals( 5 );
     DCtrl2_DSpinBox->setSingleStep( 0.01 );
     DCtrl2_DSpinBox->setRange( DCntrlMin2, DCntrlMax2 );
}

void PMWidget::Configure_PID_SpinBoxes( )
{
     P_DSpinBox->setDecimals( 1 );
     P_DSpinBox->setSingleStep( 1.0 );
     P_DSpinBox->setRange( P_Min, P_Max );

     I_DSpinBox->setDecimals( 1 );
     I_DSpinBox->setSingleStep( 1.0 );
     I_DSpinBox->setRange( I_Min, I_Max );

     D_DSpinBox->setDecimals( 1 );
     D_DSpinBox->setSingleStep( 1.0 );
     D_DSpinBox->setRange( 0.0, 10000.0 );
}

void PMWidget::On_BFiletoolButton_Clicked( )
{
    
     BootFile_LineEdit->setText( QFileDialog::getOpenFileName(this, tr("Open ADwin Boot File"), "/opt/adwin/share/btl", tr("ADwin boot file (*btl)")) );

}

void PMWidget::On_BootButton_Clicked( )
{
     int Device_Number;
     std::string File_Path;
     int32_t Error;
     QString Error_Message;

     /* When ADwin is Booted all the buttons should be disabled */
     Disable_All_Buttons( );

     /* Setup the ADwin connection */
     Device_Number = DeviceNum_LineEdit->text().toInt();
     File_Path = BootFile_LineEdit->text().toStdString();
     Error = Boot_ADwin( Device_Number, File_Path.c_str());

     /* Check for errors */
     if ( Error != 0 ){
	  ADwin_is_Booted = 0;
	  Error_Message = Get_Error_Text( Error );
	  Print_Error_Message( (int) Error, Error_Message );
	  Print_To_Log( (int) Error, Error_Message );
     } else {
	  /* ADwin is correctly booted */
	  ADwin_is_Booted = 1;
	  CntrlP_is_Loaded = 0;
	  CntrlP_is_Started = 0;
	  
	  TestP_is_Loaded = 0;
	  TestP_is_Started = 0;

	  /* Enable the Load Control and Test Process buttons */
	  LoadCtrlPButton->setEnabled( ADwin_is_Booted );
	  LoadTestPButton->setEnabled( ADwin_is_Booted );

	  /* Enable Displacement control values and limitation buttons */
	  DefValPButton->setEnabled( ADwin_is_Booted );
	  SetLimPButton->setEnabled( ADwin_is_Booted );

	  /* Enable Check boxes */
	  EnabLimCheckBox->setEnabled( ADwin_is_Booted );

	  /* Print the message to the log */
	  Error_Message = "ADwin successfully booted.";
	  Print_To_Log( (int) Error, Error_Message );
     }
	  
}

void PMWidget::On_PFiletoolButton_Clicked( )
{
    
     PFile_LineEdit->setText( QFileDialog::getOpenFileName(this, tr("Open ADwin Binary File"), "/home/ferran/workspace/ADwinGUI/PMWidget", tr("ADwin binary file (*TB1)")) );

}

void PMWidget::On_LoadCtrlPButton_Clicked( )
{

     std::string File_Path;
     int32_t Error;
     QString Error_Message;

     File_Path = PFile_LineEdit->text().toStdString();
     Load_Process_ADwin( File_Path.c_str(), &CntrlP_is_Loaded, &Error );

     if ( Error > 0 ){
	  Error_Message = Get_Error_Text( Error );
	  Print_Error_Message( (int) Error, Error_Message );
	  Print_To_Log( (int) Error, Error_Message );
     } else if ( Error == -1 ){
	  Error_Message = "The Control process is already loaded. Boot ADwin or clear the process in order to reload it";
	  Print_Error_Message( (int) Error, Error_Message );
	  Print_To_Log( (int) Error, Error_Message );
     } else if ( Error == 0 ){
	  CntrlP_is_Loaded = 1;
	  /* Set the Start Control process button to enabled */
	  StartCntrlPButton->setEnabled( CntrlP_is_Loaded );
	  /* Print the message to the log */
	  Error_Message = "Displacement control process successfully loaded.";
	  Print_To_Log( (int) Error, Error_Message );
     } else {
	  assert( Error < -1 );   /* The error can never have a value less than -1 */
     }
}

void PMWidget::On_StartCtrlPButton_Clicked( )
{
     int Success_Cntr_Val, Success_Limit, Success_PID;
     int32_t Error;
     QString Error_Message;
     QTimer *timer = new QTimer(this);

     /* Set initial values */
     Set_Initial_Control_Values( &Success_Cntr_Val );
     Set_Initial_Limits( &Success_Limit );
     Set_Initial_PID_Values( &Success_PID );

     /* Start the process */
     if ( Success_Cntr_Val && Success_Limit && Success_PID ){
	  Start_Process_ADwin( (int32_t) 1, &CntrlP_is_Started, &Error );

	  if ( Error > 0 ){
	       Error_Message = Get_Error_Text( Error );
	       Print_Error_Message( (int) Error, Error_Message );
	       Print_To_Log( (int) Error, Error_Message );
	  } else if ( Error == -1 ){
	       Error_Message = "The Control process is already started. Boot ADwin or stop the process in order to reload it";
	       Print_Error_Message( (int) Error, Error_Message );
	       Print_To_Log( (int) Error, Error_Message );
	  } else if ( Error == 0 ){
	       CntrlP_is_Started = 1;
	       /* Set the Start Control process button to enabled */
	       ADwinReadyButton->setEnabled( CntrlP_is_Started );
	       StopCntrlPButton->setEnabled( CntrlP_is_Started );

	       DCtrl1_DSpinBox->setEnabled( CntrlP_is_Started );
	       DCtrl2_DSpinBox->setEnabled( CntrlP_is_Started );
	       P_DSpinBox->setEnabled( CntrlP_is_Started );
	       I_DSpinBox->setEnabled( CntrlP_is_Started );
	       D_DSpinBox->setEnabled( CntrlP_is_Started );

	       connect( timer, SIGNAL(timeout() ), this, SLOT( Read_Measurement_Values() ) );
	       timer->start( 200 );   /* Perform this operation every 200 ms */

	       /* Print the message to the log */
	       Error_Message = "Displacement control process successfully started.";
	       Print_To_Log( (int) Error, Error_Message );
	  } else {
	       assert( Error < -1 );   /* The error can never have a value less than -1 */
	  }
     } else {
	  Error = -1;
	  Error_Message = "The operation to start the control process could not be completed";
	  Print_Error_Message( (int) Error, Error_Message );
	  Print_To_Log( (int) Error, Error_Message );
     }

}

void PMWidget::Set_Initial_Control_Values( int *const Is_Successful )
{
     float temp;
     int32_t Error;
     QString Error_Message;

     /* Assuming everything is OK */
     *Is_Successful = 1;

     Get_FPar_ADwin( DMEAS1, &temp, &Error );
     if ( Error > 0 ){
	  Error_Message = Get_Error_Text( Error );
	  Print_Error_Message( (int) Error, Error_Message );
	  *Is_Successful = 0;
	  Print_To_Log( (int) Error, Error_Message );
     } else {
	  DCtrl1_DSpinBox->setValue( (double) temp );
	  DMeas1_lineEdit->setText( QString::number( (double) temp, 'f', 5 ) );
     }

     Get_FPar_ADwin( DMEAS2, &temp, &Error );
     if ( Error > 0 ){
	  Error_Message = Get_Error_Text( Error );
	  Print_Error_Message( (int) Error, Error_Message );
	  *Is_Successful = 0;
	  Print_To_Log( (int) Error, Error_Message );
     } else {
	  DCtrl2_DSpinBox->setValue( (double) temp );
	  DMeas2_lineEdit->setText( QString::number( (double) temp, 'f', 5 ) );
     }

}

void PMWidget::Set_Initial_Limits( int *const Is_Successful )
{

     int32_t Error;
     QString Error_Message;

     /* Assuming everything is OK */
     *Is_Successful = 1;

     DMin1_lineEdit->setText( QString::number( (double) Default_DMin1, 'f', 5 ) );
     DMin2_lineEdit->setText( QString::number( (double) Default_DMin2, 'f', 5 ) );

     DMax1_lineEdit->setText( QString::number( (double) Default_DMax1, 'f', 5 ) );
     DMax2_lineEdit->setText( QString::number( (double) Default_DMax2, 'f', 5 ) );

     Set_FPar_ADwin( DMIN1, Default_DMin1, &Error );
     if ( Error > 0 ){
	  Error_Message = Get_Error_Text( Error );
	  Print_Error_Message( (int) Error, Error_Message );
	  Print_To_Log( (int) Error, Error_Message );
	  *Is_Successful = 0;
     }

     Set_FPar_ADwin( DMIN2, Default_DMin2, &Error );
     if ( Error > 0 ){
	  Error_Message = Get_Error_Text( Error );
	  Print_Error_Message( (int) Error, Error_Message );
	  Print_To_Log( (int) Error, Error_Message );
	  *Is_Successful = 0;
     }

     Set_FPar_ADwin( DMAX1, Default_DMax1, &Error );
     if ( Error > 0 ){
	  Error_Message = Get_Error_Text( Error );
	  Print_Error_Message( (int) Error, Error_Message );
	  Print_To_Log( (int) Error, Error_Message );
	  *Is_Successful = 0;
     }

     Set_FPar_ADwin( DMAX2, Default_DMax2, &Error );
     if ( Error > 0 ){
	  Error_Message = Get_Error_Text( Error );
	  Print_Error_Message( (int) Error, Error_Message );
	  Print_To_Log( (int) Error, Error_Message );
	  *Is_Successful = 0;
     }

}


void PMWidget::Set_Initial_PID_Values( int *const Is_Successful )
{

     int32_t Error;
     QString Error_Message;

     /* Assuming everything is OK */
     *Is_Successful = 1;

     P_DSpinBox->setValue( (double) Default_P );
     I_DSpinBox->setValue( (double) Default_I );
     D_DSpinBox->setValue( (double) Default_D );

     Set_FPar_ADwin( PID_P, Default_P, &Error );
     if ( Error > 0 ){
	  Error_Message = Get_Error_Text( Error );
	  Print_Error_Message( (int) Error, Error_Message );
	  Print_To_Log( (int) Error, Error_Message );
	  *Is_Successful = 0;
     }

     Set_FPar_ADwin( PID_I, Default_I, &Error );
     if ( Error > 0 ){
	  Error_Message = Get_Error_Text( Error );
	  Print_Error_Message( (int) Error, Error_Message );
	  Print_To_Log( (int) Error, Error_Message );
	  *Is_Successful = 0;
     }

     Set_FPar_ADwin( PID_D, Default_D, &Error );
     if ( Error > 0 ){
	  Error_Message = Get_Error_Text( Error );
	  Print_Error_Message( (int) Error, Error_Message );
	  Print_To_Log( (int) Error, Error_Message );
	  *Is_Successful = 0;
     }

}

void PMWidget::Read_Measurement_Values()
{
     int Message;
     int32_t Error;
     QString Error_Message;
     float temp;

     Get_FPar_ADwin( DMEAS1, &temp, &Error );
     if ( Error > 0 ){
	  Error_Message = Get_Error_Text( Error );
	  Error_Message.append( " Read_Measurement_Values(): " );
	  Print_Error_Message( (int) Error, Error_Message );
	  Print_To_Log( (int) Error, Error_Message );
     } else {
	  DMeas1_lineEdit->setText( QString::number( (double) temp, 'f', 5 ) );
     }

     Get_FPar_ADwin( DMEAS2, &temp, &Error );
     if ( Error > 0 ){
	  Error_Message = Get_Error_Text( Error );
	  Error_Message.append( " Read_Measurement_Values(): " );
	  Print_Error_Message( (int) Error, Error_Message );
	  Print_To_Log( (int) Error, Error_Message );
     } else {
	  DMeas2_lineEdit->setText( QString::number( (double) temp, 'f', 5 ) );
     }

     Get_Par_ADwin( MESSAGE_FLAG, &Message, &Error );
     if ( Error > 0 ){
	  Error_Message = Get_Error_Text( Error );
	  Error_Message.append( " Read_Measurement_Values(): " );
	  Print_Error_Message( (int) Error, Error_Message );
	  Print_To_Log( (int) Error, Error_Message );
     } else {
	  if ( Message > 0 ){
	       Error = -1;
	       switch ( Message ){
	       case OVER_DMIN :
		    Error_Message = "EMERGENCY_STOP. The measured displacement is under the minimum allowed value.";
		    Print_Error_Message( (int) Error, Error_Message );
		    Print_To_Log( (int) Error, Error_Message );
		    break;
	       case OVER_DMAX :
		    Error_Message = "EMERGENCY_STOP. The measured displacement is over the maximum allowed value.";
		    Print_Error_Message( (int) Error, Error_Message );
		    Print_To_Log( (int) Error, Error_Message );
		    break;						
	       case OVER_VMIN : 
		    Error_Message = "EMERGENCY_STOP. The measured velocity is under the minimum allowed value.";
		    Print_Error_Message( (int) Error, Error_Message );
		    Print_To_Log( (int) Error, Error_Message );
		    break;
	       case OVER_VMAX : 
		    Error_Message = "EMERGENCY_STOP. The measured velocity is over the maximum allowed value.";
		    Print_Error_Message( (int) Error, Error_Message );
		    Print_To_Log( (int) Error, Error_Message );	    
		    break;
	       case OVER_AMIN : 
		    Error_Message = "EMERGENCY_STOP. The measured acceleration is under the minimum allowed value.";
		    Print_Error_Message( (int) Error, Error_Message );
		    Print_To_Log( (int) Error, Error_Message );	    
		    break;
	       case OVER_AMAX : 
		    Error_Message = "EMERGENCY_STOP. The measured acceleration is over the maximum allowed value.";
		    Print_Error_Message( (int) Error, Error_Message );
		    Print_To_Log( (int) Error, Error_Message );	    
		    break;
	       case OVER_FMAX : 
		    Error_Message = "EMERGENCY_STOP. The measured force is over the maximum allowed value.";
		    Print_Error_Message( (int) Error, Error_Message );
		    Print_To_Log( (int) Error, Error_Message );	    
		    break;
	       }
	  }
     }
}

void PMWidget::On_StopCntrlPButton_Clicked( )
{
     int32_t Error;
     QString Error_Message;

     Stop_Process_ADwin( 1, &CntrlP_is_Started, &Error );

     if ( Error > 0 ){
	  Error_Message = Get_Error_Text( Error );
	  Print_Error_Message( (int) Error, Error_Message );
	  Print_To_Log( (int) Error, Error_Message );
     } else if ( Error == -1 ){
	  Error_Message = "The Control process is stoped. Doing nothing.";
	  Print_Error_Message( (int) Error, Error_Message );
	  Print_To_Log( (int) Error, Error_Message );
     } else if ( Error == 0 ){
	  CntrlP_is_Started = 0;
	  StartCntrlPButton->setEnabled( 1 );
	  Error_Message = "The Control process has been successfully stopped.";
	  Print_To_Log( (int) Error, Error_Message );
     } else {
	  Error_Message = "Received invalid process status.";
	  Print_Error_Message( (int) -1, Error_Message );
	  Print_To_Log( (int) Error, Error_Message );
     }
}

void PMWidget::On_TFiletoolButton_Clicked( )
{
    
     TFile_LineEdit->setText( QFileDialog::getOpenFileName(this, tr("Open ADwin Binary File"), "/home/ferran/workspace/ADwinGUI/PMWidget", tr("ADwin binary file (*TB2)")) );

}

void PMWidget::On_LoadTestPButton_Clicked( )
{

     std::string File_Path;
     int32_t Error;
     QString Error_Message;

     File_Path = TFile_LineEdit->text().toStdString();
     Load_Process_ADwin( File_Path.c_str(), &TestP_is_Loaded, &Error );

     if ( Error > 0 ){
	  Error_Message = Get_Error_Text( Error );
	  Print_Error_Message( (int) Error, Error_Message );
     } else if ( Error == -1 ){
	  Error_Message = "The Test process is already loaded. Boot ADwin or clear the process in order to reload it";
	  Print_Error_Message( (int) Error, Error_Message );
     } else if ( Error == 0 ){
	  TestP_is_Loaded = 1;
	  /* Set the Start Test button to enabled */
	  StartTestButton->setEnabled( TestP_is_Loaded );
	  /* Set the Load Test Data button to enabled */
	  LoadTDataButton->setEnabled( TestP_is_Loaded );

	  /* Print the message to the log */ 
	  Error_Message = "Testing process successfully loaded.";
	  Print_To_Log( (int) Error, Error_Message );
     } else {
	  assert( Error < -1 );   /* The error can never have a value less than -1 */
	  exit( EXIT_FAILURE );
     }
}

void PMWidget::On_StartTestButton_Clicked( )
{

     int32_t Error;
     QString Error_Message;

     Start_Process_ADwin( (int32_t) 2, &TestP_is_Started, &Error );

     if ( Error > 0 ){
	  Error_Message = Get_Error_Text( Error );
	  Print_Error_Message( (int) Error, Error_Message );
	  Print_To_Log( (int) Error, Error_Message );
     } else if ( Error == -1 ){
	  Error_Message = "The Test process is already started. Boot ADwin or stop the process in order to reload it";
	  Print_Error_Message( (int) Error, Error_Message );
	  Print_To_Log( (int) Error, Error_Message );
     } else if ( Error == 0 ){
	  TestP_is_Started = 1;   /* The test process is running */
	  StopTestButton->setEnabled( TestP_is_Started );
	  ReadDataButton->setEnabled( TestP_is_Started );

	  /* Print the message to the log */
	  Error_Message = "Test process successfully started.";
	  Print_To_Log( (int) Error, Error_Message );
     } else {
	  assert( Error < -1 );   /* The error can never have a value less than -1 */
	  exit( EXIT_FAILURE );
     }
}

void PMWidget::On_StopTestButton_Clicked( )
{
     int32_t Error;
     QString Error_Message;

     Stop_Process_ADwin( 2, &TestP_is_Started, &Error );

     if ( Error > 0 ){
	  Error_Message = Get_Error_Text( Error );
	  Print_Error_Message( (int) Error, Error_Message );
	  Print_To_Log( (int) Error, Error_Message );
     } else if ( Error == -1 ){
	  Error_Message = "The Test process is stoped. Doing nothing.";
	  Print_Error_Message( (int) Error, Error_Message );
	  Print_To_Log( (int) Error, Error_Message );
     } else if ( Error == 0 ){
	  TestP_is_Started = 0;
	  StartTestButton->setEnabled( 1 );
	  /* Print the message to the log */
	  Error_Message = "The Test process has been successfully stopped.";
	  Print_To_Log( (int) Error, Error_Message );
     } else {
	  Error_Message = "Received invalid process status.";
	  Print_Error_Message( (int) -1, Error_Message );
	  Print_To_Log( (int) Error, Error_Message );
     }
}
     
void PMWidget::On_ReadDataButton_Clicked( )
{
 
     int i;
     int32_t Error;
     QString Error_Message;

     std::ofstream OutputFile;
     float *Data_182, *Data_183, *Data_184, *Data_185, *Data_186;
     int32_t Data_Length;

     Get_DataLength_ADwin( 182, &Data_Length, &Error );
     if ( Error > 0 ){
	  Error_Message = Get_Error_Text( Error );
	  Print_Error_Message( (int) Error, Error_Message );
	  Print_To_Log( (int) Error, Error_Message );
     } else {
	  Data_182 = new float[Data_Length];
	  Get_DataFloat_ADwin( 182, Data_182, (int32_t) 1, (int32_t) Data_Length, &Error );
	  if ( Error > 0 ){
	       Error_Message = Get_Error_Text( Error );
	       Print_Error_Message( (int) Error, Error_Message );
	       Print_To_Log( (int) Error, Error_Message );
	  }

	  Data_183 = new float[Data_Length];
	  Get_DataFloat_ADwin( 183, Data_183, (int32_t) 1, (int32_t) Data_Length, &Error );
	  if ( Error > 0 ){
	       Error_Message = Get_Error_Text( Error );
	       Print_Error_Message( (int) Error, Error_Message );
	       Print_To_Log( (int) Error, Error_Message );
	  }

     	  Data_184 = new float[Data_Length];
	  Get_DataFloat_ADwin( 184, Data_184, (int32_t) 1, (int32_t) Data_Length, &Error );
	  if ( Error > 0 ){
	       Error_Message = Get_Error_Text( Error );
	       Print_Error_Message( (int) Error, Error_Message );
	       Print_To_Log( (int) Error, Error_Message );
	  }

     	  Data_185 = new float[Data_Length];
	  Get_DataFloat_ADwin( 185, Data_185, (int32_t) 1, (int32_t) Data_Length, &Error );
	  if ( Error > 0 ){
	       Error_Message = Get_Error_Text( Error );
	       Print_Error_Message( (int) Error, Error_Message );
	       Print_To_Log( (int) Error, Error_Message );
	  }

	  Data_186 = new float[Data_Length];
	  Get_DataFloat_ADwin( 186, Data_186, (int32_t) 1, (int32_t) Data_Length, &Error );
	  if ( Error > 0 ){
	       Error_Message = Get_Error_Text( Error );
	       Print_Error_Message( (int) Error, Error_Message );
	       Print_To_Log( (int) Error, Error_Message );
	  }
  
	  OutputFile.open( OFile_LineEdit->text().toStdString().c_str() );

	  if (OutputFile.is_open() ){
	       OutputFile.precision(5);
	       for (i = 0; i < Data_Length; i++ ){
		    OutputFile << std::scientific << i << "," << Data_182[i] << "," << Data_183[i] << ",";
		    OutputFile << std::scientific << Data_184[i] << ","<< Data_185[i] << "," << Data_186[i] << std::endl;
	       }
	       OutputFile.close();
	  } else {
	       Error_Message = "Could not open the output file";
	       Print_Error_Message( (int) -1, Error_Message );
	       Print_To_Log( (int) Error, Error_Message );
	  }    

	  delete[ ] Data_182;
	  delete[ ] Data_183;
	  delete[ ] Data_184;
	  delete[ ] Data_185;
	  delete[ ] Data_186;

	  Error_Message = "The Data has been transfered successfully";
	  Print_Error_Message( (int) 0, Error_Message );
	  Print_To_Log( (int) 0, Error_Message );
     }
}

void PMWidget::On_ADwinReadyButton_Clicked( )
{
     int32_t Error;
     QString Error_Message;

     /* Set the variable PAR_11 ADWIN_READY in ADwin to 1. */
     Set_Par_ADwin( ADWIN_READY, (int32_t) 1, &Error );
     
     if ( Error > 0 ){
	  Error_Message = Get_Error_Text( Error );
	  Print_Error_Message( (int) Error, Error_Message );
	  Print_To_Log( (int) 0, Error_Message );
     } else if ( Error == 0 ){
	  LowPressureButton->setEnabled( 1 );
     }
}

void PMWidget::On_LowPressureButton_Clicked( )
{
     int32_t Value_Box, Value_Pressure, Error;
     QString Error_Message;

     Value_Box = 0;
     Value_Pressure = 0;

     /* Get the value of PAR_14 CONTROLBOX_READY */
     Get_Par_ADwin( CONTROLBOX_READY, &Value_Box, &Error );

     if( Error > 0 ){
	  Error_Message = Get_Error_Text( Error );
	  Print_Error_Message( (int) Error, Error_Message );
	  Print_To_Log( (int) 0, Error_Message );
     } else if ( Error == 0 ){
	  if ( Value_Box == 1 ){
	       Get_Par_ADwin( CURRENT_PRESSURE, &Value_Pressure, &Error );
	       if( Error > 0 ){
		    Error_Message = Get_Error_Text( Error );
		    Print_Error_Message( (int) Error, Error_Message );
		    Print_To_Log( (int) 0, Error_Message );
	       } else if ( Error == 0 ){
		    if ( Value_Pressure == NO_PRESSURE || Value_Pressure == HIGH_PRESSURE ){
			 Set_Par_ADwin( NEW_PRESSURE_COMMAND, LOW_PRESSURE, &Error );
			 if( Error > 0 ){
			      Error_Message = Get_Error_Text( Error );
			      Print_Error_Message( (int) Error, Error_Message );
			      Print_To_Log( (int) 0, Error_Message );
			 } else if ( Error == 0 ){
			      HighPressureButton->setEnabled( 1 );
			      HighPressureButton->setChecked( 0 );
			      StopPressureButton->setEnabled( 1 );
			      LowPressureButton->setEnabled( 0 );
			      EmergencyButton->setEnabled( 1 );
			      /* Print the message to the log */
			      Error_Message = "Pressure set to 'Low Pressure'";
			      Print_To_Log( (int) Error, Error_Message );
			 }
		    } else if ( Value_Pressure == LOW_PRESSURE ){
			 Error_Message = "Pressure is already in 'Low Pressure'. Doing nothing.";
			 Print_Error_Message( (int) -1, Error_Message );
			 Print_To_Log( (int) -1, Error_Message );
		    } else {
			 Error_Message = "Received invalid status of hydraulic pressure.";
			 Print_Error_Message( (int) -1, Error_Message );
			 Print_To_Log( (int) -1, Error_Message );
		    }
	       }
	  } else if ( Value_Box == 0 ){
	       Error_Message = "The control box is not yet ready. Enable it before setting pressure";
	       Print_Error_Message( (int) Error, Error_Message );
	       Print_To_Log( (int) Error, Error_Message );
	  }
     }
}

void PMWidget::On_HighPressureButton_Clicked( )
{
     int32_t Value_Box, Value_Pressure, Error;
     QString Error_Message;

     Value_Box = 0;
     Value_Pressure = 0;

     /* Get the value of PAR_14 CONTROLBOX_READY */
     Get_Par_ADwin( CONTROLBOX_READY, &Value_Box, &Error );

     if( Error > 0 ){
	  Error_Message = Get_Error_Text( Error );
	  Print_Error_Message( (int) Error, Error_Message );
	  Print_To_Log( (int) Error, Error_Message );
     } else if ( Error == 0 ){
	  if ( Value_Box == 1 ){
	       Get_Par_ADwin( CURRENT_PRESSURE, &Value_Pressure, &Error );
	       if( Error > 0 ){
		    Error_Message = Get_Error_Text( Error );
		    Print_Error_Message( (int) Error, Error_Message );
		    Print_To_Log( (int) Error, Error_Message );
	       } else if ( Error == 0 ){
		    if ( Value_Pressure == LOW_PRESSURE ){
			 Set_Par_ADwin( NEW_PRESSURE_COMMAND, HIGH_PRESSURE, &Error );
			 if( Error > 0 ){
			      Error_Message = Get_Error_Text( Error );
			      Print_Error_Message( (int) Error, Error_Message );
			      Print_To_Log( (int) Error, Error_Message );
			 } else if ( Error == 0 ){
			      HighPressureButton->setEnabled( 0 );
			      LowPressureButton->setEnabled( 1 );
			      LowPressureButton->setChecked( 0 );

			      /* Print the message to the log */
			      Error_Message = "Pressure set to 'High Pressure'";
			      Print_To_Log( (int) Error, Error_Message );
			 }
		    } else if ( Value_Pressure == HIGH_PRESSURE ){
			 Error_Message = "Pressure is already in 'Low Pressure'. Doing nothing.";
			 Print_Error_Message( (int) -1, Error_Message );
			 Print_To_Log( (int) -1, Error_Message );
		    } else if ( Value_Pressure == NO_PRESSURE ){
			 Error_Message = "There is no pressure in the system. Cannot set 'high pressure' from a 'no pressure' status";
			 Print_Error_Message( (int) -1, Error_Message );
			 Print_To_Log( (int) -1, Error_Message );
		    } else {
			 Error_Message = "Received invalid status of hydraulic pressure.";
			 Print_Error_Message( (int) -1, Error_Message );
			 Print_To_Log( (int) -1, Error_Message );
		    }
	       }
	  } else if ( Value_Box == 0 ){
	       Error_Message = "The control box is not yet ready. Enable it before setting pressure";
	       Print_Error_Message( (int) Error, Error_Message );
	       Print_To_Log( (int) Error, Error_Message );
	  }
     }
}


void PMWidget::On_StopPressureButton_Clicked( )
{
     int32_t Value_Box, Value_Pressure, Error;
     QString Error_Message;

     Value_Box = 0;
     Value_Pressure = 0;

     /* Get the value of PAR_14 CONTROLBOX_READY */
     Get_Par_ADwin( CONTROLBOX_READY, &Value_Box, &Error );

     if( Error > 0 ){
	  Error_Message = Get_Error_Text( Error );
	  Print_Error_Message( (int) Error, Error_Message );
	  Print_To_Log( (int) Error, Error_Message );
     } else if ( Error == 0 ){
	  if ( Value_Box == 1 ){
	       Get_Par_ADwin( CURRENT_PRESSURE, &Value_Pressure, &Error );
	       if( Error > 0 ){
		    Error_Message = Get_Error_Text( Error );
		    Print_Error_Message( (int) Error, Error_Message );
		    Print_To_Log( (int) Error, Error_Message );
	       } else if ( Error == 0 ){
		    if ( Value_Pressure == LOW_PRESSURE || Value_Pressure == HIGH_PRESSURE ){
			 Set_Par_ADwin( NEW_PRESSURE_COMMAND, STOP_PRESSURE, &Error );
			 if( Error > 0 ){
			      Error_Message = Get_Error_Text( Error );
			      Print_Error_Message( (int) Error, Error_Message );
			      Print_To_Log( (int) Error, Error_Message );
			 } else if ( Error == 0 ){
			      StopPressureButton->setEnabled( 0 );
			      LowPressureButton->setEnabled( 1 );
			      LowPressureButton->setChecked( 0 );
			      HighPressureButton->setEnabled( 0 );
			      HighPressureButton->setChecked( 0 );

			      /* Print the message to the log */
			      Error_Message = "Pressure successfully stopped";
			      Print_To_Log( (int) Error, Error_Message );
			 }
		    } else if ( Value_Pressure == NO_PRESSURE ){
			 Error_Message = "There is no pressure in the system. There is no need to stop it.";
			 Print_Error_Message( (int) -1, Error_Message );
			 Print_To_Log( (int) -1,  Error_Message );
		    } else {
			 Error_Message = "Received invalid status of hydraulic pressure.";
			 Print_Error_Message( (int) -1, Error_Message );
			 Print_To_Log( (int) -1, Error_Message );
		    }
	       }
	  } else if ( Value_Box == 0 ){
	       Error_Message = "The control box is not yet ready. Enable it before setting pressure";
	       Print_Error_Message( (int) Error, Error_Message );
	       Print_To_Log( (int) Error, Error_Message );
	  }
     }
}

void PMWidget::On_EmergencyButton_Clicked( )
{
     int32_t Error;
     QString Error_Message;

     /* Set the variable PAR_11 ADWIN_READY in ADwin to 1. */
     Set_Par_ADwin( EMERGENCY_STOP, (int32_t) 1, &Error );
     
     if ( Error > 0 ){
	  Error_Message = Get_Error_Text( Error );
	  Print_Error_Message( (int) Error, Error_Message );
	  Print_To_Log( (int) Error, Error_Message );
     } else if ( Error == 0 ){
	  Error_Message = "Emergency stop activated. Save data and reboot ADwin" ;
	  Print_Error_Message( (int) -1, Error_Message );
	  Print_To_Log( (int) -1, Error_Message );
     }
}

void PMWidget::On_IFiletoolButton_Clicked( )
{
    
     IFile_LineEdit->setText( QFileDialog::getOpenFileName(this, tr("Input Data File"), "/home/ferran/workspace/ADwinGUI", tr("Input Data File (*)")) );

}

void PMWidget::On_OFiletoolButton_Clicked( )
{
    
     OFile_LineEdit->setText( QFileDialog::getOpenFileName(this, tr("Output Data File"), "/home/ferran/workspace/ADwinGUI", tr("Output Data File (*)")) );

}

void PMWidget::Disable_All_Buttons( )
{
     
     LoadCtrlPButton->setEnabled( 0 );
     StartCntrlPButton->setEnabled( 0 );
     StopCntrlPButton->setEnabled( 0 );

     LoadTestPButton->setEnabled( 0 );
     StartTestButton->setEnabled( 0 );
     StopTestButton->setEnabled( 0 );

     LoadTDataButton->setEnabled( 0 );

     ReadDataButton->setEnabled( 0 );

     ADwinReadyButton->setEnabled( 0 );
     EmergencyButton->setEnabled( 0 );

     LowPressureButton->setEnabled( 0 );
     LowPressureButton->setChecked( 0 );
     HighPressureButton->setEnabled( 0 );
     HighPressureButton->setChecked( 0 );
     StopPressureButton->setEnabled( 0 );

     DefValPButton->setEnabled( 0 );
     SetLimPButton->setEnabled( 0 );

     /* Disable check boxes */
     EnabLimCheckBox->setEnabled( 0 );

     /* Disable spin boxes */
     DCtrl1_DSpinBox->setEnabled( 0 );
     DCtrl2_DSpinBox->setEnabled( 0 );
     P_DSpinBox->setEnabled( 0 );
     I_DSpinBox->setEnabled( 0 );
     D_DSpinBox->setEnabled( 0 );

}

void PMWidget::On_DCtrl1_DSpinBox_valueChanged( const double New_Value_Double )
{
     float New_Value;
     int32_t Is_Testing, Error;
     QString Error_Message;

     New_Value = (float) New_Value_Double;

     if ( New_Value < DCntrlMin1 ){
	  Error = -1;
	  Error_Message = "The input value is greater than the minimum allowed value. It will be changed to the minimum value allowed.";
	  Print_Error_Message( (int) Error, Error_Message );
	  Print_To_Log( (int) Error, Error_Message );
     } else if ( New_Value > DCntrlMax1 ){
	  Error = -1;
	  Error_Message = "The input value is greater than the minimum allowed value. It will be changed to the maximum value allowed.";
	  Print_Error_Message( (int) Error, Error_Message );
	  Print_To_Log( (int) Error, Error_Message );
     } else if ( TestP_is_Started ){
	  /* Get the value of the variable TESTING in ADwin to see if it is running a test or not */
	  Get_Par_ADwin( TESTING, &Is_Testing, &Error );
	  if( Error > 0 ){
	       Error_Message = Get_Error_Text( Error );
	       Print_Error_Message( (int) Error, Error_Message );
	       Print_To_Log( (int) Error, Error_Message );
	  } else if ( Error == 0 ){
	       /* If there is no test running set the new control value */
	       if ( !Is_Testing ){
		    Set_FPar_ADwin( DCNTRL1, New_Value, &Error );
		    if ( Error > 0 ){
			 Error_Message = Get_Error_Text( Error );
			 Print_Error_Message( (int) Error, Error_Message );
			 Print_To_Log( (int) Error, Error_Message );
		    }
	       } else {
		    Error = -1;
		    Error_Message = "The control value cannot be changed during a Test. Aborting operation";
		    Print_Error_Message( (int) Error, Error_Message );
		    Print_To_Log( (int) Error, Error_Message );
	       }
	  }
     }
}

void PMWidget::On_DCtrl2_DSpinBox_valueChanged( const double New_Value_Double )
{
     float New_Value;
     int32_t Is_Testing, Error;
     QString Error_Message;

     New_Value = (float) New_Value_Double;

     if ( New_Value < DCntrlMin2 ){
	  Error = -1;
	  Error_Message = "The input value is greater than the minimum allowed value. It will be changed to the minimum value allowed.";
	  Print_Error_Message( (int) Error, Error_Message );
	  Print_To_Log( (int) Error, Error_Message );
     } else if ( New_Value > DCntrlMax2 ){
	  Error = -1;
	  Error_Message = "The input value is greater than the minimum allowed value. It will be changed to the maximum value allowed.";
	  Print_Error_Message( (int) Error, Error_Message );
	  Print_To_Log( (int) Error, Error_Message );
     } else if ( TestP_is_Started ){
	  /* Get the value of the variable TESTING in ADwin to see if it is running a test or not */
	  Get_Par_ADwin( TESTING, &Is_Testing, &Error );
	  if( Error > 0 ){
	       Error_Message = Get_Error_Text( Error );
	       Print_Error_Message( (int) Error, Error_Message );
	       Print_To_Log( (int) Error, Error_Message );
	  } else if ( Error == 0 ){
	       /* If there is no test running set the new control value */
	       if ( !Is_Testing ){
		    Set_FPar_ADwin( DCNTRL2, New_Value, &Error );
		    if ( Error > 0 ){
			 Error_Message = Get_Error_Text( Error );
			 Print_Error_Message( (int) Error, Error_Message );
			 Print_To_Log( (int) Error, Error_Message );
		    }
	       } else {
		    Error = -1;
		    Error_Message = "The control value cannot be changed during a Test. Aborting operation";
		    Print_Error_Message( (int) Error, Error_Message );
		    Print_To_Log( (int) Error, Error_Message );
	       }
	  }
     }
}

void PMWidget::On_P_DSpinBox_valueChanged( const double New_Value_Double )
{
     float New_Value;
     int32_t Is_Testing, Error;
     QString Error_Message;

     New_Value = (float) New_Value_Double;

     if ( New_Value < P_Min ){
	  Error = -1;
	  Error_Message = "The input value is greater than the minimum allowed value. It will be changed to the minimum value allowed.";
	  Print_Error_Message( (int) Error, Error_Message );
	  Print_To_Log( (int) Error, Error_Message );
     } else if ( New_Value > P_Max ){
	  Error = -1;
	  Error_Message = "The input value is greater than the minimum allowed value. It will be changed to the maximum value allowed.";
	  Print_Error_Message( (int) Error, Error_Message );
	  Print_To_Log( (int) Error, Error_Message );
     } else if ( TestP_is_Started ){
	  /* Get the value of the variable TESTING in ADwin to see if it is running a test or not */
	  Get_Par_ADwin( TESTING, &Is_Testing, &Error );
	  if( Error > 0 ){
	       Error_Message = Get_Error_Text( Error );
	       Print_Error_Message( (int) Error, Error_Message );
	       Print_To_Log( (int) Error, Error_Message );
	  } else if ( Error == 0 ){
	       /* If there is no test running set the new control value */
	       if ( !Is_Testing ){
		    Set_FPar_ADwin( PID_P, New_Value, &Error );
		    if ( Error > 0 ){
			 Error_Message = Get_Error_Text( Error );
			 Print_Error_Message( (int) Error, Error_Message );
			 Print_To_Log( (int) Error, Error_Message );
		    }
	       } else {
		    Error = -1;
		    Error_Message = "The control value cannot be changed during a Test. Aborting operation";
		    Print_Error_Message( (int) Error, Error_Message );
		    Print_To_Log( (int) Error, Error_Message );
	       }
	  }
     }
}

void PMWidget::On_I_DSpinBox_valueChanged( const double New_Value_Double )
{
     float New_Value;
     int32_t Is_Testing, Error;
     QString Error_Message;

     New_Value = (float) New_Value_Double;

     if ( New_Value < I_Min ){
	  Error = -1;
	  Error_Message = "The input value is greater than the minimum allowed value. It will be changed to the minimum value allowed.";
	  Print_Error_Message( (int) Error, Error_Message );
	  Print_To_Log( (int) Error, Error_Message );
     } else if ( New_Value > I_Max ){
	  Error = -1;
	  Error_Message = "The input value is greater than the minimum allowed value. It will be changed to the maximum value allowed.";
	  Print_Error_Message( (int) Error, Error_Message );
	  Print_To_Log( (int) Error, Error_Message );
     } else if ( TestP_is_Started ){
	  /* Get the value of the variable TESTING in ADwin to see if it is running a test or not */
	  Get_Par_ADwin( TESTING, &Is_Testing, &Error );
	  if( Error > 0 ){
	       Error_Message = Get_Error_Text( Error );
	       Print_Error_Message( (int) Error, Error_Message );
	       Print_To_Log( (int) Error, Error_Message );
	  } else if ( Error == 0 ){
	       /* If there is no test running set the new control value */
	       if ( !Is_Testing ){
		    Set_FPar_ADwin( PID_I, New_Value, &Error );
		    if ( Error > 0 ){
			 Error_Message = Get_Error_Text( Error );
			 Print_Error_Message( (int) Error, Error_Message );
			 Print_To_Log( (int) Error, Error_Message );
		    }
	       } else {
		    Error = -1;
		    Error_Message = "The control value cannot be changed during a Test. Aborting operation";
		    Print_Error_Message( (int) Error, Error_Message );
		    Print_To_Log( (int) Error, Error_Message );
	       }
	  }
     }
}

void PMWidget::On_D_DSpinBox_valueChanged( const double New_Value_Double )
{
     float New_Value;
     int32_t Is_Testing, Error;
     QString Error_Message;

     New_Value = (float) New_Value_Double;

     if ( New_Value < D_Min ){
	  Error = -1;
	  Error_Message = "The input value is greater than the minimum allowed value. It will be changed to the minimum value allowed.";
	  Print_Error_Message( (int) Error, Error_Message );
	  Print_To_Log( (int) Error, Error_Message );
     } else if ( New_Value > D_Max ){
	  Error = -1;
	  Error_Message = "The input value is greater than the minimum allowed value. It will be changed to the maximum value allowed.";
	  Print_Error_Message( (int) Error, Error_Message );
	  Print_To_Log( (int) Error, Error_Message );
     } else if ( TestP_is_Started ){
	  /* Get the value of the variable TESTING in ADwin to see if it is running a test or not */
	  Get_Par_ADwin( TESTING, &Is_Testing, &Error );
	  if( Error > 0 ){
	       Error_Message = Get_Error_Text( Error );
	       Print_Error_Message( (int) Error, Error_Message );
	       Print_To_Log( (int) Error, Error_Message );
	  } else if ( Error == 0 ){
	       /* If there is no test running set the new control value */
	       if ( !Is_Testing ){
		    Set_FPar_ADwin( PID_D, New_Value, &Error );
		    if ( Error > 0 ){
			 Error_Message = Get_Error_Text( Error );
			 Print_Error_Message( (int) Error, Error_Message );
			 Print_To_Log( (int) Error, Error_Message );
		    }
	       } else {
		    Error = -1;
		    Error_Message = "The control value cannot be changed during a Test. Aborting operation";
		    Print_Error_Message( (int) Error, Error_Message );
		    Print_To_Log( (int) Error, Error_Message );
	       }
	  }
     }
}

void PMWidget::On_SetLimCheckBox_Checked( const int state )
{

     int32_t Error;
     QString Error_Message;

     if ( state == 0 ){
	  Set_Par_ADwin( LIMIT_CHECK, (int32_t) 0, &Error );
     
	  if ( Error > 0 ){
	       Error_Message = Get_Error_Text( Error );
	       Print_Error_Message( (int) Error, Error_Message );
	       Print_To_Log( (int) Error, Error_Message );
	  }
     } else {
	  Set_Par_ADwin( LIMIT_CHECK, (int32_t) 1, &Error );
     
	  if ( Error > 0 ){
	       Error_Message = Get_Error_Text( Error );
	       Print_Error_Message( (int) Error, Error_Message );
	       Print_To_Log( (int) Error, Error_Message );
	  }
     }
}

void PMWidget::On_DefValPButton_Clicked( )
{
     int Success_Cntr_Val, Success_Limit, Success_PID;
     int Error;
     QString Error_Message;

     /* Set default values */
     Set_Initial_Control_Values( &Success_Cntr_Val );
     Set_Initial_Limits( &Success_Limit );
     Set_Initial_PID_Values( &Success_PID );

     if ( Success_Cntr_Val && Success_Limit && Success_PID ){
	  Error = 0;
	  Error_Message = "Default values successfully implemented";

     } else {
	  Error = -1;
	  Error_Message = "The operation to start the control process could not be completed";
     }

     Print_Error_Message( Error, Error_Message );
     Print_To_Log( (int) Error, Error_Message );
}

void PMWidget::Print_To_Log( int Error, QString Error_Text )
{

     QString Full_Error_Text;
     QTime Time = QTime::currentTime( );

     Full_Error_Text = Time.toString( );
     Full_Error_Text.append( " - " );
     if ( Error == 0 ){
	  LogBrowser->setTextColor( Qt::black );
     } else {
	  /* Append the error number */
	  Full_Error_Text.append(QString("Error %1").arg(Error));
	  Full_Error_Text.append(": ");
	  LogBrowser->setTextColor( Qt::red );
     }
     Full_Error_Text.append( Error_Text );
     LogBrowser->append( Full_Error_Text );
}
