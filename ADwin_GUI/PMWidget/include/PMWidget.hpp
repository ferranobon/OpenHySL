#ifndef PMWidget_HPP
#define PMWidget_HPP

#include <QWidget>

#include "ui_PMWidget.h"
#include "ADwin.hpp"

class PMWidget : public QWidget, public Ui::PMWidget, public ADwin_Class
{
     Q_OBJECT

     public:
     PMWidget( QWidget *parent = 0 );
     int ADwin_is_Booted;

     int CntrlP_is_Loaded;
     int CntrlP_is_Started;

     int TestP_is_Loaded;
     int TestP_is_Started;

     void Print_To_Log( int Error, QString Error_Text );

signals:
     void Enable_Load_Process();

private slots:
     void On_BFiletoolButton_Clicked( );
     void On_BootButton_Clicked( );

     void On_PFiletoolButton_Clicked( );
     void On_LoadCtrlPButton_Clicked( );
     void On_StartCtrlPButton_Clicked( );
     void On_StopCntrlPButton_Clicked( );

     void On_TFiletoolButton_Clicked( );
     void On_LoadTestPButton_Clicked( );
     void On_StartTestButton_Clicked( );
     void On_StopTestButton_Clicked( );

     void On_ReadDataButton_Clicked( );
     
     void On_ADwinReadyButton_Clicked( );
     void On_EmergencyButton_Clicked( );

     void On_LowPressureButton_Clicked( );
     void On_HighPressureButton_Clicked( );
     void On_StopPressureButton_Clicked( );

     void On_IFiletoolButton_Clicked( );
     void On_OFiletoolButton_Clicked( );

     void On_DCtrl1_DSpinBox_valueChanged( const double New_Value_Double );
     void On_DCtrl2_DSpinBox_valueChanged( const double New_Value_Double );

     void On_P_DSpinBox_valueChanged( const double New_Value_Double );
     void On_I_DSpinBox_valueChanged( const double New_Value_Double );
     void On_D_DSpinBox_valueChanged( const double New_Value_Double );

     void On_SetLimCheckBox_Checked( const int state );
     
     void On_DefValPButton_Clicked( );

     void Read_Measurement_Values();

private:
     void Configure_Cntrl_SpinBoxes( );
     void Configure_PID_SpinBoxes( );
     void Disable_All_Buttons( );
     void Set_Initial_Control_Values( int *const Is_Successful );
     void Set_Initial_Limits( int *const Is_Successful );
     void Set_Initial_PID_Values( int *const Is_Successful );
};

#endif
