#include "ui_ConfFileGUI.h"

#include <QProcess>

class ConfFile : public QMainWindow, public Ui::ConfFile
{
     Q_OBJECT

     QProcess *SAlgorithm;

     public:
     ConfFile( QMainWindow *parent = 0);

private slots:
     void Check_If_Zero( );
     void Match_TimeIntegration_Cnts( );
     void Refresh_Preview( );
     void Set_TimeIntegration_Cnts(const QString str);

     void On_GenCFileButton_Clicked( );
     void On_ConfButton_Clicked( );
     void On_MassButton_Clicked( );
     void On_StiffButton_Clicked( );
     void On_DampButton_Clicked( );
     void On_GainButton_Clicked( );
     void On_LVectorButton_Clicked( );
     void On_CNodesButton_Clicked( );
     void On_GMotionButton_Clicked( );
     void On_OutFileButton_Clicked( );

     void On_Run_pushButton_Clicked();
     void readStdOutput();
     void On_OpenOutput_pushButton_Clicked();
};
