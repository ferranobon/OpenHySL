#include <QApplication>
#include <QFileDialog>  // To enable the file dialog
#include <QFile>
#include <QTextStream>
#include <QProcess>
#include <QTextCodec>
#include <cassert>

#include "ConfFile.hpp"

ConfFile::ConfFile( QMainWindow *parent )
     : QMainWindow( parent )
{
     setupUi( this );

     /* Set the beta value of rayleigh damping due to limitations in the QT4-Designer program */
     RayBeta_dSpinBox->setValue(0.0000929);

     /* Hide/show the elements that are (not) required */
     Damp_label->setVisible( ReadDamp_checkBox->isChecked() );
     Damp_lineEdit->setVisible( ReadDamp_checkBox->isChecked() );
     Damp_toolButton->setVisible( ReadDamp_checkBox->isChecked() );

     Gain_label->setVisible( ReadGain_checkBox->isChecked() );
     Gain_lineEdit->setVisible( ReadGain_checkBox->isChecked() );
     Gain_toolButton->setVisible( ReadGain_checkBox->isChecked() );
     /* Enable or disable Rayleigh damping */
     RayAlpha_dSpinBox->setDisabled( ReadGain_checkBox->isChecked() );
     RayBeta_dSpinBox->setDisabled( ReadGain_checkBox->isChecked() );

     LVector_label->setVisible( ReadLVector_checkBox->isChecked() );
     LVector_lineEdit->setVisible( ReadLVector_checkBox->isChecked() );
     LVector_toolButton->setVisible( ReadLVector_checkBox->isChecked() );

     LVectorForm_label->setHidden( ReadLVector_checkBox->isChecked() );
     LVectorForm_lineEdit->setHidden( ReadLVector_checkBox->isChecked() );

     //   adjustSize();
     this->layout()->setSizeConstraint(QLayout::SetFixedSize);

     /* Initialise the contents of the preview tab */
     Refresh_Preview();

     connect( AbsValues_checkBox, SIGNAL( stateChanged(int) ), this, SLOT( Refresh_Preview() ) );
     connect( ReadSparse_checkBox, SIGNAL( stateChanged(int) ), this, SLOT( Refresh_Preview() ) );
     connect( UsePacked_checkBox, SIGNAL( stateChanged(int) ), this, SLOT( Refresh_Preview() ) );
     connect( UseSparse_checkBox, SIGNAL( stateChanged(int) ), this, SLOT( Refresh_Preview() ) );
     connect( ReadLVector_checkBox, SIGNAL( stateChanged(int) ), this, SLOT( Refresh_Preview() ) );
     connect( ReadDamp_checkBox, SIGNAL( stateChanged(int) ), this, SLOT( Refresh_Preview() ) );
     connect( ReadGain_checkBox, SIGNAL( stateChanged(int) ), this, SLOT( Refresh_Preview() ) );

     connect( Order_SpinBox, SIGNAL( valueChanged(int) ), this, SLOT( Refresh_Preview() ) );
     connect( NumSteps_SpinBox, SIGNAL( valueChanged(int) ), this, SLOT( Refresh_Preview() ) );
     connect( TimeStep_dSpinBox, SIGNAL( valueChanged(double) ), this, SLOT( Refresh_Preview() ) );

     /* Newmark */
     connect( IntCnt_ComboBox, SIGNAL( currentIndexChanged(QString)), this, SLOT( Set_TimeIntegration_Cnts(QString) ));
     connect( Beta_dSpinBox, SIGNAL( valueChanged(double) ), this, SLOT( Refresh_Preview() ) );
     connect( Beta_dSpinBox, SIGNAL( valueChanged(double) ), this, SLOT( Match_TimeIntegration_Cnts() ) );
     connect( Gamma_dSpinBox, SIGNAL( valueChanged(double) ), this, SLOT( Refresh_Preview() ) );
     connect( Gamma_dSpinBox, SIGNAL( valueChanged(double) ), this, SLOT( Match_TimeIntegration_Cnts() ) );

     /* PID */
     connect( P_dSpinBox, SIGNAL( valueChanged(double) ), this, SLOT( Refresh_Preview() ) );
     connect( I_dSpinBox, SIGNAL( valueChanged(double) ), this, SLOT( Refresh_Preview() ) );
     connect( D_dSpinBox, SIGNAL( valueChanged(double) ), this, SLOT( Refresh_Preview() ) );

     /* Rayleigh */
     connect( RayAlpha_dSpinBox, SIGNAL( valueChanged(double) ), this, SLOT( Refresh_Preview() ) );
     connect( RayBeta_dSpinBox, SIGNAL( valueChanged(double) ), this, SLOT( Refresh_Preview() ) );

     /* Substructure */
     connect( NumSub_SpinBox, SIGNAL( valueChanged(int) ), this, SLOT( Refresh_Preview() ) );
     connect( NumSub_SpinBox, SIGNAL( valueChanged(int) ), this, SLOT( Check_If_Zero() ) );
     connect( NumSubSteps_SpinBox, SIGNAL( valueChanged(int) ), this, SLOT( Refresh_Preview() ) );

     /* FileNames */
     connect( Mass_lineEdit, SIGNAL( textChanged(QString) ), this, SLOT( Refresh_Preview() ) );
     connect( Stiff_lineEdit, SIGNAL( textChanged(QString) ), this, SLOT( Refresh_Preview() ) );
     connect( Damp_lineEdit, SIGNAL( textChanged(QString) ), this, SLOT( Refresh_Preview() ) );
     connect( Gain_lineEdit, SIGNAL( textChanged(QString) ), this, SLOT( Refresh_Preview() ) );
     connect( LVector_lineEdit, SIGNAL( textChanged(QString) ), this, SLOT( Refresh_Preview() ) );
     connect( CNodes_lineEdit, SIGNAL( textChanged(QString) ), this, SLOT( Refresh_Preview() ) );
     connect( GMotion_lineEdit, SIGNAL( textChanged(QString) ), this, SLOT( Refresh_Preview() ) );
     connect( OutFile_lineEdit, SIGNAL( textChanged(QString) ), this, SLOT( Refresh_Preview() ) );

     connect( ConfFile_toolButton, SIGNAL( clicked() ), this, SLOT( On_ConfButton_Clicked()));
     connect( Mass_toolButton, SIGNAL( clicked() ), this, SLOT( On_MassButton_Clicked()));
     connect( Stiff_toolButton, SIGNAL( clicked() ), this, SLOT( On_StiffButton_Clicked()));
     connect( Damp_toolButton, SIGNAL( clicked() ), this, SLOT( On_DampButton_Clicked()));
     connect( Gain_toolButton, SIGNAL( clicked() ), this, SLOT( On_GainButton_Clicked()));
     connect( LVector_toolButton, SIGNAL( clicked() ), this, SLOT( On_LVectorButton_Clicked()));
     connect( CNodes_toolButton, SIGNAL( clicked() ), this, SLOT( On_CNodesButton_Clicked()));
     connect( GMotion_toolButton, SIGNAL( clicked() ), this, SLOT( On_GMotionButton_Clicked()));
     connect( OutFile_toolButton, SIGNAL( clicked() ), this, SLOT( On_OutFileButton_Clicked()));

     connect( GenCFile_pushButton, SIGNAL( clicked() ), this, SLOT( On_GenCFileButton_Clicked()));

     connect( Run_pushButton, SIGNAL( clicked() ), this, SLOT( On_Run_pushButton_Clicked()));
     connect( OpenOutput_pushButton, SIGNAL( clicked() ), this, SLOT( On_OpenOutput_pushButton_Clicked()));
}

void ConfFile::Set_TimeIntegration_Cnts(const QString str)
{

     if( strcmp(str.toStdString().c_str(), "Central explicit") == 0 ){
	  Beta_dSpinBox->setValue( 0.0 );
	  Gamma_dSpinBox->setValue( -0.5 );
     } else if ( strcmp(str.toStdString().c_str(), "Backward") == 0 ){
	  Beta_dSpinBox->setValue( 1.0 );
	  Gamma_dSpinBox->setValue( 3.0/2.0 );
     } else if ( strcmp(str.toStdString().c_str(), "Linear acceleration") == 0 ){
	  Beta_dSpinBox->setValue( 0.1 );
	  Gamma_dSpinBox->setValue( 0.5 );
     } else if ( strcmp(str.toStdString().c_str(), "Galerkin") == 0 ){
	  Beta_dSpinBox->setValue( 4.0/5.0 );
	  Gamma_dSpinBox->setValue( 3.0/2.0 );
     } else if ( strcmp(str.toStdString().c_str(), "Fox Goodwin") == 0 ){
	  Beta_dSpinBox->setValue( 1.0/12.0 );
	  Gamma_dSpinBox->setValue( 0.5 );
     } else if ( strcmp(str.toStdString().c_str(), "Average acceleration") == 0 ){
	  Beta_dSpinBox->setValue( 0.25 );
	  Gamma_dSpinBox->setValue( 0.5 );
     } else if ( strcmp(str.toStdString().c_str(), "Custom" ) == 0 ){
	  /* Do nothing */
     } else assert(0);
}

void ConfFile::Match_TimeIntegration_Cnts( )
{
     double Beta;
     double Gamma;

     Beta = Beta_dSpinBox->value();
     Gamma = Gamma_dSpinBox->value();

     if( Beta == 0.0 && Gamma == -0.5 ){
	  if( IntCnt_ComboBox->findText( "Central explicit" ) != -1 ){
	       IntCnt_ComboBox->setCurrentIndex( IntCnt_ComboBox->findText( "Central explicit" ));
	  } else assert(0);	       
     } else if ( Beta == 1.0 && Gamma == 3.0/2.0 ){
	  if( IntCnt_ComboBox->findText( "Backward" ) != -1 ){
	       IntCnt_ComboBox->setCurrentIndex( IntCnt_ComboBox->findText( "Backward" ));
	  } else assert(0);
     } else if ( Beta == 0.1 && Gamma == 0.5 ){
	  if( IntCnt_ComboBox->findText( "Linear acceleration" ) != -1 ){
	       IntCnt_ComboBox->setCurrentIndex( IntCnt_ComboBox->findText( "Linear acceleration" ));
	  } else assert(0);
     } else if ( Beta == 4.0/5.0 && Gamma == 3.0/2.0 ){
	  if( IntCnt_ComboBox->findText( "Galerkin" ) != -1 ){
	       IntCnt_ComboBox->setCurrentIndex( IntCnt_ComboBox->findText( "Galerkin" ));
	  } else assert(0);
     } else if ( Beta == 0.0833 && Gamma == 0.5 ){
	  if( IntCnt_ComboBox->findText( "Fox Goodwin" ) != -1 ){
	       IntCnt_ComboBox->setCurrentIndex( IntCnt_ComboBox->findText( "Fox Goodwin" ));
	  } else assert(0);
     } else if ( Beta == 0.25 && Gamma == 0.5 ){
	  if( IntCnt_ComboBox->findText( "Average acceleration" ) != -1 ){
	       IntCnt_ComboBox->setCurrentIndex( IntCnt_ComboBox->findText( "Average acceleration" ));
	  } else assert(0);
     } else {
	  if( IntCnt_ComboBox->findText( "Custom" ) != -1 ){
	       IntCnt_ComboBox->setCurrentIndex( IntCnt_ComboBox->findText( "Custom" ));
	  } else assert(0);
     }
}
     
void ConfFile::Check_If_Zero( )
{
     if( NumSub_SpinBox->value() == 0 ){
	  NumSubSteps_SpinBox->setDisabled( true );
	  ConfFile_tabWidget->setTabEnabled( 1, false );
	  CNodes_label->setVisible( false );
	  CNodes_lineEdit->setVisible( false );
	  CNodes_toolButton->setVisible( false );
     } else {
	  NumSubSteps_SpinBox->setEnabled( true );
	  ConfFile_tabWidget->setTabEnabled( 1, true );
	  CNodes_label->setVisible( true );
	  CNodes_lineEdit->setVisible( true );
	  CNodes_toolButton->setVisible( true );
     }
}

void ConfFile::Refresh_Preview( )
{

     QString str;

     /* Clear all the previous text */
     CFPreview_textBrowser->clear();

     /* General section */
     CFPreview_textBrowser->append( "[General]" );

     str = "Use_Absolute_Values = ";
     str.append( QString::number( AbsValues_checkBox->isChecked()) );
     CFPreview_textBrowser->append( str );

     str = "Read_Sparse = ";
     str.append( QString::number( ReadSparse_checkBox->isChecked()) );
     CFPreview_textBrowser->append( str );

     str = "Use_Sparse = ";
     str.append( QString::number( UseSparse_checkBox->isChecked()) );
     CFPreview_textBrowser->append( str );

     str = "Use_Packed = ";
     str.append( QString::number( UsePacked_checkBox->isChecked()) );
     CFPreview_textBrowser->append( str );

     str = "Order = ";
     str.append( QString::number( Order_SpinBox->value()) );
     CFPreview_textBrowser->append( str );

     str = "Read_CMatrix = ";
     str.append( QString::number( ReadDamp_checkBox->isChecked()) );
     CFPreview_textBrowser->append( str );

     str = "Read_GMatrix = ";
     str.append( QString::number( ReadGain_checkBox->isChecked()) );
     CFPreview_textBrowser->append( str );

     str = "Read_LVector = ";
     str.append( QString::number( ReadLVector_checkBox->isChecked()) );
     CFPreview_textBrowser->append( str );

     if( !ReadLVector_checkBox->isChecked() ){
	  str = "Excited_DOF = \"6 1 1 1 0 0 0\"";
	  CFPreview_textBrowser->append( str );
     }

     str = "Num_Steps = ";
     str.append( QString::number( NumSteps_SpinBox->value()) );
     CFPreview_textBrowser->append( str );

     str = "Delta = ";
     str.append( QString::number( TimeStep_dSpinBox->value()) );
     CFPreview_textBrowser->append( str );

     /* If the damping matrix is not read from a file, the entries specifying the Rayleigh damping are
      * required */
     if( !ReadDamp_checkBox->isChecked() ){
	  CFPreview_textBrowser->append( "\n[Rayleigh]" );
	  str = "Alpha = ";
	  str.append( QString::number( RayAlpha_dSpinBox->value()) );
	  CFPreview_textBrowser->append( str );
	  str = "Beta = ";
	  str.append( QString::number( RayBeta_dSpinBox->value()) );
	  CFPreview_textBrowser->append( str );
     }

     CFPreview_textBrowser->append( "\n[Newmark]" );
     str = "Beta = ";
     str.append( QString::number( Beta_dSpinBox->value()) );
     CFPreview_textBrowser->append( str );
     str = "Gamma = ";
     str.append( QString::number( Gamma_dSpinBox->value()) );
     CFPreview_textBrowser->append( str );

     CFPreview_textBrowser->append( "\n[PID]" );
     str = "P = ";
     str.append( QString::number( P_dSpinBox->value()) );
     CFPreview_textBrowser->append( str );
     str = "I = ";
     str.append( QString::number( I_dSpinBox->value()) );
     CFPreview_textBrowser->append( str );
     str = "D = ";
     str.append( QString::number( I_dSpinBox->value()) );
     CFPreview_textBrowser->append( str );

     CFPreview_textBrowser->append( "\n[Substructure]" );
     str = "Order = ";
     str.append( QString::number( NumSub_SpinBox->value()) );
     CFPreview_textBrowser->append( str );
     /* The entry with the number of steps is only required if there is at least one substructure */
     if( NumSub_SpinBox->value() != 0 ){
	  str = "Num_Substeps = ";
	  str.append( QString::number( NumSubSteps_SpinBox->value()) );
	  CFPreview_textBrowser->append( str );
     }

     CFPreview_textBrowser->append( "\n[FileNames]" );
     str = "Mass_Matrix = ";
     str.append( "\"" );
     str.append( Mass_lineEdit->text() );
     str.append( "\"" );
     CFPreview_textBrowser->append( str );
     str = "Stiffness_Matrix = ";
     str.append( "\"" );
     str.append( Stiff_lineEdit->text() );
     str.append( "\"" );
     CFPreview_textBrowser->append( str );

     /* The entry specifying the file with the damping matrix is only required if the appropiate checkbox is
      * checked */
     if( ReadDamp_checkBox->isChecked() ){
	  CFPreview_textBrowser->append( str );
	  str = "Damping_Matrix = ";
	  str.append( "\"" );
	  str.append( Damp_lineEdit->text() );
	  str.append( "\"" );
	  CFPreview_textBrowser->append( str );
     }
     /* The entry specifying the file with the gain matrix is only required if the appropiate checkbox is
      * checked */
     if( ReadGain_checkBox->isChecked() ){
	  str = "Gain_Matrix = ";
	  str.append( "\"" );
	  str.append( Gain_lineEdit->text() );
	  str.append( "\"" );
	  CFPreview_textBrowser->append( str );
     }
     /* The entry specifying the file with the load vector form is only required if the appropiate checkbox is
      * checked */
     if( ReadLVector_checkBox->isChecked() ){
	  str = "Load_Vector = ";
	  str.append( "\"" );
	  str.append( LVector_lineEdit->text() );
	  str.append( "\"" );
	  CFPreview_textBrowser->append( str );
     }

     /* The file specifying the coupling nodes is only required if there are substructures */
     if( NumSub_SpinBox->value() != 0 ){
	  str = "Coupling_Nodes = ";
	  str.append( "\"" );
	  str.append( CNodes_lineEdit->text() );
	  str.append( "\"" );
	  CFPreview_textBrowser->append( str );
     }
     str = "Ground_Motion = ";
     str.append( "\"" );
     str.append( GMotion_lineEdit->text() );
     str.append( "\"" );
     CFPreview_textBrowser->append( str );
     str = "OutputFile = ";
     str.append( "\"" );
     str.append( OutFile_lineEdit->text() );
     str.append( "\"" );
     CFPreview_textBrowser->append( str );
}

void ConfFile::On_ConfButton_Clicked( )
{
    
     ConfFile_lineEdit->setText( QFileDialog::getOpenFileName(this, tr("Configuration file"), "ConfFile.conf", tr("Configuration file (*)")) );

}

void ConfFile::On_MassButton_Clicked( )
{
    
     Mass_lineEdit->setText( QFileDialog::getOpenFileName(this, tr("Mass matrix"), "504M_Sp_MM.txt", tr("Mass matrix file (*)")) );

}

void ConfFile::On_StiffButton_Clicked( )
{
    
     Stiff_lineEdit->setText( QFileDialog::getOpenFileName(this, tr("Stiffness matrix"), "504K_Sp_MM.txt", tr("Stiffness matrix file (*)")) );

}

void ConfFile::On_DampButton_Clicked( )
{
    
     Damp_lineEdit->setText( QFileDialog::getOpenFileName(this, tr("Damping matrix"), "504C_Sp_MM.txt", tr("Damping matrix file (*)")) );

}

void ConfFile::On_GainButton_Clicked( )
{
    
     Gain_lineEdit->setText( QFileDialog::getOpenFileName(this, tr("Gain matrix"), "504G_MM.txt", tr("Gain matrix file (*)")) );

}

void ConfFile::On_LVectorButton_Clicked( )
{
    
     LVector_lineEdit->setText( QFileDialog::getOpenFileName(this, tr("Load vector"), "504LV_MM.txt", tr("Load vector file (*)")) );

}

void ConfFile::On_CNodesButton_Clicked( )
{
    
     CNodes_lineEdit->setText( QFileDialog::getOpenFileName(this, tr("Coupling nodes"), "Couple_Nodes.txt", tr("Coupling nodes file (*)")) );

}

void ConfFile::On_GMotionButton_Clicked( )
{
    
     GMotion_lineEdit->setText( QFileDialog::getOpenFileName(this, tr("Ground motion"), "GroundMovement.txt", tr("Ground motion file (*)")) );

}

void ConfFile::On_OutFileButton_Clicked( )
{
    
     OutFile_lineEdit->setText( QFileDialog::getOpenFileName(this, tr("Mass matrix"), "Out-Sparse.h5", tr("Output File (*h5 *txt)")) );

}

void ConfFile::On_GenCFileButton_Clicked( )
{
     QFile OutFile( ConfFile_lineEdit->text()) ;

     if( !OutFile.open( QIODevice::WriteOnly | QIODevice::Text ) ){
	  
     }

     QTextStream out( &OutFile );
     out << CFPreview_textBrowser->document()->toPlainText();
     OutFile.close();
}

void ConfFile::On_Run_pushButton_Clicked()
{

     SAlgorithm = new QProcess();

     SAlgorithm->setProcessChannelMode( QProcess::MergedChannels );

     SAlgorithm->start( Exec_lineEdit->text() );
     connect( SAlgorithm, SIGNAL(readyReadStandardOutput()), this, SLOT(readStdOutput()));
     SAlgorithm->waitForFinished();
}

void ConfFile::readStdOutput()
{

     SAlg_textBrowser->append(SAlgorithm->readAllStandardOutput());
}     


void ConfFile::On_OpenOutput_pushButton_Clicked()
{

     QProcess View;

     View.startDetached( "hdfview", QStringList() << OutFile_lineEdit->text() );

}
