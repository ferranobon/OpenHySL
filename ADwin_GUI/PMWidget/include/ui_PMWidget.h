/********************************************************************************
** Form generated from reading UI file 'PMWidget.ui'
**
** Created: Tue Mar 27 10:47:12 2012
**      by: Qt User Interface Compiler version 4.7.4
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_PMWIDGET_H
#define UI_PMWIDGET_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QDoubleSpinBox>
#include <QtGui/QGridLayout>
#include <QtGui/QGroupBox>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QPushButton>
#include <QtGui/QSpacerItem>
#include <QtGui/QTextBrowser>
#include <QtGui/QToolButton>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_PMWidget
{
public:
    QGridLayout *gridLayout_5;
    QGroupBox *groupBox_5;
    QGridLayout *gridLayout_8;
    QGroupBox *groupBox;
    QVBoxLayout *verticalLayout;
    QHBoxLayout *horizontalLayout_3;
    QLabel *label_18;
    QLineEdit *DeviceNum_LineEdit;
    QPushButton *BootButton;
    QSpacerItem *verticalSpacer;
    QPushButton *LoadCtrlPButton;
    QPushButton *StartCntrlPButton;
    QPushButton *ADwinReadyButton;
    QSpacerItem *verticalSpacer_2;
    QPushButton *LoadTDataButton;
    QPushButton *StopCntrlPButton;
    QSpacerItem *verticalSpacer_3;
    QGridLayout *gridLayout_3;
    QGroupBox *groupBox_2;
    QGridLayout *gridLayout;
    QPushButton *LowPressureButton;
    QPushButton *HighPressureButton;
    QPushButton *StopPressureButton;
    QGroupBox *groupBox_3;
    QVBoxLayout *verticalLayout_3;
    QPushButton *LoadTestPButton;
    QPushButton *StartTestButton;
    QPushButton *StopTestButton;
    QPushButton *ReadDataButton;
    QPushButton *EmergencyButton;
    QGroupBox *groupBox_4;
    QGridLayout *gridLayout_2;
    QLineEdit *PFile_LineEdit;
    QToolButton *PFiletoolButton;
    QLabel *label_2;
    QLineEdit *TFile_LineEdit;
    QToolButton *TFiletoolButton;
    QLabel *label_3;
    QLineEdit *IFile_LineEdit;
    QToolButton *IFiletoolButton;
    QLabel *label_4;
    QLineEdit *OFile_LineEdit;
    QToolButton *OFiletoolButton;
    QLabel *label;
    QLabel *label_19;
    QLineEdit *BootFile_LineEdit;
    QToolButton *BFiletoolButton;
    QGroupBox *groupBox_6;
    QVBoxLayout *verticalLayout_2;
    QPushButton *DefValPButton;
    QHBoxLayout *horizontalLayout_2;
    QLabel *label_9;
    QDoubleSpinBox *P_DSpinBox;
    QSpacerItem *horizontalSpacer;
    QLabel *label_10;
    QDoubleSpinBox *I_DSpinBox;
    QSpacerItem *horizontalSpacer_2;
    QLabel *label_11;
    QDoubleSpinBox *D_DSpinBox;
    QGroupBox *groupBox_8;
    QGridLayout *gridLayout_4;
    QLabel *label_5;
    QDoubleSpinBox *DCtrl1_DSpinBox;
    QLabel *label_7;
    QLineEdit *DMeas1_lineEdit;
    QLabel *label_6;
    QLabel *label_8;
    QLineEdit *DMeas2_lineEdit;
    QDoubleSpinBox *DCtrl2_DSpinBox;
    QGroupBox *groupBox_9;
    QGridLayout *gridLayout_7;
    QPushButton *SetLimPButton;
    QCheckBox *EnabLimCheckBox;
    QLabel *label_14;
    QLineEdit *DMin1_lineEdit;
    QLabel *label_17;
    QLineEdit *DMax1_lineEdit;
    QLabel *label_15;
    QLineEdit *DMin2_lineEdit;
    QLabel *label_16;
    QLineEdit *DMax2_lineEdit;
    QGroupBox *groupBox_7;
    QGridLayout *gridLayout_6;
    QTextBrowser *LogBrowser;

    void setupUi(QWidget *PMWidget)
    {
        if (PMWidget->objectName().isEmpty())
            PMWidget->setObjectName(QString::fromUtf8("PMWidget"));
        PMWidget->resize(848, 536);
        gridLayout_5 = new QGridLayout(PMWidget);
        gridLayout_5->setObjectName(QString::fromUtf8("gridLayout_5"));
        groupBox_5 = new QGroupBox(PMWidget);
        groupBox_5->setObjectName(QString::fromUtf8("groupBox_5"));
        gridLayout_8 = new QGridLayout(groupBox_5);
        gridLayout_8->setObjectName(QString::fromUtf8("gridLayout_8"));
        groupBox = new QGroupBox(groupBox_5);
        groupBox->setObjectName(QString::fromUtf8("groupBox"));
        verticalLayout = new QVBoxLayout(groupBox);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        horizontalLayout_3 = new QHBoxLayout();
        horizontalLayout_3->setObjectName(QString::fromUtf8("horizontalLayout_3"));
        label_18 = new QLabel(groupBox);
        label_18->setObjectName(QString::fromUtf8("label_18"));

        horizontalLayout_3->addWidget(label_18);

        DeviceNum_LineEdit = new QLineEdit(groupBox);
        DeviceNum_LineEdit->setObjectName(QString::fromUtf8("DeviceNum_LineEdit"));
        DeviceNum_LineEdit->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        horizontalLayout_3->addWidget(DeviceNum_LineEdit);


        verticalLayout->addLayout(horizontalLayout_3);

        BootButton = new QPushButton(groupBox);
        BootButton->setObjectName(QString::fromUtf8("BootButton"));
        BootButton->setCheckable(false);

        verticalLayout->addWidget(BootButton);

        verticalSpacer = new QSpacerItem(20, 18, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer);

        LoadCtrlPButton = new QPushButton(groupBox);
        LoadCtrlPButton->setObjectName(QString::fromUtf8("LoadCtrlPButton"));
        LoadCtrlPButton->setEnabled(false);

        verticalLayout->addWidget(LoadCtrlPButton);

        StartCntrlPButton = new QPushButton(groupBox);
        StartCntrlPButton->setObjectName(QString::fromUtf8("StartCntrlPButton"));
        StartCntrlPButton->setEnabled(false);

        verticalLayout->addWidget(StartCntrlPButton);

        ADwinReadyButton = new QPushButton(groupBox);
        ADwinReadyButton->setObjectName(QString::fromUtf8("ADwinReadyButton"));
        ADwinReadyButton->setEnabled(false);

        verticalLayout->addWidget(ADwinReadyButton);

        verticalSpacer_2 = new QSpacerItem(20, 18, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer_2);

        LoadTDataButton = new QPushButton(groupBox);
        LoadTDataButton->setObjectName(QString::fromUtf8("LoadTDataButton"));
        LoadTDataButton->setEnabled(false);

        verticalLayout->addWidget(LoadTDataButton);

        StopCntrlPButton = new QPushButton(groupBox);
        StopCntrlPButton->setObjectName(QString::fromUtf8("StopCntrlPButton"));
        StopCntrlPButton->setEnabled(false);

        verticalLayout->addWidget(StopCntrlPButton);

        verticalSpacer_3 = new QSpacerItem(20, 18, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer_3);


        gridLayout_8->addWidget(groupBox, 0, 0, 1, 1);

        gridLayout_3 = new QGridLayout();
        gridLayout_3->setObjectName(QString::fromUtf8("gridLayout_3"));
        groupBox_2 = new QGroupBox(groupBox_5);
        groupBox_2->setObjectName(QString::fromUtf8("groupBox_2"));
        gridLayout = new QGridLayout(groupBox_2);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        LowPressureButton = new QPushButton(groupBox_2);
        LowPressureButton->setObjectName(QString::fromUtf8("LowPressureButton"));
        LowPressureButton->setEnabled(false);
        LowPressureButton->setMinimumSize(QSize(0, 31));
        LowPressureButton->setCheckable(true);

        gridLayout->addWidget(LowPressureButton, 0, 0, 1, 1);

        HighPressureButton = new QPushButton(groupBox_2);
        HighPressureButton->setObjectName(QString::fromUtf8("HighPressureButton"));
        HighPressureButton->setEnabled(false);
        HighPressureButton->setMinimumSize(QSize(0, 31));
        HighPressureButton->setCheckable(true);

        gridLayout->addWidget(HighPressureButton, 0, 1, 1, 1);

        StopPressureButton = new QPushButton(groupBox_2);
        StopPressureButton->setObjectName(QString::fromUtf8("StopPressureButton"));
        StopPressureButton->setEnabled(false);
        StopPressureButton->setMinimumSize(QSize(0, 71));

        gridLayout->addWidget(StopPressureButton, 1, 0, 1, 2);


        gridLayout_3->addWidget(groupBox_2, 0, 0, 1, 1);

        groupBox_3 = new QGroupBox(groupBox_5);
        groupBox_3->setObjectName(QString::fromUtf8("groupBox_3"));
        verticalLayout_3 = new QVBoxLayout(groupBox_3);
        verticalLayout_3->setObjectName(QString::fromUtf8("verticalLayout_3"));
        LoadTestPButton = new QPushButton(groupBox_3);
        LoadTestPButton->setObjectName(QString::fromUtf8("LoadTestPButton"));
        LoadTestPButton->setEnabled(false);

        verticalLayout_3->addWidget(LoadTestPButton);

        StartTestButton = new QPushButton(groupBox_3);
        StartTestButton->setObjectName(QString::fromUtf8("StartTestButton"));
        StartTestButton->setEnabled(false);

        verticalLayout_3->addWidget(StartTestButton);

        StopTestButton = new QPushButton(groupBox_3);
        StopTestButton->setObjectName(QString::fromUtf8("StopTestButton"));
        StopTestButton->setEnabled(false);

        verticalLayout_3->addWidget(StopTestButton);

        ReadDataButton = new QPushButton(groupBox_3);
        ReadDataButton->setObjectName(QString::fromUtf8("ReadDataButton"));
        ReadDataButton->setEnabled(false);

        verticalLayout_3->addWidget(ReadDataButton);


        gridLayout_3->addWidget(groupBox_3, 0, 1, 1, 1);

        EmergencyButton = new QPushButton(groupBox_5);
        EmergencyButton->setObjectName(QString::fromUtf8("EmergencyButton"));
        EmergencyButton->setEnabled(false);
        EmergencyButton->setMinimumSize(QSize(0, 181));

        gridLayout_3->addWidget(EmergencyButton, 1, 0, 1, 2);


        gridLayout_8->addLayout(gridLayout_3, 0, 1, 1, 1);

        groupBox_4 = new QGroupBox(groupBox_5);
        groupBox_4->setObjectName(QString::fromUtf8("groupBox_4"));
        gridLayout_2 = new QGridLayout(groupBox_4);
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        PFile_LineEdit = new QLineEdit(groupBox_4);
        PFile_LineEdit->setObjectName(QString::fromUtf8("PFile_LineEdit"));

        gridLayout_2->addWidget(PFile_LineEdit, 1, 1, 1, 1);

        PFiletoolButton = new QToolButton(groupBox_4);
        PFiletoolButton->setObjectName(QString::fromUtf8("PFiletoolButton"));

        gridLayout_2->addWidget(PFiletoolButton, 1, 2, 1, 1);

        label_2 = new QLabel(groupBox_4);
        label_2->setObjectName(QString::fromUtf8("label_2"));

        gridLayout_2->addWidget(label_2, 2, 0, 1, 1);

        TFile_LineEdit = new QLineEdit(groupBox_4);
        TFile_LineEdit->setObjectName(QString::fromUtf8("TFile_LineEdit"));

        gridLayout_2->addWidget(TFile_LineEdit, 2, 1, 1, 1);

        TFiletoolButton = new QToolButton(groupBox_4);
        TFiletoolButton->setObjectName(QString::fromUtf8("TFiletoolButton"));

        gridLayout_2->addWidget(TFiletoolButton, 2, 2, 1, 1);

        label_3 = new QLabel(groupBox_4);
        label_3->setObjectName(QString::fromUtf8("label_3"));

        gridLayout_2->addWidget(label_3, 3, 0, 1, 1);

        IFile_LineEdit = new QLineEdit(groupBox_4);
        IFile_LineEdit->setObjectName(QString::fromUtf8("IFile_LineEdit"));

        gridLayout_2->addWidget(IFile_LineEdit, 3, 1, 1, 1);

        IFiletoolButton = new QToolButton(groupBox_4);
        IFiletoolButton->setObjectName(QString::fromUtf8("IFiletoolButton"));

        gridLayout_2->addWidget(IFiletoolButton, 3, 2, 1, 1);

        label_4 = new QLabel(groupBox_4);
        label_4->setObjectName(QString::fromUtf8("label_4"));

        gridLayout_2->addWidget(label_4, 4, 0, 1, 1);

        OFile_LineEdit = new QLineEdit(groupBox_4);
        OFile_LineEdit->setObjectName(QString::fromUtf8("OFile_LineEdit"));

        gridLayout_2->addWidget(OFile_LineEdit, 4, 1, 1, 1);

        OFiletoolButton = new QToolButton(groupBox_4);
        OFiletoolButton->setObjectName(QString::fromUtf8("OFiletoolButton"));

        gridLayout_2->addWidget(OFiletoolButton, 4, 2, 1, 1);

        label = new QLabel(groupBox_4);
        label->setObjectName(QString::fromUtf8("label"));

        gridLayout_2->addWidget(label, 1, 0, 1, 1);

        label_19 = new QLabel(groupBox_4);
        label_19->setObjectName(QString::fromUtf8("label_19"));

        gridLayout_2->addWidget(label_19, 0, 0, 1, 1);

        BootFile_LineEdit = new QLineEdit(groupBox_4);
        BootFile_LineEdit->setObjectName(QString::fromUtf8("BootFile_LineEdit"));

        gridLayout_2->addWidget(BootFile_LineEdit, 0, 1, 1, 1);

        BFiletoolButton = new QToolButton(groupBox_4);
        BFiletoolButton->setObjectName(QString::fromUtf8("BFiletoolButton"));

        gridLayout_2->addWidget(BFiletoolButton, 0, 2, 1, 1);


        gridLayout_8->addWidget(groupBox_4, 1, 0, 1, 2);


        gridLayout_5->addWidget(groupBox_5, 0, 0, 2, 1);

        groupBox_6 = new QGroupBox(PMWidget);
        groupBox_6->setObjectName(QString::fromUtf8("groupBox_6"));
        verticalLayout_2 = new QVBoxLayout(groupBox_6);
        verticalLayout_2->setObjectName(QString::fromUtf8("verticalLayout_2"));
        DefValPButton = new QPushButton(groupBox_6);
        DefValPButton->setObjectName(QString::fromUtf8("DefValPButton"));
        DefValPButton->setEnabled(false);

        verticalLayout_2->addWidget(DefValPButton);

        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        label_9 = new QLabel(groupBox_6);
        label_9->setObjectName(QString::fromUtf8("label_9"));

        horizontalLayout_2->addWidget(label_9);

        P_DSpinBox = new QDoubleSpinBox(groupBox_6);
        P_DSpinBox->setObjectName(QString::fromUtf8("P_DSpinBox"));
        P_DSpinBox->setEnabled(false);

        horizontalLayout_2->addWidget(P_DSpinBox);

        horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_2->addItem(horizontalSpacer);

        label_10 = new QLabel(groupBox_6);
        label_10->setObjectName(QString::fromUtf8("label_10"));

        horizontalLayout_2->addWidget(label_10);

        I_DSpinBox = new QDoubleSpinBox(groupBox_6);
        I_DSpinBox->setObjectName(QString::fromUtf8("I_DSpinBox"));
        I_DSpinBox->setEnabled(false);

        horizontalLayout_2->addWidget(I_DSpinBox);

        horizontalSpacer_2 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_2->addItem(horizontalSpacer_2);

        label_11 = new QLabel(groupBox_6);
        label_11->setObjectName(QString::fromUtf8("label_11"));

        horizontalLayout_2->addWidget(label_11);

        D_DSpinBox = new QDoubleSpinBox(groupBox_6);
        D_DSpinBox->setObjectName(QString::fromUtf8("D_DSpinBox"));
        D_DSpinBox->setEnabled(false);

        horizontalLayout_2->addWidget(D_DSpinBox);


        verticalLayout_2->addLayout(horizontalLayout_2);

        groupBox_8 = new QGroupBox(groupBox_6);
        groupBox_8->setObjectName(QString::fromUtf8("groupBox_8"));
        gridLayout_4 = new QGridLayout(groupBox_8);
        gridLayout_4->setObjectName(QString::fromUtf8("gridLayout_4"));
        label_5 = new QLabel(groupBox_8);
        label_5->setObjectName(QString::fromUtf8("label_5"));

        gridLayout_4->addWidget(label_5, 0, 0, 1, 1);

        DCtrl1_DSpinBox = new QDoubleSpinBox(groupBox_8);
        DCtrl1_DSpinBox->setObjectName(QString::fromUtf8("DCtrl1_DSpinBox"));
        DCtrl1_DSpinBox->setEnabled(false);

        gridLayout_4->addWidget(DCtrl1_DSpinBox, 0, 1, 1, 1);

        label_7 = new QLabel(groupBox_8);
        label_7->setObjectName(QString::fromUtf8("label_7"));

        gridLayout_4->addWidget(label_7, 0, 2, 1, 1);

        DMeas1_lineEdit = new QLineEdit(groupBox_8);
        DMeas1_lineEdit->setObjectName(QString::fromUtf8("DMeas1_lineEdit"));
        DMeas1_lineEdit->setReadOnly(true);

        gridLayout_4->addWidget(DMeas1_lineEdit, 0, 3, 1, 1);

        label_6 = new QLabel(groupBox_8);
        label_6->setObjectName(QString::fromUtf8("label_6"));

        gridLayout_4->addWidget(label_6, 1, 0, 1, 1);

        label_8 = new QLabel(groupBox_8);
        label_8->setObjectName(QString::fromUtf8("label_8"));

        gridLayout_4->addWidget(label_8, 1, 2, 1, 1);

        DMeas2_lineEdit = new QLineEdit(groupBox_8);
        DMeas2_lineEdit->setObjectName(QString::fromUtf8("DMeas2_lineEdit"));
        DMeas2_lineEdit->setReadOnly(true);

        gridLayout_4->addWidget(DMeas2_lineEdit, 1, 3, 1, 1);

        DCtrl2_DSpinBox = new QDoubleSpinBox(groupBox_8);
        DCtrl2_DSpinBox->setObjectName(QString::fromUtf8("DCtrl2_DSpinBox"));
        DCtrl2_DSpinBox->setEnabled(false);

        gridLayout_4->addWidget(DCtrl2_DSpinBox, 1, 1, 1, 1);


        verticalLayout_2->addWidget(groupBox_8);

        groupBox_9 = new QGroupBox(groupBox_6);
        groupBox_9->setObjectName(QString::fromUtf8("groupBox_9"));
        gridLayout_7 = new QGridLayout(groupBox_9);
        gridLayout_7->setObjectName(QString::fromUtf8("gridLayout_7"));
        SetLimPButton = new QPushButton(groupBox_9);
        SetLimPButton->setObjectName(QString::fromUtf8("SetLimPButton"));
        SetLimPButton->setEnabled(false);

        gridLayout_7->addWidget(SetLimPButton, 0, 0, 1, 2);

        EnabLimCheckBox = new QCheckBox(groupBox_9);
        EnabLimCheckBox->setObjectName(QString::fromUtf8("EnabLimCheckBox"));
        EnabLimCheckBox->setEnabled(false);
        EnabLimCheckBox->setChecked(true);

        gridLayout_7->addWidget(EnabLimCheckBox, 0, 2, 1, 2);

        label_14 = new QLabel(groupBox_9);
        label_14->setObjectName(QString::fromUtf8("label_14"));

        gridLayout_7->addWidget(label_14, 1, 0, 1, 1);

        DMin1_lineEdit = new QLineEdit(groupBox_9);
        DMin1_lineEdit->setObjectName(QString::fromUtf8("DMin1_lineEdit"));

        gridLayout_7->addWidget(DMin1_lineEdit, 1, 1, 1, 1);

        label_17 = new QLabel(groupBox_9);
        label_17->setObjectName(QString::fromUtf8("label_17"));

        gridLayout_7->addWidget(label_17, 1, 2, 1, 1);

        DMax1_lineEdit = new QLineEdit(groupBox_9);
        DMax1_lineEdit->setObjectName(QString::fromUtf8("DMax1_lineEdit"));

        gridLayout_7->addWidget(DMax1_lineEdit, 1, 3, 1, 1);

        label_15 = new QLabel(groupBox_9);
        label_15->setObjectName(QString::fromUtf8("label_15"));

        gridLayout_7->addWidget(label_15, 2, 0, 1, 1);

        DMin2_lineEdit = new QLineEdit(groupBox_9);
        DMin2_lineEdit->setObjectName(QString::fromUtf8("DMin2_lineEdit"));

        gridLayout_7->addWidget(DMin2_lineEdit, 2, 1, 1, 1);

        label_16 = new QLabel(groupBox_9);
        label_16->setObjectName(QString::fromUtf8("label_16"));

        gridLayout_7->addWidget(label_16, 2, 2, 1, 1);

        DMax2_lineEdit = new QLineEdit(groupBox_9);
        DMax2_lineEdit->setObjectName(QString::fromUtf8("DMax2_lineEdit"));

        gridLayout_7->addWidget(DMax2_lineEdit, 2, 3, 1, 1);


        verticalLayout_2->addWidget(groupBox_9);

        groupBox_8->raise();
        groupBox_9->raise();

        gridLayout_5->addWidget(groupBox_6, 0, 1, 1, 1);

        groupBox_7 = new QGroupBox(PMWidget);
        groupBox_7->setObjectName(QString::fromUtf8("groupBox_7"));
        gridLayout_6 = new QGridLayout(groupBox_7);
        gridLayout_6->setObjectName(QString::fromUtf8("gridLayout_6"));
        LogBrowser = new QTextBrowser(groupBox_7);
        LogBrowser->setObjectName(QString::fromUtf8("LogBrowser"));

        gridLayout_6->addWidget(LogBrowser, 0, 0, 1, 1);


        gridLayout_5->addWidget(groupBox_7, 1, 1, 1, 1);


        retranslateUi(PMWidget);

        QMetaObject::connectSlotsByName(PMWidget);
    } // setupUi

    void retranslateUi(QWidget *PMWidget)
    {
        PMWidget->setWindowTitle(QApplication::translate("PMWidget", "Process and Pressure Management", 0, QApplication::UnicodeUTF8));
        groupBox_5->setTitle(QApplication::translate("PMWidget", "ADwin and Process Management", 0, QApplication::UnicodeUTF8));
        groupBox->setTitle(QApplication::translate("PMWidget", "ADwin and Control Process", 0, QApplication::UnicodeUTF8));
        label_18->setText(QApplication::translate("PMWidget", "Device Number:", 0, QApplication::UnicodeUTF8));
        DeviceNum_LineEdit->setText(QApplication::translate("PMWidget", "336", 0, QApplication::UnicodeUTF8));
        BootButton->setText(QApplication::translate("PMWidget", "Boot ADwin", 0, QApplication::UnicodeUTF8));
        LoadCtrlPButton->setText(QApplication::translate("PMWidget", "Load Control Process", 0, QApplication::UnicodeUTF8));
        StartCntrlPButton->setText(QApplication::translate("PMWidget", "Start Control Process", 0, QApplication::UnicodeUTF8));
        ADwinReadyButton->setText(QApplication::translate("PMWidget", "ADwin Ready", 0, QApplication::UnicodeUTF8));
        LoadTDataButton->setText(QApplication::translate("PMWidget", "Load Test Data", 0, QApplication::UnicodeUTF8));
        StopCntrlPButton->setText(QApplication::translate("PMWidget", "Stop Control Process", 0, QApplication::UnicodeUTF8));
        groupBox_2->setTitle(QApplication::translate("PMWidget", "Pressure control", 0, QApplication::UnicodeUTF8));
        LowPressureButton->setText(QApplication::translate("PMWidget", "Low Pressure", 0, QApplication::UnicodeUTF8));
        HighPressureButton->setText(QApplication::translate("PMWidget", "High Pressure", 0, QApplication::UnicodeUTF8));
        StopPressureButton->setText(QApplication::translate("PMWidget", "Stop Pressure", 0, QApplication::UnicodeUTF8));
        groupBox_3->setTitle(QApplication::translate("PMWidget", "Test Process", 0, QApplication::UnicodeUTF8));
        LoadTestPButton->setText(QApplication::translate("PMWidget", "Load Test Process", 0, QApplication::UnicodeUTF8));
        StartTestButton->setText(QApplication::translate("PMWidget", "Start Test Process", 0, QApplication::UnicodeUTF8));
        StopTestButton->setText(QApplication::translate("PMWidget", "Sto&p Test Process", 0, QApplication::UnicodeUTF8));
        ReadDataButton->setText(QApplication::translate("PMWidget", "Read and &Save Data", 0, QApplication::UnicodeUTF8));
        EmergencyButton->setText(QApplication::translate("PMWidget", "EMERGENCY STOP", 0, QApplication::UnicodeUTF8));
        groupBox_4->setTitle(QApplication::translate("PMWidget", "Files", 0, QApplication::UnicodeUTF8));
        PFile_LineEdit->setText(QApplication::translate("PMWidget", "/home/ferran/workspace/ADwinGUI/PMWidget/ADbasic2_Pr1.TB1", 0, QApplication::UnicodeUTF8));
        PFiletoolButton->setText(QApplication::translate("PMWidget", "...", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("PMWidget", "Test Process File:", 0, QApplication::UnicodeUTF8));
        TFile_LineEdit->setText(QApplication::translate("PMWidget", "/home/ferran/workspace/ADwinGUI/PMWidget/ADbasic1_Pr2.TB2", 0, QApplication::UnicodeUTF8));
        TFiletoolButton->setText(QApplication::translate("PMWidget", "...", 0, QApplication::UnicodeUTF8));
        label_3->setText(QApplication::translate("PMWidget", "Input Data File:", 0, QApplication::UnicodeUTF8));
        IFiletoolButton->setText(QApplication::translate("PMWidget", "...", 0, QApplication::UnicodeUTF8));
        label_4->setText(QApplication::translate("PMWidget", "Output Data File:", 0, QApplication::UnicodeUTF8));
        OFiletoolButton->setText(QApplication::translate("PMWidget", "...", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("PMWidget", "Control Process File:", 0, QApplication::UnicodeUTF8));
        label_19->setText(QApplication::translate("PMWidget", "ADwin Boot File:", 0, QApplication::UnicodeUTF8));
        BootFile_LineEdit->setText(QApplication::translate("PMWidget", "/opt/adwin/share/btl/ADwin11.btl", 0, QApplication::UnicodeUTF8));
        BFiletoolButton->setText(QApplication::translate("PMWidget", "...", 0, QApplication::UnicodeUTF8));
        groupBox_6->setTitle(QApplication::translate("PMWidget", "Displacement Control", 0, QApplication::UnicodeUTF8));
        DefValPButton->setText(QApplication::translate("PMWidget", "Reset to default values", 0, QApplication::UnicodeUTF8));
        label_9->setText(QApplication::translate("PMWidget", "P:", 0, QApplication::UnicodeUTF8));
        label_10->setText(QApplication::translate("PMWidget", "I:", 0, QApplication::UnicodeUTF8));
        label_11->setText(QApplication::translate("PMWidget", "D:", 0, QApplication::UnicodeUTF8));
        groupBox_8->setTitle(QApplication::translate("PMWidget", "Current Data", 0, QApplication::UnicodeUTF8));
        label_5->setText(QApplication::translate("PMWidget", "DCtrl1", 0, QApplication::UnicodeUTF8));
        label_7->setText(QApplication::translate("PMWidget", "DMea1", 0, QApplication::UnicodeUTF8));
        label_6->setText(QApplication::translate("PMWidget", "DCtrl2", 0, QApplication::UnicodeUTF8));
        label_8->setText(QApplication::translate("PMWidget", "DMea2", 0, QApplication::UnicodeUTF8));
        groupBox_9->setTitle(QApplication::translate("PMWidget", "Limitation", 0, QApplication::UnicodeUTF8));
        SetLimPButton->setText(QApplication::translate("PMWidget", "Set Limitation", 0, QApplication::UnicodeUTF8));
        EnabLimCheckBox->setText(QApplication::translate("PMWidget", "Enable Limits", 0, QApplication::UnicodeUTF8));
        label_14->setText(QApplication::translate("PMWidget", "Dmin1:", 0, QApplication::UnicodeUTF8));
        label_17->setText(QApplication::translate("PMWidget", "Dmax1:", 0, QApplication::UnicodeUTF8));
        label_15->setText(QApplication::translate("PMWidget", "Dmin2:", 0, QApplication::UnicodeUTF8));
        label_16->setText(QApplication::translate("PMWidget", "Dmax2:", 0, QApplication::UnicodeUTF8));
        groupBox_7->setTitle(QApplication::translate("PMWidget", "Log", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class PMWidget: public Ui_PMWidget {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_PMWIDGET_H
