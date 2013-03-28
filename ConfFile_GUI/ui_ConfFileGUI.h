/********************************************************************************
** Form generated from reading UI file 'ConfFileGUI.ui'
**
** Created: Thu Mar 28 11:56:23 2013
**      by: Qt User Interface Compiler version 4.8.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_CONFFILEGUI_H
#define UI_CONFFILEGUI_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QComboBox>
#include <QtGui/QDoubleSpinBox>
#include <QtGui/QGridLayout>
#include <QtGui/QGroupBox>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QListWidget>
#include <QtGui/QMainWindow>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QPushButton>
#include <QtGui/QSpacerItem>
#include <QtGui/QSpinBox>
#include <QtGui/QStatusBar>
#include <QtGui/QTabWidget>
#include <QtGui/QTableView>
#include <QtGui/QTextBrowser>
#include <QtGui/QToolButton>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_ConfFile
{
public:
    QAction *actionSave;
    QAction *actionOpen;
    QAction *actionSave_As;
    QAction *actionExit;
    QWidget *centralwidget;
    QGridLayout *gridLayout_11;
    QTabWidget *ConfFile_tabWidget;
    QWidget *ConfFile_tab;
    QGridLayout *gridLayout_12;
    QGroupBox *groupBox_4;
    QGridLayout *gridLayout_4;
    QComboBox *IntCnt_ComboBox;
    QLabel *label_10;
    QDoubleSpinBox *Beta_dSpinBox;
    QLabel *label_11;
    QDoubleSpinBox *Gamma_dSpinBox;
    QGroupBox *groupBox_3;
    QGridLayout *gridLayout_3;
    QLabel *label_6;
    QDoubleSpinBox *P_dSpinBox;
    QLabel *label_7;
    QDoubleSpinBox *I_dSpinBox;
    QLabel *label_8;
    QDoubleSpinBox *D_dSpinBox;
    QVBoxLayout *verticalLayout_3;
    QTextBrowser *CFPreview_textBrowser;
    QPushButton *GenCFile_pushButton;
    QGroupBox *groupBox_5;
    QGridLayout *gridLayout;
    QLabel *label_14;
    QSpinBox *NumSub_SpinBox;
    QLabel *label_9;
    QSpinBox *NumSubSteps_SpinBox;
    QGroupBox *groupBox_2;
    QGridLayout *gridLayout_2;
    QLabel *label_4;
    QDoubleSpinBox *RayAlpha_dSpinBox;
    QLabel *label_5;
    QDoubleSpinBox *RayBeta_dSpinBox;
    QGroupBox *groupBox_6;
    QHBoxLayout *horizontalLayout;
    QGridLayout *gridLayout_6;
    QLabel *ConfFile_label;
    QLineEdit *ConfFile_lineEdit;
    QToolButton *ConfFile_toolButton;
    QLabel *Mass_label;
    QLineEdit *Mass_lineEdit;
    QToolButton *Mass_toolButton;
    QLabel *Stiffness_label;
    QLineEdit *Stiff_lineEdit;
    QToolButton *Stiff_toolButton;
    QLabel *Damp_label;
    QLineEdit *Damp_lineEdit;
    QToolButton *Damp_toolButton;
    QLabel *Gain_label;
    QLineEdit *Gain_lineEdit;
    QToolButton *Gain_toolButton;
    QLabel *LVector_label;
    QLineEdit *LVector_lineEdit;
    QToolButton *LVector_toolButton;
    QLabel *LVectorForm_label;
    QLineEdit *LVectorForm_lineEdit;
    QLabel *CNodes_label;
    QLineEdit *CNodes_lineEdit;
    QToolButton *CNodes_toolButton;
    QLabel *GMotion_label;
    QLineEdit *GMotion_lineEdit;
    QToolButton *GMotion_toolButton;
    QLabel *OutFile_label;
    QLineEdit *OutFile_lineEdit;
    QToolButton *OutFile_toolButton;
    QGroupBox *groupBox;
    QVBoxLayout *verticalLayout_5;
    QGridLayout *gridLayout_5;
    QCheckBox *AbsValues_checkBox;
    QCheckBox *ReadGain_checkBox;
    QCheckBox *ReadSparse_checkBox;
    QLabel *label;
    QSpinBox *Order_SpinBox;
    QCheckBox *UseSparse_checkBox;
    QLabel *label_3;
    QDoubleSpinBox *TimeStep_dSpinBox;
    QCheckBox *UsePacked_checkBox;
    QLabel *label_2;
    QSpinBox *NumSteps_SpinBox;
    QVBoxLayout *verticalLayout_4;
    QCheckBox *ReadLVector_checkBox;
    QCheckBox *ReadDamp_checkBox;
    QWidget *CNodes_tab;
    QGridLayout *gridLayout_8;
    QLabel *label_23;
    QSpacerItem *horizontalSpacer;
    QLabel *label_24;
    QSpacerItem *horizontalSpacer_6;
    QListWidget *listWidget;
    QSpacerItem *verticalSpacer_2;
    QTableView *tableView;
    QPushButton *pushButton;
    QPushButton *pushButton_2;
    QSpacerItem *verticalSpacer_4;
    QWidget *Algorithm_tab;
    QGridLayout *gridLayout_7;
    QHBoxLayout *horizontalLayout_2;
    QLabel *Exec_label;
    QLineEdit *Exec_lineEdit;
    QToolButton *Exec_toolButton;
    QVBoxLayout *verticalLayout;
    QPushButton *Run_pushButton;
    QPushButton *OpenOutput_pushButton;
    QSpacerItem *verticalSpacer;
    QPushButton *pushButton_4;
    QTextBrowser *SAlg_textBrowser;
    QMenuBar *menubar;
    QMenu *menuFile;
    QStatusBar *statusbar;

    void setupUi(QMainWindow *ConfFile)
    {
        if (ConfFile->objectName().isEmpty())
            ConfFile->setObjectName(QString::fromUtf8("ConfFile"));
        ConfFile->resize(962, 593);
        QSizePolicy sizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(ConfFile->sizePolicy().hasHeightForWidth());
        ConfFile->setSizePolicy(sizePolicy);
        actionSave = new QAction(ConfFile);
        actionSave->setObjectName(QString::fromUtf8("actionSave"));
        actionOpen = new QAction(ConfFile);
        actionOpen->setObjectName(QString::fromUtf8("actionOpen"));
        actionSave_As = new QAction(ConfFile);
        actionSave_As->setObjectName(QString::fromUtf8("actionSave_As"));
        actionExit = new QAction(ConfFile);
        actionExit->setObjectName(QString::fromUtf8("actionExit"));
        centralwidget = new QWidget(ConfFile);
        centralwidget->setObjectName(QString::fromUtf8("centralwidget"));
        gridLayout_11 = new QGridLayout(centralwidget);
        gridLayout_11->setObjectName(QString::fromUtf8("gridLayout_11"));
        ConfFile_tabWidget = new QTabWidget(centralwidget);
        ConfFile_tabWidget->setObjectName(QString::fromUtf8("ConfFile_tabWidget"));
        ConfFile_tab = new QWidget();
        ConfFile_tab->setObjectName(QString::fromUtf8("ConfFile_tab"));
        gridLayout_12 = new QGridLayout(ConfFile_tab);
        gridLayout_12->setObjectName(QString::fromUtf8("gridLayout_12"));
        groupBox_4 = new QGroupBox(ConfFile_tab);
        groupBox_4->setObjectName(QString::fromUtf8("groupBox_4"));
        QSizePolicy sizePolicy1(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(groupBox_4->sizePolicy().hasHeightForWidth());
        groupBox_4->setSizePolicy(sizePolicy1);
        gridLayout_4 = new QGridLayout(groupBox_4);
        gridLayout_4->setObjectName(QString::fromUtf8("gridLayout_4"));
        IntCnt_ComboBox = new QComboBox(groupBox_4);
        IntCnt_ComboBox->setObjectName(QString::fromUtf8("IntCnt_ComboBox"));

        gridLayout_4->addWidget(IntCnt_ComboBox, 0, 0, 1, 2);

        label_10 = new QLabel(groupBox_4);
        label_10->setObjectName(QString::fromUtf8("label_10"));

        gridLayout_4->addWidget(label_10, 1, 0, 1, 1);

        Beta_dSpinBox = new QDoubleSpinBox(groupBox_4);
        Beta_dSpinBox->setObjectName(QString::fromUtf8("Beta_dSpinBox"));
        Beta_dSpinBox->setDecimals(4);
        Beta_dSpinBox->setMinimum(-5);
        Beta_dSpinBox->setSingleStep(0.5);
        Beta_dSpinBox->setValue(0.25);

        gridLayout_4->addWidget(Beta_dSpinBox, 1, 1, 1, 1);

        label_11 = new QLabel(groupBox_4);
        label_11->setObjectName(QString::fromUtf8("label_11"));

        gridLayout_4->addWidget(label_11, 2, 0, 1, 1);

        Gamma_dSpinBox = new QDoubleSpinBox(groupBox_4);
        Gamma_dSpinBox->setObjectName(QString::fromUtf8("Gamma_dSpinBox"));
        Gamma_dSpinBox->setDecimals(4);
        Gamma_dSpinBox->setMinimum(-5);
        Gamma_dSpinBox->setMaximum(5);
        Gamma_dSpinBox->setSingleStep(0.5);
        Gamma_dSpinBox->setValue(0.5);

        gridLayout_4->addWidget(Gamma_dSpinBox, 2, 1, 1, 1);


        gridLayout_12->addWidget(groupBox_4, 0, 1, 1, 1);

        groupBox_3 = new QGroupBox(ConfFile_tab);
        groupBox_3->setObjectName(QString::fromUtf8("groupBox_3"));
        sizePolicy1.setHeightForWidth(groupBox_3->sizePolicy().hasHeightForWidth());
        groupBox_3->setSizePolicy(sizePolicy1);
        gridLayout_3 = new QGridLayout(groupBox_3);
        gridLayout_3->setObjectName(QString::fromUtf8("gridLayout_3"));
        label_6 = new QLabel(groupBox_3);
        label_6->setObjectName(QString::fromUtf8("label_6"));

        gridLayout_3->addWidget(label_6, 0, 0, 1, 1);

        P_dSpinBox = new QDoubleSpinBox(groupBox_3);
        P_dSpinBox->setObjectName(QString::fromUtf8("P_dSpinBox"));
        P_dSpinBox->setDecimals(2);
        P_dSpinBox->setSingleStep(0.5);
        P_dSpinBox->setValue(0.95);

        gridLayout_3->addWidget(P_dSpinBox, 0, 1, 1, 1);

        label_7 = new QLabel(groupBox_3);
        label_7->setObjectName(QString::fromUtf8("label_7"));

        gridLayout_3->addWidget(label_7, 1, 0, 1, 1);

        I_dSpinBox = new QDoubleSpinBox(groupBox_3);
        I_dSpinBox->setObjectName(QString::fromUtf8("I_dSpinBox"));
        I_dSpinBox->setDecimals(2);
        I_dSpinBox->setSingleStep(0.5);

        gridLayout_3->addWidget(I_dSpinBox, 1, 1, 1, 1);

        label_8 = new QLabel(groupBox_3);
        label_8->setObjectName(QString::fromUtf8("label_8"));

        gridLayout_3->addWidget(label_8, 2, 0, 1, 1);

        D_dSpinBox = new QDoubleSpinBox(groupBox_3);
        D_dSpinBox->setObjectName(QString::fromUtf8("D_dSpinBox"));
        D_dSpinBox->setDecimals(2);
        D_dSpinBox->setSingleStep(0.5);

        gridLayout_3->addWidget(D_dSpinBox, 2, 1, 1, 1);


        gridLayout_12->addWidget(groupBox_3, 0, 2, 1, 1);

        verticalLayout_3 = new QVBoxLayout();
        verticalLayout_3->setObjectName(QString::fromUtf8("verticalLayout_3"));
        CFPreview_textBrowser = new QTextBrowser(ConfFile_tab);
        CFPreview_textBrowser->setObjectName(QString::fromUtf8("CFPreview_textBrowser"));
        CFPreview_textBrowser->setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded);
        CFPreview_textBrowser->setAutoFormatting(QTextEdit::AutoNone);
        CFPreview_textBrowser->setLineWrapMode(QTextEdit::NoWrap);

        verticalLayout_3->addWidget(CFPreview_textBrowser);

        GenCFile_pushButton = new QPushButton(ConfFile_tab);
        GenCFile_pushButton->setObjectName(QString::fromUtf8("GenCFile_pushButton"));

        verticalLayout_3->addWidget(GenCFile_pushButton);


        gridLayout_12->addLayout(verticalLayout_3, 0, 3, 3, 1);

        groupBox_5 = new QGroupBox(ConfFile_tab);
        groupBox_5->setObjectName(QString::fromUtf8("groupBox_5"));
        sizePolicy1.setHeightForWidth(groupBox_5->sizePolicy().hasHeightForWidth());
        groupBox_5->setSizePolicy(sizePolicy1);
        gridLayout = new QGridLayout(groupBox_5);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        label_14 = new QLabel(groupBox_5);
        label_14->setObjectName(QString::fromUtf8("label_14"));

        gridLayout->addWidget(label_14, 0, 0, 1, 1);

        NumSub_SpinBox = new QSpinBox(groupBox_5);
        NumSub_SpinBox->setObjectName(QString::fromUtf8("NumSub_SpinBox"));
        NumSub_SpinBox->setMaximum(9999);
        NumSub_SpinBox->setValue(1);

        gridLayout->addWidget(NumSub_SpinBox, 0, 1, 1, 1);

        label_9 = new QLabel(groupBox_5);
        label_9->setObjectName(QString::fromUtf8("label_9"));

        gridLayout->addWidget(label_9, 1, 0, 1, 1);

        NumSubSteps_SpinBox = new QSpinBox(groupBox_5);
        NumSubSteps_SpinBox->setObjectName(QString::fromUtf8("NumSubSteps_SpinBox"));
        NumSubSteps_SpinBox->setMaximum(9999);
        NumSubSteps_SpinBox->setValue(4);

        gridLayout->addWidget(NumSubSteps_SpinBox, 1, 1, 1, 1);


        gridLayout_12->addWidget(groupBox_5, 1, 1, 1, 1);

        groupBox_2 = new QGroupBox(ConfFile_tab);
        groupBox_2->setObjectName(QString::fromUtf8("groupBox_2"));
        sizePolicy1.setHeightForWidth(groupBox_2->sizePolicy().hasHeightForWidth());
        groupBox_2->setSizePolicy(sizePolicy1);
        gridLayout_2 = new QGridLayout(groupBox_2);
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        label_4 = new QLabel(groupBox_2);
        label_4->setObjectName(QString::fromUtf8("label_4"));

        gridLayout_2->addWidget(label_4, 0, 0, 1, 1);

        RayAlpha_dSpinBox = new QDoubleSpinBox(groupBox_2);
        RayAlpha_dSpinBox->setObjectName(QString::fromUtf8("RayAlpha_dSpinBox"));
        RayAlpha_dSpinBox->setDecimals(7);
        RayAlpha_dSpinBox->setSingleStep(0.5);
        RayAlpha_dSpinBox->setValue(0.197);

        gridLayout_2->addWidget(RayAlpha_dSpinBox, 0, 1, 1, 1);

        label_5 = new QLabel(groupBox_2);
        label_5->setObjectName(QString::fromUtf8("label_5"));

        gridLayout_2->addWidget(label_5, 1, 0, 1, 1);

        RayBeta_dSpinBox = new QDoubleSpinBox(groupBox_2);
        RayBeta_dSpinBox->setObjectName(QString::fromUtf8("RayBeta_dSpinBox"));
        RayBeta_dSpinBox->setDecimals(7);
        RayBeta_dSpinBox->setSingleStep(0.5);
        RayBeta_dSpinBox->setValue(9.2e-05);

        gridLayout_2->addWidget(RayBeta_dSpinBox, 1, 1, 1, 1);


        gridLayout_12->addWidget(groupBox_2, 1, 2, 1, 1);

        groupBox_6 = new QGroupBox(ConfFile_tab);
        groupBox_6->setObjectName(QString::fromUtf8("groupBox_6"));
        horizontalLayout = new QHBoxLayout(groupBox_6);
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        gridLayout_6 = new QGridLayout();
        gridLayout_6->setObjectName(QString::fromUtf8("gridLayout_6"));
        ConfFile_label = new QLabel(groupBox_6);
        ConfFile_label->setObjectName(QString::fromUtf8("ConfFile_label"));

        gridLayout_6->addWidget(ConfFile_label, 0, 0, 1, 1);

        ConfFile_lineEdit = new QLineEdit(groupBox_6);
        ConfFile_lineEdit->setObjectName(QString::fromUtf8("ConfFile_lineEdit"));

        gridLayout_6->addWidget(ConfFile_lineEdit, 0, 1, 1, 1);

        ConfFile_toolButton = new QToolButton(groupBox_6);
        ConfFile_toolButton->setObjectName(QString::fromUtf8("ConfFile_toolButton"));

        gridLayout_6->addWidget(ConfFile_toolButton, 0, 2, 1, 1);

        Mass_label = new QLabel(groupBox_6);
        Mass_label->setObjectName(QString::fromUtf8("Mass_label"));

        gridLayout_6->addWidget(Mass_label, 1, 0, 1, 1);

        Mass_lineEdit = new QLineEdit(groupBox_6);
        Mass_lineEdit->setObjectName(QString::fromUtf8("Mass_lineEdit"));

        gridLayout_6->addWidget(Mass_lineEdit, 1, 1, 1, 1);

        Mass_toolButton = new QToolButton(groupBox_6);
        Mass_toolButton->setObjectName(QString::fromUtf8("Mass_toolButton"));

        gridLayout_6->addWidget(Mass_toolButton, 1, 2, 1, 1);

        Stiffness_label = new QLabel(groupBox_6);
        Stiffness_label->setObjectName(QString::fromUtf8("Stiffness_label"));

        gridLayout_6->addWidget(Stiffness_label, 2, 0, 1, 1);

        Stiff_lineEdit = new QLineEdit(groupBox_6);
        Stiff_lineEdit->setObjectName(QString::fromUtf8("Stiff_lineEdit"));

        gridLayout_6->addWidget(Stiff_lineEdit, 2, 1, 1, 1);

        Stiff_toolButton = new QToolButton(groupBox_6);
        Stiff_toolButton->setObjectName(QString::fromUtf8("Stiff_toolButton"));

        gridLayout_6->addWidget(Stiff_toolButton, 2, 2, 1, 1);

        Damp_label = new QLabel(groupBox_6);
        Damp_label->setObjectName(QString::fromUtf8("Damp_label"));

        gridLayout_6->addWidget(Damp_label, 3, 0, 1, 1);

        Damp_lineEdit = new QLineEdit(groupBox_6);
        Damp_lineEdit->setObjectName(QString::fromUtf8("Damp_lineEdit"));

        gridLayout_6->addWidget(Damp_lineEdit, 3, 1, 1, 1);

        Damp_toolButton = new QToolButton(groupBox_6);
        Damp_toolButton->setObjectName(QString::fromUtf8("Damp_toolButton"));

        gridLayout_6->addWidget(Damp_toolButton, 3, 2, 1, 1);

        Gain_label = new QLabel(groupBox_6);
        Gain_label->setObjectName(QString::fromUtf8("Gain_label"));

        gridLayout_6->addWidget(Gain_label, 4, 0, 1, 1);

        Gain_lineEdit = new QLineEdit(groupBox_6);
        Gain_lineEdit->setObjectName(QString::fromUtf8("Gain_lineEdit"));

        gridLayout_6->addWidget(Gain_lineEdit, 4, 1, 1, 1);

        Gain_toolButton = new QToolButton(groupBox_6);
        Gain_toolButton->setObjectName(QString::fromUtf8("Gain_toolButton"));

        gridLayout_6->addWidget(Gain_toolButton, 4, 2, 1, 1);

        LVector_label = new QLabel(groupBox_6);
        LVector_label->setObjectName(QString::fromUtf8("LVector_label"));

        gridLayout_6->addWidget(LVector_label, 5, 0, 1, 1);

        LVector_lineEdit = new QLineEdit(groupBox_6);
        LVector_lineEdit->setObjectName(QString::fromUtf8("LVector_lineEdit"));

        gridLayout_6->addWidget(LVector_lineEdit, 5, 1, 1, 1);

        LVector_toolButton = new QToolButton(groupBox_6);
        LVector_toolButton->setObjectName(QString::fromUtf8("LVector_toolButton"));

        gridLayout_6->addWidget(LVector_toolButton, 5, 2, 1, 1);

        LVectorForm_label = new QLabel(groupBox_6);
        LVectorForm_label->setObjectName(QString::fromUtf8("LVectorForm_label"));

        gridLayout_6->addWidget(LVectorForm_label, 6, 0, 1, 1);

        LVectorForm_lineEdit = new QLineEdit(groupBox_6);
        LVectorForm_lineEdit->setObjectName(QString::fromUtf8("LVectorForm_lineEdit"));

        gridLayout_6->addWidget(LVectorForm_lineEdit, 6, 1, 1, 1);

        CNodes_label = new QLabel(groupBox_6);
        CNodes_label->setObjectName(QString::fromUtf8("CNodes_label"));

        gridLayout_6->addWidget(CNodes_label, 7, 0, 1, 1);

        CNodes_lineEdit = new QLineEdit(groupBox_6);
        CNodes_lineEdit->setObjectName(QString::fromUtf8("CNodes_lineEdit"));

        gridLayout_6->addWidget(CNodes_lineEdit, 7, 1, 1, 1);

        CNodes_toolButton = new QToolButton(groupBox_6);
        CNodes_toolButton->setObjectName(QString::fromUtf8("CNodes_toolButton"));

        gridLayout_6->addWidget(CNodes_toolButton, 7, 2, 1, 1);

        GMotion_label = new QLabel(groupBox_6);
        GMotion_label->setObjectName(QString::fromUtf8("GMotion_label"));

        gridLayout_6->addWidget(GMotion_label, 8, 0, 1, 1);

        GMotion_lineEdit = new QLineEdit(groupBox_6);
        GMotion_lineEdit->setObjectName(QString::fromUtf8("GMotion_lineEdit"));

        gridLayout_6->addWidget(GMotion_lineEdit, 8, 1, 1, 1);

        GMotion_toolButton = new QToolButton(groupBox_6);
        GMotion_toolButton->setObjectName(QString::fromUtf8("GMotion_toolButton"));

        gridLayout_6->addWidget(GMotion_toolButton, 8, 2, 1, 1);

        OutFile_label = new QLabel(groupBox_6);
        OutFile_label->setObjectName(QString::fromUtf8("OutFile_label"));

        gridLayout_6->addWidget(OutFile_label, 9, 0, 1, 1);

        OutFile_lineEdit = new QLineEdit(groupBox_6);
        OutFile_lineEdit->setObjectName(QString::fromUtf8("OutFile_lineEdit"));

        gridLayout_6->addWidget(OutFile_lineEdit, 9, 1, 1, 1);

        OutFile_toolButton = new QToolButton(groupBox_6);
        OutFile_toolButton->setObjectName(QString::fromUtf8("OutFile_toolButton"));

        gridLayout_6->addWidget(OutFile_toolButton, 9, 2, 1, 1);


        horizontalLayout->addLayout(gridLayout_6);


        gridLayout_12->addWidget(groupBox_6, 2, 0, 1, 3);

        groupBox = new QGroupBox(ConfFile_tab);
        groupBox->setObjectName(QString::fromUtf8("groupBox"));
        sizePolicy.setHeightForWidth(groupBox->sizePolicy().hasHeightForWidth());
        groupBox->setSizePolicy(sizePolicy);
        verticalLayout_5 = new QVBoxLayout(groupBox);
        verticalLayout_5->setObjectName(QString::fromUtf8("verticalLayout_5"));
        gridLayout_5 = new QGridLayout();
        gridLayout_5->setObjectName(QString::fromUtf8("gridLayout_5"));
        AbsValues_checkBox = new QCheckBox(groupBox);
        AbsValues_checkBox->setObjectName(QString::fromUtf8("AbsValues_checkBox"));
        AbsValues_checkBox->setChecked(true);

        gridLayout_5->addWidget(AbsValues_checkBox, 0, 0, 1, 1);

        ReadGain_checkBox = new QCheckBox(groupBox);
        ReadGain_checkBox->setObjectName(QString::fromUtf8("ReadGain_checkBox"));
        ReadGain_checkBox->setCheckable(true);
        ReadGain_checkBox->setChecked(false);

        gridLayout_5->addWidget(ReadGain_checkBox, 0, 1, 1, 3);

        ReadSparse_checkBox = new QCheckBox(groupBox);
        ReadSparse_checkBox->setObjectName(QString::fromUtf8("ReadSparse_checkBox"));
        ReadSparse_checkBox->setChecked(true);

        gridLayout_5->addWidget(ReadSparse_checkBox, 1, 0, 1, 1);

        label = new QLabel(groupBox);
        label->setObjectName(QString::fromUtf8("label"));

        gridLayout_5->addWidget(label, 1, 1, 1, 1);

        Order_SpinBox = new QSpinBox(groupBox);
        Order_SpinBox->setObjectName(QString::fromUtf8("Order_SpinBox"));
        Order_SpinBox->setMaximum(9999999);
        Order_SpinBox->setValue(504);

        gridLayout_5->addWidget(Order_SpinBox, 1, 2, 1, 2);

        UseSparse_checkBox = new QCheckBox(groupBox);
        UseSparse_checkBox->setObjectName(QString::fromUtf8("UseSparse_checkBox"));
        UseSparse_checkBox->setChecked(true);

        gridLayout_5->addWidget(UseSparse_checkBox, 2, 0, 1, 1);

        label_3 = new QLabel(groupBox);
        label_3->setObjectName(QString::fromUtf8("label_3"));

        gridLayout_5->addWidget(label_3, 2, 1, 1, 2);

        TimeStep_dSpinBox = new QDoubleSpinBox(groupBox);
        TimeStep_dSpinBox->setObjectName(QString::fromUtf8("TimeStep_dSpinBox"));
        TimeStep_dSpinBox->setDecimals(3);
        TimeStep_dSpinBox->setValue(0.01);

        gridLayout_5->addWidget(TimeStep_dSpinBox, 2, 3, 1, 1);

        UsePacked_checkBox = new QCheckBox(groupBox);
        UsePacked_checkBox->setObjectName(QString::fromUtf8("UsePacked_checkBox"));
        UsePacked_checkBox->setChecked(true);

        gridLayout_5->addWidget(UsePacked_checkBox, 3, 0, 1, 1);

        label_2 = new QLabel(groupBox);
        label_2->setObjectName(QString::fromUtf8("label_2"));

        gridLayout_5->addWidget(label_2, 3, 1, 1, 2);

        NumSteps_SpinBox = new QSpinBox(groupBox);
        NumSteps_SpinBox->setObjectName(QString::fromUtf8("NumSteps_SpinBox"));
        NumSteps_SpinBox->setMaximum(9999);
        NumSteps_SpinBox->setValue(4096);

        gridLayout_5->addWidget(NumSteps_SpinBox, 3, 3, 1, 1);


        verticalLayout_5->addLayout(gridLayout_5);

        verticalLayout_4 = new QVBoxLayout();
        verticalLayout_4->setObjectName(QString::fromUtf8("verticalLayout_4"));
        ReadLVector_checkBox = new QCheckBox(groupBox);
        ReadLVector_checkBox->setObjectName(QString::fromUtf8("ReadLVector_checkBox"));
        ReadLVector_checkBox->setChecked(true);

        verticalLayout_4->addWidget(ReadLVector_checkBox);

        ReadDamp_checkBox = new QCheckBox(groupBox);
        ReadDamp_checkBox->setObjectName(QString::fromUtf8("ReadDamp_checkBox"));
        ReadDamp_checkBox->setCheckable(true);
        ReadDamp_checkBox->setChecked(false);

        verticalLayout_4->addWidget(ReadDamp_checkBox);


        verticalLayout_5->addLayout(verticalLayout_4);


        gridLayout_12->addWidget(groupBox, 0, 0, 2, 1);

        ConfFile_tabWidget->addTab(ConfFile_tab, QString());
        CNodes_tab = new QWidget();
        CNodes_tab->setObjectName(QString::fromUtf8("CNodes_tab"));
        gridLayout_8 = new QGridLayout(CNodes_tab);
        gridLayout_8->setObjectName(QString::fromUtf8("gridLayout_8"));
        label_23 = new QLabel(CNodes_tab);
        label_23->setObjectName(QString::fromUtf8("label_23"));

        gridLayout_8->addWidget(label_23, 0, 0, 1, 1);

        horizontalSpacer = new QSpacerItem(119, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_8->addItem(horizontalSpacer, 0, 2, 1, 1);

        label_24 = new QLabel(CNodes_tab);
        label_24->setObjectName(QString::fromUtf8("label_24"));

        gridLayout_8->addWidget(label_24, 0, 3, 1, 1);

        horizontalSpacer_6 = new QSpacerItem(119, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_8->addItem(horizontalSpacer_6, 0, 4, 1, 1);

        listWidget = new QListWidget(CNodes_tab);
        new QListWidgetItem(listWidget);
        new QListWidgetItem(listWidget);
        new QListWidgetItem(listWidget);
        new QListWidgetItem(listWidget);
        new QListWidgetItem(listWidget);
        new QListWidgetItem(listWidget);
        new QListWidgetItem(listWidget);
        new QListWidgetItem(listWidget);
        listWidget->setObjectName(QString::fromUtf8("listWidget"));
        QSizePolicy sizePolicy2(QSizePolicy::Fixed, QSizePolicy::Expanding);
        sizePolicy2.setHorizontalStretch(0);
        sizePolicy2.setVerticalStretch(0);
        sizePolicy2.setHeightForWidth(listWidget->sizePolicy().hasHeightForWidth());
        listWidget->setSizePolicy(sizePolicy2);
        listWidget->setMaximumSize(QSize(181, 16777215));

        gridLayout_8->addWidget(listWidget, 1, 0, 4, 1);

        verticalSpacer_2 = new QSpacerItem(20, 157, QSizePolicy::Minimum, QSizePolicy::Expanding);

        gridLayout_8->addItem(verticalSpacer_2, 1, 1, 1, 1);

        tableView = new QTableView(CNodes_tab);
        tableView->setObjectName(QString::fromUtf8("tableView"));

        gridLayout_8->addWidget(tableView, 1, 2, 4, 3);

        pushButton = new QPushButton(CNodes_tab);
        pushButton->setObjectName(QString::fromUtf8("pushButton"));
        sizePolicy1.setHeightForWidth(pushButton->sizePolicy().hasHeightForWidth());
        pushButton->setSizePolicy(sizePolicy1);

        gridLayout_8->addWidget(pushButton, 2, 1, 1, 1);

        pushButton_2 = new QPushButton(CNodes_tab);
        pushButton_2->setObjectName(QString::fromUtf8("pushButton_2"));
        sizePolicy1.setHeightForWidth(pushButton_2->sizePolicy().hasHeightForWidth());
        pushButton_2->setSizePolicy(sizePolicy1);

        gridLayout_8->addWidget(pushButton_2, 3, 1, 1, 1);

        verticalSpacer_4 = new QSpacerItem(20, 198, QSizePolicy::Minimum, QSizePolicy::Expanding);

        gridLayout_8->addItem(verticalSpacer_4, 4, 1, 1, 1);

        ConfFile_tabWidget->addTab(CNodes_tab, QString());
        Algorithm_tab = new QWidget();
        Algorithm_tab->setObjectName(QString::fromUtf8("Algorithm_tab"));
        gridLayout_7 = new QGridLayout(Algorithm_tab);
        gridLayout_7->setObjectName(QString::fromUtf8("gridLayout_7"));
        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        Exec_label = new QLabel(Algorithm_tab);
        Exec_label->setObjectName(QString::fromUtf8("Exec_label"));

        horizontalLayout_2->addWidget(Exec_label);

        Exec_lineEdit = new QLineEdit(Algorithm_tab);
        Exec_lineEdit->setObjectName(QString::fromUtf8("Exec_lineEdit"));

        horizontalLayout_2->addWidget(Exec_lineEdit);

        Exec_toolButton = new QToolButton(Algorithm_tab);
        Exec_toolButton->setObjectName(QString::fromUtf8("Exec_toolButton"));

        horizontalLayout_2->addWidget(Exec_toolButton);


        gridLayout_7->addLayout(horizontalLayout_2, 0, 0, 1, 2);

        verticalLayout = new QVBoxLayout();
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        Run_pushButton = new QPushButton(Algorithm_tab);
        Run_pushButton->setObjectName(QString::fromUtf8("Run_pushButton"));
        Run_pushButton->setMinimumSize(QSize(0, 50));

        verticalLayout->addWidget(Run_pushButton);

        OpenOutput_pushButton = new QPushButton(Algorithm_tab);
        OpenOutput_pushButton->setObjectName(QString::fromUtf8("OpenOutput_pushButton"));
        OpenOutput_pushButton->setMinimumSize(QSize(0, 50));

        verticalLayout->addWidget(OpenOutput_pushButton);

        verticalSpacer = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer);

        pushButton_4 = new QPushButton(Algorithm_tab);
        pushButton_4->setObjectName(QString::fromUtf8("pushButton_4"));
        pushButton_4->setMinimumSize(QSize(0, 0));

        verticalLayout->addWidget(pushButton_4);


        gridLayout_7->addLayout(verticalLayout, 1, 0, 1, 1);

        SAlg_textBrowser = new QTextBrowser(Algorithm_tab);
        SAlg_textBrowser->setObjectName(QString::fromUtf8("SAlg_textBrowser"));
        QPalette palette;
        QBrush brush(QColor(255, 255, 255, 255));
        brush.setStyle(Qt::SolidPattern);
        palette.setBrush(QPalette::Active, QPalette::Text, brush);
        QBrush brush1(QColor(0, 0, 0, 255));
        brush1.setStyle(Qt::SolidPattern);
        palette.setBrush(QPalette::Active, QPalette::Base, brush1);
        palette.setBrush(QPalette::Inactive, QPalette::Text, brush);
        palette.setBrush(QPalette::Inactive, QPalette::Base, brush1);
        QBrush brush2(QColor(169, 167, 167, 255));
        brush2.setStyle(Qt::SolidPattern);
        palette.setBrush(QPalette::Disabled, QPalette::Text, brush2);
        QBrush brush3(QColor(244, 244, 244, 255));
        brush3.setStyle(Qt::SolidPattern);
        palette.setBrush(QPalette::Disabled, QPalette::Base, brush3);
        SAlg_textBrowser->setPalette(palette);

        gridLayout_7->addWidget(SAlg_textBrowser, 1, 1, 1, 1);

        ConfFile_tabWidget->addTab(Algorithm_tab, QString());

        gridLayout_11->addWidget(ConfFile_tabWidget, 0, 0, 1, 1);

        ConfFile->setCentralWidget(centralwidget);
        menubar = new QMenuBar(ConfFile);
        menubar->setObjectName(QString::fromUtf8("menubar"));
        menubar->setGeometry(QRect(0, 0, 962, 20));
        menuFile = new QMenu(menubar);
        menuFile->setObjectName(QString::fromUtf8("menuFile"));
        ConfFile->setMenuBar(menubar);
        statusbar = new QStatusBar(ConfFile);
        statusbar->setObjectName(QString::fromUtf8("statusbar"));
        ConfFile->setStatusBar(statusbar);

        menubar->addAction(menuFile->menuAction());
        menuFile->addAction(actionOpen);
        menuFile->addSeparator();
        menuFile->addAction(actionSave);
        menuFile->addAction(actionSave_As);
        menuFile->addSeparator();
        menuFile->addAction(actionExit);

        retranslateUi(ConfFile);
        QObject::connect(ReadGain_checkBox, SIGNAL(toggled(bool)), Gain_label, SLOT(setVisible(bool)));
        QObject::connect(ReadGain_checkBox, SIGNAL(toggled(bool)), Gain_lineEdit, SLOT(setVisible(bool)));
        QObject::connect(ReadGain_checkBox, SIGNAL(toggled(bool)), Gain_toolButton, SLOT(setVisible(bool)));
        QObject::connect(ReadDamp_checkBox, SIGNAL(toggled(bool)), Damp_label, SLOT(setVisible(bool)));
        QObject::connect(ReadDamp_checkBox, SIGNAL(toggled(bool)), Damp_lineEdit, SLOT(setVisible(bool)));
        QObject::connect(ReadDamp_checkBox, SIGNAL(toggled(bool)), Damp_toolButton, SLOT(setVisible(bool)));
        QObject::connect(ReadLVector_checkBox, SIGNAL(toggled(bool)), LVector_label, SLOT(setVisible(bool)));
        QObject::connect(ReadLVector_checkBox, SIGNAL(toggled(bool)), LVector_lineEdit, SLOT(setVisible(bool)));
        QObject::connect(ReadLVector_checkBox, SIGNAL(toggled(bool)), LVector_toolButton, SLOT(setVisible(bool)));
        QObject::connect(ReadLVector_checkBox, SIGNAL(toggled(bool)), LVectorForm_label, SLOT(setHidden(bool)));
        QObject::connect(ReadLVector_checkBox, SIGNAL(toggled(bool)), LVectorForm_lineEdit, SLOT(setHidden(bool)));
        QObject::connect(ReadDamp_checkBox, SIGNAL(toggled(bool)), RayAlpha_dSpinBox, SLOT(setDisabled(bool)));
        QObject::connect(ReadDamp_checkBox, SIGNAL(toggled(bool)), RayBeta_dSpinBox, SLOT(setDisabled(bool)));

        ConfFile_tabWidget->setCurrentIndex(2);
        IntCnt_ComboBox->setCurrentIndex(5);


        QMetaObject::connectSlotsByName(ConfFile);
    } // setupUi

    void retranslateUi(QMainWindow *ConfFile)
    {
        ConfFile->setWindowTitle(QApplication::translate("ConfFile", "SAlgorithm", 0, QApplication::UnicodeUTF8));
        actionSave->setText(QApplication::translate("ConfFile", "Save", 0, QApplication::UnicodeUTF8));
        actionSave->setShortcut(QApplication::translate("ConfFile", "Ctrl+S", 0, QApplication::UnicodeUTF8));
        actionOpen->setText(QApplication::translate("ConfFile", "Open...", 0, QApplication::UnicodeUTF8));
        actionOpen->setShortcut(QApplication::translate("ConfFile", "Ctrl+O", 0, QApplication::UnicodeUTF8));
        actionSave_As->setText(QApplication::translate("ConfFile", "Save As...", 0, QApplication::UnicodeUTF8));
        actionExit->setText(QApplication::translate("ConfFile", "Exit", 0, QApplication::UnicodeUTF8));
        actionExit->setShortcut(QApplication::translate("ConfFile", "Ctrl+Q", 0, QApplication::UnicodeUTF8));
        groupBox_4->setTitle(QApplication::translate("ConfFile", "Integration process", 0, QApplication::UnicodeUTF8));
        IntCnt_ComboBox->clear();
        IntCnt_ComboBox->insertItems(0, QStringList()
         << QApplication::translate("ConfFile", "Central explicit", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("ConfFile", "Backward", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("ConfFile", "Linear acceleration", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("ConfFile", "Galerkin", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("ConfFile", "Fox Goodwin", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("ConfFile", "Average acceleration", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("ConfFile", "Custom", 0, QApplication::UnicodeUTF8)
        );
        label_10->setText(QApplication::translate("ConfFile", "Beta:", 0, QApplication::UnicodeUTF8));
        label_11->setText(QApplication::translate("ConfFile", "Gamma:", 0, QApplication::UnicodeUTF8));
        groupBox_3->setTitle(QApplication::translate("ConfFile", "PID", 0, QApplication::UnicodeUTF8));
        label_6->setText(QApplication::translate("ConfFile", "Proportional:", 0, QApplication::UnicodeUTF8));
        label_7->setText(QApplication::translate("ConfFile", "Integral:", 0, QApplication::UnicodeUTF8));
        label_8->setText(QApplication::translate("ConfFile", "Derivative:", 0, QApplication::UnicodeUTF8));
        GenCFile_pushButton->setText(QApplication::translate("ConfFile", "Generate", 0, QApplication::UnicodeUTF8));
        groupBox_5->setTitle(QApplication::translate("ConfFile", "Substructure", 0, QApplication::UnicodeUTF8));
        label_14->setText(QApplication::translate("ConfFile", "N. Substructures:", 0, QApplication::UnicodeUTF8));
        label_9->setText(QApplication::translate("ConfFile", "N. Sub-steps:", 0, QApplication::UnicodeUTF8));
        groupBox_2->setTitle(QApplication::translate("ConfFile", "Rayleigh", 0, QApplication::UnicodeUTF8));
        label_4->setText(QApplication::translate("ConfFile", "Alpha:", 0, QApplication::UnicodeUTF8));
        label_5->setText(QApplication::translate("ConfFile", "Beta:", 0, QApplication::UnicodeUTF8));
        groupBox_6->setTitle(QApplication::translate("ConfFile", "Input/Output Files", 0, QApplication::UnicodeUTF8));
        ConfFile_label->setText(QApplication::translate("ConfFile", "Conf. File:", 0, QApplication::UnicodeUTF8));
        ConfFile_lineEdit->setText(QApplication::translate("ConfFile", "ConfFile.conf", 0, QApplication::UnicodeUTF8));
        ConfFile_toolButton->setText(QApplication::translate("ConfFile", "...", 0, QApplication::UnicodeUTF8));
        Mass_label->setText(QApplication::translate("ConfFile", "Mass Matrix:", 0, QApplication::UnicodeUTF8));
        Mass_lineEdit->setText(QApplication::translate("ConfFile", "504M_Sp_MM.txt", 0, QApplication::UnicodeUTF8));
        Mass_toolButton->setText(QApplication::translate("ConfFile", "...", 0, QApplication::UnicodeUTF8));
        Stiffness_label->setText(QApplication::translate("ConfFile", "Stiffness Matrix:", 0, QApplication::UnicodeUTF8));
        Stiff_lineEdit->setText(QApplication::translate("ConfFile", "504K_Sp_MM.txt", 0, QApplication::UnicodeUTF8));
        Stiff_toolButton->setText(QApplication::translate("ConfFile", "...", 0, QApplication::UnicodeUTF8));
        Damp_label->setText(QApplication::translate("ConfFile", "Damping Matrix:", 0, QApplication::UnicodeUTF8));
        Damp_lineEdit->setText(QApplication::translate("ConfFile", "504C_Sp_MM.txt", 0, QApplication::UnicodeUTF8));
        Damp_toolButton->setText(QApplication::translate("ConfFile", "...", 0, QApplication::UnicodeUTF8));
        Gain_label->setText(QApplication::translate("ConfFile", "Gain Matrix:", 0, QApplication::UnicodeUTF8));
        Gain_lineEdit->setText(QApplication::translate("ConfFile", "504G_MM.txt", 0, QApplication::UnicodeUTF8));
        Gain_toolButton->setText(QApplication::translate("ConfFile", "...", 0, QApplication::UnicodeUTF8));
        LVector_label->setText(QApplication::translate("ConfFile", "Load Vector:", 0, QApplication::UnicodeUTF8));
        LVector_lineEdit->setText(QApplication::translate("ConfFile", "504LV_MM.txt", 0, QApplication::UnicodeUTF8));
        LVector_toolButton->setText(QApplication::translate("ConfFile", "...", 0, QApplication::UnicodeUTF8));
        LVectorForm_label->setText(QApplication::translate("ConfFile", "Excited DOFs:", 0, QApplication::UnicodeUTF8));
        LVectorForm_lineEdit->setText(QApplication::translate("ConfFile", "6 1 1 1 0 0 0", 0, QApplication::UnicodeUTF8));
        CNodes_label->setText(QApplication::translate("ConfFile", "Coupling Nodes:", 0, QApplication::UnicodeUTF8));
        CNodes_lineEdit->setText(QApplication::translate("ConfFile", "Couple_Nodes.txt", 0, QApplication::UnicodeUTF8));
        CNodes_toolButton->setText(QApplication::translate("ConfFile", "...", 0, QApplication::UnicodeUTF8));
        GMotion_label->setText(QApplication::translate("ConfFile", "Ground Motion:", 0, QApplication::UnicodeUTF8));
        GMotion_lineEdit->setText(QApplication::translate("ConfFile", "GroundMovement.txt", 0, QApplication::UnicodeUTF8));
        GMotion_toolButton->setText(QApplication::translate("ConfFile", "...", 0, QApplication::UnicodeUTF8));
        OutFile_label->setText(QApplication::translate("ConfFile", "Output File:", 0, QApplication::UnicodeUTF8));
        OutFile_lineEdit->setText(QApplication::translate("ConfFile", "Out-Sparse.h5", 0, QApplication::UnicodeUTF8));
        OutFile_toolButton->setText(QApplication::translate("ConfFile", "...", 0, QApplication::UnicodeUTF8));
        groupBox->setTitle(QApplication::translate("ConfFile", "General", 0, QApplication::UnicodeUTF8));
        AbsValues_checkBox->setText(QApplication::translate("ConfFile", "Absolute Values", 0, QApplication::UnicodeUTF8));
        ReadGain_checkBox->setText(QApplication::translate("ConfFile", "Read Gain Matrix", 0, QApplication::UnicodeUTF8));
        ReadSparse_checkBox->setText(QApplication::translate("ConfFile", "Read Sparse", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("ConfFile", "Order:", 0, QApplication::UnicodeUTF8));
        UseSparse_checkBox->setText(QApplication::translate("ConfFile", "Use Sparse", 0, QApplication::UnicodeUTF8));
        label_3->setText(QApplication::translate("ConfFile", "Time step:", 0, QApplication::UnicodeUTF8));
        UsePacked_checkBox->setText(QApplication::translate("ConfFile", "Use Packed Storage", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("ConfFile", "N. Steps:", 0, QApplication::UnicodeUTF8));
        ReadLVector_checkBox->setText(QApplication::translate("ConfFile", "Read Load Vector", 0, QApplication::UnicodeUTF8));
        ReadDamp_checkBox->setText(QApplication::translate("ConfFile", "Read Damping Matrix", 0, QApplication::UnicodeUTF8));
        ConfFile_tabWidget->setTabText(ConfFile_tabWidget->indexOf(ConfFile_tab), QApplication::translate("ConfFile", "Configuration File", 0, QApplication::UnicodeUTF8));
        label_23->setText(QApplication::translate("ConfFile", "Available Substructure Types:", 0, QApplication::UnicodeUTF8));
        label_24->setText(QApplication::translate("ConfFile", "Defined substructures", 0, QApplication::UnicodeUTF8));

        const bool __sortingEnabled = listWidget->isSortingEnabled();
        listWidget->setSortingEnabled(false);
        QListWidgetItem *___qlistwidgetitem = listWidget->item(0);
        ___qlistwidgetitem->setText(QApplication::translate("ConfFile", "Exact Integration", 0, QApplication::UnicodeUTF8));
        QListWidgetItem *___qlistwidgetitem1 = listWidget->item(1);
        ___qlistwidgetitem1->setText(QApplication::translate("ConfFile", "Experimental (ADwin)", 0, QApplication::UnicodeUTF8));
        QListWidgetItem *___qlistwidgetitem2 = listWidget->item(2);
        ___qlistwidgetitem2->setText(QApplication::translate("ConfFile", "Measured (File)", 0, QApplication::UnicodeUTF8));
        QListWidgetItem *___qlistwidgetitem3 = listWidget->item(3);
        ___qlistwidgetitem3->setText(QApplication::translate("ConfFile", "Remote - NSEP", 0, QApplication::UnicodeUTF8));
        QListWidgetItem *___qlistwidgetitem4 = listWidget->item(4);
        ___qlistwidgetitem4->setText(QApplication::translate("ConfFile", "Remote - OpenFresco", 0, QApplication::UnicodeUTF8));
        QListWidgetItem *___qlistwidgetitem5 = listWidget->item(5);
        ___qlistwidgetitem5->setText(QApplication::translate("ConfFile", "Remote - TCP", 0, QApplication::UnicodeUTF8));
        QListWidgetItem *___qlistwidgetitem6 = listWidget->item(6);
        ___qlistwidgetitem6->setText(QApplication::translate("ConfFile", "Remote - UDP", 0, QApplication::UnicodeUTF8));
        QListWidgetItem *___qlistwidgetitem7 = listWidget->item(7);
        ___qlistwidgetitem7->setText(QApplication::translate("ConfFile", "UHYDE-fbr", 0, QApplication::UnicodeUTF8));
        listWidget->setSortingEnabled(__sortingEnabled);

        pushButton->setText(QApplication::translate("ConfFile", ">>", 0, QApplication::UnicodeUTF8));
        pushButton_2->setText(QApplication::translate("ConfFile", "<<", 0, QApplication::UnicodeUTF8));
        ConfFile_tabWidget->setTabText(ConfFile_tabWidget->indexOf(CNodes_tab), QApplication::translate("ConfFile", "Coupling Nodes", 0, QApplication::UnicodeUTF8));
        Exec_label->setText(QApplication::translate("ConfFile", "Executable:", 0, QApplication::UnicodeUTF8));
        Exec_lineEdit->setText(QApplication::translate("ConfFile", "./SAlgorithm", 0, QApplication::UnicodeUTF8));
        Exec_toolButton->setText(QApplication::translate("ConfFile", "...", 0, QApplication::UnicodeUTF8));
        Run_pushButton->setText(QApplication::translate("ConfFile", "Run", 0, QApplication::UnicodeUTF8));
        OpenOutput_pushButton->setText(QApplication::translate("ConfFile", "Open Output File", 0, QApplication::UnicodeUTF8));
        pushButton_4->setText(QApplication::translate("ConfFile", "Save Log", 0, QApplication::UnicodeUTF8));
        ConfFile_tabWidget->setTabText(ConfFile_tabWidget->indexOf(Algorithm_tab), QApplication::translate("ConfFile", "Algorithm", 0, QApplication::UnicodeUTF8));
        menuFile->setTitle(QApplication::translate("ConfFile", "File", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class ConfFile: public Ui_ConfFile {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_CONFFILEGUI_H
