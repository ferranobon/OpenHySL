/****************************************************************************
** Meta object code from reading C++ file 'PMWidget.hpp'
**
** Created: Tue Apr 24 13:43:03 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.4)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "include/PMWidget.hpp"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'PMWidget.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.4. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_PMWidget[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
      27,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: signature, parameters, type, tag, flags
      10,    9,    9,    9, 0x05,

 // slots: signature, parameters, type, tag, flags
      32,    9,    9,    9, 0x08,
      61,    9,    9,    9, 0x08,
      85,    9,    9,    9, 0x08,
     114,    9,    9,    9, 0x08,
     143,    9,    9,    9, 0x08,
     173,    9,    9,    9, 0x08,
     203,    9,    9,    9, 0x08,
     232,    9,    9,    9, 0x08,
     261,    9,    9,    9, 0x08,
     290,    9,    9,    9, 0x08,
     318,    9,    9,    9, 0x08,
     346,    9,    9,    9, 0x08,
     376,    9,    9,    9, 0x08,
     405,    9,    9,    9, 0x08,
     436,    9,    9,    9, 0x08,
     468,    9,    9,    9, 0x08,
     500,    9,    9,    9, 0x08,
     529,    9,    9,    9, 0x08,
     575,  558,    9,    9, 0x08,
     615,  558,    9,    9, 0x08,
     655,  558,    9,    9, 0x08,
     690,  558,    9,    9, 0x08,
     725,  558,    9,    9, 0x08,
     766,  760,    9,    9, 0x08,
     797,    9,    9,    9, 0x08,
     824,    9,    9,    9, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_PMWidget[] = {
    "PMWidget\0\0Enable_Load_Process()\0"
    "On_BFiletoolButton_Clicked()\0"
    "On_BootButton_Clicked()\0"
    "On_PFiletoolButton_Clicked()\0"
    "On_LoadCtrlPButton_Clicked()\0"
    "On_StartCtrlPButton_Clicked()\0"
    "On_StopCntrlPButton_Clicked()\0"
    "On_TFiletoolButton_Clicked()\0"
    "On_LoadTestPButton_Clicked()\0"
    "On_StartTestButton_Clicked()\0"
    "On_StopTestButton_Clicked()\0"
    "On_ReadDataButton_Clicked()\0"
    "On_ADwinReadyButton_Clicked()\0"
    "On_EmergencyButton_Clicked()\0"
    "On_LowPressureButton_Clicked()\0"
    "On_HighPressureButton_Clicked()\0"
    "On_StopPressureButton_Clicked()\0"
    "On_IFiletoolButton_Clicked()\0"
    "On_OFiletoolButton_Clicked()\0"
    "New_Value_Double\0"
    "On_DCtrl1_DSpinBox_valueChanged(double)\0"
    "On_DCtrl2_DSpinBox_valueChanged(double)\0"
    "On_P_DSpinBox_valueChanged(double)\0"
    "On_I_DSpinBox_valueChanged(double)\0"
    "On_D_DSpinBox_valueChanged(double)\0"
    "state\0On_SetLimCheckBox_Checked(int)\0"
    "On_DefValPButton_Clicked()\0"
    "Read_Measurement_Values()\0"
};

const QMetaObject PMWidget::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_PMWidget,
      qt_meta_data_PMWidget, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &PMWidget::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *PMWidget::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *PMWidget::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_PMWidget))
        return static_cast<void*>(const_cast< PMWidget*>(this));
    if (!strcmp(_clname, "Ui::PMWidget"))
        return static_cast< Ui::PMWidget*>(const_cast< PMWidget*>(this));
    if (!strcmp(_clname, "ADwin_Class"))
        return static_cast< ADwin_Class*>(const_cast< PMWidget*>(this));
    return QWidget::qt_metacast(_clname);
}

int PMWidget::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: Enable_Load_Process(); break;
        case 1: On_BFiletoolButton_Clicked(); break;
        case 2: On_BootButton_Clicked(); break;
        case 3: On_PFiletoolButton_Clicked(); break;
        case 4: On_LoadCtrlPButton_Clicked(); break;
        case 5: On_StartCtrlPButton_Clicked(); break;
        case 6: On_StopCntrlPButton_Clicked(); break;
        case 7: On_TFiletoolButton_Clicked(); break;
        case 8: On_LoadTestPButton_Clicked(); break;
        case 9: On_StartTestButton_Clicked(); break;
        case 10: On_StopTestButton_Clicked(); break;
        case 11: On_ReadDataButton_Clicked(); break;
        case 12: On_ADwinReadyButton_Clicked(); break;
        case 13: On_EmergencyButton_Clicked(); break;
        case 14: On_LowPressureButton_Clicked(); break;
        case 15: On_HighPressureButton_Clicked(); break;
        case 16: On_StopPressureButton_Clicked(); break;
        case 17: On_IFiletoolButton_Clicked(); break;
        case 18: On_OFiletoolButton_Clicked(); break;
        case 19: On_DCtrl1_DSpinBox_valueChanged((*reinterpret_cast< const double(*)>(_a[1]))); break;
        case 20: On_DCtrl2_DSpinBox_valueChanged((*reinterpret_cast< const double(*)>(_a[1]))); break;
        case 21: On_P_DSpinBox_valueChanged((*reinterpret_cast< const double(*)>(_a[1]))); break;
        case 22: On_I_DSpinBox_valueChanged((*reinterpret_cast< const double(*)>(_a[1]))); break;
        case 23: On_D_DSpinBox_valueChanged((*reinterpret_cast< const double(*)>(_a[1]))); break;
        case 24: On_SetLimCheckBox_Checked((*reinterpret_cast< const int(*)>(_a[1]))); break;
        case 25: On_DefValPButton_Clicked(); break;
        case 26: Read_Measurement_Values(); break;
        default: ;
        }
        _id -= 27;
    }
    return _id;
}

// SIGNAL 0
void PMWidget::Enable_Load_Process()
{
    QMetaObject::activate(this, &staticMetaObject, 0, 0);
}
QT_END_MOC_NAMESPACE
