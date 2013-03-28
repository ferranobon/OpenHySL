/****************************************************************************
** Meta object code from reading C++ file 'ConfFile.hpp'
**
** Created: Thu Mar 28 11:56:27 2013
**      by: The Qt Meta Object Compiler version 63 (Qt 4.8.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "include/ConfFile.hpp"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'ConfFile.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_ConfFile[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
      17,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      10,    9,    9,    9, 0x08,
      26,    9,    9,    9, 0x08,
      55,    9,    9,    9, 0x08,
      77,   73,    9,    9, 0x08,
     111,    9,    9,    9, 0x08,
     139,    9,    9,    9, 0x08,
     163,    9,    9,    9, 0x08,
     187,    9,    9,    9, 0x08,
     212,    9,    9,    9, 0x08,
     236,    9,    9,    9, 0x08,
     260,    9,    9,    9, 0x08,
     287,    9,    9,    9, 0x08,
     313,    9,    9,    9, 0x08,
     340,    9,    9,    9, 0x08,
     367,    9,    9,    9, 0x08,
     395,    9,    9,    9, 0x08,
     411,    9,    9,    9, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_ConfFile[] = {
    "ConfFile\0\0Check_If_Zero()\0"
    "Match_TimeIntegration_Cnts()\0"
    "Refresh_Preview()\0str\0"
    "Set_TimeIntegration_Cnts(QString)\0"
    "On_GenCFileButton_Clicked()\0"
    "On_ConfButton_Clicked()\0On_MassButton_Clicked()\0"
    "On_StiffButton_Clicked()\0"
    "On_DampButton_Clicked()\0On_GainButton_Clicked()\0"
    "On_LVectorButton_Clicked()\0"
    "On_CNodesButton_Clicked()\0"
    "On_GMotionButton_Clicked()\0"
    "On_OutFileButton_Clicked()\0"
    "On_Run_pushButton_Clicked()\0readStdOutput()\0"
    "On_OpenOutput_pushButton_Clicked()\0"
};

void ConfFile::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        ConfFile *_t = static_cast<ConfFile *>(_o);
        switch (_id) {
        case 0: _t->Check_If_Zero(); break;
        case 1: _t->Match_TimeIntegration_Cnts(); break;
        case 2: _t->Refresh_Preview(); break;
        case 3: _t->Set_TimeIntegration_Cnts((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 4: _t->On_GenCFileButton_Clicked(); break;
        case 5: _t->On_ConfButton_Clicked(); break;
        case 6: _t->On_MassButton_Clicked(); break;
        case 7: _t->On_StiffButton_Clicked(); break;
        case 8: _t->On_DampButton_Clicked(); break;
        case 9: _t->On_GainButton_Clicked(); break;
        case 10: _t->On_LVectorButton_Clicked(); break;
        case 11: _t->On_CNodesButton_Clicked(); break;
        case 12: _t->On_GMotionButton_Clicked(); break;
        case 13: _t->On_OutFileButton_Clicked(); break;
        case 14: _t->On_Run_pushButton_Clicked(); break;
        case 15: _t->readStdOutput(); break;
        case 16: _t->On_OpenOutput_pushButton_Clicked(); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData ConfFile::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject ConfFile::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_ConfFile,
      qt_meta_data_ConfFile, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &ConfFile::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *ConfFile::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *ConfFile::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_ConfFile))
        return static_cast<void*>(const_cast< ConfFile*>(this));
    if (!strcmp(_clname, "Ui::ConfFile"))
        return static_cast< Ui::ConfFile*>(const_cast< ConfFile*>(this));
    return QMainWindow::qt_metacast(_clname);
}

int ConfFile::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 17)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 17;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
