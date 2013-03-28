#include <QApplication>

#include "ConfFile.hpp"

int main( int argc, char **argv )
{

     QApplication app( argc, argv );

     ConfFile *mainWin = new ConfFile;
     mainWin->show( );

     return app.exec( );
}
