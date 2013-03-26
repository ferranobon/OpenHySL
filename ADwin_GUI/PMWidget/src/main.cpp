#include <QApplication>
#include <QWidget>

#include "PMWidget.hpp"

int main( int argc, char **argv )
{

     QApplication app( argc, argv );

     PMWidget *widget = new PMWidget;
     widget->show( );

     return app.exec( );
}
