#include <QString>
#include <QErrorMessage>

void Print_Error_Message( const int Error, QString Error_Text )
{
     QString Full_Error_Text;
     QErrorMessage Error_Message;

     Full_Error_Text = "Error ";
     /* Append the error number */
     Full_Error_Text.append(QString("%1").arg(Error));
     Full_Error_Text.append(": ");
     Full_Error_Text.append( Error_Text );

     Error_Message.showMessage( Full_Error_Text );
     Error_Message.exec();
}
