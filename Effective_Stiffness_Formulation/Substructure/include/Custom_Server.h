#ifndef CUSTOM_SERVER_H
#define CUSTOM_SERVER_H

int Init_TCP_Server_Socket( const unsigned short int Server_Port );
int Accept_TCP_Client_Connection( int Server_Socket );

#endif /* CUSTOM_SERVER_H */
