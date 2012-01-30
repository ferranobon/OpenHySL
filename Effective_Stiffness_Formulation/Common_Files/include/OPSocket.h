/*
 * OPSocket.h
 *
 *  Created on: 12/07/2011
 *      Author: ferran
 */

#ifndef OPSOCKET_H_
#define OPSOCKET_H_

#include <sys/socket.h>
#include <sys/types.h>

#define MAX_UDP_DATAGRAM 9126
#define MAX_INET_ADDR 28

#ifdef _WIN32
typedef SOCKET socket_type;
typedef int socklen_type;
static int numSockets = 0;
#define bzero(s,n) memset((s),0,(n))
#define bcmp(s1,s2,n) memcmp((s1),(s2),(n))
#else
typedef int socket_type;
typedef socklen_t socklen_type;
#endif

typedef struct socketConnection {
    unsigned int port;
    char *machineInetAddr;
    int socketID;
    socket_type sockfd;
    struct socketConnection *next;
} SocketConnection;

static SocketConnection *theSockets = NULL;
static int socketIDs = 0;

void startupsockets(int *ierr);
void cleanupsockets();
void setupconnectionserver(unsigned int *port, int *socketID);
void setupconnectionclient(unsigned int *other_Port, const char other_InetAddr[], int *lengthInet, int *socketID);
void closeconnection(int *socketID, int *ierr);
void senddata(int *socketID, int *dataTypeSize, char data[], int *lenData, int *ierr);
void sendnbdata(int *socketID, int *dataTypeSize, char data[], int *lenData, int *ierr);
void recvdata(int *socketID, int *dataTypeSize, char data[], int *lenData, int *ierr);
void recvnbdata(int *socketID, int *dataTypeSize, char data[], int *lenData, int *ierr);
void getsocketid(unsigned int *port, const char machineInetAddr[], int *lengthInet, int *socketID);

#endif /* OPSOCKET_H_ */
