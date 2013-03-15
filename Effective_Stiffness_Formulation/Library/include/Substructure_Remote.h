#ifndef SUBSTRUCTURE_REMOTE_H_
#define SUBSTRUCTURE_REMOTE_H_

/**
 * \brief Structure to handle the Server IP and the communication port.
 *
 * This structure is used to store the IP address of the Server (IPv4) where the sub-stepping part
 * of the algorithm is performed. The port that will be used throughout the TCP/IP communication is also stored.
 */
typedef struct Remote{
     char *IP;               /*!< IP address (IPv4) of the server. */
     char *Port;             /*!< Port that will be used for TCP/IP communication. */
     char *Account_Name;     /*!< Account information */
     char *Account_Password; /*!< Account password */
     int   Socket;
     int   NSub;
     int  *DOFs;
     char *Description;
} Remote_t;

void Substructure_Remote_Connect( Remote_t *const RemoteNode, const int Type );
void Substructure_Remote_Init( const char *IPAddress, const char *Port, const int NSub, const int *const DOF,
			       const char *Description, Remote_t *const Sub );
void Substructure_Remote_Destroy( Remote_t *const Sub );


/**
 * \brief Frees the memory allocated during GetServerInformation.
 *
 * The memory allocated during the GetServerInformation() routine by the strdup()
 * function is deallocated.
 * 
 * \param[out] Server Struct storing the IP address, port and login information.
 *
 */
void Delete_NetworkInformation( Remote_t *const Server );


#endif /* SUBSTRUCTURE_REMOTE_H_ */
