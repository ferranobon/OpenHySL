#include <stdio.h>
#include <stdlib.h>

#include "Substructure_Remote.h"
#include "Substructure_Remote_OpenFresco.h"

void Substructure_Remote_OpenFresco( const int Socket, const int WhatToDo, const unsigned int Size, const hysl_float_t *const Data_To_Send, hysl_float_t *const Data_To_Receive )
{

     /* Local Variables */
     static unsigned int i, iData[OF_INFO_DATA_LENGTH];
     static unsigned int DataSize = OF_DATA_SIZE;

     /* OpenFresco uses type HYSL_FLOAT to send/receive messages */
     static double *sData, *rData;

     /* Data communication variables */
     size_t Datatype_Size;

     /* Setup the connection with the OpenFresco's SimAppServer */
     if( WhatToDo == OF_REMOTE_SETUP ){

	  /* Allocate memory to send and receive vectors */
	  DataSize = ( 2*Size > DataSize ) ? 2*Size : DataSize;
	  DataSize = ( Size*Size > DataSize ) ? Size*Size : DataSize;
     
	  sData = (double *) calloc( (size_t) DataSize, sizeof (double) );
	  rData = (double *) calloc( (size_t) DataSize, sizeof (double) );

	  /* Set the Data size for the experimental element */
	  /* SizeCtrl */
	  iData[0] = Size;   /* Displacement */
	  iData[1] = 0;      /* Velocity */
	  iData[2] = 0;      /* Acceleration */
	  iData[3] = 0;      /* Force */
	  iData[4] = 0;      /* Time */
	  /* sizeDaq */
	  iData[5] = Size;   /* Displacement */
	  iData[6] = Size;   /* Velocity. This is done so that we can send a vector of size 3 */
	  iData[7] = 0;      /* Acceleration */
	  iData[8] = 2*Size; /* Force */
	  iData[9] = 0;      /* Time */
	  /* DataSize */
	  iData[10] = DataSize;

	  Datatype_Size = sizeof( int );
	  Substructure_Remote_Send( Socket, OF_INFO_DATA_LENGTH, Datatype_Size, (char *const) iData );
     } else if ( WhatToDo == OF_REMOTE_SET_TRIAL_RESPONSE ) {
	  /* Send Trial response to the experimental site */
	  sData[0] = (double) OF_REMOTE_SET_TRIAL_RESPONSE;
	  for ( i = 0; i < Size; i++ ){
	       /* ADwin can only handle float correctly. Therefore we must stick with it and make the
		* conversion from float to double so that the messages can be sent correctly using
		* OpenFresco. */
	       sData[i + 1] = (double) Data_To_Send[i];
	  }
	  Substructure_Remote_Send( Socket, DataSize, sizeof(double), (char *const) sData );
     } else if ( WhatToDo == OF_REMOTE_GET_DAQ_RESPONSE ) {
	  /* Ask for DAQ values.*/
	  sData[0] = (double) OF_REMOTE_GET_DAQ_RESPONSE;
	  Substructure_Remote_Send( Socket, DataSize, sizeof(double), (char *const) sData );

	  /* The DAQ values are received. The data received is bounded to a particular order in order to work
	   * properly: the new displacement in the first position, followed by the coupling force from the
	   * previous sub-step and the coupling force at the last sub-step. */
	  Substructure_Remote_Receive( Socket, DataSize, sizeof(double), (char *const) rData );

	  for ( i = 0; i < Size; i++ ){
	       Data_To_Receive[i]= (hysl_float_t) rData[i];
	       Data_To_Receive[i + Size]= (hysl_float_t) rData[i + Size];
	       Data_To_Receive[i + 2*Size]= (hysl_float_t) rData[i + 2*Size];
	  }
	  
     } else if ( WhatToDo == OF_REMOTE_DIE ){
	  /* Disconnect from the experimental site */
	  sData[0] = (double) OF_REMOTE_DIE;
	  Substructure_Remote_Send( Socket, DataSize, sizeof(double), (char *const) sData );

	  /* Release dinamically allocated memory */
	  free( sData );
	  free( rData );
     }
}
