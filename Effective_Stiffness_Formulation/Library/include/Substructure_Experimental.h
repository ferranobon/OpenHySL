#ifndef SUBSTRUCTURE_EXPERIMENTAL_H_
#define SUBSTRUCTURE_EXPERIMENTAL_H_

#include "Substructure.h"

#define NUM_CHANNELS 25  /*!< Number of channels of recorded data in ADwin */

typedef struct ExpSub{
     char *Description;
} ExpSub_t;

void Substructure_Experimental_Init( const char *Description, ExpSub_t *const Sub );
void Substructure_Experimental_Destroy( ExpSub_t *const Sub );

#endif /* SUBSTRUCTURE_EXACT_H_ */
