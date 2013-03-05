#ifndef SUBSTRUCTURE_EXPERIMENTAL_H_
#define SUBSTRUCTURE_EXPERIMENTAL_H_

#include "Substructure.h"

#define EXACT_NUMPARAM_INIT 3  /*!< Subber of required parameters in order to initialise a substructure of
			        * type Exact fbr */

typedef struct ExpSub{
     char *Description;
} ExpSub_t;

void Substructure_Experimental_Init( const char *Description, ExpSub_t *const Sub );
void Substructure_Experimental_Destroy( ExpSub_t *const Sub );

#endif /* SUBSTRUCTURE_EXACT_H_ */
