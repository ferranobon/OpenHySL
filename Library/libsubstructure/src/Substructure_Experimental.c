#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Print_Messages.h"
#include "Substructure_Experimental.h"

void Substructure_Experimental_Init (const char *Description, ExpSub_t *const Sub) {
    Sub->Description = strdup(Description);
}

void Substructure_Experimental_Destroy (ExpSub_t *const Sub) {
    free(Sub->Description);
}
