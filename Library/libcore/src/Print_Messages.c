#include <stdio.h>  /* For fprintf() and stderr */
#include <stdlib.h> /* For exit() */
#include <assert.h> /* For assert() */

#include "Print_Messages.h"  /* Function prototypes and definition of ansi colors */

void Print_Header (const int32_t Mess_Type) {

    FILE *stream; /* File stream */

    if (Mess_Type == SUCCESS) {
        stream = stdout;
        fprintf(stream, "[" GREEN "  OK  " RESET "] ");
    } else if (Mess_Type == INFO) {
        stream = stdout;
        fprintf(stream, "[" MAGENTA " INFO " RESET "] ");
    } else if (Mess_Type == WARNING) {
        stream = stderr;
        fprintf(stream, "[" YELLOW " WARN " RESET "] ");
    } else if (Mess_Type == ERROR) {
        stream = stderr;
        fprintf(stream, "[" RED "FAILED" RESET "] ");
    } else {
        assert((Mess_Type >= ERROR) && (Mess_Type <= WARNING));
    }
}
