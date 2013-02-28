#include <stdio.h>  /* For fprintf() and stderr */
#include <stdlib.h> /* For exit() */
#include <assert.h> /* For assert() */
#include <stdarg.h> /* For va_list, va_start(), va_arg() and va_end() */

#include "Print_Messages.h"  /* Function prototypes and definition of ansi colors */

void Print_Message( const int Mess_Type, const int Num_Pairs, ... )
{

     int i;        /* A counter */
     int type;     /* To be used to identify the type of argument in va_arg() */
     FILE *stream; /* File stream */

     va_list ap;   /* List of arguments */

     if( Mess_Type == SUCCESS ){
	  stream = stdout;
	  fprintf( stream, "[" GREEN "  OK  " RESET "]" );
     } else if ( Mess_Type == INFO ){
	  stream = stdout;
	  fprintf( stream, "[ INFO ]" );
     } else if ( Mess_Type == WARNING ){
	  stream = stderr;
	  fprintf( stream, "[" YELLOW " WARN " RESET "]" );
     } else if ( Mess_Type == ERROR ){
	  stream = stderr;
	  fprintf( stream, "[" RED "FAILED" RESET "]" );
     } else assert( Mess_Type >= ERROR && Mess_Type <= WARNING );

     va_start( ap, Num_Pairs );
     for( i=0 ; i < Num_Pairs; i++ ){
	  type = va_arg( ap, enum MyTypes );
	  switch (type){
	  case INT:
	       fprintf( stream,"%d", va_arg(ap,int) );
	       break;
	  case FLOAT:
	       /* float is automatically promoted to double when passed to va_arg */
	  case DOUBLE:
	       fprintf( stream,"%lf", va_arg(ap,double) );
	       break;
	  case STRING:
	       fprintf( stream,"%s", va_arg(ap,char *) );
	       break;
	  default: /* unknown type */
	       Print_Message( ERROR, 2, STRING, "PrintMessage: Unknown type." );
	       break;
	  }
     }
     va_end( ap );

     fprintf( stream, "\n" );
}
