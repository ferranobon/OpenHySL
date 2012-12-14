#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <ctype.h>  /* For isspace */

#include "Conf_Parser.h"

#define MAX_LENGTH 256

ConfFile* ConfFile_Create( const size_t Size )
{
     ConfFile *Config;

     Config = (ConfFile *) calloc( 1, sizeof(ConfFile) );
     if( Config == NULL ){
	  fprintf( stderr, "Error allocating memory\n" );
	  return NULL;
     }

     Config->Length = Size;
     Config->Keys = (char **) calloc( Size, sizeof(char*));
     if( Config->Keys == NULL ){
	  free( Config );
	  fprintf( stderr, "Error allocating memory\n" );
	  return NULL;
     }

     Config->Values = (char **) calloc(Size, sizeof(char*));
     if( Config->Keys == NULL ){
	  free( Config->Keys );
	  free( Config );
	  fprintf( stderr, "Error allocating memory\n" );
	  return NULL;
     }
     return Config;
}
     
int ConfFile_ReadFile( ConfFile *const CFile, const char *FileName )
{
     unsigned int Line_Num;
     short int Line_Type;
     FILE *InFile;
     size_t StrLength;
     int i;
     char Line[MAX_LENGTH];
     char Section[MAX_LENGTH];
     char Key[MAX_LENGTH];
     char Value[MAX_LENGTH];
     char EntryLine[MAX_LENGTH];

     InFile = fopen( FileName, "r" );

     if( InFile == NULL ){
	  fprintf( stderr, "Could not open %s\n.", FileName );
	  return EXIT_FAILURE;
     }

     memset( Line, 0, MAX_LENGTH );
     memset( Section, 0, MAX_LENGTH );
     memset( Key, 0, MAX_LENGTH );
     memset( Value, 0, MAX_LENGTH );
     memset( EntryLine, 0, MAX_LENGTH );

     Line_Num = 0;
     /* Read the file until the end */
     while( fgets( Line, MAX_LENGTH, InFile ) != NULL ){
	  Line_Num = Line_Num + 1;
	  StrLength = strlen( Line ) - 1;

	  if ( StrLength == 0 ){
	       /* This is an empty line. Ignore */
	       continue;
	  }

	  /* Check if the line was too long */
	  if( Line[StrLength] != '\n' && !feof(InFile) ){
	       fprintf( stderr, "Line %d was too long in file %s.\n", Line_Num, FileName );
	       return EXIT_FAILURE;
	  }

	  /* Remove the end spaces and the new line mark */
	  i = (int) StrLength;
	  while( i >= 0 && (Line[i] == '\n' || isspace( Line[i] ))){
	       Line[i] = 0;
	       i = i - 1;
	  }
	  /* Evaluate the contens of a line */
	  Line_Type = ConfFile_EvalLine( Line );

	  switch (Line_Type){
	  case CONFFILE_LEMPTY:
	       /* Do nothing */
	       break;
	  case CONFFILE_LCOMMENT:
	       /* Do nothing */
	       break;
	  case CONFFILE_LSECTION:
//	       memset( Section, 0, MAX_LENGTH );

	       sscanf( Line, "[%[^]]", Section );
	       memset( Line, 0, MAX_LENGTH );
	       /* Do nothing else since empty sections are not going to be processed */
	       break;
	  case CONFFILE_LVALUE:
	       ConfFile_ProcessLine( Line, Key, Value );
	       sprintf( EntryLine, "%s:%s", Section, Key );
	       ConfFile_SetEntry( CFile, EntryLine, Value );
	       memset( EntryLine, 0, MAX_LENGTH );
	       memset( Line, 0, MAX_LENGTH );
	       memset( Key, 0, MAX_LENGTH );
	       memset( Value, 0, MAX_LENGTH );
	       break;
	  case CONFFILE_ERROR:
	       return EXIT_FAILURE;
	       break;
	  default:
	       /* Do nothing */
	       break;
	  }
     }
     /* Close the file */
     fclose( InFile );
     return 0;
}

short int ConfFile_EvalLine( const char *Line )
{
     size_t StrLength;
     short int TypeLine;

     StrLength = strlen( Line ) - 1;
     TypeLine = CONFFILE_ERROR;

     if ( StrLength < 1 ){
	  /* This is an empty line */
	  TypeLine = CONFFILE_LEMPTY;
     } else if ( Line[0] == '#' || Line[0] == ';' ){
	  /* This line is a comment and it should be ignored */
	  TypeLine = CONFFILE_LCOMMENT;
     } else if ( Line[0] == '[' && Line[StrLength] == ']' ){
	  /* This line contains a section identifier */
	  TypeLine = CONFFILE_LSECTION;
     } else {
	  TypeLine = CONFFILE_LVALUE;
     }
     return TypeLine;
}

int ConfFile_ProcessLine( const char *Line, char *Key, char *Value )
{
     int StrLength;

     if( sscanf( Line, "%[^=] = \"%[^\"]\"", Key, Value) == 2 ||
	 sscanf( Line, "%[^=] = '%[^\']\'" , Key, Value) == 2 ||
	 sscanf( Line, "%[^=] = %[^;#]",     Key, Value) == 2 ){

	  /* Remove white spaces at the end of the line */
	  StrLength = (int) strlen(Key) - 1;
	  while( StrLength >= 0 && isspace( Key[StrLength] )){
	       Key[StrLength] = 0;
	       StrLength = StrLength - 1;
	  }
	  StrLength = (int) strlen(Value) - 1;
	  while( StrLength >= 0 && isspace( Value[StrLength] )){
	       Value[StrLength] = 0;
	       StrLength = StrLength - 1;
	  }
	  return 1;
     } else if ( sscanf(Line, "%[^=] = %[;#]", Key, Value) == 2 ||
		 sscanf(Line, "%[^=] %[=]",    Key, Value) == 2) {
	  /* There is no value */
	  Value[0] = 0;
	  return 1;
     } else return 0;
}	  

void ConfFile_SetEntry( ConfFile *const CFile, const char *Entry, const char *Value )
{
     if( CFile->NumEntries < CFile->Length ){
	  CFile->Keys[CFile->NumEntries] = strdup( Entry );
	  CFile->Values[CFile->NumEntries] = strdup( Value );

	  CFile->NumEntries = CFile->NumEntries + 1;
     } else {
	  fprintf( stderr, "Maximum number of entries reached in Configuration File.\n" );
	  exit( EXIT_FAILURE );
     }
}

char* ConfFile_GetString( const ConfFile *const CFile, const char *Key )
{
     size_t Position;

     Position = ConfFile_FindPosition( CFile, Key );
     return CFile->Values[Position];
}

int ConfFile_GetInt( const ConfFile *const CFile, const char *Key )
{
     size_t Position;

     Position = ConfFile_FindPosition( CFile, Key );
     return atoi( CFile->Values[Position] );
}

float ConfFile_GetFloat( const ConfFile *const CFile, const char *Key )
{
     size_t Position;
     
     Position = ConfFile_FindPosition( CFile, Key );
     return strtof( CFile->Values[Position], NULL );
}

double ConfFile_GetDouble( const ConfFile *const CFile, const char *Key )
{
     size_t Position;

     Position = ConfFile_FindPosition( CFile, Key );
     return atof( CFile->Values[Position] );
}

size_t ConfFile_FindPosition( const ConfFile *const CFile, const char *Key )
{
     size_t i;
     bool Found;

     Found = false;
     i = 0;
     while( i < CFile->NumEntries && !Found ){
	  if( strcmp( CFile->Keys[i], Key ) == 0 ){
	       Found = true;
	  } else {
	       i = i + 1;
	  }
     }

     if ( !Found ){
	  fprintf( stderr, "The entry \"%s\" was not found.\n", Key );
	  exit( EXIT_FAILURE );
     } else {
	  return i;
     }
}
     

void ConfFile_Free( ConfFile *const CFile )
{
     size_t i;

     /* Free the allocated memory
	Deallocate first the memory in Keys or Values */
     for ( i = 0; i < CFile->NumEntries; i++ ){
	  free( CFile->Keys[i] );
	  free( CFile->Values[i] );
     }

     free( CFile->Keys );
     free( CFile->Values );
     free( CFile );
}
