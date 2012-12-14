#ifndef _CONF_PARSER_H
#define _CONF_PARSER_H

#include <string.h> /* For size_t */

typedef struct ConfigFile {
     size_t NumEntries;
     size_t Length;
     char **Keys;
     char **Values;
} ConfFile;

#define CONFFILE_LEMPTY   0
#define CONFFILE_LCOMMENT 1
#define CONFFILE_LSECTION 2
#define CONFFILE_LVALUE   3
#define CONFFILE_ERROR   -1

ConfFile* ConfFile_Create( const size_t Size );
int ConfFile_ReadFile( ConfFile *const CFile, const char *FileName );
int ConfFile_ProcessLine( const char *Line, char *Key, char *Value );
void ConfFile_SetEntry( ConfFile *const CFile, const char *Entry, const char *Value );
short int ConfFile_EvalLine( const char *Line );

char* ConfFile_GetString( const ConfFile *const CFile, const char *Key );
int ConfFile_GetInt( const ConfFile *const CFile, const char *Key );
float ConfFile_GetFloat( const ConfFile *const CFile, const char *Key );
double ConfFile_GetDouble( const ConfFile *const CFile, const char *Key );
size_t ConfFile_FindPosition( const ConfFile *const CFile, const char *Key );

void ConfFile_Free( ConfFile *const CFile );

#endif
