#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <stdarg.h>  /* For va_arg( ), ... */  

#include "Algorithm_Aux.h"
#include "Conf_Parser.h"
#include "Print_Messages.h"

#include "Definitions.h"

void Algorithm_Init( const char *FileName, AlgConst_t *const InitConst )
{

     ConfFile_t *Config;
     bool Error = false;

     Config = ConfFile_Create( 70 );

     ConfFile_ReadFile( FileName, Config );

     /* Use Relative or absolute values */
     InitConst->Use_Absolute_Values = ConfFile_GetInt( Config, "General:Use_Absolute_Values" );
     if ( !Valid_Value( InitConst->Use_Absolute_Values ) ){
	  Error = true;
	  Print_Header( ERROR );
	  fprintf( stderr, "Algorithm_Init(): Invalid option for Use_Absolute_Values.\n" );
     }

     InitConst->Read_Sparse = ConfFile_GetInt( Config, "General:Read_Sparse" );
     if ( !Valid_Value( InitConst->Read_Sparse) ){
	  Error = true;
	  Print_Header( ERROR );
	  fprintf( stderr, "Algorithm_Init(): Invalid option for Read_Sparse.\n" );
     }

     InitConst->Use_Sparse = ConfFile_GetInt( Config, "General:Use_Sparse" );
#if _SPARSE_
     if ( !Valid_Value( InitConst->Use_Sparse ) ){
	  Error = true;
	  Print_Header( ERROR );
	  fprintf( stderr, "Algorithm_Init(): Invalid option for Use_Sparse.\n" );
     }
#else
     if ( InitConst->Use_Sparse != 0 ){
	  Error = true;
	  Print_Header( ERROR );
	  fprintf( stderr, "Algorithm_Init(): The algorithm has been compiled without sparse support. ");
	  fprintf( stderr, "Please, set Use_Sparse variable to 0 in the configuration file.\n" );
     }
#endif

     InitConst->Use_Packed = ConfFile_GetInt( Config, "General:Use_Packed" );
     if ( !Valid_Value( InitConst->Use_Packed) ){
	  Error = true;
	  Print_Header( ERROR );
	  fprintf( stderr, "Algorithm_Init(): Invalid option for Use_Packed.\n" );
     }

     InitConst->Read_LVector = ConfFile_GetInt( Config, "General:Read_LVector" );
     if ( !Valid_Value( InitConst->Read_LVector )){
	  Error = true;
	  Print_Header( ERROR );
	  fprintf( stderr, "Algorithm_Init(): Invalid option for Read_LVector.\n" );
     }

     /* Order of the matrices */
     InitConst->Order = ConfFile_GetInt( Config, "General:Order" );
     if ( InitConst->Order <= 0 ){
	  Error = true;
	  Print_Header( ERROR );
	  fprintf( stderr, "Algorithm_Init(): Invalid order of the matrices.\n" );
     }

     /* Do nothing if there is also an error */
     if( !InitConst->Read_LVector && !Error ){
	  InitConst->ExcitedDOF = Algorithm_GetExcitedDOF( Config, "General:Excited_DOF" );
     } else {
	  InitConst->ExcitedDOF = NULL;
     }

     /* Read the damping matrix or use rayleigh */
     InitConst->Read_CMatrix = ConfFile_GetInt( Config, "General:Read_CMatrix" );
     if ( !Valid_Value( InitConst->Read_CMatrix )){
	  Error = true;
	  Print_Header( ERROR );
	  fprintf( stderr, "Algorithm_Init(): Invalid option for Read_CMatrix.\n" );
     }

     /* Number of steps and Time step */
     InitConst->NStep = (unsigned int) ConfFile_GetInt( Config, "General:Num_Steps" );
     if ( InitConst->NStep <= 0 ){
	  Error = true;
	  Print_Header( ERROR );
	  fprintf( stderr, "Algorithm_Init(): Invalid number of steps.\n" );
     }

#if _FLOAT_
     InitConst->Delta_t = ConfFile_GetFloat( Config, "General:Delta" );
#else
     InitConst->Delta_t = ConfFile_GetDouble( Config, "General:Delta" );
#endif
     if ( InitConst->Delta_t <= 0.0 ){
	  Error = true;
	  Print_Header( ERROR );
	  fprintf( stderr, "Algorithm_Init(): Invalid time step.\n" );
     }

#if _FLOAT_
     InitConst->Scale_Factor = ConfFile_GetFloat( Config, "General:Scale_Factor" );
#else
     InitConst->Scale_Factor = ConfFile_GetDouble( Config, "General:Scale_Factor" );
#endif

     /* Rayleigh values */
     if( !InitConst->Read_CMatrix ){
#if _FLOAT_
	  InitConst->Rayleigh.Alpha = ConfFile_GetFloat( Config, "Rayleigh:Alpha" );
	  InitConst->Rayleigh.Beta = ConfFile_GetFloat( Config, "Rayleigh:Beta" );
#else
	  InitConst->Rayleigh.Alpha = ConfFile_GetDouble( Config, "Rayleigh:Alpha" );
	  InitConst->Rayleigh.Beta = ConfFile_GetDouble( Config, "Rayleigh:Beta" );
#endif
     }
	  

     /* Newmark integration constants */
#if _FLOAT_
     InitConst->Newmark.Gamma = ConfFile_GetFloat( Config, "Newmark:Gamma" );
     InitConst->Newmark.Beta = ConfFile_GetFloat( Config, "Newmark:Beta" );
#else
     InitConst->Newmark.Gamma = ConfFile_GetDouble( Config, "Newmark:Gamma" );
     InitConst->Newmark.Beta = ConfFile_GetDouble( Config, "Newmark:Beta" );
#endif

     /* PID Constants */
#if _FLOAT_
     InitConst->PID.P = ConfFile_GetFloat( Config, "PID:P" );
     InitConst->PID.I = ConfFile_GetFloat( Config, "PID:I" );
     InitConst->PID.D = ConfFile_GetFloat( Config, "PID:D" );
#else
     InitConst->PID.P = ConfFile_GetDouble( Config, "PID:P" );
     InitConst->PID.I = ConfFile_GetDouble( Config, "PID:I" );
     InitConst->PID.D = ConfFile_GetDouble( Config, "PID:D" );
#endif

     /* Constants for Ending Step */
#if _FLOAT_
     InitConst->a0 = 1.0f/(InitConst->Newmark.Beta*InitConst->Delta_t*InitConst->Delta_t);
     InitConst->a2 = 1.0f/(InitConst->Newmark.Beta*InitConst->Delta_t);
     InitConst->a3 = 1.0f/(2.0f*InitConst->Newmark.Beta) - 1.0f;
     InitConst->a4 = InitConst->Newmark.Gamma/InitConst->Newmark.Beta - 1.0f;
     InitConst->a5 = (InitConst->Delta_t/2.0f)*(InitConst->Newmark.Gamma/InitConst->Newmark.Beta - 2.0f);
     InitConst->a6 = (1.0f - InitConst->Newmark.Gamma)*InitConst->Delta_t;
     InitConst->a10 = (0.5f - InitConst->Newmark.Beta)*InitConst->Delta_t*InitConst->Delta_t;
     InitConst->a16 = (1.0f - 2.0f*InitConst->Newmark.Gamma)*InitConst->Delta_t;
     InitConst->a17 = (0.5f -2.0f*InitConst->Newmark.Beta + InitConst->Newmark.Gamma)*InitConst->Delta_t*InitConst->Delta_t;
     InitConst->a18 = (0.5f + InitConst->Newmark.Beta - InitConst->Newmark.Gamma)*InitConst->Delta_t*InitConst->Delta_t;
#else
     InitConst->a0 = 1.0/(InitConst->Newmark.Beta*InitConst->Delta_t*InitConst->Delta_t);
     InitConst->a2 = 1.0/(InitConst->Newmark.Beta*InitConst->Delta_t);
     InitConst->a3 = 1.0/(2.0*InitConst->Newmark.Beta) - 1.0;
     InitConst->a4 = InitConst->Newmark.Gamma/InitConst->Newmark.Beta - 1.0;
     InitConst->a5 = (InitConst->Delta_t/2.0)*(InitConst->Newmark.Gamma/InitConst->Newmark.Beta - 2.0);
     InitConst->a6 = (1.0 - InitConst->Newmark.Gamma)*InitConst->Delta_t;
     InitConst->a10 = (0.5 - InitConst->Newmark.Beta)*InitConst->Delta_t*InitConst->Delta_t;
     InitConst->a16 = (1.0 - 2.0*InitConst->Newmark.Gamma)*InitConst->Delta_t;
     InitConst->a17 = (0.5 -2.0*InitConst->Newmark.Beta + InitConst->Newmark.Gamma)*InitConst->Delta_t*InitConst->Delta_t;
     InitConst->a18 = (0.5 + InitConst->Newmark.Beta - InitConst->Newmark.Gamma)*InitConst->Delta_t*InitConst->Delta_t;
#endif

     InitConst->a1 = InitConst->Newmark.Gamma/(InitConst->Newmark.Beta*InitConst->Delta_t);
     InitConst->a7 = InitConst->Newmark.Gamma*InitConst->Delta_t;
     InitConst->a8 = InitConst->Newmark.Beta*InitConst->Delta_t*InitConst->Delta_t;
     InitConst->a9 = InitConst->Delta_t;

     /* File Names */
     InitConst->FileM = strdup( ConfFile_GetString( Config, "FileNames:Mass_Matrix" ) );
     if( !Valid_File( InitConst->FileM ) ){
	  Error = true;
	  Print_Header( ERROR );
	  fprintf( stderr, "%s: No such file or directory.\n", InitConst->FileM );
     }

     InitConst->FileK = strdup( ConfFile_GetString( Config, "FileNames:Stiffness_Matrix" ) );
     if( !Valid_File( InitConst->FileK ) ){
	  Error = true;
	  Print_Header( ERROR );
	  fprintf( stderr, "%s: No such file or directory.\n", InitConst->FileK );
     }

     if( InitConst->Read_CMatrix ){
	  InitConst->FileC = strdup( ConfFile_GetString( Config, "FileNames:Damping_Matrix" ) );
	  if( !Valid_File( InitConst->FileC ) ){
	       Error = true;
	       Print_Header( ERROR );
	       fprintf( stderr, "%s: No such file or directory.\n", InitConst->FileC );
	  }
     } else {
	  InitConst->FileC = NULL;
     }

     if( InitConst->Read_LVector ){
	  InitConst->FileLV1 = strdup( ConfFile_GetString( Config, "FileNames:Load_Vector1" ) );
	  if( !Valid_File( InitConst->FileLV1 ) ){
	       Error = true;
	       Print_Header( ERROR );
	       fprintf( stderr, "%s: No such file or directory.\n", InitConst->FileLV1 );
	  }
	  InitConst->FileLV2 = strdup( ConfFile_GetString( Config, "FileNames:Load_Vector2" ) );
	  if( !Valid_File( InitConst->FileLV2 ) ){
	       Error = true;
	       Print_Header( ERROR );
	       fprintf( stderr, "%s: No such file or directory.\n", InitConst->FileLV2 );
	  }

	  InitConst->FileLV3 = strdup( ConfFile_GetString( Config, "FileNames:Load_Vector3" ) );
	  if( !Valid_File( InitConst->FileLV3 ) ){
	       Error = true;
	       Print_Header( ERROR );
	       fprintf( stderr, "%s: No such file or directory.\n", InitConst->FileLV3 );
	  }
     } else {
	  InitConst->FileLV1 = NULL;
	  InitConst->FileLV2 = NULL;
	  InitConst->FileLV3 = NULL;
     }

     InitConst->FileCNodes = strdup( ConfFile_GetString( Config, "FileNames:Coupling_Nodes" ) );
     if( !Valid_File( InitConst->FileCNodes ) ){
	  Error = true;
	  Print_Header( ERROR );
	  fprintf( stderr, "%s: No such file or directory.\n", InitConst->FileCNodes);
     }

     InitConst->FileData1 = strdup( ConfFile_GetString( Config, "FileNames:Ground_Motion1" ) );
     if( !Valid_File( InitConst->FileData1 ) ){
	  Error = true;
	  Print_Header( ERROR );
	  fprintf( stderr, "%s: No such file or directory.\n", InitConst->FileData1 );
     }

     InitConst->FileData2 = strdup( ConfFile_GetString( Config, "FileNames:Ground_Motion2" ) );
     if( !Valid_File( InitConst->FileData2 ) ){
	  Error = true;
	  Print_Header( ERROR );
	  fprintf( stderr, "%s: No such file or directory.\n", InitConst->FileData2 );
     }

     InitConst->FileData3 = strdup( ConfFile_GetString( Config, "FileNames:Ground_Motion3" ) );
     if( !Valid_File( InitConst->FileData3 ) ){
	  Error = true;
	  Print_Header( ERROR );
	  fprintf( stderr, "%s: No such file or directory.\n", InitConst->FileData3 );
     }

     /* Since this is a write operation, a warning should be issued and the filename should be changed so that
      * it does not get overwriten. */
     if( Valid_File( Concatenate_Strings( 2, ConfFile_GetString( Config, "FileNames:OutputFile" ), ".h5") ) ){
	  Print_Header( WARNING );
	  fprintf( stderr, "Output data file %s would have been overwritten. ", Concatenate_Strings( 2, ConfFile_GetString( Config, "FileNames:OutputFile" ), ".h5") );
	  
	  InitConst->FileOutput = Change_Filename( ConfFile_GetString( Config, "FileNames:OutputFile" ) );

	  fprintf( stderr, "Renaming it to: %s\n", Concatenate_Strings( 2, InitConst->FileOutput, ".h5") );
     } else {
	  InitConst->FileOutput = strdup( ConfFile_GetString( Config, "FileNames:OutputFile" ) );
     }

     /* Read the information regarding the numerical sub-structures */

     /* Number of substructures */
     InitConst->OrderSub = ConfFile_GetInt( Config, "Substructure:Order" );
     if ( InitConst->OrderSub < 0 ){
	  Error = true;
	  Print_Header( ERROR );
	  fprintf( stderr, "Algorithm_Init(): Invalid option for the number of sub-structures.\n" );
     }
     
     /* Number of substructures */
     InitConst->NSubstep = (unsigned int) ConfFile_GetInt( Config, "Substructure:Num_Substeps" );

     InitConst->DeltaT_Sub = InitConst->Delta_t/(HYSL_FLOAT) InitConst->NSubstep;

     ConfFile_Destroy( Config );

     if( Error ){
	  Print_Header( ERROR );
	  fprintf( stderr, "Algorithm_Init(): Initialisation errors. Aborting.\n" );
	  exit( EXIT_FAILURE );
     } else {
	  Print_Header( SUCCESS );
	  printf( "Algorithm_Init(): Initialisation succcessfully completed.\n" );
     }
}

bool Valid_Value( const int Value )
{
     if ( Value != 0 && Value != 1 ){
	  return false;
     } else {
	  return true;
     }
}

bool Valid_File( const char *Filename )
{
     FILE *TheFile = NULL;

     TheFile = fopen( Filename, "r" );

     if ( TheFile == NULL ){
	  return false;
     } else {
	  fclose( TheFile );
	  return true;
     }
}

char* Change_Filename( const char *const Name )
{

     char *NewName = NULL;
     char HelpChar[5];
     short unsigned int i;

     i = 1;
     do {
	  if( NewName != NULL ){
	       free( NewName );
	  }
	  sprintf( HelpChar, "_%hu", i );
	  i = i + (short unsigned int) 1;
	  NewName = Concatenate_Strings( 2, Name, HelpChar );
     } while( Valid_File( Concatenate_Strings( 2, NewName, ".h5")));

     return NewName;
}

void Algorithm_Destroy( AlgConst_t *const InitConst )
{

     if( InitConst->ExcitedDOF != NULL ){
	  free( InitConst->ExcitedDOF );
     }
     free( InitConst->FileM );
     free( InitConst->FileK );
     
     if( InitConst->FileC != NULL ){
	  free( InitConst->FileC );
     }
     if( InitConst->FileLV1 != NULL ){
	  free( InitConst->FileLV1 );
     }

     if( InitConst->FileLV2 != NULL ){
	  free( InitConst->FileLV2 );
     }

     if( InitConst->FileLV3 != NULL ){
	  free( InitConst->FileLV3 );
     }
     
     free( InitConst->FileCNodes );

     free( InitConst->FileData1 );
     free( InitConst->FileData2 );
     free( InitConst->FileData3 );

     free( InitConst->FileOutput );
}

int* Algorithm_GetExcitedDOF( const ConfFile_t *const Config, const char *Expression )
{
     unsigned int i, j;
     int *DOF_Table;  
     char *FullString;
     char Temp[1];

     FullString = strdup( ConfFile_GetString( Config, Expression ) );

     /* The first position contains the number of degrees of Freedom per node present in the
      * structure */
     strncpy( Temp, &FullString[0], (size_t) 1 );

     DOF_Table = (int *) calloc( (size_t) (atoi(Temp)+1), sizeof(int) );
     DOF_Table[0] = atoi( Temp );

     j = 1;
     for( i = 1; i < strlen( FullString ); i++ ){
	  if ( FullString[i] != ' ' ){
	       strncpy( Temp, &FullString[i], (size_t) 1 );
	       DOF_Table[j] = atoi( Temp );
	       j = j + 1;
	  }
     }

     free( FullString );
     return DOF_Table;
}

void Algorithm_ReadDataEarthquake( const unsigned int NumSteps, const char *Filename,
				   const HYSL_FLOAT Scale_Factor, HYSL_FLOAT *const Acceleration,
				   HYSL_FLOAT *const Velocity, HYSL_FLOAT *const Displacement )
{

     unsigned int i;		    /* A counter */
     HYSL_FLOAT unnecessary;	    /* Variable to store unnecessary data */
     HYSL_FLOAT temp1, temp2, temp3;
     FILE *InFile;

     InFile = fopen( Filename, "r" );

     if ( InFile == NULL ){
	  Print_Header( ERROR );
	  fprintf( stderr, "The earthquake data cannot be read because it was not possible to open %s.\n",
		   Filename );
	  exit( EXIT_FAILURE );
     }

     for ( i = 0; i < NumSteps; i++ ){
#if _FLOAT_
	  fscanf( InFile, "%E %E %E %E", &unnecessary, &temp1, &temp2, &temp3 );
#else
	  fscanf( InFile, "%lE %lE %lE %lE", &unnecessary, &temp1, &temp2, &temp3 );
#endif
	  Acceleration[i] = temp1*Scale_Factor;
	  Velocity[i] = temp2*Scale_Factor;
	  Displacement[i] = temp3*Scale_Factor;
     }

     /* Close File */
     fclose( InFile );
}

void Algorithm_PrintHelp( const char *Program_Name )
{

     fprintf( stderr, "Usage: %s [-h] -c <ConfFile>\n", Program_Name );
     fprintf( stderr,
	      "  -h  --help           This help text.\n"
	      "  -c  --config-file    The name of the configuration file. Default value: ConfFile.conf\n" );
}

char* Concatenate_Strings( int count, ... )
{
     va_list ap;
     int i;
     char *merged, *s;
     size_t len = 1, null_pos; /* Space for NULL */

     /* Find required length to store merged string */
     va_start( ap, count );

     for( i = 0; i < count; i++ )
	  len += strlen( va_arg( ap, char* ) );
     va_end( ap );

     /* Allocate memory to concat strings */
     merged = calloc( sizeof(char), len );
     null_pos = 0;

     /* Actually concatenate strings */
     va_start( ap, count );

     for( i = 0; i < count; i++ )
     {
	  s = va_arg( ap, char* );
	  strcpy( merged + null_pos, s );
	  null_pos += strlen( s );
     }
     va_end( ap );

     return merged;
}
