#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include "Algorithm_Aux.h"
#include "Conf_Parser.h"
#include "Print_Messages.h"

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

     /* Number of steps and Time step */
     InitConst->NStep = (unsigned int) ConfFile_GetInt( Config, "General:Num_Steps" );
     if ( InitConst->NStep <= 0 ){
	  Error = true;
	  Print_Header( ERROR );
	  fprintf( stderr, "Algorithm_Init(): Invalid number of steps.\n" );
     }

     InitConst->Delta_t = ConfFile_GetDouble( Config, "General:Delta" );
     if ( InitConst->Delta_t <= 0.0 ){
	  Error = true;
	  Print_Header( ERROR );
	  fprintf( stderr, "Algorithm_Init(): Invalid time step.\n" );
     }

     InitConst->Scale_Factor = ConfFile_GetDouble( Config, "General:Scale_Factor" );

     /* Rayleigh values */
     InitConst->Rayleigh.Alpha = ConfFile_GetDouble( Config, "Rayleigh:Alpha" );
     InitConst->Rayleigh.Beta = ConfFile_GetDouble( Config, "Rayleigh:Beta" );

     /* Newmark integration constants */
     InitConst->Newmark.Gamma = ConfFile_GetDouble( Config, "Newmark:Gamma" );
     InitConst->Newmark.Beta = ConfFile_GetDouble( Config, "Newmark:Beta" );

     /* PID Constants */
     InitConst->PID.P = ConfFile_GetDouble( Config, "PID:P" );
     InitConst->PID.I = ConfFile_GetDouble( Config, "PID:I" );
     InitConst->PID.D = ConfFile_GetDouble( Config, "PID:D" );

     /* Constants for Ending Step */
     InitConst->a0 = 1.0/(InitConst->Newmark.Beta*InitConst->Delta_t*InitConst->Delta_t);
     InitConst->a1 = InitConst->Newmark.Gamma/(InitConst->Newmark.Beta*InitConst->Delta_t);
     InitConst->a2 = 1.0/(InitConst->Newmark.Beta*InitConst->Delta_t);
     InitConst->a3 = 1.0/(2.0*InitConst->Newmark.Beta) - 1.0;
     InitConst->a4 = InitConst->Newmark.Gamma/InitConst->Newmark.Beta - 1.0;
     InitConst->a5 = (InitConst->Delta_t/2.0)*(InitConst->Newmark.Gamma/InitConst->Newmark.Beta - 2.0);
     InitConst->a6 = (1.0 - InitConst->Newmark.Gamma)*InitConst->Delta_t;
     InitConst->a7 = InitConst->Newmark.Gamma*InitConst->Delta_t;
     InitConst->a8 = InitConst->Newmark.Beta*InitConst->Delta_t*InitConst->Delta_t;
     InitConst->a9 = InitConst->Delta_t;
     InitConst->a10 = (0.5 - InitConst->Newmark.Beta)*InitConst->Delta_t*InitConst->Delta_t;

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

     InitConst->FileC = strdup( ConfFile_GetString( Config, "FileNames:Damping_Matrix" ) );

     if( InitConst->Read_LVector ){
	  InitConst->FileLV = strdup( ConfFile_GetString( Config, "FileNames:Load_Vector" ) );
	  if( !Valid_File( InitConst->FileLV ) ){
	       Error = true;
	       Print_Header( ERROR );
	       fprintf( stderr, "%s: No such file or directory.\n", InitConst->FileLV );
	  }
     } else {
	  InitConst->FileLV = NULL;
     }
     InitConst->FileCNodes = strdup( ConfFile_GetString( Config, "FileNames:Coupling_Nodes" ) );
     if( !Valid_File( InitConst->FileCNodes ) ){
	  Error = true;
	  Print_Header( ERROR );
	  fprintf( stderr, "%s: No such file or directory.\n", InitConst->FileCNodes);
     }

     InitConst->FileData = strdup( ConfFile_GetString( Config, "FileNames:Ground_Motion" ) );
     if( !Valid_File( InitConst->FileData ) ){
	  Error = true;
	  Print_Header( ERROR );
	  fprintf( stderr, "%s: No such file or directory.\n", InitConst->FileData );
     }

     /* Since this is a write operation, a warning should be issued and the filename should be changed so that
      * it does not overwrite. */
     InitConst->FileOutput = strdup( ConfFile_GetString( Config, "FileNames:OutputFile" ) );
     if( Valid_File( InitConst->FileOutput ) ){
	  Print_Header( WARNING );
	  fprintf( stderr, "Output data file %s would have been overwritten. ", InitConst->FileOutput );
	  
	  Change_Filename( InitConst->FileOutput );

	  fprintf( stderr, "Renaming it to: %s\n", InitConst->FileOutput );

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

     InitConst->DeltaT_Sub = InitConst->Delta_t/(double) InitConst->NSubstep;

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

void Change_Filename( char *Name )
{

     char NewName[80];
     char HelpChar[4];
     char Extension[6];
     unsigned int i;
     size_t Pos;

     Pos = strcspn( Name, "." );

     i = 0;
     while( Name[Pos + i] != '\0' && i < 6){
	  Extension[i] = Name[Pos + i];
	  i = i + 1;
     }

     strcpy( NewName, Name );
     i = 1;
     while( Valid_File( NewName )){
	  /* Reset the name */
	  memset(NewName, 0, sizeof(NewName));
	  sprintf( HelpChar, "%hu", i );
	  strncpy( NewName, Name, Pos );
	  strcat( NewName, "_" );
	  strcat( NewName, HelpChar );	  
	  strcat( NewName, Extension );
	  i = i + 1;
     }

     free( Name );
     Name =  strdup( NewName );
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
     if( InitConst->FileLV != NULL ){
	  free( InitConst->FileLV );
     }
     free( InitConst->FileCNodes );
     free( InitConst->FileData );
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

void Algorithm_ReadDataEarthquake_AbsValues( const unsigned int NumSteps, const char *Filename,
					     const double Scale_Factor, double *const Velocity,
					     double *const Displacement )
{

     unsigned int i;		    /* A counter */
     double unnecessary;	    /* Variable to store unnecessary data */
     double temp1, temp2, temp3;
     FILE *InFile;

     InFile = fopen( Filename, "r" );

     if ( InFile == NULL ){
	  Print_Header( ERROR );
	  fprintf( stderr, "The earthquake data cannot be read because it was not possible to open %s.\n",
		   Filename );
	  exit( EXIT_FAILURE );
     }

     for ( i = 0; i < NumSteps; i++ ){
	  fscanf( InFile, "%lE %lE %lE %lE", &unnecessary, &temp1, &temp2, &temp3 );
	  Velocity[i] = temp2*Scale_Factor;
	  Displacement[i] = temp3*Scale_Factor;
     }

     /* Close File */
     fclose( InFile );
}

void Algorithm_ReadDataEarthquake_RelValues( const unsigned int NumSteps, const char *Filename, 
					     const double Scale_Factor, double *const Acceleration )
{

     unsigned int i;		    /* A counter */
     double unnecessary;	    /* Variable to store unnecessary data */
     double temp1, temp2, temp3;
     FILE *InFile;

     InFile = fopen( Filename, "r" );

     if ( InFile == NULL ){
	  Print_Header( ERROR );
	  fprintf( stderr, "The earthquake data cannot be read because it was not possible to open %s.\n",
		   Filename );
	  exit( EXIT_FAILURE );
     }

     for ( i = 0; i < NumSteps; i++ ){
//	  fscanf( InFile, "%lE %lE", &unnecessary, &temp1 );
	  fscanf( InFile, "%lE %lE %lE %lE", &unnecessary, &temp1, &temp2, &temp3 );
	  Acceleration[i] = temp1*Scale_Factor;
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
