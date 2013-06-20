#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Algorithm_Aux.h"
#include "Conf_Parser.h"
#include "Print_Messages.h"

void Algorithm_Init( const char *FileName, AlgConst_t *const InitConst )
{

     ConfFile_t *Config;
     
     Config = ConfFile_Create( 70 );

     ConfFile_ReadFile( FileName, Config );

     /* Use Relative or absolute values */
     InitConst->Use_Absolute_Values = ConfFile_GetInt( Config, "General:Use_Absolute_Values" );
     if ( InitConst->Use_Absolute_Values != 0 && InitConst->Use_Absolute_Values != 1 ){
	  Print_Header( ERROR );
	  fprintf( stderr, "Invalid option for Use_Absolute_Values.\n" );
	  exit( EXIT_FAILURE );
     }

     InitConst->Read_Sparse = ConfFile_GetInt( Config, "General:Read_Sparse" );
     if ( InitConst->Read_Sparse != 0 && InitConst->Read_Sparse != 1 ){
	  Print_Header( ERROR );
	  fprintf( stderr, "Invalid option for Read_Sparse.\n" );
     }

     InitConst->Use_Sparse = ConfFile_GetInt( Config, "General:Use_Sparse" );
     if ( InitConst->Use_Sparse != 0 && InitConst->Use_Sparse != 1 ){
	  Print_Header( ERROR );
	  fprintf( stderr, "Invalid option for Use_Sparse.\n" );
     }

     InitConst->Use_Packed = ConfFile_GetInt( Config, "General:Use_Packed" );
     if ( InitConst->Use_Packed != 0 && InitConst->Use_Packed != 1 ){
	  Print_Header( ERROR );
	  fprintf( stderr, "Invalid option for Use_Packed.\n" );
     }

     InitConst->Read_LVector = ConfFile_GetInt( Config, "General:Read_LVector" );
     if ( InitConst->Read_LVector != 0 && InitConst->Read_LVector != 1 ){
	  Print_Header( ERROR );
	  fprintf( stderr, "Invalid option for Read_LVector.\n" );
     }

     /* Order of the matrices */
     InitConst->Order = ConfFile_GetInt( Config, "General:Order" );
     if ( InitConst->Order <= 0 ){
	  Print_Header( ERROR );
	  fprintf( stderr, "Invalid option for the order of the matrices.\n" );
	  exit( EXIT_FAILURE );
     }

     if( !InitConst->Read_LVector ){
	  InitConst->ExcitedDOF = Algorithm_GetExcitedDOF( Config, "General:Excited_DOF" );
     } else {
	  InitConst->ExcitedDOF = NULL;
     }

     /* Number of steps and Time step */
     InitConst->NStep = (unsigned int) ConfFile_GetInt( Config, "General:Num_Steps" );
     if ( InitConst->NStep <= 0 ){
	  Print_Header( ERROR );
	  fprintf( stderr, "Invalid number of steps.\n" );
	  exit( EXIT_FAILURE );
     }

     InitConst->Delta_t = ConfFile_GetDouble( Config, "General:Delta" );
     if ( InitConst->Delta_t <= 0.0 ){
	  Print_Header( ERROR );
	  fprintf( stderr, "Invalid time step.\n" );
	  exit( EXIT_FAILURE );
     }

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

     /* Several constants to multiply the vectors */
     InitConst->Const1 = InitConst->Newmark.Beta*InitConst->Delta_t*InitConst->Delta_t;
     InitConst->Const2 = (0.5 - 2.0*InitConst->Newmark.Beta + InitConst->Newmark.Gamma)*InitConst->Delta_t*InitConst->Delta_t;
     InitConst->Const3 = (0.5 + InitConst->Newmark.Beta - InitConst->Newmark.Gamma)*InitConst->Delta_t*InitConst->Delta_t;

     /* Constants for Ending Step */
     InitConst->a0 = 1.0/(InitConst->Newmark.Beta*InitConst->Delta_t*InitConst->Delta_t);
     InitConst->a1 = InitConst->Newmark.Gamma/(InitConst->Newmark.Beta*InitConst->Delta_t);
     InitConst->a2 = 1.0/(InitConst->Newmark.Beta*InitConst->Delta_t);
     InitConst->a3 = 1.0/(2.0*InitConst->Newmark.Beta) - 1.0;
     InitConst->a4 = InitConst->Newmark.Gamma/InitConst->Newmark.Beta - 1.0;
     InitConst->a5 = (InitConst->Delta_t/2.0)*(InitConst->Newmark.Gamma/InitConst->Newmark.Beta - 2.0);
     InitConst->a6 = (1.0 - InitConst->Newmark.Gamma)*InitConst->Delta_t;
     InitConst->a7 = InitConst->Newmark.Gamma*InitConst->Delta_t;

     /* File Names */
     InitConst->FileM = strdup( ConfFile_GetString( Config, "FileNames:Mass_Matrix" ) );
     InitConst->FileK = strdup( ConfFile_GetString( Config, "FileNames:Stiffness_Matrix" ) );
     InitConst->FileC = strdup( ConfFile_GetString( Config, "FileNames:Damping_Matrix" ) );
     if( InitConst->Read_LVector ){
	  InitConst->FileLV = strdup( ConfFile_GetString( Config, "FileNames:Load_Vector" ) );

     } else {
	  InitConst->FileLV = NULL;
     }
     InitConst->FileCNodes = strdup( ConfFile_GetString( Config, "FileNames:Coupling_Nodes" ) );
     InitConst->FileData = strdup( ConfFile_GetString( Config, "FileNames:Ground_Motion" ) );
     InitConst->FileOutput = strdup( ConfFile_GetString( Config, "FileNames:OutputFile" ) );

     /* Read the information regarding the numerical sub-structures */

     /* Number of substructures */
     InitConst->OrderSub = ConfFile_GetInt( Config, "Substructure:Order" );
     if ( InitConst->OrderSub < 0 ){
	  Print_Header( ERROR );
	  fprintf( stderr, "Invalid option for the number of sub-structures of the matrices.\n" );
	  exit( EXIT_FAILURE );
     }
     
     /* Number of substructures */
     InitConst->NSubstep = (unsigned int) ConfFile_GetInt( Config, "Substructure:Num_Substeps" );

     InitConst->DeltaT_Sub = InitConst->Delta_t/(double) InitConst->NSubstep;

     ConfFile_Destroy( Config );

     Print_Header( SUCCESS );
     printf( "Initialisation succcessfully completed.\n" );
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
					     double *Velocity, double *Displacement )
{

     unsigned int i;					/* A counter */
     double unnecessary;		/* Variable to store unnecessary data */
     double temp1, temp2, temp3;
     FILE *InFile;

     InFile = fopen( Filename, "r" );

     if ( InFile == NULL ){
	  Print_Header( ERROR );
	  fprintf( stderr, "The earthquake data cannot be read because it was not possible to open %s.\n",
		   Filename );
     }

     for ( i = 0; i < NumSteps; i++ ){
	  fscanf( InFile, "%lE %lE %lE %lE", &unnecessary, &temp1, &temp2, &temp3 );
	  Velocity[i] = temp2/1000.0;
	  Displacement[i] = temp3/1000.0;
     }

     /* Close File */
     fclose( InFile );
}

void Algorithm_ReadDataEarthquake_RelValues( const unsigned int NumSteps, const char *Filename,
					     double *Acceleration )
{

     unsigned int i;					/* A counter */
     double unnecessary;		/* Variable to store unnecessary data */
     double temp1, temp2, temp3;
     FILE *InFile;

     InFile = fopen( Filename, "r" );

     if ( InFile == NULL ){
	  Print_Header( ERROR );
	  fprintf( stderr, "The earthquake data cannot be read because it was not possible to open %s.\n", Filename );
     }

     for ( i = 0; i < NumSteps; i++ ){
	  fscanf( InFile, "%lE %lE", &unnecessary, &temp1 );
//	  fscanf( InFile, "%lE %lE %lE %lE", &unnecessary, &temp1, &temp2, &temp3 );
	  Acceleration[i] = temp1/1000.0;
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
