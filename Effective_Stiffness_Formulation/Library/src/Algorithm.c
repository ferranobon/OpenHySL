#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Algorithm.h"

void Init_Algorithm( AlgConst *const InitConst, const char *FileName )
{

     ConfFile *Config;
     
     Config = ConfFile_Create( 35 );

     ConfFile_ReadFile( Config, FileName );

     /* Use Relative or absolute values */
     (*InitConst).Use_Absolute_Values = ConfFile_GetInt( Config, "General:Use_Absolute_Values" );
     if ( InitConst->Use_Absolute_Values != 0 && InitConst->Use_Absolute_Values != 1 ){
	  PrintErrorAndExit( "Invalid option for Use_Absolute_Values" );
     }

     /* Order of the matrices */
     (*InitConst).Order = ConfFile_GetInt( Config, "General:Order" );
     if ( InitConst->Order <= 0 ){
	  PrintErrorAndExit( "Invalid option for the order of the matrices" );
     }

     /* Number of steps and Time step */
     (*InitConst).Nstep = (unsigned int) ConfFile_GetInt( Config, "General:Num_Steps" );
     if ( InitConst->Nstep <= 0 ){
	  PrintErrorAndExit( "Invalid number of steps" );
     }

     (*InitConst).Delta_t = ConfFile_GetFloat( Config, "General:Delta" );
     if ( InitConst->Delta_t <= 0.0f ){
	  PrintErrorAndExit( "Invalid time step" );
     }

     /* Rayleigh values */
     (*InitConst).Rayleigh.Alpha = ConfFile_GetFloat( Config, "Rayleigh:Alpha" );
     (*InitConst).Rayleigh.Beta = ConfFile_GetFloat( Config, "Rayleigh:Beta" );

     /* Newmark integration constants */
     (*InitConst).Newmark.Gamma = ConfFile_GetFloat( Config, "Newmark:Gamma" );
     (*InitConst).Newmark.Beta = ConfFile_GetFloat( Config, "Newmark:Beta" );

     /* PID Constants */
     (*InitConst).PID.P = ConfFile_GetFloat( Config, "PID:P" );
     (*InitConst).PID.I = ConfFile_GetFloat( Config, "PID:I" );
     (*InitConst).PID.D = ConfFile_GetFloat( Config, "PID:D" );

     /* Several constants to multiply the vectors */
     (*InitConst).Const1 = (*InitConst).Newmark.Beta*(*InitConst).Delta_t*(*InitConst).Delta_t;
     (*InitConst).Const2 = (0.5f - 2.0f*(*InitConst).Newmark.Beta + (*InitConst).Newmark.Gamma)*(*InitConst).Delta_t*(*InitConst).Delta_t;
     (*InitConst).Const3 = (0.5f + (*InitConst).Newmark.Beta - (*InitConst).Newmark.Gamma)*(*InitConst).Delta_t*(*InitConst).Delta_t;

     /* Constants for Ending Step */
     (*InitConst).a0 = 1.0f/((*InitConst).Newmark.Beta*(*InitConst).Delta_t*(*InitConst).Delta_t);
     (*InitConst).a1 = (*InitConst).Newmark.Gamma/((*InitConst).Newmark.Beta*(*InitConst).Delta_t);
     (*InitConst).a2 = 1.0f/((*InitConst).Newmark.Beta*(*InitConst).Delta_t);
     (*InitConst).a3 = 1.0f/(2.0f*(*InitConst).Newmark.Beta) - 1.0f;
     (*InitConst).a4 = (*InitConst).Newmark.Gamma/(*InitConst).Newmark.Beta - 1.0f;
     (*InitConst).a5 = ((*InitConst).Delta_t/2.0f)*((*InitConst).Newmark.Gamma/(*InitConst).Newmark.Beta - 2.0f);
     (*InitConst).a6 = (1.0f - (*InitConst).Newmark.Gamma)*(*InitConst).Delta_t;
     (*InitConst).a7 = (*InitConst).Newmark.Gamma*(*InitConst).Delta_t;

     /* File Names */
/*EFAST*/
     (*InitConst).FileM = strdup( ConfFile_GetString( Config, "FileNames:Mass_Matrix" ) );
     (*InitConst).FileK = strdup( ConfFile_GetString( Config, "FileNames:Stiffness_Matrix" ) );
     (*InitConst).FileC = strdup( ConfFile_GetString( Config, "FileNames:Damping_Matrix" ) );
     (*InitConst).FileLVector = strdup( ConfFile_GetString( Config, "FileNames:FileLVector" ) );
     (*InitConst).FileCNodes = strdup( ConfFile_GetString( Config, "FileNames:Coupling_Nodes" ) );
     (*InitConst).FileData = strdup( ConfFile_GetString( Config, "FileNames:Ground_Motion" ) );
     (*InitConst).FileOutput = strdup( ConfFile_GetString( Config, "FileNames:OutputFile" ) );

     /* Get the communication protocol to be used */
     GetServerInformation( &InitConst->Remote, Config );
     ConfFile_Free( Config );
}

void Delete_Algorithm( AlgConst *const InitConst )
{

     free( InitConst->FileM );
     free( InitConst->FileK );
     if( InitConst->FileC != NULL ){
	  free( InitConst->FileC );
     }
     free( InitConst->FileLVector );
     free( InitConst->FileCNodes );
     free( InitConst->FileData );
     free( InitConst->FileOutput );

     Delete_ServerInformation( &InitConst->Remote );
}
