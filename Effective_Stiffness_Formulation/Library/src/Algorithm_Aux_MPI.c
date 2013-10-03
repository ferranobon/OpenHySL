#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include "mpi.h"

#include "Algorithm_Aux.h"
#include "Conf_Parser.h"
#include "Print_Messages.h"

void Algorithm_Init_MPI( const char *FileName, AlgConst_t *const InitConst )
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
	  fprintf( stderr, "Algorithm_Init_MPI(): Invalid option for Use_Absolute_Values.\n" );
     }

     InitConst->Read_Sparse = ConfFile_GetInt( Config, "General:Read_Sparse" );
     if ( !Valid_Value( InitConst->Read_Sparse) ){
	  Error = true;
	  Print_Header( ERROR );
	  fprintf( stderr, "Algorithm_Init_MPI(): Invalid option for Read_Sparse.\n" );
     }

     InitConst->Use_Sparse = ConfFile_GetInt( Config, "General:Use_Sparse" );
     if ( InitConst->Use_Sparse != 0 ){
	  Error = true;
	  Print_Header( ERROR );
	  fprintf( stderr, "Algorithm_Init_MPI(): Sparse storage is not supported in the MPI version.\n" );
     }

     InitConst->Use_Packed = ConfFile_GetInt( Config, "General:Use_Packed" );
     if ( InitConst->Use_Packed != 0 ){
	  Error = true;
	  Print_Header( ERROR );
	  fprintf( stderr, "Algorithm_Init_MPI(): Packed storage is not supported in the MPI version.\n" );
     }

     InitConst->Read_LVector = ConfFile_GetInt( Config, "General:Read_LVector" );
     if ( !Valid_Value( InitConst->Read_LVector )){
	  Error = true;
	  Print_Header( ERROR );
	  fprintf( stderr, "Algorithm_Init_MPI(): Invalid option for Read_LVector.\n" );
     }

     /* Order of the matrices */
     InitConst->Order = ConfFile_GetInt( Config, "General:Order" );
     if ( InitConst->Order <= 0 ){
	  Error = true;
	  Print_Header( ERROR );
	  fprintf( stderr, "Algorithm_Init_MPI(): Invalid order of the matrices.\n" );
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
	  fprintf( stderr, "Algorithm_Init_MPI(): Invalid number of steps.\n" );
     }

     InitConst->Delta_t = ConfFile_GetDouble( Config, "General:Delta" );
     if ( InitConst->Delta_t <= 0.0 ){
	  Error = true;
	  Print_Header( ERROR );
	  fprintf( stderr, "Algorithm_Init_MPI(): Invalid time step.\n" );
     }

     InitConst->Scale_Factor = ConfFile_GetDouble( Config, "General:Scale_Factor" );

     /* Grid information */
     InitConst->ProcessGrid.Rows = ConfFile_GetInt( Config, "Grid:Rows" );
     InitConst->ProcessGrid.Cols = ConfFile_GetInt( Config, "Grid:Cols" );
     if ( InitConst->ProcessGrid.Rows <= 0 || InitConst->ProcessGrid.Cols <= 0 ){
	  Error = true;
	  Print_Header( ERROR );
	  fprintf( stderr, "Algorithm_Init_MPI(): Invalid grid specification (%d,%d).\n", InitConst->ProcessGrid.Rows,
		   InitConst->ProcessGrid.Cols );
     }

     /* Block information */
     InitConst->BlockSize.Rows = ConfFile_GetInt( Config, "Block_Size:Rows" );
     InitConst->BlockSize.Cols = ConfFile_GetInt( Config, "Block_Size:Cols" );
     if ( InitConst->ProcessGrid.Rows <= 0 || InitConst->ProcessGrid.Cols <= 0 ){
	  Error = true;
	  Print_Header( ERROR );
	  fprintf( stderr, "Algorithm_Init_MPI(): Invalid block size specification (%d,%d).\n", InitConst->BlockSize.Rows,
		   InitConst->BlockSize.Cols );
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
	  fprintf( stderr, "Algorithm_Init_MPI(): Invalid option for the number of sub-structures.\n" );
     }
     
     /* Number of substructures */
     InitConst->NSubstep = (unsigned int) ConfFile_GetInt( Config, "Substructure:Num_Substeps" );

     InitConst->DeltaT_Sub = InitConst->Delta_t/(double) InitConst->NSubstep;

     ConfFile_Destroy( Config );

     if( Error ){
	  Print_Header( ERROR );
	  fprintf( stderr, "Algorithm_Init_MPI(): Initialisation errors. Aborting.\n" );
	  exit( EXIT_FAILURE );
     } else {
	  Print_Header( SUCCESS );
	  printf( "Algorithm_Init_MPI(): Initialisation succcessfully completed.\n" );
     }
}

void Algorithm_BroadcastConfFile( AlgConst_t *const InitConst )
{

     /* MPI Variables */
     int rank;
     
     int i;     /* A counter */
     
     /* Setup three blocks */
     int          blockcounts[2] = {13, 21};
     MPI_Datatype types[3];
     MPI_Aint     displs[3];
     MPI_Datatype InfoFile;

     MPI_Comm_rank( MPI_COMM_WORLD, &rank );
     
     InitConst->ExcitedDOF = Algorithm_BroadcastExcitedDOF( InitConst->ExcitedDOF );

     InitConst->FileM = Algorithm_BroadcastString( InitConst->FileM );
     InitConst->FileK = Algorithm_BroadcastString( InitConst->FileK );
     InitConst->FileC = Algorithm_BroadcastString( InitConst->FileC );

     InitConst->FileLV = Algorithm_BroadcastString( InitConst->FileLV );
     InitConst->FileCNodes = Algorithm_BroadcastString( InitConst->FileCNodes );
     InitConst->FileData = Algorithm_BroadcastString( InitConst->FileData );

     InitConst->FileOutput = Algorithm_BroadcastString( InitConst->FileOutput );

     /* Initialize types and displs with addresses anof items */
     MPI_Address( &InitConst->ProcessGrid.Rows, &displs[0] );
     MPI_Address( &InitConst->Delta_t, &displs[1] );

     types[0] = MPI_INT;
     types[1] = MPI_DOUBLE;

     /* Adjust the displacement array so that the displacements are offsets from the beginning of the
      * structure */
     for (i = 1; i >=0; i--){
	  displs[i] -= displs[0];
     }

     MPI_Type_create_struct( 2, blockcounts, displs, types, &InfoFile );
     MPI_Type_commit( &InfoFile );

     MPI_Bcast( &(*InitConst), 1, InfoFile, 0, MPI_COMM_WORLD );

     Print_Header( SUCCESS );
     printf( "Successfully broadcasted the contents of the configuration file to process %d.\n", rank );
}

int* Algorithm_BroadcastExcitedDOF( int *ExcitedDOF )
{
     int rank;
     int ArrayLength;

     MPI_Comm_rank( MPI_COMM_WORLD, &rank );

     /* Retrieve the array length from process 0 */
     if( rank == 0 ){
	  if( ExcitedDOF == NULL ){
	       ArrayLength = 0;
	  } else {
	       ArrayLength = ExcitedDOF[0] + 1;
	  }
     }

     /* Broadcast the length of the array to the rest of the processes */
     MPI_Bcast( &ArrayLength, 1, MPI_INT, 0, MPI_COMM_WORLD );

     if ( ArrayLength == 0 ){
	  return NULL; /* Exit the routine if this is not initialised in rank 0 */
     } else {
	  /* Allocate the memory in the other processes */
	  if (rank != 0 ){
	       ExcitedDOF = (int *) calloc( (size_t) ArrayLength, sizeof(int) );
	  }

	  /* Broadcast the contents of the string to the rest of the processes */
	  MPI_Bcast( ExcitedDOF, ArrayLength, MPI_INT, 0, MPI_COMM_WORLD );

	  return ExcitedDOF;
     }
}

char* Algorithm_BroadcastString( char *String )
{
     int rank;
     int ArrayLength;

     MPI_Comm_rank( MPI_COMM_WORLD, &rank );

     /* Retrieve the array length from process 0 */
     if( rank == 0 ){
	  ArrayLength = (int) strlen( String ) + 1;
     }

     /* Broadcast the length of the array to the rest of the processes */
     MPI_Bcast( &ArrayLength, 1, MPI_INT, 0, MPI_COMM_WORLD );

     /* Allocate the memory in the other processes */
     if (rank != 0 ){
	  String = (char *) calloc( (size_t) ArrayLength, sizeof(char) );
     }

     /* Broadcast the contents of the string to the rest of the processes */
     MPI_Bcast( String, ArrayLength, MPI_CHAR, 0, MPI_COMM_WORLD );

     return String;
}
