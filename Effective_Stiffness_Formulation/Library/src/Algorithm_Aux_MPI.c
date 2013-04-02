#include <stdio.h>
#include <stdlib.h>

#include "mpi.h"

#include "Algorithm_Aux.h"
#include "Conf_Parser.h"
#include "Print_Messages.h"


void Algorithm_Init_MPI( const char *FileName, AlgConst_t *const InitConst )
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

     /* Grid information */
     InitConst->ProcessGrid.Rows = ConfFile_GetInt( Config, "Grid:Rows" );
     InitConst->ProcessGrid.Cols = ConfFile_GetInt( Config, "Grid:Cols" );
     if ( InitConst->ProcessGrid.Rows <= 0 || InitConst->ProcessGrid.Cols <= 0 ){
	  Print_Header( ERROR );
	  fprintf( stderr, "Invalid grid specification (%d,%d).\n", InitConst->ProcessGrid.Rows,
		   InitConst->ProcessGrid.Cols );
	  exit( EXIT_FAILURE );
     }

     /* Block information */
     InitConst->BlockSize.Rows = ConfFile_GetInt( Config, "Block_Size:Rows" );
     InitConst->BlockSize.Cols = ConfFile_GetInt( Config, "Block_Size:Cols" );
     if ( InitConst->ProcessGrid.Rows <= 0 || InitConst->ProcessGrid.Cols <= 0 ){
	  Print_Header( ERROR );
	  fprintf( stderr, "Invalid block size specification (%d,%d).\n", InitConst->BlockSize.Rows,
		   InitConst->BlockSize.Cols );
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
	  fprintf( stderr, "Invalid option for the number of sub-structuresr of the matrices.\n" );
	  exit( EXIT_FAILURE );
     }
     
     /* Number of substructures */
     InitConst->NSubstep = (unsigned int) ConfFile_GetInt( Config, "Substructure:Num_Substeps" );

     InitConst->DeltaT_Sub = InitConst->Delta_t/(double) InitConst->NSubstep;

     ConfFile_Destroy( Config );

     Print_Header( SUCCESS );
     printf( "Initialisation succcessfully completed.\n" );
}


void Algorithm_BroadcastConfFile( AlgConst_t *const InitConst )
{

     /* MPI Variables */
     int rank;
     
     size_t LengthArrays;
     int i;     /* A counter */
     
     /* Setup three blocks */
     int          blockcounts[3] = {14, 19, 0};
     MPI_Datatype types[3];
     MPI_Aint     displs[3];
     MPI_Datatype InfoFile;

     MPI_Comm_rank( MPI_COMM_WORLD, &rank );
     
     if ( rank == 0 ){
	  LengthArrays = strlen( InitConst->FileM ) + 1;
	  LengthArrays = LengthArrays + strlen( InitConst->FileK ) + 1;
	  LengthArrays = LengthArrays + strlen( InitConst->FileC ) + 1;
	  LengthArrays = LengthArrays + strlen( InitConst->FileLV ) + 1;
	  LengthArrays = LengthArrays + strlen( InitConst->FileCNodes ) + 1;
	  LengthArrays = LengthArrays + strlen( InitConst->FileData ) + 1;
	  LengthArrays = LengthArrays + strlen( InitConst->FileOutput ) + 1;
     }

     MPI_Bcast( &LengthArrays, 1, MPI_INT, 0, MPI_COMM_WORLD );

     blockcounts[2] = (int) LengthArrays;

     /* Initialize types and displs with addresses anof items */
     MPI_Address( &InitConst->ProcessGrid, &displs[0] );
     MPI_Address( &InitConst->Delta_t,   &displs[1] );
     MPI_Address( &InitConst->FileM, &displs[2] );

     types[0] = MPI_INT;
     types[1] = MPI_DOUBLE;
     types[2] = MPI_CHAR;

     /* Adjust the displacement array so that the displacements are offsets from the beginning of the
      * structure */
     for (i = 2; i >=0; i--){
	  displs[i] -= displs[0];
     }

     MPI_Type_create_struct( 3, blockcounts, displs, types, &InfoFile );
     MPI_Type_commit( &InfoFile );

     MPI_Bcast( &(*InitConst), 1, InfoFile, 0, MPI_COMM_WORLD );

     Print_Header( SUCCESS );
     printf( "Successfully broadcasted the contents of the configuration file to the rest of the MPI processes.\n" );
}
