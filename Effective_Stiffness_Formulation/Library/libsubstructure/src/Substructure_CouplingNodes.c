#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include "MatrixVector.h"
#include "Algorithm_Aux.h"
#include "Print_Messages.h"
#include "Substructure_Exact.h"
#include "Substructure_Remote.h"
#include "Substructure_UHYDEfbr.h"
#include "Substructure_SimMeasured.h"
#include "Substructure_Experimental.h"
#include "Substructure_CouplingNodes.h"


const char *Substructure_Type[] = {"Sim_Exact_MDOF",
				   "Sim_Exact_SDOF",
				   "Sim_Exact_ESP",
				   "Sim_UHYDEfbr",
				   "Sim_Measured",
				   "Exp_ADwin",
				   "Remote" };


void Substructure_ReadCouplingNodes( const AlgConst_t *const InitCnt, CouplingNode_t *const CNodes )
{
     FILE *InFile;
     int Count_Type;
     int i, j, k;
     int itemp, ndof;
     double *ftemp, *rayleigh;
     char *ctemp;
     MatrixVector_t mass, stiff;
     char Type[MAX_SUBTYPE], Description[MAX_DESCRIPTION], FileMeas[MAX_FILENAME];
     char RemoteType[MAX_SUBTYPE];
     char InLine[MAX_LINE];
     char IPAddress[20], Port[20];
     double DeltaTSub;

     DeltaTSub = InitCnt->Delta_t/(double) InitCnt->NSubstep;

     InFile = fopen( InitCnt->FileCNodes, "r" );

     if( InFile == NULL ){
	  Print_Header( ERROR );
	  fprintf( stderr, "Substructure_ReadCouplingNodes: could not open file %s.\n", InitCnt->FileCNodes );
	  exit( EXIT_FAILURE );
     }

     /* The first value should be the number of Coupling nodes */
     fscanf( InFile, "%i", &CNodes->Order );
     fgets( InLine, MAX_LINE, InFile );

     if( CNodes->Order != InitCnt->OrderSub ){
	  fclose( InFile );
	  Print_Header( ERROR );
	  fprintf( stderr, "Substructure_ReadCouplingNodes: Invalid number of substructures.\n" );
	  exit( EXIT_FAILURE );
     }
	  
     /* Allocate the necessary memory */
     CNodes->Array = (int *) calloc( (size_t) CNodes->Order, sizeof(int) );
     CNodes->Sub = (Substructure_t *) malloc( (size_t) CNodes->Order*sizeof(Substructure_t) );
     CNodes->VecTdT0_c0 = (double *) calloc( (size_t) CNodes->Order, sizeof(double) );

     /* Read the contents of the file */
     i = 0;
     while( i < CNodes->Order ){

	  /* Read until the coma */
	  fscanf( InFile, "%[^,], %d", Type, &Count_Type );

	  /* Check if the number of sub-structures is still valid */
	  if( (i +  Count_Type) > CNodes->Order ){
	       fclose( InFile );
	       Print_Header( ERROR );
	       fprintf( stderr, "Substructure_ReadCouplingNodes: Number of substructures exceeded.\n" );
	       Print_Header( ERROR );
	       fprintf( stderr, "Substructure_ReadCouplingNodes: Specified %d but read %d.\n", CNodes->Order,
			i + Count_Type );
	       exit( EXIT_FAILURE );
	  }
	  Substructure_Identify( Type, &CNodes->Sub[i].Type );

	  for( j = 0; j < Count_Type; j++ ){
	       fscanf( InFile, "%d", &CNodes->Array[i + j] );
	       if( j > 0 ){
		    /* Only copy the values if j > 0 */
		    CNodes->Sub[i + j].Type = CNodes->Sub[i].Type;
	       }
	       for( k = 0; k < (i + j); k++ ){
		    if ( CNodes->Array[i + j] == CNodes->Array[k] ){
			 Print_Header( ERROR );
			 fprintf( stderr, "Substructure_ReadCouplingNodes: There is already a substructure assigned to coupling node %d.\n",
				  CNodes->Array[k]);
			 Print_Header( ERROR );
			 fprintf( stderr, "Substructre_ReadCouplingNodes: Error when reading line %d.\n",
				  i + 1 );
			 exit( EXIT_FAILURE );
		    }
	       }
	  }
	  switch (CNodes->Sub[i].Type) {
	  case SIM_EXACT_MDOF:
	       /* Ignore coma */
	       fscanf( InFile, "%*[,] %i", &itemp );
	       if ( itemp != EXACTMDOF_NUMPARAM_INIT ){
		    Print_Header( ERROR );
		    fprintf( stderr, "Wrong number of parameters for the substructue number %i of type Exact_MDOF.\n", i );
		    fprintf( stderr, "The number of init parameters should be %i\n", EXACTMDOF_NUMPARAM_INIT );
		    exit( EXIT_FAILURE );
	       } else {
		    ftemp = NULL; ctemp = NULL;
		    ftemp = (double *) calloc( (size_t) EXACTMDOF_NUMPARAM_INIT, sizeof( double ) );
		    ctemp = (char *) calloc( (size_t) 20, sizeof( char ) );
		    rayleigh = (double *) calloc ( 2, sizeof(double) );

		    /* Read the number of degrees of freedom */
		    fscanf( InFile, "%i", &ndof );
		    MatrixVector_Create( ndof, ndof, &mass );
		    MatrixVector_Create( ndof, ndof, &stiff );

		    /* Read the matrices from a MM file */
		    fscanf( InFile, "%s", ctemp );
		    MatrixVector_FromFile_MM( ctemp, &mass );
		    fscanf( InFile, "%s", ctemp );
		    MatrixVector_FromFile_MM( ctemp, &stiff );
		    fscanf( InFile, "%lf %lf", &rayleigh[0], &rayleigh[1] );

		    /* Read the optional description */
		    Substructure_GetDescription( InFile, i, Description );
		    
		    for( j = 0; j < Count_Type; j++ ){
			 CNodes->Sub[i + j].SimStruct = (void *) malloc( sizeof(ExactSim_t) );

/*			 mass.Array[0] = 600.0; mass.Array[1] = 0.0; mass.Array[2] = 0.0; mass.Array[3] = 285.0;
			 stiff.Array[0] = 150000; stiff.Array[1] = -68000.0; stiff.Array[2] = -68000.0; stiff.Array[3] = 68000.0;
			 ray[0] = 0.900; ray[1] = 0.000015;*/
			 Substructure_ExactSolutionMDOF_Init( mass.Array, stiff.Array, ndof, InitCnt->Rayleigh.Alpha, InitCnt->Rayleigh.Beta,
							      rayleigh[0], rayleigh[1], InitCnt->a0, InitCnt->a2, InitCnt->a3, InitCnt->a6,
							      InitCnt->a7, Description, (ExactSim_t *) CNodes->Sub[i + j].SimStruct );
			 
			 /*Substructure_ExactSolutionSDOF_Init( ftemp[0], ftemp[1], ftemp[2], InitCnt->a0, InitCnt->a2, InitCnt->a3, InitCnt->a6, InitCnt->a7, Description,
			   (ExactSim_t *) CNodes->Sub[i + j].SimStruct );*/
			 Print_Header( INFO );
			 printf( "Simulating the substructure in the coupling node %d as an exact integration method (MDOF).\n",
				 CNodes->Array[i + j] );
		    }
		    MatrixVector_Destroy( &mass );
		    MatrixVector_Destroy( &stiff );
		    free( rayleigh );
		    free( ftemp );
		    free( ctemp );
	       }
	       break;
	  case SIM_EXACT_SDOF:
	       /* Ignore coma */
	       fscanf( InFile, "%*[,] %i", &itemp );
	       if ( itemp != EXACTSDOF_NUMPARAM_INIT ){
		    Print_Header( ERROR );
		    fprintf( stderr, "Wrong number of parameters for the substructue number %i of type Exact_SDOF.\n", i );
		    fprintf( stderr, "The number of init parameters should be %i\n", EXACTSDOF_NUMPARAM_INIT );
		    exit( EXIT_FAILURE );
	       } else {
		    ftemp = NULL;
		    ftemp = (double *) calloc( (size_t) EXACTSDOF_NUMPARAM_INIT, sizeof( double ) );

		    /* Read the input parameters */
		    for( j = 0; j < EXACTSDOF_NUMPARAM_INIT; j++ ){
			 fscanf( InFile, "%lf", &ftemp[j] );
		    }

		    /* Read the optional description */
		    Substructure_GetDescription( InFile, i, Description );
		    
		    for( j = 0; j < Count_Type; j++ ){
			 CNodes->Sub[i + j].SimStruct = (void *) malloc( sizeof(ExactSim_t) );
			 
			 Substructure_ExactSolutionSDOF_Init( ftemp[0], ftemp[1], ftemp[2],
							      InitCnt->a0, InitCnt->a2, InitCnt->a3, InitCnt->a6,
							      InitCnt->a7, Description,
							      (ExactSim_t *) CNodes->Sub[i + j].SimStruct );

			 Print_Header( INFO );
			 printf( "Simulating the substructure in the coupling node %d as an exact integration method (SDOF).\n",
				 CNodes->Array[i + j] );
		    }
		    free( ftemp );
	       }
	       break;
	  case SIM_EXACT_ESP:
	       /* Ignore coma */
	       fscanf( InFile, "%*[,] %i", &itemp );
	       if ( itemp != EXACTSDOF_NUMPARAM_INIT ){
		    Print_Header( ERROR );
		    fprintf( stderr, "Wrong number of parameters for the substructue number %i of type Exact_ESP.\n", i );
		    fprintf( stderr, "The number of init parameters should be %i\n", EXACTSDOF_NUMPARAM_INIT );
		    exit( EXIT_FAILURE );
	       } else {
		    ftemp = NULL;
		    ftemp = (double *) calloc( (size_t) EXACTSDOF_NUMPARAM_INIT, sizeof( double ) );

		    /* Read the input parameters */
		    for( j = 0; j < EXACTSDOF_NUMPARAM_INIT; j++ ){
			 fscanf( InFile, "%lf", &ftemp[j] );
		    }

		    /* Read the optional description */
		    Substructure_GetDescription( InFile, i, Description );
		    
		    for( j = 0; j < Count_Type; j++ ){
			 CNodes->Sub[i + j].SimStruct = (void *) malloc( sizeof(ExactSim_t) );
			 
			 Substructure_ExactSolutionESP_Init( ftemp[0], ftemp[1], ftemp[2], DeltaTSub, Description, (ExactSimESP_t *) CNodes->Sub[i + j].SimStruct );

			 Print_Header( INFO );
			 printf( "Simulating the substructure in the coupling node %d as an exact integration method (ESP).\n",
				 CNodes->Array[i + j] );
		    }
		    free( ftemp );
	       }
	       break;
	  case SIM_UHYDE:
	       /* Ignore coma */
	       fscanf( InFile, "%*[,] %i", &itemp );
	       if ( itemp != UHYDE_NUMPARAM_INIT ){
		    Print_Header( ERROR );
		    fprintf( stderr, "Wrong number of parameters for the substructue number %i of type UHYDE.\n", i );
		    fprintf( stderr, "The number of init parameters should be %i\n", UHYDE_NUMPARAM_INIT );
		    exit( EXIT_FAILURE );
	       } else {
		    ftemp = NULL;
		    ftemp = (double *) calloc( (size_t) UHYDE_NUMPARAM_INIT, sizeof( double ) );

		    for( j = 0; j < UHYDE_NUMPARAM_INIT; j++ ){
			 fscanf( InFile, "%lf", &ftemp[j] );
		    }

		    /* Read the optional description */
		    Substructure_GetDescription( InFile, i, Description );

		    for( j = 0; j < Count_Type; j++ ){			
			 CNodes->Sub[i + j].SimStruct = (void *) malloc( sizeof(UHYDEfbrSim_t) );
			 Substructure_SimUHYDE_1D_Init( ftemp[0], ftemp[1], ftemp[2], Description,
							(UHYDEfbrSim_t *) CNodes->Sub[i + j].SimStruct );
			 Print_Header( INFO );
			 printf( "Simulating the substructure in the coupling node %d as a UHYDE-fbr device.\n", CNodes->Array[i + j] );
		    }
		    free( ftemp );
	       }
	       break;
	  case SIM_MEASURED:
	       fscanf( InFile, "%*[,] %[^,]", FileMeas );

	       /* Read the optional description */
	       Substructure_GetDescription( InFile, i, Description );

	       for( j = 0; j <  Count_Type; j++ ){    
		    CNodes->Sub[i + j].SimStruct = (void *) malloc( sizeof(MeasuredSim_t) );
		    Substructure_SimMeasured_Init( FileMeas, InitCnt->NStep, InitCnt->NSubstep, Description,
						   (MeasuredSim_t *) CNodes->Sub[i + j].SimStruct );
		    Print_Header( INFO );
		    printf( "Simulating the substructure in the coupling node %d using time history measured forces.\n", CNodes->Array[i + j] );
	       }
	       break;
	  case EXP_ADWIN:
	       /* Read the optional description */
	       Substructure_GetDescription( InFile, i, Description );

	       for( j = 0; j <  Count_Type; j++ ){
		    CNodes->Sub[i + j].SimStruct = (void *) malloc( sizeof(ExpSub_t) );
		    Substructure_Experimental_Init( Description, (ExpSub_t *) CNodes->Sub[i + j].SimStruct );
		    Print_Header( INFO );
		    printf( "The substructure in the coupling node %d is computed in ADwin.\n",
			    CNodes->Array[i + j] );
	       }
	       break;
	  case REMOTE:
	       /* Read IP Address and Port */
	       fscanf( InFile, "%*[,] %s %s %[^,]", RemoteType, IPAddress, Port );

	       /* Read the optional description */
	       Substructure_GetDescription( InFile, i, Description );
	       
	       for( j = 0; j <  Count_Type; j++ ){
		    CNodes->Sub[i + j].SimStruct = (void *) malloc( sizeof(Remote_t) );
		    Substructure_Remote_Init( RemoteType, IPAddress, Port, Count_Type, &CNodes->Array[i], Description,
					      (Remote_t *) CNodes->Sub[i + j].SimStruct );
		    Print_Header( INFO );
		    printf( "The substructure in the coupling node %d is computed is computed at %s:%s using %s protocol.\n", CNodes->Array[i + j],
			    IPAddress, Port, RemoteType );
	       }
	       break;
	  }
	  i = i + Count_Type;
     }

     /* Close the file */
     fclose( InFile );
    
     /* Sort the coupling nodes in ascending order */
     Substructure_SortCouplingNodes( CNodes );
}

void Substructure_SortCouplingNodes( CouplingNode_t *const CNodes )
{

     int i, j, Tmp;
     Substructure_t TmpSub;

     for( i = 0; i < CNodes->Order; i++ ){
	  for( j = 0; j < CNodes->Order - 1 ; j++ ){
	       if( CNodes->Array[j] > CNodes->Array[j+1] ){
		    Tmp = CNodes->Array[j+1];
		    TmpSub = CNodes->Sub[j+1];

		    CNodes->Array[j+1] = CNodes->Array[j];
		    CNodes->Sub[j+1] = CNodes->Sub[j];

		    CNodes->Array[j] = Tmp;
		    CNodes->Sub[j] = TmpSub;
	       }
	  }
     }
}

void Substructure_Identify( char *const Type, int *const Identity_Num )
{
     int ID; /* A counter */
     bool Found = false;

     ID = 0;
     /* Identify with substructure is in Type. Exit when found */
     while( ID < NUM_TYPE_SUB && !Found ){
	  if ( strcmp( Substructure_Type[ID], Type ) == 0 ){
	       Found = true;
	  } else {
	       ID = ID + 1;
	  }
     }

     /* The substructure in Type is not supported.*/
     if ( !Found ){
	  Print_Header( ERROR );
	  fprintf( stderr, "Substructure_Identify: The substructure type '%s' is not supported. Valid substructures are:\n", Type );
	  for( ID = 0; ID < NUM_TYPE_SUB; ID++ ){
	       fprintf( stderr, "[......] %d) %s.\n", ID+1, Substructure_Type[ID] );
	  }
	  exit( EXIT_FAILURE );
     } else {
	  /* Assign the identity */
	  (*Identity_Num) = ID;
     }
}

void Substructure_GetDescription( FILE *const InFile, const int LineNum, char *const Description )
{

     fscanf( InFile, "%*[,]" );
     fgets( Description, MAX_DESCRIPTION, InFile );
     
     if( Description[strlen(Description) - 1] != '\n'  && !feof(InFile) ){
	  Print_Header( ERROR );
	  fprintf( stderr, "Substructure_ReadCouplingNodes: Maximum description length (%d) exceeded in line %d.\n", MAX_DESCRIPTION, LineNum + 1 );
	  exit( EXIT_FAILURE );
     }

     if( Description[strlen(Description) - 2] != ';' && Description[strlen(Description) - 1] != ';' ){
	  Print_Header( ERROR );
	  fprintf( stderr, "Substructure_ReadCouplingNodes: Line number %d should terminate with ';'.\n",
		   LineNum + 1 );
	  exit( EXIT_FAILURE );
     }
}

void Substructure_DeleteCouplingNodes( CouplingNode_t *CNodes )
{
     int i;

     for( i = 0; i < CNodes->Order; i++ ){
	  switch( CNodes->Sub[i].Type ){
	  case SIM_EXACT_MDOF:
	       Substructure_ExactSolutionMDOF_Destroy( (ExactSim_t *) CNodes->Sub[i].SimStruct );
	       break;
	  case SIM_EXACT_SDOF:
	       Substructure_ExactSolutionSDOF_Destroy( (ExactSim_t *) CNodes->Sub[i].SimStruct );
	       break;
	  case SIM_EXACT_ESP:
	       Substructure_ExactSolutionESP_Destroy( (ExactSimESP_t *) CNodes->Sub[i].SimStruct );
	       break;
	  case SIM_UHYDE:
	       Substructure_SimUHYDE_Destroy( (UHYDEfbrSim_t *) CNodes->Sub[i].SimStruct );
	       break;
	  case SIM_MEASURED:
	       Substructure_SimMeasured_Destroy( (MeasuredSim_t *) CNodes->Sub[i].SimStruct );
	       break;
	  case EXP_ADWIN:
	       Substructure_Experimental_Destroy( (ExpSub_t *) CNodes->Sub[i].SimStruct );
	       break;
	  case REMOTE:
	       Substructure_Remote_Destroy( (Remote_t *) CNodes->Sub[i].SimStruct, CNodes->Order );
	       break;
	  }
	  free( CNodes->Sub[i].SimStruct );
     }

     CNodes->Order = 0;
     free( CNodes->Sub );
     free( CNodes->Array );
     free( CNodes->VecTdT0_c0 );
}
