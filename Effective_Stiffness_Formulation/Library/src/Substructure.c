#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>  /* For bool, false, true */
#include <string.h>   /* For strlen() */

#include "Substructure.h"
#include "Substructure_Exact.h"
#include "Substructure_UHYDEfbr.h"
#include "Substructure_SimMeasured.h"

#include "Auxiliary_Math.h"

#include "Print_Messages.h"
#include "Netlib.h"

const char *Substructure_Type[] = {"Sim_Exact",
				   "Sim_UHYDEfbr",
				   "Sim_Measured",
				   "Exp_ADwin",
				   "Remote_TCP",
				   "Remote_UDP",
				   "Remote_NSEP",
				   "Remote_OF" };

void Join_NonCouplingPart( MatrixVector_t *const VecTdT_m, const MatrixVector_t *const Gain_m,
			   const MatrixVector_t *const fcprevsub, const CouplingNode_t *const CNodes,
			   MatrixVector_t *const VecTdT )			  
{
     static int icoup;                 /* Counter for the coupling nodes */
     static int incx, incy;            /* Stride in the vectors */
     static double Alpha, Beta;        /* Constants for the BLAS routines */
     static char trans;                /* Use or not the transpose */
     static int Rows, Cols;            /* Number of Rows and columns */
     static int lda;                   /* Leading dimension */
     static int Length, PosX, PosXm;   /* Length and position counters */
     
     incx = 1; incy = 1;
     trans = 'N';
     Alpha = 1.0; Beta = 1.0;
     Rows = Gain_m->Rows;
     Cols = Gain_m->Cols;
     lda = Max( 1, Gain_m->Rows);

     /* Update the VecTdT_m displacments to include the effects of the coupling force */
     /* BLAS: VecTdT_m = Gain_m*fcprevsub */
     dgemv_( &trans, &Rows, &Cols, &Alpha, Gain_m->Array, &lda,
	     fcprevsub->Array, &incx, &Beta, VecTdT_m->Array, &incy );

     /* Copy the updated values into the complete displacement vector */
     PosX = 0; PosXm = 0;
     for ( icoup = 0; icoup < CNodes->Order; icoup++ ){
	  Length = CNodes->Array[icoup] - PosX -1;
	  dcopy( &Length, &VecTdT_m->Array[PosXm], &incx, &VecTdT->Array[PosX], &incy );
	  PosX = CNodes->Array[icoup];
	  PosXm = PosXm + Length;
     }

     /* Add the elements between the final coupling node and the final element
      * of the complete displacement vector */
     Length = VecTdT->Rows - CNodes->Array[CNodes->Order -1];
     dcopy( &Length, &VecTdT_m->Array[PosXm], &incx, &VecTdT->Array[PosX], &incy );	
}


void Substructure_ReadCouplingNodes( CouplingNode_t *const CNodes, const unsigned int NSteps, const unsigned int NSubsteps,
				     const int OrderSub, const double DeltaTSub, const char *Filename )
{
     FILE *InFile;
     int i, j;
     int itemp;
     double *ftemp;
     char Type[8], Description[MAX_DESCRIPTION], FileMeas[MAX_FILENAME];

     ExpSub_t *Experimental;     

     InFile = fopen( Filename, "r" );

     if( InFile == NULL ){
	  Print_Header( ERROR );
	  fprintf( stderr, "Substructure_ReadCouplingNodes: could not open file %s.\n", Filename );
	  exit( EXIT_FAILURE );
     }

     /* The first value should be the number of Coupling nodes */
     fscanf( InFile, "%i", &CNodes->Order );

     if( CNodes->Order != OrderSub ){
	  fclose( InFile );
	  Print_Header( ERROR );
	  fprintf( stderr, "Substructure_ReadCouplingNodes: Invalid number of substructures.\n" );
	  exit( EXIT_FAILURE );
     }
	  
     /* Allocate the necessary memory */
     CNodes->Array = (int *) calloc( (size_t) CNodes->Order, sizeof(int) );
     CNodes->Sub = (Substructure_t *) malloc( (size_t) CNodes->Order*sizeof(Substructure_t) );
     CNodes->u0c0 = (double *) calloc( (size_t) CNodes->Order, sizeof(double) );

     /* Read the contents of the file */
     for( i = 0; i < CNodes->Order; i++ ){
	  fscanf( InFile, "%s %i", Type, &CNodes->Array[i] );
	  
	  Substructure_Identify( Type, &CNodes->Sub[i].Type );

	  switch (CNodes->Sub[i].Type) {
	  case SIM_EXACT:
	       fscanf( InFile, "%i", &itemp );
	       if ( itemp != UHYDE_NUMPARAM_INIT ){
		    Print_Header( ERROR );
		    fprintf( stderr, "Wrong number of parameters for the substructue number %i of type TMD.\n", i );
		    fprintf( stderr, "The number of init parameters should be %i\n", EXACT_NUMPARAM_INIT );
		    exit( EXIT_FAILURE );
	       } else {
		    Print_Header( INFO );
		    printf( "Simulating the sub-structure in the coupling node %d as an exact integration method.\n", CNodes->Array[i] );
		    CNodes->Sub[i].SimStruct = (void *) malloc( sizeof(ExactSim_t) );
		    ftemp = NULL;
		    ftemp = (double *) calloc( (size_t) EXACT_NUMPARAM_INIT, sizeof( double ) );

		    for( j = 0; j < EXACT_NUMPARAM_INIT; j++ ){
			 fscanf( InFile, "%lf", &ftemp[j] );
		    }
			 
		    Substructure_ExactSolution_Init( ftemp[0], ftemp[1], ftemp[2], DeltaTSub, (ExactSim_t *) CNodes->Sub[i].SimStruct );
		    free( ftemp );
	       }
	       break;
	  case SIM_UHYDE:
	       fscanf( InFile, "%i", &itemp );
	       if ( itemp != UHYDE_NUMPARAM_INIT ){
		    Print_Header( ERROR );
		    fprintf( stderr, "Wrong number of parameters for the substructue number %i of type UHYDE.\n", i );
		    fprintf( stderr, "The number of init parameters should be %i\n", UHYDE_NUMPARAM_INIT );
		    exit( EXIT_FAILURE );
	       } else {
		    Print_Header( INFO );
		    printf( "Simulating the sub-structure in the coupling node %d as a UHYDE-fbr device.\n", CNodes->Array[i] );
		    CNodes->Sub[i].SimStruct = (void *) malloc( sizeof(UHYDEfbrSim_t) );
		    ftemp = NULL;
		    ftemp = (double *) calloc( (size_t) UHYDE_NUMPARAM_INIT, sizeof( double ) );

		    for( j = 0; j < UHYDE_NUMPARAM_INIT; j++ ){
			 fscanf( InFile, "%lf", &ftemp[j] );
		    }
			 
		    Substructure_SimUHYDE_1D_Init( ftemp[0], ftemp[1], ftemp[2], (UHYDEfbrSim_t *) CNodes->Sub[i].SimStruct );
		    free( ftemp );
	       }
	       break;
	  case SIM_MEASURED:
	       fscanf( InFile, "%s", FileMeas );
	       fgets( Description, MAX_DESCRIPTION, InFile );

	       if( Description[strlen(Description) - 1] != '\n' && !feof(InFile) ){
		    Print_Header( ERROR );
		    fprintf( stderr, "Substructure_ReadCouplingNodes: The desciption of the substructure %d is too long. Maximum number of characters %d.",
			     i + 1, MAX_DESCRIPTION );
		    exit( EXIT_FAILURE );
	       } else {
		    Print_Header( INFO );
		    printf( "Simulating the sub-structure in the coupling node %d using time history measured forces.\n", CNodes->Array[i] );
		    CNodes->Sub[i].SimStruct = (void *) malloc( sizeof(MeasuredSim_t) );
		    Substructure_SimMeasured_Init( FileMeas, NSteps, NSubsteps, (MeasuredSim_t *) CNodes->Sub[i].SimStruct );
	       }
	       break;
	  case EXP_ADWIN:
	       CNodes->Sub[i].SimStruct = (void *) malloc( sizeof(ExpSub_t) );
	       Experimental = CNodes->Sub[i].SimStruct;
	       /* Dynamic string input. Reads everything between " " */
	       fgets( Experimental->Description, MAX_DESCRIPTION, InFile );

	       if( Experimental->Description[strlen(Experimental->Description) - 1] != '\n' && !feof(InFile) ){
		    Print_Header( ERROR );
		    fprintf( stderr, "Substructure_ReadCouplingNodes: The desciption of the substructure %d is too long. Maximum characters %d.",
			     i + 1, MAX_DESCRIPTION );
		    exit( EXIT_FAILURE );
	       }
	       break;
	  case REMOTE_TCP:
	       break;
	  case REMOTE_UDP:
	       break;
	  case REMOTE_NSEP:
	       break;
	  case REMOTE_OF:
	       break;
	  }
     }
     /* Close the file */
     fclose( InFile );
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
	  fprintf( stderr, "Substructure_Identify: The substructure type '%s' is not supported. Valid operations are:\n", Type );
	  for( ID = 0; ID < NUM_TYPE_SUB; ID++ ){
	       fprintf( stderr, "[......] %d) %s.\n", ID+1, Substructure_Type[ID] );
	  }
	  exit( EXIT_FAILURE );
     } else {
	  /* Assign the identity */
	  (*Identity_Num) = ID;
     }
}

void Substructure_DeleteCouplingNodes( CouplingNode_t *CNodes )
{
     int i;

     for( i = 0; i < CNodes->Order; i++ ){
	  free( CNodes->Sub[i].SimStruct );
     }

     CNodes->Order = 0;
     free( CNodes->Sub );
     free( CNodes->Array );
     free( CNodes->u0c0 );
}


#if _1_
void Substructure_Substepping( double *const IGain, double *const DispTdT0_c, double *const DispTdT, double *const fcprevsub, double *const fc, const double Time, CouplingNode_t *const CNodes, const int NSubstep, const double DeltaT_Sub )
{

     unsigned int i;
     bool Called_Sub = false;
     double *Recv = NULL;

     Recv = (double *) calloc( (size_t) 3*(size_t)CNodes->Order, sizeof(double) );

     for( i = 0; i < CNodes->Order; i++ ){
	  switch ( CNodes->Sub[i].Type ){
	  case SIM_EXACT:
	       /* This is the same case as SIM_MEASURED. All the simulated substructures are treated together in
		* the same routine.*/
	  case SIM_UHYDE:
	       /* This is the same case as SIM_MEASURED. All the simulated substructures are treated together in
		* the same routine.*/
	  case SIM_MEASURED:
	       /* Call the Simulate_Substructures() function only once. All the simulated substructures are handled
		* together in this routine */
	       if( !Called_Sub ){
		    Simulate_Substructures( CNodes, IGain, DispTdT0_c, &Recv[0], &Recv[CNodes->Order], &Recv[2*CNodes->Order], NSubstep, DeltaT_Sub );
		    Called_Sub = true;
	       }
	       break;
#if _ADWIN_
	  case EXP_ADWIN:
	       /* Tell ADwin to perform the substepping process */
	       ADWIN_Substep( DispTdT0_c, &Recv[0], &Recv[1], &Recv[2], CNodes->Order );
	       break;
#endif
	  case REMOTE_TCP:
	       /* Using TCP communication protocol */
	       Send_Data( DispTdT0_c, CNodes->Order, Socket );

	       Receive_Data( Recv, 3*CNodes->Order, Socket );
	       break;
	  case REMOTE_UDP:
	       /* Using UDP communication protocol */

	       Send_Data( DispTdT0_c, CNodes->Order, Socket );
	       if ( recv( Socket, Recv, sizeof(double)*3*(size_t) CNodes->Order,0) != (int) sizeof(double)*3*CNodes->Order ){    /* sizeof returns an unsigned integer ? */
		    PrintErrorAndExit( "recv() failed in connected UDP mode" );
	       }
	       break;
	  case REMOTE_NSEP:
	       /* Using NSEP Protocol */
	       Communicate_With_PNSE( 1, Time, DispTdT0_c, Recv, CNodes->Order );
	       /* Receive the force from the PNSE server. WhatToDo = 2 */
	       Communicate_With_PNSE( 2, Time, DispTdT0_c, Recv, 3*CNodes->Order );
	       break;
	  case REMOTE_OF:
	       /* Using OpenFresco */
	       Communicate_With_OpenFresco( DispTdT0_c, Recv, CNodes->Order, 3 ); 
	       break;
	  }
     }

#pragma omp parallel for
     for ( i = 0; i < CNodes->Order; i++ ){
#if _MPI_

	  DispTdT[i] = Recv[i];
	  fcprevsub[i] = Recv[CNodes->Order + i];
	  fc[i] = Recv[2*CNodes->Order + i];
#else
	  DispTdT[CNodes->Array[i] - 1] = Recv[i];
	  fcprevsub[i] = Recv[CNodes->Order + i];
	  fc[CNodes->Array[i] - 1] = Recv[2*CNodes->Order + i];
#endif
     }

     free( Recv );
}


void Substructure_Simulate( CouplingNode_t *const CNodes, double *IGain, double *const u0c, double *const uc, double *const fcprev, double *const fc, const unsigned int NSubstep, const double DeltaT_Sub )
{

     unsigned int i, Substep;
     double ramp0, ramp;
     int incx = 1, incy = 1;
     double One;
     int Length;
     char uplo = 'L';
     ExactSim_t *Exact;
     UHYDEfbrSim_t *UHYDE;  
     MeasuredSim_t *Measured;

     Length = (int) CNodes->Order;
     One = 1.0;

     for ( Substep = 1; Substep <= NSubstep; Substep++ ){

	  /* Backup data so that fcprev contains always the last coupling force */
	  dcopy( &Length, fc, &incx, fcprev, &incy );
	       
	  ramp = (double) Substep / (double) NSubstep;

	  ramp0 = 1.0 - ramp;   

	  if ( CNodes->Order > 1 ){
	       dcopy( &Length, CNodes->u0c0, &incx, uc, &incy );
	       dscal( &Length, &ramp0, uc, &incx );
	       daxpy( &Length, &ramp, u0c, &incx, uc, &incy );
	       dsymv( &uplo, &Length, &One, IGain, &Length, fc, &incx, &One, uc, &incy ); 
	  } else {
	       uc[0] = ramp0*CNodes->u0c0[0] + ramp*u0c[0] + IGain[0]*fc[0];
	  }
	  
	  /* Compute the new fc */
	  for( i = 0; i < CNodes->Order; i ++ ){
	       switch( CNodes->Sub[i].Type )
	       if( CNodes->Sub[i].Type == SIM_EXACT ){
		    Exact = (ExactSim_t *) CNodes->Sub[i].SimStruct;
		    ExactSolution_SDOF( u0c[i], DeltaT_Sub, Exact, &fc[i] );
	       } else if ( CNodes->Sub[i].Type == SIM_UHYDE ){
		    UHYDE = (UHYDEfbrSim_t *) CNodes->Sub[i].SimStruct;
		    Simulate_UHYDE_1D( u0c[i], DeltaT_Sub, UHYDE, &fc[i] );
	       } else if ( CNodes->Sub[i].Type == SIM_MEASURED ){
		    Measured = (MeasuredSim_t *) CNodes->Sub[i].SimStruct;
		    Simulate_Measured( Measured, &fc[i] );
	       }
	  }
	  
     }

     /* Backup u0c */
     dcopy( &Length, u0c, &incx, CNodes->u0c0, &incy );
}
#endif
