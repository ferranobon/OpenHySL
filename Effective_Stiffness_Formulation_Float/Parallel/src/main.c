#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#include "PMatrixVector.h"
#include "Initiation.h"
#include "Precalculations.h"
#include "ComputeU0.h"
#include "EndingStep.h"
#include "ErrorHandling.h"
#include "Netlib.h"
#include "Send_Receive_Data.h"

int main( int argc, char **argv )
{

     /* Output file */
     FILE *OutputFile;

     /* MPI Variables */
     int rank, size;
     MPI_Status Status;
     /* BLACS Variables */
     int icontxt, bhandle;
     int nprow, npcol, myrow, mycol;

     /* NETLIB Variables */
     int ione = 1, i, istep, aux;
     Scalars Constants;

     /* Variables required by the Substructure Algorithm */
     AlgConst InitCnt;

     float *AccAll, *VelAll, *DispAll;

     /* Matrices and Vectors */
     PMatrixVector M, C, K;               /* Mass, Damping and Stiffness matrices */
     PMatrixVector Keinv;

     PMatrixVector Keinv_m;

     PMatrixVector EffT;

     PMatrixVector DispT, DispTdT0, DispTdT;
     PMatrixVector DispTdT0_m;
     PMatrixVector tempvec;

     PMatrixVector VelT, VelTdT;
     PMatrixVector AccT, AccTdT;

     PMatrixVector LoadVectorForm, LoadTdT, LoadTdT1;

     PMatrixVector fc, fcprevsub;
     PMatrixVector fu;

     PMatrixVector Disp, Vel, Acc;

     int *LRowIndex_Coupling, *LColIndex_Coupling, *RowProcess_Coupling, *ColProcess_Coupling;  /* Variables to store the information regarding where the coupling position is in the fc vector and */
     int *LRowIndex_Coupling_fcp, *LColIndex_Coupling_fcp, *RowProcess_Coupling_fcp, *ColProcess_Coupling_fcp;  /* Variables to store the information regarding where the coupling position is in the fc vector and */

     float *Keinv_c;
     float *DispTdT0_c;
     float *Recv;

     double TimeElapsed;
     double TimeElapsedEnd;

     Coupling_Node CNodes;

     time_t clock;

     /* TCP socket connection variables */
     int Socket;

     MPI_Init( &argc, &argv );
     MPI_Comm_size( MPI_COMM_WORLD, &size );
     MPI_Comm_rank( MPI_COMM_WORLD, &rank );

     if ( rank == 0 ){
	  /* TODO: Review */
	  InitConstants( &InitCnt, size );
	  
	  /* Read the coupling nodes from a file */
	  Read_Coupling_Nodes( &CNodes, InitCnt.FileCNodes );
     }

     /* Send the information read in process 0 to all other processes using MPI_Bcast */
     BroadcastConfFile( &InitCnt );
     BroadCast_Coupling_Nodes( &CNodes );

     if( InitCnt.Use_Absolute_Values ){
	  AccAll == NULL;
	  VelAll = (float *) calloc( (size_t) InitCnt.Nstep, sizeof(float) );
	  DispAll = (float *) calloc( (size_t) InitCnt.Nstep, sizeof(float) );
     } else {
	  AccAll = (float *) calloc( (size_t) InitCnt.Nstep, sizeof(float) );
	  VelAll = NULL;
	  DispAll = NULL;
     }

     LRowIndex_Coupling = (int *) calloc( (size_t) CNodes.Order, sizeof(int) );
     LColIndex_Coupling = (int *) calloc( (size_t) CNodes.Order, sizeof(int) );
     RowProcess_Coupling = (int *) calloc( (size_t) CNodes.Order, sizeof(int) );
     ColProcess_Coupling = (int *) calloc( (size_t) CNodes.Order, sizeof(int) );

     LRowIndex_Coupling_fcp = (int *) calloc( (size_t) CNodes.Order, sizeof(int) );
     LColIndex_Coupling_fcp = (int *) calloc( (size_t) CNodes.Order, sizeof(int) );
     RowProcess_Coupling_fcp = (int *) calloc( (size_t) CNodes.Order, sizeof(int) );
     ColProcess_Coupling_fcp = (int *) calloc( (size_t) CNodes.Order, sizeof(int) );

     /* BLACS: Initialise BLACS */
     bhandle = Csys2blacs_handle( MPI_COMM_WORLD );
     icontxt = bhandle;

     /* BLACS: Initialise the process grid */
     Cblacs_gridinit( &icontxt, (char *) "Row", InitCnt.ProcessGrid.Rows, InitCnt.ProcessGrid.Cols );
     Cblacs_gridinfo( icontxt, &nprow, &npcol, &myrow, &mycol );

     CreateDistMatrix( icontxt, &M, InitCnt.Order, InitCnt.Order,
		       InitCnt.BlockSize.Rows, InitCnt.BlockSize.Cols );
     CreateDistMatrix( icontxt, &C, InitCnt.Order, InitCnt.Order,
		       InitCnt.BlockSize.Rows, InitCnt.BlockSize.Cols );
     CreateDistMatrix( icontxt, &K, InitCnt.Order, InitCnt.Order,
		       InitCnt.BlockSize.Rows, InitCnt.BlockSize.Cols );

     CreateDistMatrix( icontxt, &Keinv, InitCnt.Order, InitCnt.Order,
		       InitCnt.BlockSize.Rows, InitCnt.BlockSize.Cols );

     CreateDistMatrix( icontxt, &Keinv_m, InitCnt.Order - CNodes.Order,
		       CNodes.Order, InitCnt.BlockSize.Rows,
		       InitCnt.BlockSize.Cols );
     
     CreateDistMatrix( icontxt, &tempvec, InitCnt.Order, 1,
		       InitCnt.BlockSize.Rows, 1 );

     CreateDistMatrix( icontxt, &EffT, InitCnt.Order, 1,
		       InitCnt.BlockSize.Rows, 1 );

     CreateDistMatrix( icontxt, &DispT, InitCnt.Order, 1,
		       InitCnt.BlockSize.Rows, 1 );
     CreateDistMatrix( icontxt, &DispTdT0, InitCnt.Order, 1,
		       InitCnt.BlockSize.Rows, 1 );
     CreateDistMatrix( icontxt, &DispTdT, InitCnt.Order, 1,
		       InitCnt.BlockSize.Rows, 1 );

     CreateDistMatrix( icontxt, &DispTdT0_m, InitCnt.Order - CNodes.Order,
		       1, InitCnt.BlockSize.Rows, 1 );

     CreateDistMatrix( icontxt, &VelT, InitCnt.Order, 1,
		       InitCnt.BlockSize.Rows, 1 );
     CreateDistMatrix( icontxt, &VelTdT, InitCnt.Order, 1,
		       InitCnt.BlockSize.Rows, 1 );

     CreateDistMatrix( icontxt, &AccT, InitCnt.Order, 1,
		       InitCnt.BlockSize.Rows, 1 );
     CreateDistMatrix( icontxt, &AccTdT, InitCnt.Order, 1,
		       InitCnt.BlockSize.Rows, 1 );

     CreateDistMatrix( icontxt, &LoadTdT, InitCnt.Order, 1,
		       InitCnt.BlockSize.Rows, 1 );
     CreateDistMatrix( icontxt, &LoadTdT1, InitCnt.Order, 1,
		       InitCnt.BlockSize.Rows, 1 );
     CreateDistMatrix( icontxt, &LoadVectorForm, InitCnt.Order, 1,
		       InitCnt.BlockSize.Rows, 1 );

     CreateDistMatrix( icontxt, &fc, InitCnt.Order, 1,
		       InitCnt.BlockSize.Rows, 1 );
     CreateDistMatrix( icontxt, &fcprevsub, CNodes.Order, 1,
		       InitCnt.BlockSize.Rows, 1 );
     CreateDistMatrix( icontxt, &fu, InitCnt.Order, 1,
		       InitCnt.BlockSize.Rows, 1 );

     CreateDistMatrix( icontxt, &Disp, InitCnt.Order, 1,
		       InitCnt.BlockSize.Rows, 1 );
     CreateDistMatrix( icontxt, &Vel, InitCnt.Order, 1,
		       InitCnt.BlockSize.Rows, 1 );
     CreateDistMatrix( icontxt, &Acc, InitCnt.Order, 1,
		       InitCnt.BlockSize.Rows, 1 );

     if ( rank == 0 ){
	  Keinv_c = (float *) calloc( (size_t) CNodes.Order*(size_t) CNodes.Order, sizeof(float) );
	  DispTdT0_c = (float *) calloc( (size_t) CNodes.Order, sizeof(float) );
	  Recv = (float *) calloc( (size_t) (3*CNodes.Order), sizeof(float) );
     } else {
	  Keinv_c = NULL;
	  DispTdT0_c = NULL;
	  Recv = NULL;
     }

     /***********************************************************
      **************   BEGIN of INITIATION PHASE    *************
      ************** Function codes can be found in *************
      **************         Initiation.cpp         *************
      **********************************************************/

     /* Read Matrices M and K */
     DistMatrixFromFile( &M, InitCnt.FileM );
     DistMatrixFromFile( &K, InitCnt.FileK );
     DistMatrixFromFile( &LoadVectorForm, InitCnt.FileLVector ); 

     /* Calculate the Damping matrix */
     CalculateMatrixC( &M, &K, &C, InitCnt.Rayleigh );
   
     /* Calculate Matrix Keinv = [K + a0*M + a1*C]^(-1) */
     Constants.Alpha = 1.0f;
     Constants.Beta = InitCnt.a0;
     Constants.Gamma = InitCnt.a1;
     CalculateMatrixKeinv( &Keinv, &M, &C, &K, Constants );

     /* TODO: Check Order Column/Row */
     BuildMatrixXc( MPI_COMM_WORLD, &Keinv, Keinv_c, &CNodes );
     BuildMatrixXcm( MPI_COMM_WORLD, &Keinv, &Keinv_m, &CNodes );

     if ( rank == 0 ){
	  Send_Effective_Matrix( Keinv_c, (unsigned int) CNodes.Order, &Socket, InitCnt.Remote );
     }

     /*
      * Get the information about the position of the coupling force in the vector fc and its
      * coordinates in the process grid. Used when receiving the information from the server
      */
     for( i = 0; i < CNodes.Order; i++ ){
	  aux = i + 1;
	  infog2l_( &CNodes.Array[i], &ione, fc.Desc, &nprow, &npcol, &myrow, &mycol, &LRowIndex_Coupling[i], &LColIndex_Coupling[i], &RowProcess_Coupling[i], &ColProcess_Coupling[i] );
	  infog2l_( &aux, &ione, fcprevsub.Desc, &nprow, &npcol, &myrow, &mycol, &LRowIndex_Coupling_fcp[i], &LColIndex_Coupling_fcp[i], &RowProcess_Coupling_fcp[i], &ColProcess_Coupling_fcp[i] );
     }

     /***********************************************************
      **************    END of INITIATION PHASE     *************
      **********************************************************/

     /***********************************************************
      ************** BEGIN of PRECALCULATIONS PHASE *************
      ************** Function codes can be found in *************
      **************       Precalculations.cpp      *************
      **********************************************************/

     /* Read the earthquake data from a file and distribute it to all processes */
     if ( rank == 0 ){
	  if( InitCnt.Use_Absolute_Values ){
	       ReadDataEarthquake_AbsValues( VelAll, DispAll, (unsigned int ) InitCnt.Nstep, InitCnt.FileData );
	  } else {
	       ReadDataEarthquake_RelValues( AccAll, (unsigned int) InitCnt.Nstep, InitCnt.FileData );
	  }
     }

     if ( InitCnt.Use_Absolute_Values ){
	  MPI_Bcast( DispAll, InitCnt.Nstep, MPI_FLOAT, 0, MPI_COMM_WORLD );
	  MPI_Bcast( VelAll, InitCnt.Nstep, MPI_FLOAT, 0, MPI_COMM_WORLD );
     } else {
	  MPI_Bcast( AccAll, InitCnt.Nstep, MPI_FLOAT, 0, MPI_COMM_WORLD );
     }
     
     
     if ( rank == Cblacs_pnum ( fc.Desc[1], RowProcess_Coupling[CNodes.Order-1], ColProcess_Coupling[CNodes.Order-1] ) ){
	  OutputFile = fopen( InitCnt.FileOutput, "w" );
	  if ( OutputFile == NULL ){
	       PrintErrorAndExit( "Cannot proceed because the file Out.txt could not be opened" );
	  } else {
	       clock = time (NULL);	  	       
	       fprintf( OutputFile, "Test started at %s", ctime( &clock ) );
	       fprintf( OutputFile, "Number of DOF: %d, ", InitCnt.Order );
	       fprintf( OutputFile, "Number of Steps: %d, Time step: %f, Number of substeps: %d, P (PID): %f\n",
			InitCnt.Nstep, InitCnt.Delta_t, 4, InitCnt.PID.P );
	       fprintf( OutputFile, "li\t ai1(m/s^2)\t ai(m/s^2)\t vi1 (m/s)\t vi (m/s)\t ui1 (m)\t ui (m)\t fc (N)\t fu(N)\n" );

	  }

     }

     istep = 1;

     /* Calculate the input load */ 
     if( InitCnt.Use_Absolute_Values ){
	  Apply_LoadVectorForm( &Disp, &LoadVectorForm, DispAll[istep - 1] );
	  Apply_LoadVectorForm( &Vel, &LoadVectorForm, VelAll[istep - 1] );
	  Calc_Input_Load_AbsValues( &LoadTdT, &K, &C, &Disp, &Vel );
     } else {
	  Apply_LoadVectorForm( &Acc, &LoadVectorForm, AccAll[istep -1] );
	  Calc_Input_Load_RelValues( &LoadTdT, &M, &Acc );
     }

     if ( rank == 0 ){
	  printf( "Starting stepping process.\n" );
     }

     TimeElapsed = MPI_Wtime( );

     while ( istep <= InitCnt.Nstep ){

	  /* Calculate the effective force vector
	     Fe = M*(a0*u + a2*v + a3*a) + C*(a1*u + a4*v + a5*a) */
	  EffK_Calc_Effective_Force( &M, &C, &DispT, &VelT, &AccT, &tempvec,
				     InitCnt.a0, InitCnt.a1, InitCnt.a2, InitCnt.a3, InitCnt.a4, InitCnt.a5,
				     &EffT );

	  /* Compute the new Displacement u0 */
	  EffK_ComputeU0( &EffT, &LoadTdT, &fu, InitCnt.PID.P, &Keinv, &tempvec, &DispTdT0 );

	  /* Split DispTdT into coupling and non-coupling part */
	  CreateVectorXm( MPI_COMM_WORLD, &DispTdT0, &DispTdT0_m, &CNodes );
	  CreateVectorXc( MPI_COMM_WORLD, &DispTdT0, DispTdT0_c, &CNodes );

	  if ( rank == 0 ){
	       /* Perform substepping */
	       Do_Substepping( DispTdT0_c, &Recv[0], &Recv[CNodes.Order], &Recv[CNodes.Order*2], InitCnt.Remote.Type,
			       InitCnt.Delta_t*(float) istep, Socket, (unsigned int) CNodes.Order, (unsigned int *) CNodes.Array  );
	  }

	  /* Send the received information from the server to the respective vectors */
	  for( i = 0; i < CNodes.Order; i++ ){
	       if ( rank == 0 && Cblacs_pnum( icontxt, RowProcess_Coupling[i], ColProcess_Coupling[i] ) == 0 ){
		    DispTdT.Array[(LRowIndex_Coupling[i] - 1) + DispTdT.LocalSize.Row*(LColIndex_Coupling[i] -1)] = Recv[i];    
		    fc.Array[(LRowIndex_Coupling[i] - 1 ) + fc.LocalSize.Row*(LColIndex_Coupling[i] - 1)] = Recv[2*CNodes.Order + i];
	       } else {
		    if ( rank == 0 ){
			 MPI_Send( &Recv[i], CNodes.Order, MPI_FLOAT,
				   Cblacs_pnum( DispTdT.Desc[1], RowProcess_Coupling[i], ColProcess_Coupling[i] ),
				   2, MPI_COMM_WORLD );
			 MPI_Send( &Recv[2*CNodes.Order + i], CNodes.Order, MPI_FLOAT,
				   Cblacs_pnum( fc.Desc[1], RowProcess_Coupling[i], ColProcess_Coupling[i] ),
				   11, MPI_COMM_WORLD );
		    } else if ( rank == Cblacs_pnum( icontxt, RowProcess_Coupling[i], ColProcess_Coupling[i] ) ){
			 MPI_Recv( &DispTdT.Array[(LRowIndex_Coupling[i] - 1) + DispTdT.LocalSize.Row*(LColIndex_Coupling[i] -1)],
				   CNodes.Order, MPI_FLOAT, 0, 2, MPI_COMM_WORLD, &Status );
			 MPI_Recv( &fc.Array[(LRowIndex_Coupling[i] - 1 ) + fc.LocalSize.Row*(LColIndex_Coupling[i] - 1)],
				   CNodes.Order, MPI_FLOAT, 0, 11, MPI_COMM_WORLD, &Status );
		    }
	       }

	       if ( rank == 0 && Cblacs_pnum( icontxt, RowProcess_Coupling_fcp[i], ColProcess_Coupling_fcp[i] ) == 0 ){

		    fcprevsub.Array[(LRowIndex_Coupling_fcp[i] - 1) + fcprevsub.LocalSize.Row*(LColIndex_Coupling_fcp[i] -1)] = Recv[CNodes.Order +i];

	       } else {
		    if( rank == 0 ){
			 MPI_Send( &Recv[CNodes.Order + i], CNodes.Order, MPI_FLOAT,
				   Cblacs_pnum( fcprevsub.Desc[1], RowProcess_Coupling_fcp[i], ColProcess_Coupling_fcp[i] ),
				   10, MPI_COMM_WORLD );
		    } else if (rank == Cblacs_pnum( icontxt, RowProcess_Coupling_fcp[i], ColProcess_Coupling_fcp[i] ) ){
			 MPI_Recv( &fcprevsub.Array[(LRowIndex_Coupling_fcp[i] - 1) + fcprevsub.LocalSize.Row*(LColIndex_Coupling_fcp[i] -1)],
				   CNodes.Order, MPI_FLOAT, 0, 10, MPI_COMM_WORLD, &Status );
		    }
	       }
	  }

	  if ( istep < InitCnt.Nstep ){
	       /* Calculate the input load for the next step during the
		  sub-stepping process */ 
	       if( InitCnt.Use_Absolute_Values ){
		    Apply_LoadVectorForm( &Disp, &LoadVectorForm, DispAll[istep] );
		    Apply_LoadVectorForm( &Vel, &LoadVectorForm, VelAll[istep] );
		    Calc_Input_Load_AbsValues( &LoadTdT1, &K, &C, &Disp, &Vel );	       
	       } else {
		    Apply_LoadVectorForm( &Acc, &LoadVectorForm, AccAll[istep] );
		    Calc_Input_Load_RelValues( &LoadTdT1, &M, &Acc );
	       }
	  }
	  /* Join the non-coupling part. DispTdT_m = Keinv_m*fc + DispTdT0_m. Although DispTdT0_m is what has been received from the other computer,
	     it has the name of DispTdT_m to avoid further operations if using the NETLIB libraries. */
	  JoinNonCouplingPart( &DispTdT0_m, &Keinv_m, &fcprevsub, &DispTdT, &CNodes );

	  /* Compute acceleration ai1 = a0*(ui1 -ui) - a2*vi -a3*ai */
	  Compute_Acceleration( &DispTdT, &DispT, &VelT, &AccT, InitCnt.a0, InitCnt.a2, InitCnt.a3, &AccTdT );

	  /* Compute Velocity. vi = vi + a6*ai +a7*ai */
	  Compute_Velocity( &VelT, &AccT, &AccTdT, InitCnt.a6, InitCnt.a7, &VelTdT );

	  /* Error Compensation. fu = LoatTdT + fc -(Mass*AccTdT + Damp*VelTdT + Stiff*DispTdT) */
	  Compute_Force_Error( &M, &C, &K, &AccTdT, &VelTdT, &DispTdT, &fc, &LoadTdT, &fu );

	  if ( rank == Cblacs_pnum ( fc.Desc[1], RowProcess_Coupling[CNodes.Order-1], ColProcess_Coupling[CNodes.Order-1] ) ){

	       fprintf( OutputFile, "%E\t", LoadTdT.Array[(LRowIndex_Coupling[CNodes.Order-1] - 1 ) + DispTdT.LocalSize.Row*(LColIndex_Coupling[CNodes.Order-1] - 1)]);
	       fprintf( OutputFile, "%E\t", AccTdT.Array[(LRowIndex_Coupling[CNodes.Order-1] - 1 ) + DispTdT.LocalSize.Row*(LColIndex_Coupling[CNodes.Order-1] - 1)]);
	       fprintf( OutputFile, "%E\t", AccT.Array[(LRowIndex_Coupling[CNodes.Order-1] - 1 ) + DispTdT.LocalSize.Row*(LColIndex_Coupling[CNodes.Order-1] - 1)]);
	       fprintf( OutputFile, "%E\t",VelTdT.Array[(LRowIndex_Coupling[CNodes.Order-1] - 1 ) + DispTdT.LocalSize.Row*(LColIndex_Coupling[CNodes.Order-1] - 1)]);
	       fprintf( OutputFile, "%E\t",VelT.Array[(LRowIndex_Coupling[CNodes.Order-1] - 1 ) + DispTdT.LocalSize.Row*(LColIndex_Coupling[CNodes.Order-1] - 1)]);
	       fprintf( OutputFile, "%E\t", DispTdT.Array[(LRowIndex_Coupling[CNodes.Order-1] - 1 ) + DispTdT.LocalSize.Row*(LColIndex_Coupling[CNodes.Order-1] - 1)]);
	       fprintf( OutputFile, "%E\t", DispT.Array[(LRowIndex_Coupling[CNodes.Order-1] - 1 ) + DispTdT.LocalSize.Row*(LColIndex_Coupling[CNodes.Order-1] - 1)]);
	       fprintf( OutputFile, "%E\t", fc.Array[(LRowIndex_Coupling[CNodes.Order-1] - 1 ) + fc.LocalSize.Row*(LColIndex_Coupling[CNodes.Order-1] - 1)] );
	       fprintf( OutputFile, "%E\n", fu.Array[ (LRowIndex_Coupling[CNodes.Order-1] - 1 ) + fu.LocalSize.Row*(LColIndex_Coupling[CNodes.Order-1] - 1)] );

	  }
	  pscopy_( &LoadTdT.GlobalSize.Row, LoadTdT1.Array, &ione, &ione, LoadTdT1.Desc, &ione, LoadTdT.Array, &ione, &ione, LoadTdT.Desc, &ione ); /* ui = ui1 */
	  /* Backup vectors */
	  pscopy_( &DispTdT.GlobalSize.Row, DispTdT.Array, &ione, &ione, DispTdT.Desc, &ione, DispT.Array, &ione, &ione, DispT.Desc, &ione ); /* ui = ui1 */
	  pscopy_( &VelTdT.GlobalSize.Row, VelTdT.Array, &ione, &ione, VelTdT.Desc, &ione, VelT.Array, &ione, &ione, VelT.Desc, &ione ); /* vi = vi1 */
	  pscopy_( &AccTdT.GlobalSize.Row, AccTdT.Array, &ione, &ione, AccTdT.Desc, &ione, AccT.Array, &ione, &ione, AccT.Desc, &ione ); /* ai = ai1 */

	  if(rank == 0){
	       printf("Step %d\n", istep );
	  }
	  istep = istep + 1; 
     }
     MPI_Barrier( MPI_COMM_WORLD );
     if ( rank == Cblacs_pnum ( fc.Desc[1], RowProcess_Coupling[CNodes.Order-1], ColProcess_Coupling[CNodes.Order-1] ) ){
	  fclose( OutputFile );
     }

     TimeElapsedEnd = MPI_Wtime( );
     if ( rank == 0 ){
	  printf( "Time Integration process finished successfully in %f ms\n", TimeElapsedEnd - TimeElapsed );
	  /* Close the Connection */
	  Close_Connection( &Socket, InitCnt.Remote.Type, (unsigned int) CNodes.Order, (unsigned int ) InitCnt.Nstep, 4 );
     }

     if( InitCnt.Use_Absolute_Values ){
	  free( VelAll );
	  free( DispAll );
     } else {
	  free( AccAll );
     }

     /* Destroy the data structures */
     DestroyDistMatrix( &M );
     DestroyDistMatrix( &C );
     DestroyDistMatrix( &K );

     DestroyDistMatrix( &Keinv );
     DestroyDistMatrix( &Keinv_m );

     DestroyDistMatrix( &tempvec );

     DestroyDistMatrix( &DispT );
     DestroyDistMatrix( &DispTdT0 );

     DestroyDistMatrix( &DispTdT0_m );
     DestroyDistMatrix( &DispTdT );

     DestroyDistMatrix( &VelT );
     DestroyDistMatrix( &VelTdT );

     DestroyDistMatrix( &AccT );
     DestroyDistMatrix( &AccTdT );

     DestroyDistMatrix( &LoadTdT );

     DestroyDistMatrix( &LoadTdT1 );
     DestroyDistMatrix( &LoadVectorForm );

     DestroyDistMatrix( &EffT );

     DestroyDistMatrix( &fc );
     DestroyDistMatrix( &fcprevsub );
     DestroyDistMatrix( &fu );

     DestroyDistMatrix( &Disp );
     DestroyDistMatrix( &Vel );
     DestroyDistMatrix( &Acc );


     if ( rank == 0 ){
	  free( Keinv_c );
	  free( DispTdT0_c );
	  free( Recv );
     }

     free( CNodes.Array );
     free( LRowIndex_Coupling );
     free( LColIndex_Coupling );
     free( RowProcess_Coupling );
     free( ColProcess_Coupling );

     free( LRowIndex_Coupling_fcp );
     free( LColIndex_Coupling_fcp );
     free( RowProcess_Coupling_fcp );
     free( ColProcess_Coupling_fcp );

     /* BLACS: Exit Grid and release bhandle */
     Cblacs_gridexit( icontxt );
     Cfree_blacs_system_handle( bhandle );

     MPI_Finalize();
     return 0;
}
