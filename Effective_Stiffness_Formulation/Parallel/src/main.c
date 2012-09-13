#include <stdio.h>
#include <stdlib.h>

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
     int ione = 1;
     Scalars Constants;

     /* Variables required by the Substructure Algorithm */
     unsigned int i, istep;       /* A counter */
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

     PMatrixVector ErrCompT, ErrCompTdT;

     PMatrixVector Disp, Vel, Acc;

     int LRowIndex_fc, LColIndex_fc, RowProcess_fc, ColProcess_fc;  /* Variables to store the information regarding where the coupling position is in the fc vector and */

     float *Keinv_c;
     float *DispTdT0_c;
     float *Recv;

     float TimeElapsed;
     float TimeElapsedEnd;

     Coupling_Node CNodes;

     /* TCP socket connection variables */
     int Socket;

     MPI_Init( &argc, &argv );
     MPI_Comm_size( MPI_COMM_WORLD, &size );
     MPI_Comm_rank( MPI_COMM_WORLD, &rank );

     if ( rank == 0 ){
	  /* TODO: Review */
	  InitConstants( &InitCnt, size );
     }

     /* Send the information read in process 0 to all other processes using MPI_Bcast */
     BroadcastConfFile( &InitCnt );

     AccAll = (float *) calloc( (size_t) InitCnt.Nstep, sizeof(float) );
     if( InitCnt.Use_Absolute_Values ){
	  VelAll = (float *) calloc( (size_t) InitCnt.Nstep, sizeof(float) );
	  DispAll = (float *) calloc( (size_t) InitCnt.Nstep, sizeof(float) );
     } else {
	  VellAll = NULL;
	  DispAll = NULL;
     }

     /* BLACS: Initialise BLACS */
     bhandle = Csys2blacs_handle( MPI_COMM_WORLD );
     icontxt = bhandle;

     /* BLACS: Initialise the process grid */
     Cblacs_gridinit( &icontxt, (char *)"Row", InitCnt.ProcessGrid.Rows, InitCnt.ProcessGrid.Cols );
     Cblacs_gridinfo( icontxt, &nprow, &npcol, &myrow, &mycol );

     CreateDistMatrix( icontxt, &M, InitCnt.Order, InitCnt.Order,
		       InitCnt.BlockSize.Rows, InitCnt.BlockSize.Cols );
     CreateDistMatrix( icontxt, &C, InitCnt.Order, InitCnt.Order,
		       InitCnt.BlockSize.Rows, InitCnt.BlockSize.Cols );
     CreateDistMatrix( icontxt, &K, InitCnt.Order, InitCnt.Order,
		       InitCnt.BlockSize.Rows, InitCnt.BlockSize.Cols );

     CreateDistMatrix( icontxt, &Keinv, InitCnt.Order, InitCnt.Order,
		       InitCnt.BlockSize.Rows, InitCnt.BlockSize.Cols );

     CreateDistMatrix( icontxt, &Keinv_m, InitCnt.Order - InitCnt.OrderC,
		       InitCnt.OrderC, InitCnt.BlockSize.Rows,
		       InitCnt.BlockSize.Cols );

     CreateDistMatrix( icontxt, &tempvec, InitCnt.Order, 1,
		       InitCnt.BlockSize.Rows, 1 );

     CreateDistMatrix( icontxt, &DispT, InitCnt.Order, 1,
		       InitCnt.BlockSize.Rows, 1 );
     CreateDistMatrix( icontxt, &DispTdT0, InitCnt.Order, 1,
		       InitCnt.BlockSize.Rows, 1 );
     CreateDistMatrix( icontxt, &DispTdT, InitCnt.Order, 1,
		       InitCnt.BlockSize.Rows, 1 );

     CreateDistMatrix( icontxt, &DispTdT0_m, InitCnt.Order - InitCnt.OrderC,
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

     CreateDistMatrix( icontxt, &fc, InitCnt.Order, 1,
		       InitCnt.BlockSize.Rows, 1 );
     CreateDistMatrix( icontxt, &fcprevsub, InitCnt.Order, 1,
		       InitCnt.BlockSize.Rows, 1 );
     CreateDistMatrix( icontxt, &fu, InitCnt.Order, 1,
		       InitCnt.BlockSize.Rows, 1 );

     CreateDistMatrix( icontxt, &ErrCompT, InitCnt.Order, 1,
		       InitCnt.BlockSize.Rows, 1 );

     CreateDistMatrix( icontxt, &Disp, InitCnt.Order, 1,
		       InitCnt.BlockSize.Rows, 1 );
     CreateDistMatrix( icontxt, &Vel, InitCnt.Order, 1,
		       InitCnt.BlockSize.Rows, 1 );
     CreateDistMatrix( icontxt, &Acc, InitCnt.Order, 1,
		       InitCnt.BlockSize.Rows, 1 );

     if ( rank == 0 ){
	  Keinv_c = (float *) calloc( (size_t) InitCnt.OrderC*InitCnt.OrderC, sizeof(float) );
	  for ( i = 0; i < InitCnt.OrderC*InitCnt.OrderC; i++ ){
	       Keinv_c[i] = 0.0f;
	  }

	  DispTdT0_c = (float *) calloc( (size_t) InitCnt.OrderC, sizeof(float) );
	  for ( i = 0; i < InitCnt.OrderC; i++ ){
	       DispTdT0_c[i] = 0.0f;
	  }
	  Recv = (float *) calloc( (size_t) 3*InitCnt.OrderC, sizeof(float) );
	  for ( i = 0; i < 3*InitCnt.OrderC; i++ ){
	       Recv[i] = 0.0f;
	  }

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
     BuildMatrixXc( MPI_COMM_WORLD, &Keinv, Keinv_c, InitCnt.PosCouple, InitCnt.OrderC );
     BuildMatrixXcm( MPI_COMM_WORLD, &Keinv, &Keinv_m, InitCnt.PosCouple, InitCnt.OrderC );

     if ( rank == 0 ){
	  Send_Effective_Matrix( Keinv_c, InitCnt.Type_Protocol, InitCnt.OrderC, &Socket );
     }

     /*
      * Get the information about the position of the coupling force in the vector fc and its
      * coordinates in the process grid. Used when receiving the information from the server
      */
     infog2l_( &InitCnt.PosCouple, &ione, fc.Desc, &nprow, &npcol, &myrow, &mycol, &LRowIndex_fc, &LColIndex_fc, &RowProcess_fc, &ColProcess_fc );

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
	       ReadDataEarthquake_AbsValues( AccAll, VelAll, DispAll, InitCnt.Nstep, InitCnt.FileData );
	  } else {
	       ReadDataEarthquake_RelValues( AccAll, InitCnt.Nstep, InitCnt.FileData );
	  }
     }

     if ( InitCnt.Use_Absolute_Values ){
	  MPI_Bcast( DispAll, InitCnt.Nstep, MPI_FLOAT, 0, MPI_COMM_WORLD );
	  MPI_Bcast( VelAll, InitCnt.Nstep, MPI_FLOAT, 0, MPI_COMM_WORLD );
	  MPI_Bcast( AccAll, InitCnt.Nstep, MPI_FLOAT, 0, MPI_COMM_WORLD );
     } else {
	  MPI_Bcast( AccAll, InitCnt.Nstep, MPI_FLOAT, 0, MPI_COMM_WORLD );
     }


     if ( rank == Cblacs_pnum ( fc.Desc[1], RowProcess_fc, ColProcess_fc ) ){
	  OutputFile = fopen( Icnt.FileOutput, "w" );
	  if ( OutputFile == NULL ){
	       PrintErrorAndExit( "Cannot proceed because the file Out.txt could not be opened" );
	  } else {
	       clock = time (NULL);	  
	       fprintf( OutputFile, "Test started at %s", ctime( &clock ) );
	  }
     }

     istep = 1;

     /* Calculate the input load */ 
     Set2Value( &Disp, DispAll[istep - 1] );
     Set2Value( &Vel, VelAll[istep - 1] );
     Set2Value( &Acc, AccAll[istep - 1] );
     Calc_Input_Load( &LoadTdT, &K, &C, &M, &DiagM, &Disp, &Vel, &Acc );

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
	  CreateVectorXm( &DispTdT0, &DispTdT0_m, InitCnt.PosCouple, InitCnt.OrderC );
	  CreateVectorXc( MPI_COMM_WORLD, &DispTdT0, DispTdT0_c, InitCnt.PosCouple, InitCnt.OrderC );

	  if ( rank == 0 ){
	       /* Perform substepping */
	       Do_Substepping( DispTdT0_c, &Recv[0], &Recv[InitCnt.OrderC], &Recv[2*InitCnt.OrderC], InitCnt.Type_Protocol,
			       InitCnt.Delta_t*istep, InitCnt.Delta_t, 4, Socket, InitCnt.OrderC, InitCnt.PosCouple );
	  }

	  /* Send the received information from the server to the respective vectors */
	  if ( rank == 0 && Cblacs_pnum( icontxt, RowProcess_fc, ColProcess_fc ) == 0 ){
	       DispTdT.Array[(LRowIndex_fc - 1) + DispTdT.LocalSize.Row*(LColIndex_fc -1)] = Recv[0];
	       fcprevsub.Array[(LRowIndex_fc - 1) + fc.LocalSize.Row*(LColIndex_fc -1)] = Recv[InitCnt.OrderC];
	       fc.Array[(LRowIndex_fc - 1 ) + fc.LocalSize.Row*(LColIndex_fc - 1)] = Recv[2*InitCnt.OrderC];
	  } else {
	       if ( rank == 0 ){
		    MPI_Send( &Recv[0], InitCnt.OrderC, MPI_FLOAT, Cblacs_pnum( fc.Desc[1], RowProcess_fc, ColProcess_fc ), 2, MPI_COMM_WORLD );
		    MPI_Send( &Recv[InitCnt.OrderC], InitCnt.OrderC, MPI_FLOAT, Cblacs_pnum( fc.Desc[1], RowProcess_fc, ColProcess_fc ), 10, MPI_COMM_WORLD );
		    MPI_Send( &Recv[2*InitCnt.OrderC], InitCnt.OrderC, MPI_FLOAT, Cblacs_pnum( fc.Desc[1], RowProcess_fc, ColProcess_fc ), 11, MPI_COMM_WORLD );
	       } else if ( rank == Cblacs_pnum( icontxt, RowProcess_fc, ColProcess_fc ) ){
		    MPI_Recv( &DispTdT.Array[(LRowIndex_fc - 1) + DispTdT.LocalSize.Row*(LColIndex_fc -1)], InitCnt.OrderC, MPI_FLOAT, 0, 2, MPI_COMM_WORLD, &Status );
		    MPI_Recv( &fcprevsub.Array[(LRowIndex_fc - 1) + fc.LocalSize.Row*(LColIndex_fc -1)], InitCnt.OrderC, MPI_FLOAT, 0, 10, MPI_COMM_WORLD, &Status );
		    MPI_Recv( &fc.Array[(LRowIndex_fc - 1 ) + fc.LocalSize.Row*(LColIndex_fc - 1)], InitCnt.OrderC, MPI_FLOAT, 0, 11, MPI_COMM_WORLD, &Status );
	       }
	  }

	  if ( istep < InitCnt.Nstep ){
	       Set2Value( &Disp, DispAll[istep] );
	       Set2Value( &Vel, VelAll[istep] );
	       Set2Value( &Acc, AccAll[istep] );
	       Calc_Input_Load( &LoadTdT1, &K, &C, &M, &DiagM, &Disp, &Vel, &Acc );
	  }

	  /* Join the non-coupling part. DispTdT_m = Keinv_m*fc + DispTdT0_m. Although DispTdT0_m is what has been received from the other computer,
	     it has the name of DispTdT_m to avoid further operations if using the NETLIB libraries. */
	  JoinNonCouplingPart( &DispTdT0_m, &Keinv_m, &fcprevsub, &DispTdT, InitCnt.PosCouple, InitCnt.OrderC );

	  /* Compute acceleration ai1 = a0*(ui1 -ui) - a2*vi -a3*ai */
	  Compute_Acceleration( &DispTdT, &DispT, &VelT, &AccT, InitCnt.a0, InitCnt.a2, InitCnt.a3, &AccTdT );

	  /* Compute Velocity. vi = vi + a6*ai +a7*ai */
	  Compute_Velocity( &VelT, &AccT, &AccTdT, InitCnt.a6, InitCnt.a7, &VelTdT );

	  /* Error Compensation. fu = LoatTdT + fc -(Mass*AccTdT + Damp*VelTdT + Stiff*DispTdT) */
	  Compute_Force_Error( &M, &C, &K, &AccTdT, &VelTdT, &DispTdT, &fc, &LoadTdT, &fu );

	  /* Backup vectors */
	  pscopy_( &DispTdT.GlobalSize.Row, DispTdT.Array, &ione, &ione, DispTdT.Desc, &ione, DispT.Array, &ione, &ione, DispT.Desc, &ione ); /* ui = ui1 */
	  pscopy_( &VelTdT.GlobalSize.Row, VelTdT.Array, &ione, &ione, VelTdT.Desc, &ione, VelT.Array, &ione, &ione, VelT.Desc, &ione ); /* vi = vi1 */
	  pscopy_( &AccTdT.GlobalSize.Row, AccTdT.Array, &ione, &ione, AccTdT.Desc, &ione, AccT.Array, &ione, &ione, AccT.Desc, &ione ); /* ai = ai1 */

	  if ( rank == Cblacs_pnum ( fc.Desc[1], RowProcess_fc, ColProcess_fc ) ){
	       fprintf( OutputFile, "%e\t%e\t%e\n", DispTdT.Array[(LRowIndex_fc - 1 ) + DispTdT.LocalSize.Row*(LColIndex_fc - 1)], fc.Array[(LRowIndex_fc - 1 ) + fc.LocalSize.Row*(LColIndex_fc - 1)], fu.Array[ (LRowIndex_fc - 1 ) + fu.LocalSize.Row*(LColIndex_fc - 1)] );
     }

	  istep = istep + 1;

	  
     }

     TimeElapsedEnd = MPI_Wtime( );
     if ( rank == 0 ){
	  printf( "Time Integration process finished sDispTdT0_ccessfully in %f ms\n", TimeElapsedEnd - TimeElapsed );
	  /* Close the Connection */
	  Close_Connection( &Socket, InitCnt.Type_Protocol, InitCnt.Nstep, 4 );
     }
     
     if ( rank == Cblacs_pnum ( fc.Desc[1], RowProcess_fc, ColProcess_fc ) ){
	  fclose( OutputFile );
     }

     free( AccAll );
     if( InitCnt.Use_Absolute_Values ){
	  free( VelAll );
	  free( DispAll );
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
     DestroyDistMatrix( &EffT );
     DestroyDistMatrix( &DiagM );

     DestroyDistMatrix( &fc );
     DestroyDistMatrix( &fcprevsub );
     DestroyDistMatrix( &fu );

     DestroyDistMatrix( &ErrCompT );
     DestroyDistMatrix( &ErrCompTdT );

     DestroyDistMatrix( &Disp );
     DestroyDistMatrix( &Vel );
     DestroyDistMatrix( &Acc );

     if ( rank == 0 ){
	  free( Keinv_c );
	  free( DispTdT0_c );
	  free( Recv );
     }

     /* BLACS: Exit Grid and release bhandle */
     Cblacs_gridexit( icontxt );
     Cfree_blacs_system_handle( bhandle );

     MPI_Finalize();
     return 0;
}
