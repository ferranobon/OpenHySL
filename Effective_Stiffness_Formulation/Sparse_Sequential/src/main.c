#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "MatrixVector.h"
#include "Initiation.h"
#include "Precalculations.h"
#include "ErrorHandling.h"
#include "Netlib.h"
#include "ComputeU0.h"
#include "EndingStep.h"
#include "Send_Receive_Data.h"

#if REAL_TIME_
#include <sched.h>    /* For sched_setscheduler( ) */
#include <malloc.h>
#include <sys/mman.h> /* For mlockall( ) */
#include <unistd.h>   /* For sysconf( ) */
#define SOMESIZE (150*1024*1024) // 150MB
#endif

int main( int argc, char **argv )
{
     
     /* Output file */
     FILE *OutputFile;

     int i, istep;		/* Counters */
     AlgConst InitCnt;

     /* NETLIB Variables */
     int incx, incy;
     Scalars Constants;
     
     float *AccAll, *VelAll, *DispAll;
     /* Variables to store the result we desire, so that no disk i/o is done during the test */
     float *TimeHistoryli, *TimeHistoryai1, *TimeHistoryvi1, *TimeHistoryui1, *TimeHistoryai, *TimeHistoryvi, *TimeHistoryui, *TimeHistoryfc, *TimeHistoryfu;  

     MatrixVector M, C, K;               /* Mass, Damping and Stiffness matrices */
     MatrixVector Keinv;
     MatrixVector Keinv_c, Keinv_m;

     MatrixVector DiagM;
     MatrixVector EffT;

     MatrixVector DispT, DispTdT0, DispTdT;
     MatrixVector DispTdT0_c, DispTdT0_m;
     MatrixVector tempvec;

     MatrixVector VelT, VelTdT;
     MatrixVector AccT, AccTdT;

     MatrixVector LoadTdT, LoadTdT1;

     MatrixVector fc, fcprevsub;
     MatrixVector fu;

     MatrixVector ErrCompT, ErrCompTdT;

     MatrixVector Disp, Vel, Acc;

     /* TCP socket connection Variables */
     int Socket;

#if REAL_TIME_
     struct sched_param sp;
     int policy;

     if((policy = sched_getscheduler(0) == -1)) {
	  fprintf(stderr, "Error getting scheduler" );
     }
     if( policy != SCHED_FIFO ) {
	  printf("Setting priority: SCHED_FIFO.\n");
	  sp.sched_priority = sched_get_priority_max(SCHED_FIFO);
	  if(sched_setscheduler(0, SCHED_FIFO, &sp) == -1 ){
	    PrintErrorAndExit("Filed real time\n");
	  }
     }

     // Allocate some memory
       int page_size;
       char* buffer;       

       // Now lock all current and future pages from preventing of being paged
       if (mlockall(MCL_CURRENT | MCL_FUTURE ))
       {
           perror("mlockall failed:");
       }
       

       // Turn off malloc trimming.
       mallopt (M_TRIM_THRESHOLD, -1);
       

       // Turn off mmap usage.
       mallopt (M_MMAP_MAX, 0);
       

       page_size = sysconf(_SC_PAGESIZE);
       buffer = malloc(SOMESIZE);
       

       // Touch each page in this piece of memory to get it mapped into RAM
       for (i=0; i < SOMESIZE; i+=page_size){
           // Each write to this buffer will generate a pagefault.
           // Once the pagefault is handled a page will be locked in memory and never
           // given back to the system.
           buffer[i] = 0;
       }
       free(buffer);
       // buffer is now released. As glibc is configured such that it never gives back memory to
       // the kernel, the memory allocated above is locked for this process. All malloc() and new()
       // calls come from the memory pool reserved and locked above. Issuing free() and delete()
       // does NOT make this locking undone. So, with this locking mechanism we can build C++ applications
       // that will never run into a major/minor pagefault, even with swapping enabled.
       

       //<do your RT-thing>
#endif


     /* Constants definitions. */
     InitConstants( &InitCnt );

     /* Allocate memory for saving the acceleration, displacement and velocity (input files) that will
      * be used during the test */
     AccAll = calloc( InitCnt.Nstep, sizeof(float) );
     VelAll = calloc( InitCnt.Nstep, sizeof(float) );
     DispAll = calloc( InitCnt.Nstep, sizeof(float) );

     /* Allocate the memory for the variables to store. The results will be saved each step */
     TimeHistoryli = calloc( InitCnt.Nstep, sizeof(float) );
     TimeHistoryui1 = calloc( InitCnt.Nstep, sizeof(float) );
     TimeHistoryvi1 = calloc( InitCnt.Nstep, sizeof(float) );
     TimeHistoryai1 = calloc( InitCnt.Nstep, sizeof(float) );
     TimeHistoryui = calloc( InitCnt.Nstep, sizeof(float) );
     TimeHistoryvi = calloc( InitCnt.Nstep, sizeof(float) );
     TimeHistoryai = calloc( InitCnt.Nstep, sizeof(float) );
     TimeHistoryfc = calloc( InitCnt.Nstep, sizeof(float) );
     TimeHistoryfu = calloc( InitCnt.Nstep, sizeof(float) );

     /* Initialise the matrices and vectors that will be used in the Time Integration process */
     Init_MatrixVector( &M, InitCnt.Order, InitCnt.Order );
     Init_MatrixVector( &K, InitCnt.Order, InitCnt.Order );
     Init_MatrixVector( &C, InitCnt.Order, InitCnt.Order );
     Init_MatrixVector( &Keinv, InitCnt.Order, InitCnt.Order );
 
     Init_MatrixVector( &Keinv_c, InitCnt.OrderC, InitCnt.OrderC );
     Init_MatrixVector( &Keinv_m, InitCnt.Order - InitCnt.OrderC, InitCnt.OrderC );

     Init_MatrixVector( &DiagM, InitCnt.Order, 1 );
     Init_MatrixVector( &tempvec, InitCnt.Order, 1 );

     Init_MatrixVector( &DispT, InitCnt.Order, 1 );
     Init_MatrixVector( &DispTdT0, InitCnt.Order, 1 );
     Init_MatrixVector( &DispTdT, InitCnt.Order, 1 );

     Init_MatrixVector( &DispTdT0_c, InitCnt.OrderC, 1 );
     Init_MatrixVector( &DispTdT0_m, InitCnt.Order-InitCnt.OrderC, 1 );

     Init_MatrixVector( &VelT, InitCnt.Order, 1 );
     Init_MatrixVector( &VelTdT, InitCnt.Order, 1 );

     Init_MatrixVector( &AccT, InitCnt.Order, 1 );
     Init_MatrixVector( &AccTdT, InitCnt.Order, 1 );

     Init_MatrixVector( &EffT, InitCnt.Order, 1 );

     Init_MatrixVector( &LoadTdT, InitCnt.Order, 1 );
     Init_MatrixVector( &LoadTdT1, InitCnt.Order, 1 );

     Init_MatrixVector( &fc, InitCnt.Order, 1 );
     Init_MatrixVector( &fcprevsub, InitCnt.Order, 1 );
     Init_MatrixVector( &fu, InitCnt.Order, 1 );

     Init_MatrixVector( &ErrCompT, InitCnt.Order, 1 );
     Init_MatrixVector( &ErrCompTdT, InitCnt.Order, 1 );

     Init_MatrixVector( &Disp, InitCnt.Order, 1 );
     Init_MatrixVector( &Vel, InitCnt.Order, 1 );
     Init_MatrixVector( &Acc, InitCnt.Order, 1 );

     /* Read the matrices from a file */
     MatrixVector_From_File( &M, InitCnt.FileM );
     MatrixVector_From_File( &K, InitCnt.FileK );

     //CalculateMatrixC( &M, &K, &C, &InitCnt.Rayleigh );
     C.Array[0] = 0.104*2*sqrtf(K.Array[0]*M.Array[0] ); /* EFAST */

     /* Calculate Matrix Keinv = [K + a0*M + a1*C]^(-1) */
     Constants.Alpha = 1.0;
     Constants.Beta = InitCnt.a0;
     Constants.Gamma = InitCnt.a1;
     CalculateMatrixKeinv( &Keinv, &M, &C, &K, Constants );

     BuildMatrixXc( &Keinv, Keinv_c.Array, InitCnt.PosCouple, InitCnt.OrderC );
     BuildMatrixXcm( &Keinv, &Keinv_m, InitCnt.PosCouple, InitCnt.OrderC );

     /* Send the coupling part of the effective matrix */
     Send_Effective_Matrix( Keinv_c.Array, InitCnt.Type_Protocol, InitCnt.OrderC, &Socket );

     /* Read the earthquake data from a file */
     ReadDataEarthquake( AccAll, VelAll, DispAll, InitCnt.Nstep, InitCnt.FileData );

     /* Open Output file. If the file cannot be opened, the program will exit, since the results cannot be stored. */
     OutputFile = fopen( "Out.txt", "w" );
     if ( OutputFile == NULL ){
	  PrintErrorAndExit( "Cannot proceed because the file Out.txt could not be opened" );
     } else {
	  fprintf( OutputFile, "Test with %d DOF\n", InitCnt.Order );
	  fprintf( OutputFile, "u0\t ai1(m/s^2)\t ai(m/s^2)\t vi1 (m/s)\t vi (m/s)\t ui1 (m)\t ui (m)\t fc (N)\t fu(N)\n" );
     }

     istep = 1;

     /* Copy the diagonal elements of M */
     CopyDiagonalValues( &M, &DiagM );

     /* Calculate the input load */
     Set2Value( &Disp, DispAll[istep - 1] );
     Set2Value( &Vel, VelAll[istep - 1] );
     Set2Value( &Acc, AccAll[istep - 1] );	  
     Calc_Input_Load( &LoadTdT, &K, &C, &M, &DiagM, &Disp, &Vel, &Acc );

     incx = 1; incy = 1;
     printf( "Starting stepping process\n" );
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
	  CreateVectorXc( &DispTdT0, DispTdT0_c.Array, InitCnt.PosCouple, InitCnt.OrderC );

	  /* Perform substepping */
	  Do_Substepping( DispTdT0_c.Array, DispTdT.Array, fcprevsub.Array, fc.Array, InitCnt.Type_Protocol,
			  InitCnt.Delta_t*istep, InitCnt.Delta_t, 4, Socket, InitCnt.OrderC, InitCnt.PosCouple );
	  
	  if ( istep < InitCnt.Nstep ){
	       /* Calculate the input load for the next step during the
		  sub-stepping process. */
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

	  /* Output variables */
	  TimeHistoryli[istep - 1] = LoadTdT.Array[InitCnt.PosCouple - 1];
	  TimeHistoryai1[istep - 1] = AccTdT.Array[InitCnt.PosCouple - 1];
	  TimeHistoryai[istep - 1] = AccT.Array[InitCnt.PosCouple - 1];
	  TimeHistoryvi1[istep - 1] = VelTdT.Array[InitCnt.PosCouple - 1];
	  TimeHistoryvi[istep - 1] = VelT.Array[InitCnt.PosCouple - 1];
	  TimeHistoryui1[istep - 1] = DispTdT.Array[InitCnt.PosCouple - 1];
	  TimeHistoryui[istep - 1] = DispT.Array[InitCnt.PosCouple - 1];
	  TimeHistoryfc[istep - 1] = fc.Array[InitCnt.PosCouple - 1];
	  TimeHistoryfu[istep - 1] = fu.Array[InitCnt.PosCouple - 1];

	  /* Backup vectors */
	  scopy_( &LoadTdT1.Rows, LoadTdT1.Array, &incx, LoadTdT.Array, &incy ); /* li = li1 */
	  scopy_( &DispTdT.Rows, DispTdT.Array, &incx, DispT.Array, &incy ); /* ui = ui1 */
	  scopy_( &VelTdT.Rows, VelTdT.Array, &incx, VelT.Array, &incy ); /* vi = vi1 */
	  scopy_( &AccTdT.Rows, AccTdT.Array, &incx, AccT.Array, &incy ); /* ai = ai1 */
	  istep = istep + 1;
     }

     /* Save the results into a file */
     for ( i = 0; i < InitCnt.Nstep; i++ ){
	  fprintf( OutputFile, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", TimeHistoryli[i], TimeHistoryai1[i], TimeHistoryai[i], TimeHistoryvi1[i], TimeHistoryvi[i], TimeHistoryui1[i], TimeHistoryui[i], TimeHistoryfc[i], TimeHistoryfu[i] );
     }
     /* Close the output file */
     fclose( OutputFile );

     /* Close the Connection */
     Close_Connection( &Socket, InitCnt.Type_Protocol, InitCnt.Nstep, 4 );
   
     /* Free the memory stored in TimeHistory variables */
     free( TimeHistoryli );
     free( TimeHistoryvi1 );
     free( TimeHistoryai1 );
     free( TimeHistoryui );
     free( TimeHistoryvi );
     free( TimeHistoryai );

     free( TimeHistoryui1 );
     free( TimeHistoryfc );
     free( TimeHistoryfu );

     /* Free the memory */
     free( AccAll );
     free( VelAll );
     free( DispAll );

     /* Destroy the data structures */
     Destroy_MatrixVector( &M );
     Destroy_MatrixVector( &C );
     Destroy_MatrixVector( &K );

     Destroy_MatrixVector( &Keinv );
     Destroy_MatrixVector( &Keinv_c );
     Destroy_MatrixVector( &Keinv_m );

     Destroy_MatrixVector( &DiagM );
     Destroy_MatrixVector( &tempvec );

     Destroy_MatrixVector( &DispT );
     Destroy_MatrixVector( &DispTdT0 );
     Destroy_MatrixVector( &DispTdT0_c );
     Destroy_MatrixVector( &DispTdT0_m );
     Destroy_MatrixVector( &DispTdT );

     Destroy_MatrixVector( &VelT );
     Destroy_MatrixVector( &VelTdT );

     Destroy_MatrixVector( &AccT );
     Destroy_MatrixVector( &AccTdT );

     Destroy_MatrixVector( &LoadTdT );
     Destroy_MatrixVector( &LoadTdT1 );

     Destroy_MatrixVector( &EffT );

     Destroy_MatrixVector( &fc );
     Destroy_MatrixVector( &fcprevsub );
     Destroy_MatrixVector( &fu );

     Destroy_MatrixVector( &ErrCompT );
     Destroy_MatrixVector( &ErrCompTdT );

     Destroy_MatrixVector( &Disp );
     Destroy_MatrixVector( &Vel );
     Destroy_MatrixVector( &Acc );

     return 0;
}