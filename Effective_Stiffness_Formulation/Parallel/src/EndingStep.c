/**
 * \file EndingStep.c
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 10th of August 2011
 *
 * \brief Source code of the functions used during the Compute U0 phase.
 *
 * This file contains the source code of the functions that are called during the Ending Step phase of the substructure
 * algorithm. For the moment, only the joining of the non-coupling and the coupling part routine has been considered.
 *
 * \todo Evaluate the performance impact (overhead) of writing routines to calculate the acceleration, velocity and displacement in this
 * phase as a separate routines.
 */

#include "EndingStep.h"
#include "PMatrixVector.h"
#include "Netlib.h"

void JoinNonCouplingPart( PMatrixVector *const VecXm, PMatrixVector *const Keinv_m, PMatrixVector *const fcprevsub,
			  PMatrixVector *const Vec, int PosCouple, const int OrderC )
{

     static int ione = 1;
     static float Alpha = 1.0, Beta = 0.0;
     static char trans = 'N';
     static int Rows, Cols, TempSize, Position;

     Rows = Keinv_m->GlobalSize.Row;
     Cols = Keinv_m->GlobalSize.Col;

     Position = PosCouple;
     psgemv_( &trans, &Rows, &Cols, &Alpha, Keinv_m->Array, &ione, &ione, Keinv_m->Desc, fcprevsub->Array, &PosCouple, &ione, fcprevsub->Desc, &ione, &Beta, VecXm->Array, &ione, &ione, VecXm->Desc, &ione );

     /* Copy the first elements */
     TempSize = PosCouple - 1;
     pscopy_( &TempSize, VecXm->Array, &ione, &ione, VecXm->Desc, &ione, Vec->Array, &ione, &ione, Vec->Desc, &ione );

     /* Join the part after the coupling position */
     TempSize = Vec->GlobalSize.Row - ( PosCouple + OrderC - 1 );
     Position = PosCouple + OrderC;
     pscopy_( &TempSize, VecXm->Array, &PosCouple, &ione, VecXm->Desc, &ione, Vec->Array, &Position, &ione, Vec->Desc, &ione );

}

void Compute_Acceleration( PMatrixVector *const DispTdT, PMatrixVector *const DispT, PMatrixVector *const VelT,
			   PMatrixVector *const AccT, const float a0, const float a2, const float a3,
			   PMatrixVector *const AccTdT )
{

     static int ione = 1;
     static int incx = 1, incy = 1;
     static float Alpha;

     /* PBLAS: AccTdT = DispTdT */
     pscopy_( &AccTdT->GlobalSize.Row, DispTdT->Array, &ione, &ione, DispTdT->Desc, &ione, AccTdT->Array, &ione, &ione, AccTdT->Desc, &ione );
     /* PBLAS: AccTdT = DispTdT - DispT = AccTdT - DispT */
     Alpha = -1.0;
     psaxpy_( &AccTdT->GlobalSize.Row, &Alpha, DispT->Array, &ione, &ione, DispT->Desc, &incx, AccTdT->Array, &ione, &ione, AccTdT->Desc, &incy );
     /* PBLAS: AccTdT = a0*(DispTdT - DispT) = a0*AccTdT */
     Alpha = a0;
     psscal_( &AccTdT->GlobalSize.Row, &Alpha, AccTdT->Array, &ione, &ione, AccTdT->Desc, &incx );
     /* PBLAS: AccTdT = a0*(DispTdT - DispT) - a2*VelT = AccTdT - a2*VelT */
     Alpha = -a2;
     psaxpy_( &AccTdT->GlobalSize.Row, &Alpha, VelT->Array, &ione, &ione, VelT->Desc, &incx, AccTdT->Array, &ione, &ione, AccTdT->Desc, &incy );
     /* OBLAS: AccTdT = a0*(DispTdT - DispT) - a2*VelT - a3*AccT = AccTdT - a3*AccT */
     Alpha = -a3;
     psaxpy_( &AccTdT->GlobalSize.Row, &Alpha, AccT->Array, &ione, &ione, AccT->Desc, &incx, AccTdT->Array, &ione, &ione, AccTdT->Desc, &incy );
}

void Compute_Velocity( PMatrixVector *const VelT, PMatrixVector *const AccT, PMatrixVector *const AccTdT,
		       const float a6, const float a7, PMatrixVector *const VelTdT )
{
     static int ione = 1;
     static int incx = 1, incy= 1;
     static float Alpha;

     /* BLAS: VelTdT = VelT */
     pscopy_( &VelTdT->GlobalSize.Row, VelT->Array, &ione, &ione, VelT->Desc, &ione, VelTdT->Array, &ione, &ione, VelTdT->Desc, &ione );
     /* BLAS: VelTdT = VelT + a6*AccT = VelTdT + a6*AccT */
     Alpha = a6;
     psaxpy_( &VelTdT->GlobalSize.Row, &Alpha, AccT->Array, &ione, &ione, AccT->Desc, &incx, VelTdT->Array, &ione, &ione, VelTdT->Desc, &incy );
     /* BLAS: VelTdT = VelT + a6*AccT + a7*AccTdT = VelTdT + a7*AccTdT */
     Alpha = a7;
     psaxpy_( &VelTdT->GlobalSize.Row, &Alpha, AccTdT->Array, &ione, &ione, AccTdT->Desc, &incx, VelTdT->Array, &ione, &ione, VelTdT->Desc, &incy );
}

void Compute_Force_Error( PMatrixVector *const Mass, PMatrixVector *const Damp, PMatrixVector *Stiff,
			  PMatrixVector *const AccTdT, PMatrixVector *const VelTdT, PMatrixVector *const DispTdT,
			  PMatrixVector *const fc, PMatrixVector *const LoadTdT, PMatrixVector *const fu )
{

     static int ione = 1;
     static int incx = 1, incy = 1;
     static float Alpha, Beta;
     static char uplo = 'L';

     /* PBLAS: fu = Mass*AccTdT */
     Alpha = 1.0; Beta = 0.0;
     pssymv_( &uplo, &fu->GlobalSize.Row, &Alpha, Mass->Array, &ione, &ione, Mass->Desc, AccTdT->Array, &ione, &ione, AccTdT->Desc, &incx, &Beta, fu->Array, &ione, &ione, fu->Desc, &incy );
     /* PBLAS: fu = Mass*AccTdT + Damp*VelTdT = fu + Damp*VelTdT */
     Beta = 1.0;
     pssymv_( &uplo, &fu->GlobalSize.Row, &Alpha, Damp->Array, &ione, &ione, Damp->Desc, VelTdT->Array, &ione, &ione, VelTdT->Desc, &incx, &Beta, fu->Array, &ione, &ione, fu->Desc, &incy );
     /* BLAS: fu = Mass*AccTdT + Damp*VelTdT + Stiff*DispTdT = fu + Stiff*DispTdT */
     pssymv_( &uplo, &fu->GlobalSize.Row, &Alpha, Stiff->Array, &ione, &ione, Stiff->Desc, DispTdT->Array, &ione, &ione, DispTdT->Desc, &incx, &Beta, fu->Array, &ione, &ione, fu->Desc, &incy );
     /* BLAS: fu = -(Mass*AccTdT + Damp*VelTdT + Stiff*DispTdT) = -fu */
     Alpha = -1.0;
     psscal_( &fu->GlobalSize.Row, &Alpha, fu->Array, &ione, &ione, fu->Desc, &incx );
     /* BLAS: fu = fc -(Mass*AccTdT + Damp*VelTdT + Stiff*DispTdT) = fc -fu */
     Alpha = 1.0;
     psaxpy_( &fu->GlobalSize.Row, &Alpha, fc->Array, &ione, &ione, fc->Desc, &incx, fu->Array, &ione, &ione, fu->Desc, &incy );
     /* BLAS: fu = LoatTdT + fc -(Mass*AccTdT + Damp*VelTdT + Stiff*DispTdT) = LoadTdT -fu */
     psaxpy_( &fu->GlobalSize.Row, &Alpha, LoadTdT->Array, &ione, &ione, LoadTdT->Desc, &incx, fu->Array, &ione, &ione, fu->Desc, &incy );

}
