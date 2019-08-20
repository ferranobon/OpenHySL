#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <mpi.h>

#include "Substructure.h"
#include "Substructure_Exact.h"
#include "Substructure_Remote.h"
#include "Substructure_Remote_OpenFresco.h"
#include "Substructure_Remote_NSEP.h"
#include "Substructure_UHYDEfbr.h"
#include "Substructure_SimMeasured.h"
#include "Substructure_Experimental.h"
#include "Substructure_CouplingNodes.h"
#include "Substructure_Auxiliary.h"

#include "Definitions.h"

#include "Print_Messages.h"

#if _ADWIN_
#include "ADwin_Routines.h"
#endif

#if _MKL_
#include "Cblacs.h"
#else
#include "Netlib.h"
#include "Cblacs.h"
#endif

void Substructure_Substepping_MPI (const hysl_float_t *const IGain, const hysl_float_t *const VecTdT0_c, const hysl_float_t Time, const hysl_float_t GAcc, const unsigned int NSubstep,
        const hysl_float_t DeltaT_Sub, const MPI_Comm Comm, const InfoLocation_t *const ILoc_VecTdT, const InfoLocation_t *const ILoc_CoupForcePrev, const InfoLocation_t *const ILoc_CoupForce,
        CouplingNode_t *const CNodes, PMatrixVector_t *const VecTdT, PMatrixVector_t *const CoupForcePrev, PMatrixVector_t *const CoupForce) {

    int i, rank;
    bool Called_Sub = false;
    hysl_float_t *Recv = NULL;
    hysl_float_t *Send = NULL;

    Remote_t *Remote;
    MPI_Status Status;

    MPI_Comm_rank(Comm, &rank);

    if (rank == 0) {
        Recv = (hysl_float_t*) calloc((size_t) 3 * (size_t) CNodes->Order, sizeof(hysl_float_t));
    } else {
        Recv == NULL;
    }

    /* Copy the older coupling force. This is necessary for simulations */
    for (i = 0; i < CNodes->Order; i++) {
        if (rank == 0 && Cblacs_pnum(CoupForce->Desc[1], ILoc_CoupForce->RowProcess[i], ILoc_CoupForce->ColProcess[i]) == 0) {
            Recv[2 * CNodes->Order + i] = CoupForce->Array[ILoc_CoupForce->LRowIndex[i] - 1];
        } else {
            if (rank == Cblacs_pnum(CoupForce->Desc[1], ILoc_CoupForce->RowProcess[i], ILoc_CoupForce->ColProcess[i])) {
                MPI_Send(&CoupForce->Array[ILoc_CoupForce->LRowIndex[i] - 1], 1, MPI_HYSL_FLOAT, 0, MPI_FORCE, Comm);
            } else if (rank == 0) {
                MPI_Recv(&Recv[2 * CNodes->Order + i], 1, MPI_HYSL_FLOAT, Cblacs_pnum(CoupForce->Desc[1], ILoc_CoupForce->RowProcess[i], ILoc_CoupForce->ColProcess[i]), MPI_FORCE, Comm, &Status);
            }
        }
    }

    if (rank == 0) {
        for (i = 0; i < CNodes->Order; i++) {
            switch (CNodes->Sub[i].Type) {
            case SIM_EXACT_MDOF:
                /* This is the same case as SIM_MEASURED. All the simulated substructures are treated together
                 * in the same routine.*/
            case SIM_EXACT_SDOF:
                /* This is the same case as SIM_MEASURED. All the simulated substructures are treated together
                 * in the same routine.*/
            case SIM_EXACT_ESP:
            case SIM_UHYDE:
                /* This is the same case as SIM_MEASURED. All the simulated substructures are treated together
                 * in the same routine.*/
            case SIM_MEASURED:
                /* Call the Simulate_Substructures() function only once. All the simulated substructures are
                 * handled together in this routine */
                if (!Called_Sub) {
                    Substructure_Simulate(IGain, VecTdT0_c, GAcc, NSubstep, DeltaT_Sub, CNodes, &Recv[0], &Recv[CNodes->Order], &Recv[2 * CNodes->Order]);
                    Called_Sub = true;
                }
                break;
#if _ADWIN_
	       case EXP_ADWIN:
		    /* Tell ADwin to perform the substepping process */
		    ADwin_Substep( VecTdT0_c, (unsigned int) CNodes->Order, 0.75, &Recv[0], &Recv[1], &Recv[2] );
		    break;
#endif
            case REMOTE:
                Remote = (Remote_t*) CNodes->Sub[i].SimStruct;

                if( Remote->Type == REMOTE_TCP || Remote->Type == REMOTE_UDP || Remote->Type == REMOTE_CELESTINA ) {
                    Send = (hysl_float_t *) calloc( (size_t) 1+(size_t)CNodes->Order, sizeof(hysl_float_t) );
                    for( i = 0; i < CNodes->Order; i++ ) {
                        Send[i] = VecTdT0_c[i];
                    }
                    Send[CNodes->Order] = GAcc;

                    Substructure_Remote_Send( Remote->Socket, (unsigned int) CNodes->Order + 1, sizeof(hysl_float_t), (char *const) Send );

                    Substructure_Remote_Receive( Remote->Socket, 3*(unsigned int) CNodes->Order, sizeof(hysl_float_t), (char *const) Recv );
                    free( Send );
                } else if( Remote->Type == REMOTE_NSEP ) {
                    /* Using NSEP Protocol */
                    Substructure_Remote_NSEP( Remote, NSEP_SEND_CMD, Time, Remote->NSub, Send, Recv );
                    /* Receive the force from the PNSE server. WhatToDo = 2 */
                    Substructure_Remote_NSEP( Remote, NSEP_REQUEST_CSIG, Time, 3*Remote->NSub, Send, Recv );
                } else if( Remote->Type == REMOTE_OF ) {
                    /* Using OpenFresco */
                    /* Send the trial response */
                    Substructure_Remote_OpenFresco( Remote->Socket, OF_REMOTE_SET_TRIAL_RESPONSE, Remote->NSub, Send, NULL );
                    /* Get the DAQ response */
                    Substructure_Remote_OpenFresco( Remote->Socket, OF_REMOTE_GET_DAQ_RESPONSE, Remote->NSub, NULL, Recv );
                } else assert( Remote->Type >= 0 || Remote->Type < NUM_REMOTE_TYPE );
                break;
            }
        }
    }

    for (i = 0; i < CNodes->Order; i++) {
        if (rank == 0 && Cblacs_pnum(VecTdT->Desc[1], ILoc_VecTdT->RowProcess[i], ILoc_VecTdT->ColProcess[i]) == 0) {
            VecTdT->Array[ILoc_VecTdT->LRowIndex[i] - 1] = Recv[i];
        } else {
            if (rank == 0) {
                MPI_Send(&Recv[i], 1, MPI_HYSL_FLOAT, Cblacs_pnum(VecTdT->Desc[1], ILoc_VecTdT->RowProcess[i], ILoc_VecTdT->ColProcess[i]), MPI_NEW_STATE, Comm);
            } else if (rank == Cblacs_pnum(VecTdT->Desc[1], ILoc_VecTdT->RowProcess[i], ILoc_VecTdT->ColProcess[i])) {
                MPI_Recv(&VecTdT->Array[ILoc_VecTdT->LRowIndex[i] - 1], 1, MPI_HYSL_FLOAT, 0, MPI_NEW_STATE, Comm, &Status);
            }
        }

        if (rank == 0 && Cblacs_pnum(CoupForcePrev->Desc[1], ILoc_CoupForcePrev->RowProcess[i], ILoc_CoupForcePrev->ColProcess[i]) == 0) {
            CoupForcePrev->Array[ILoc_CoupForcePrev->LRowIndex[i] - 1] = Recv[CNodes->Order + i];
        } else {
            if (rank == 0) {
                MPI_Send(&Recv[1 + i], 1, MPI_HYSL_FLOAT, Cblacs_pnum(CoupForcePrev->Desc[1], i, i), MPI_PREV_FORCE, Comm);
            } else if (rank == Cblacs_pnum(CoupForcePrev->Desc[1], ILoc_CoupForcePrev->RowProcess[i], ILoc_CoupForcePrev->ColProcess[i])) {
                MPI_Recv(&CoupForcePrev->Array[ILoc_CoupForce->LRowIndex[i] - 1], 1, MPI_HYSL_FLOAT, 0, MPI_PREV_FORCE, Comm, &Status);
            }
        }

        if (rank == 0 && Cblacs_pnum(CoupForce->Desc[1], ILoc_CoupForce->RowProcess[i], ILoc_CoupForce->ColProcess[i]) == 0) {
            CoupForce->Array[ILoc_CoupForce->LRowIndex[i] - 1] = Recv[2 * CNodes->Order + i];
        } else {
            if (rank == 0) {
                MPI_Send(&Recv[2 * 1 + i], 1, MPI_HYSL_FLOAT, Cblacs_pnum(CoupForce->Desc[1], ILoc_CoupForce->RowProcess[i], ILoc_CoupForce->ColProcess[i]), MPI_FORCE, Comm);
            } else if (rank == Cblacs_pnum(CoupForce->Desc[1], ILoc_CoupForce->RowProcess[i], ILoc_CoupForce->ColProcess[i])) {
                MPI_Recv(&CoupForce->Array[ILoc_CoupForce->LRowIndex[i] - 1], 1, MPI_HYSL_FLOAT, 0, MPI_FORCE, Comm, &Status);
            }
        }
    }

    free(Recv);
}
