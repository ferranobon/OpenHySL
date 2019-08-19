#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <assert.h>

#include "Print_Messages.h"

#include "Substructure.h"
#include "Substructure_BoucWen.h"
#include "Substructure_Exact.h"
#include "Substructure_Newmark.h"
#include "Substructure_Remote.h"
#include "Substructure_Remote_OpenFresco.h"
#include "Substructure_Remote_NSEP.h"
#include "Substructure_UHYDEfbr.h"
#include "Substructure_SimMeasured.h"
#include "Substructure_Experimental.h"
#include "Substructure_CouplingNodes.h"
#include "Substructure_StoneDrums.h"

#include "Definitions.h"

#if _ADWIN_
#include "ADwin_Routines.h"
#endif

#if _MKL_
#include <mkl_blas.h>
#else
#include "Netlib.h"
#endif

void Substructure_SendGainMatrix (const hysl_float_t *const Gain, uint32_t Order, const Substructure_t *const Substructure) {
    Remote_t *Remote = NULL;
    hysl_float_t *Send, *Recv; /* Only used in case of of substructures of remote type REMOTE_NSEP */

    if (Substructure->Type == EXP_ADWIN) {
#if _ADWIN_

        /* Send matrix Gc to ADwin. */
        /* In the code of ADWIN, the variable G is stored in DATA_70 */
        ADwin_SendArray( 70, Gain, Order*Order );

        Print_Header( SUCCESS );
        printf("Gain Matrix successfully sent to ADwin system.\n" );
#else
        Print_Header( ERROR);
        fprintf( stderr, "The support for ADwin was disabled at runtime.\n");
        exit( EXIT_FAILURE);
#endif
    } else if (Substructure->Type == REMOTE) {
        Remote = (Remote_t*) Substructure->SimStruct;

        if ((Remote->Type == REMOTE_TCP) || (Remote->Type == REMOTE_UDP) || (Remote->Type == REMOTE_CELESTINA)) {
            Substructure_Remote_Send(Remote->Socket, Order * Order, sizeof(hysl_float_t), (const char* const ) Gain);
            Print_Header( SUCCESS);
            printf("Gain Matrix successfully sent to Remote site %s:%s (%s protocol).\n", Remote->IP, Remote->Port, Substructure_Remote_Type[Remote->Type]);
        } else if (Remote->Type == REMOTE_NSEP) {
            /*
             * Send the matrix Gc to the PNSE server in order to reach the FCM.
             *
             * Note that the CGM can only send a NSEP_CMD message with the size of order, therefore it is needed
             * to send the matrix Gc per rows to fullfill this requisite. This also means, that the FCM must
             * send as many NSEP_CSIG packets as the number of rows to keep everything synchronised
             */
            Send = (hysl_float_t*) calloc((size_t) Order, sizeof(hysl_float_t));
            Recv = (hysl_float_t*) calloc((size_t) Order, sizeof(hysl_float_t));
            for (uint32_t idx = 0; idx < Order; idx++) {
                for (uint32_t jdx = 0; jdx < Order; jdx++) {
                    Send[jdx] = Gain[(idx * Order) + jdx];
                }
                /* Send the matrix Keinv_c to PNSE Server */
                Substructure_Remote_NSEP(Remote, NSEP_SEND_CMD, 0.0, Remote->NSub, Send, Recv);
                /* This is done so that PNSE do not overtake the first step */
                Substructure_Remote_NSEP(Remote, NSEP_REQUEST_CSIG, 0.0, Remote->NSub, Send, Recv);
            }

            free(Send);
            free(Recv);
        } else if (Remote->Type == REMOTE_OF) {
            /* TODO Implement Send the Matrix G in OpenFresco. Wait for the answer from Andreas. What
             * follows is an ugly hack. */
            Substructure_Remote_OpenFresco(Remote->Socket, OF_REMOTE_SET_TRIAL_RESPONSE, Remote->NSub * Remote->NSub, Gain, NULL);
        } else {
            assert((Remote->Type >= 0) || (Remote->Type < NUM_REMOTE_TYPE));
        }
    } else {
        assert((Substructure->Type == EXP_ADWIN) || (Substructure->Type == REMOTE));
    }
}

void Substructure_Substepping (const hysl_float_t *const IGain, const hysl_float_t *const VecTdT0_c, const hysl_float_t Time, const hysl_float_t GAcc, const uint32_t NSubstep,
        const hysl_float_t DeltaT_Sub, CouplingNode_t *const CNodes, hysl_float_t *const VecTdT, hysl_float_t *const CoupForcePrev,
        hysl_float_t *const CoupForce) {

    bool Called_Sub = false, Called_ADwin = false;
    hysl_float_t *Recv = NULL, *Recv_ADwin = NULL, *VecTdT0_c_ADwin = NULL;
    hysl_float_t *Send = NULL;
    bool MultipleTypes = true;

    Remote_t *Remote;

    Recv = (hysl_float_t*) calloc((size_t) 3 * (size_t) CNodes->Order, sizeof(hysl_float_t));
    if ((CNodes->OrderADwin >= 1) && MultipleTypes) {
        Recv_ADwin = (hysl_float_t*) calloc((size_t) 3 * (size_t) CNodes->OrderADwin, sizeof(hysl_float_t));
        VecTdT0_c_ADwin = (hysl_float_t*) calloc((size_t) CNodes->OrderADwin, sizeof(hysl_float_t));
    }

    /* Copy the older coupling force. This is necessary for simulations */
    int32_t pos = 0;
    for (int32_t idx = 0; idx < CNodes->Order; idx++) {
        Recv[(2 * CNodes->Order) + idx] = CoupForce[CNodes->Array[idx] - 1];
        if ((CNodes->Sub[idx].Type == EXP_ADWIN) && MultipleTypes) {
            Recv_ADwin[(2 * CNodes->OrderADwin) + pos] = CoupForce[CNodes->Array[idx] - 1];
            VecTdT0_c_ADwin[pos] = VecTdT0_c[idx];
            pos = pos + 1;
        }
    }

    for (int32_t idx = 0; idx < CNodes->Order; idx++) {
        switch (CNodes->Sub[idx].Type) {
        case SIM_EXACT_MDOF:
            /* This is the same case as SIM_MEASURED. All the simulated substructures are treated together
             * in the same routine.*/
        case SIM_EXACT_SDOF:
            /* This is the same case as SIM_MEASURED. All the simulated substructures are treated together
             * in the same routine.*/
        case SIM_EXACT_ESP:
            /* This is the same case as SIM_MEASURED. All the simulated substructures are treated together
             * in the same routine.*/
        case SIM_NEWMARK:
            /* This is the same case as SIM_MEASURED. All the simulated substructures are treated together
             * in the same routine.*/
        case SIM_BOUCWEN:
            /* This is the same case as SIM_MEASURED. All the simulated substructures are treated together
             * in the same routine.*/
        case SIM_UHYDE:
            /* This is the same case as SIM_MEASURED. All the simulated substructures are treated together
             * in the same routine.*/
        case SIM_STONEDRUMS:
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
            if( !Called_ADwin ){
                ADwin_Substep_Pre( VecTdT0_c_ADwin, (uint32_t) CNodes->OrderADwin );
                Called_ADwin = true;
            }
            break;
#endif
        case REMOTE:
            Remote = (Remote_t*) CNodes->Sub[idx].SimStruct;

            if ((Remote->Type == REMOTE_TCP) || (Remote->Type == REMOTE_UDP) || (Remote->Type == REMOTE_CELESTINA)) {
                Send = (hysl_float_t*) calloc(1u + (size_t) CNodes->Order, sizeof(hysl_float_t));
                for (int32_t jdx = 0; jdx < CNodes->Order; jdx++) {
                    Send[jdx] = VecTdT0_c[jdx];
                }
                Send[CNodes->Order] = GAcc;

                Substructure_Remote_Send(Remote->Socket, (uint32_t) CNodes->Order + 1u, sizeof(hysl_float_t), (char* const ) Send);

                Substructure_Remote_Receive(Remote->Socket, 3u * (uint32_t) CNodes->Order, sizeof(hysl_float_t), (char* const ) Recv);
                free(Send);
            } else if (Remote->Type == REMOTE_NSEP) {
                /* Using NSEP Protocol */
                Substructure_Remote_NSEP(Remote, NSEP_SEND_CMD, Time, Remote->NSub, Send, Recv);
                /* Receive the force from the PNSE server. WhatToDo = 2 */
                Substructure_Remote_NSEP(Remote, NSEP_REQUEST_CSIG, Time, 3 * Remote->NSub, Send, Recv);
            } else if (Remote->Type == REMOTE_OF) {
                /* Using OpenFresco */
                /* Send the trial response */
                Substructure_Remote_OpenFresco(Remote->Socket, OF_REMOTE_SET_TRIAL_RESPONSE, Remote->NSub, Send, NULL);
                /* Get the DAQ response */
                Substructure_Remote_OpenFresco(Remote->Socket, OF_REMOTE_GET_DAQ_RESPONSE, Remote->NSub, NULL, Recv);
            } else {
                assert((Remote->Type >= 0) || (Remote->Type < NUM_REMOTE_TYPE));
            }
            break;
        }
    }

#if _ADWIN_
    if( CNodes->OrderADwin > 0 ){
        ADwin_Substep_Post( (uint32_t) CNodes->OrderADwin, 75.0, &Recv_ADwin[0], &Recv_ADwin[CNodes->OrderADwin], &Recv_ADwin[2*CNodes->OrderADwin] );
    }
#endif

    //#pragma omp parallel for
    pos = 0;
    for (int32_t idx = 0; idx < CNodes->Order; idx++) {
        if ((CNodes->Sub[idx].Type == EXP_ADWIN) && MultipleTypes) {
            VecTdT[CNodes->Array[idx] - 1] = Recv_ADwin[pos];
            CoupForcePrev[idx] = Recv_ADwin[CNodes->OrderADwin + pos];
            CoupForce[CNodes->Array[idx] - 1] = Recv_ADwin[(2 * CNodes->OrderADwin) + pos];
            pos = pos + 1;
        } else {
            VecTdT[CNodes->Array[idx] - 1] = Recv[idx];
            CoupForcePrev[idx] = Recv[CNodes->Order + idx];
            CoupForce[CNodes->Array[idx] - 1] = Recv[(2 * CNodes->Order) + idx];
        }
    }

    free(Recv);
    if ((CNodes->OrderADwin >= 1) && MultipleTypes) {
        free(Recv_ADwin);
        free(VecTdT0_c_ADwin);
    }
}

void Substructure_Simulate (const hysl_float_t *IGain, const hysl_float_t *const VecTdT0_c, const hysl_float_t GAcc, const uint32_t NSubstep, const hysl_float_t DeltaT_Sub,
        CouplingNode_t *const CNodes,
        hysl_float_t *const VecTdT_c, hysl_float_t *const CoupForcePrev_c, hysl_float_t *const CoupForce_c) {

    int32_t incx = 1;
    int32_t incy = 1;

    ExactSim_t *Exact;
    ExactSimESP_t *ExactEsp;
    NewmarkSim_t *Newmark;
    BoucWen_t *BoucWen;
    UHYDEfbrSim_t *UHYDE;
    MeasuredSim_t *Measured;
    StoneDrums_t *StoneDrums, *StoneDrums_Prev;

    int32_t Length = CNodes->Order;
    hysl_float_t One = 1.0;

    for (uint32_t Substep = 1u; Substep <= NSubstep; Substep++) {
        /* Backup data so that CoupForcePrev_c contains always the last coupling force */
        hysl_copy(&Length, CoupForce_c, &incx, CoupForcePrev_c, &incy);

        hysl_float_t ramp = (hysl_float_t) Substep / (hysl_float_t) NSubstep;

#if _FLOAT_
        hysl_float_t ramp0 = 1.0f - ramp;
#else
        hysl_float_t ramp0 = 1.0 - ramp;
#endif

        if ((CNodes->Order - CNodes->OrderADwin) > 1) {
            hysl_copy(&Length, CNodes->VecTdT0_c0, &incx, VecTdT_c, &incy);
            hysl_scal(&Length, &ramp0, VecTdT_c, &incx);
            hysl_axpy(&Length, &ramp, VecTdT0_c, &incx, VecTdT_c, &incy);

            char uplo = 'L';
            hysl_symv(&uplo, &Length, &One, IGain, &Length, CoupForce_c, &incx, &One, VecTdT_c, &incy);
        } else {
            VecTdT_c[0] = (ramp0 * CNodes->VecTdT0_c0[0]) + (ramp * VecTdT0_c[0]) + (IGain[0] * CoupForce_c[0]);
        }

        /* Compute the new CoupForce_c */
        for (uint32_t idx = 0; idx < (uint32_t) CNodes->Order; idx++) {
            switch (CNodes->Sub[idx].Type) {
            case SIM_EXACT_MDOF:
                Exact = (ExactSim_t*) CNodes->Sub[idx].SimStruct;
                Substructure_ExactSolutionMDOF(VecTdT_c[idx], ramp, GAcc, DeltaT_Sub, Exact, &CoupForce_c[idx]);
                break;
            case SIM_EXACT_SDOF:
                Exact = (ExactSim_t*) CNodes->Sub[idx].SimStruct;
                Substructure_ExactSolutionSDOF(VecTdT_c[idx], ramp, GAcc, DeltaT_Sub, Exact, &CoupForce_c[idx]);
                break;
            case SIM_EXACT_ESP:
                ExactEsp = (ExactSimESP_t*) CNodes->Sub[idx].SimStruct;
                Substructure_ExactSolutionESP_SDOF(VecTdT_c[idx], ramp, GAcc, DeltaT_Sub, ExactEsp, &CoupForce_c[idx]);
                break;
            case SIM_NEWMARK:
                Newmark = (NewmarkSim_t*) CNodes->Sub[idx].SimStruct;
                Substructure_Newmark_SDOF(VecTdT_c[idx], ramp, GAcc, Newmark, &CoupForce_c[idx]);
                break;
            case SIM_BOUCWEN:
                BoucWen = (BoucWen_t*) CNodes->Sub[idx].SimStruct;
                Substructure_BoucWen(VecTdT_c[idx], BoucWen, &CoupForce_c[idx]);
                break;
            case SIM_UHYDE:
                UHYDE = (UHYDEfbrSim_t*) CNodes->Sub[idx].SimStruct;
                Substructure_SimUHYDE_1D(VecTdT_c[idx], DeltaT_Sub, UHYDE, &CoupForce_c[idx]);
                break;
            case SIM_MEASURED:
                Measured = (MeasuredSim_t*) CNodes->Sub[idx].SimStruct;
                Substructure_SimMeasured(Measured, &CoupForce_c[idx]);
                break;
            case SIM_STONEDRUMS:
                StoneDrums = (StoneDrums_t*) CNodes->Sub[idx].SimStruct;

                // printf("%d %lE %lE %lE %lE\n", Substep, StoneDrums->AccTdT, StoneDrums->VelTdT, VecTdT_c[i], CoupForce_c[i]);
                break;
            }
        }
    }

    /* Backup VecTdT0_c */
    hysl_copy(&Length, VecTdT0_c, &incx, CNodes->VecTdT0_c0, &incy);
    for (uint32_t idx = 0u; idx < (uint32_t) CNodes->Order; idx++) {
        switch (CNodes->Sub[idx].Type) {
        case SIM_EXACT_MDOF:
            /* Same as SIM_EXACT_SDOF */
        case SIM_EXACT_SDOF:
            Exact->Acc0 = Exact->AccT;
            Exact->AccT = Exact->AccTdT;
            Exact->Vel0 = Exact->VelT;
            Exact->VelT = Exact->VelTdT;
            Exact->Disp0 = Exact->DispT;
            Exact->DispT = VecTdT_c[idx];
            break;
        case SIM_EXACT_ESP:
            ExactEsp->Acc0 = ExactEsp->AccT;
            ExactEsp->AccT = ExactEsp->AccTdT;
            ExactEsp->Vel0 = ExactEsp->VelT;
            ExactEsp->VelT = ExactEsp->VelTdT;
            ExactEsp->Disp0 = ExactEsp->DispT;
            ExactEsp->DispT = VecTdT_c[idx];
            break;
        case SIM_NEWMARK:
            Newmark->Acc0 = Newmark->AccT;
            Newmark->AccT = Newmark->AccTdT;
            Newmark->Vel0 = Newmark->VelT;
            Newmark->VelT = Newmark->VelTdT;
            Newmark->Disp0 = Newmark->DispT;
            Newmark->DispT = VecTdT_c[idx];
            break;
        case SIM_BOUCWEN:
            BoucWen->DispT = VecTdT_c[idx];
            break;
        case SIM_UHYDE:
            break;
        case SIM_MEASURED:
            break;
        }
    }
}
