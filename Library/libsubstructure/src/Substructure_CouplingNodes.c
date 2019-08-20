#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>

#include "MatrixVector.h"
#include "Algorithm_Aux.h"
#include "Print_Messages.h"
#include "Substructure_BoucWen.h"
#include "Substructure_Exact.h"
#include "Substructure_Newmark.h"
#include "Substructure_Remote.h"
#include "Substructure_UHYDEfbr.h"
#include "Substructure_SimMeasured.h"
#include "Substructure_Experimental.h"
#include "Substructure_CouplingNodes.h"
#include "Substructure_StoneDrums.h"

#include "Definitions.h"

const char *Substructure_Type[] = {
        "Sim_Exact_MDOF",
        "Sim_Exact_SDOF",
        "Sim_Exact_ESP",
        "Sim_Newmark",
        "Sim_BoucWen",
        "Sim_BoucWen_Surface",
        "Sim_UHYDEfbr",
        "Sim_Measured",
        "Sim_StoneDrums",
        "Exp_ADwin",
        "Remote"
};

void Substructure_ReadCouplingNodes (const AlgConst_t *const InitCnt, CouplingNode_t *const CNodes) {
    FILE *InFile;
    int Count_Type;
    int i, j, k;
    int itemp, ndof;
    int intArray[3];
    hysl_float_t *ftemp, *rayleigh;
    char *ctemp;
    MatrixVector_t mass, stiff;
    char Type[MAX_SUBTYPE], Description[MAX_DESCRIPTION], FileMeas[MAX_FILENAME];
    char RemoteType[MAX_SUBTYPE];
    char InLine[MAX_LINE];
    char IPAddress[20], Port[20];
    char Account_Name[20], Account_Password[20];
    hysl_float_t DeltaTSub;

    DeltaTSub = InitCnt->Delta_t / (hysl_float_t) InitCnt->NSubstep;

    InFile = fopen(InitCnt->FileCNodes, "r");

    if (InFile == NULL) {
        Print_Header( ERROR);
        fprintf( stderr, "Substructure_ReadCouplingNodes: could not open file %s.\n", InitCnt->FileCNodes);
        exit( EXIT_FAILURE);
    }

    /* Number of ADwin coupling nodes is always 0 at the beginning */
    CNodes->OrderADwin = 0;

    /* The first value should be the number of Coupling nodes */
    fscanf(InFile, "%i", &CNodes->Order);
    fgets(InLine, MAX_LINE, InFile);

    if (CNodes->Order != InitCnt->OrderSub) {
        fclose(InFile);
        Print_Header( ERROR);
        fprintf( stderr, "Substructure_ReadCouplingNodes: Invalid number of substructures.\n");
        exit( EXIT_FAILURE);
    }

    /* Allocate the necessary memory */
    CNodes->Array = (int*) calloc((size_t) CNodes->Order, sizeof(int));
    CNodes->Sub = (Substructure_t*) malloc((size_t) CNodes->Order * sizeof(Substructure_t));
    CNodes->VecTdT0_c0 = (hysl_float_t*) calloc((size_t) CNodes->Order, sizeof(hysl_float_t));

    /* Read the contents of the file */
    i = 0;
    while (i < CNodes->Order) {

        /* Read until the coma */
        fscanf(InFile, "%[^,], %d", Type, &Count_Type);

        /* Check if the number of sub-structures is still valid */
        if ((i + Count_Type) > CNodes->Order) {
            fclose(InFile);
            Print_Header( ERROR);
            fprintf( stderr, "Substructure_ReadCouplingNodes: Number of substructures exceeded.\n");
            Print_Header( ERROR);
            fprintf( stderr, "Substructure_ReadCouplingNodes: Specified %d but read %d.\n", CNodes->Order, i + Count_Type);
            exit( EXIT_FAILURE);
        }
        Substructure_Identify(Type, &CNodes->Sub[i].Type);

        for (j = 0; j < Count_Type; j++) {
            fscanf(InFile, "%d", &CNodes->Array[i + j]);
            if (j > 0) {
                /* Only copy the values if j > 0 */
                CNodes->Sub[i + j].Type = CNodes->Sub[i].Type;
            }
            for (k = 0; k < (i + j); k++) {
                if (CNodes->Array[i + j] == CNodes->Array[k]) {
                    Print_Header( ERROR);
                    fprintf( stderr, "Substructure_ReadCouplingNodes: There is already a substructure assigned to coupling node %d.\n", CNodes->Array[k]);
                    Print_Header( ERROR);
                    fprintf( stderr, "Substructre_ReadCouplingNodes: Error when reading line %d.\n", i + 1);
                    exit( EXIT_FAILURE);
                }
            }
        }
        switch (CNodes->Sub[i].Type) {
        case SIM_EXACT_MDOF:
            /* Ignore coma */
            fscanf(InFile, "%*[,] %i", &itemp);
            if (itemp != EXACTMDOF_NUMPARAM_INIT) {
                Print_Header( ERROR);
                fprintf( stderr, "Wrong number of parameters for the substructue number %i of type Exact_MDOF.\n", i);
                fprintf( stderr, "The number of init parameters should be %i\n", EXACTMDOF_NUMPARAM_INIT);
                exit( EXIT_FAILURE);
            } else {
                ftemp = NULL;
                ctemp = NULL;
                ftemp = (hysl_float_t*) calloc((size_t) EXACTMDOF_NUMPARAM_INIT, sizeof( hysl_float_t ));
                ctemp = (char*) calloc((size_t) 20, sizeof(char));
                rayleigh = (hysl_float_t*) calloc(2, sizeof(hysl_float_t));

                /* Read the number of degrees of freedom */
                fscanf(InFile, "%i", &ndof);
                MatrixVector_Create(ndof, ndof, &mass);
                MatrixVector_Create(ndof, ndof, &stiff);

                /* Read the matrices from a MM file */
                fscanf(InFile, "%s", ctemp);
                MatrixVector_FromFile_MM(ctemp, &mass);
                fscanf(InFile, "%s", ctemp);
                MatrixVector_FromFile_MM(ctemp, &stiff);
#if _FLOAT_
		    fscanf( InFile, "%f %f", &rayleigh[0], &rayleigh[1] );
#else
                fscanf(InFile, "%lf %lf", &rayleigh[0], &rayleigh[1]);
#endif

                /* Read the optional description */
                Substructure_GetDescription(InFile, i, Description);

                for (j = 0; j < Count_Type; j++) {
                    CNodes->Sub[i + j].SimStruct = (void*) malloc(sizeof(ExactSim_t));

                    /*			 mass.Array[0] = 600.0; mass.Array[1] = 0.0; mass.Array[2] = 0.0; mass.Array[3] = 285.0;
                     stiff.Array[0] = 150000; stiff.Array[1] = -68000.0; stiff.Array[2] = -68000.0; stiff.Array[3] = 68000.0;
                     ray[0] = 0.900; ray[1] = 0.000015;*/
                    Substructure_ExactSolutionMDOF_Init(mass.Array, stiff.Array, ndof, InitCnt->Rayleigh.Alpha, InitCnt->Rayleigh.Beta, rayleigh[0], rayleigh[1], InitCnt->a0, InitCnt->a2, InitCnt->a3,
                            InitCnt->a6, InitCnt->a7, Description, (ExactSim_t*) CNodes->Sub[i + j].SimStruct);

                    /*Substructure_ExactSolutionSDOF_Init( ftemp[0], ftemp[1], ftemp[2], InitCnt->a0, InitCnt->a2, InitCnt->a3, InitCnt->a6, InitCnt->a7, Description,
                     (ExactSim_t *) CNodes->Sub[i + j].SimStruct );*/
                    Print_Header( INFO);
                    printf("Simulating the substructure in the coupling node %d as an exact integration method (MDOF).\n", CNodes->Array[i + j]);
                }
                MatrixVector_Destroy(&mass);
                MatrixVector_Destroy(&stiff);
                free(rayleigh);
                free(ftemp);
                free(ctemp);
            }
            break;
        case SIM_EXACT_SDOF:
            /* Ignore coma */
            fscanf(InFile, "%*[,] %i", &itemp);
            if (itemp != EXACTSDOF_NUMPARAM_INIT) {
                Print_Header( ERROR);
                fprintf( stderr, "Wrong number of parameters for the substructue number %i of type Exact_SDOF.\n", i);
                fprintf( stderr, "The number of init parameters should be %i\n", EXACTSDOF_NUMPARAM_INIT);
                exit( EXIT_FAILURE);
            } else {
                ftemp = NULL;
                ftemp = (hysl_float_t*) calloc((size_t) EXACTSDOF_NUMPARAM_INIT, sizeof( hysl_float_t ));

                /* Read the input parameters */
                for (j = 0; j < EXACTSDOF_NUMPARAM_INIT; j++) {
#if _FLOAT_
			 fscanf( InFile, "%f", &ftemp[j] );
#else
                    fscanf(InFile, "%lf", &ftemp[j]);
#endif
                }

                /* Read the optional description */
                Substructure_GetDescription(InFile, i, Description);

                for (j = 0; j < Count_Type; j++) {
                    CNodes->Sub[i + j].SimStruct = (void*) malloc(sizeof(ExactSim_t));

                    Substructure_ExactSolutionSDOF_Init(ftemp[0], ftemp[1], ftemp[2], InitCnt->a0, InitCnt->a2, InitCnt->a3, InitCnt->a6, InitCnt->a7, Description,
                            (ExactSim_t*) CNodes->Sub[i + j].SimStruct);

                    Print_Header( INFO);
                    printf("Simulating the substructure in the coupling node %d as an exact integration method (SDOF).\n", CNodes->Array[i + j]);
                }
                free(ftemp);
            }
            break;
        case SIM_EXACT_ESP:
            /* Ignore coma */
            fscanf(InFile, "%*[,] %i", &itemp);
            if (itemp != EXACTSDOF_NUMPARAM_INIT) {
                Print_Header( ERROR);
                fprintf( stderr, "Wrong number of parameters for the substructue number %i of type Exact_ESP.\n", i);
                fprintf( stderr, "The number of init parameters should be %i\n", EXACTSDOF_NUMPARAM_INIT);
                exit( EXIT_FAILURE);
            } else {
                ftemp = NULL;
                ftemp = (hysl_float_t*) calloc((size_t) EXACTSDOF_NUMPARAM_INIT, sizeof( hysl_float_t ));

                /* Read the input parameters */
                for (j = 0; j < EXACTSDOF_NUMPARAM_INIT; j++) {
#if _FLOAT_
			 fscanf( InFile, "%f", &ftemp[j] );
#else
                    fscanf(InFile, "%lf", &ftemp[j]);
#endif
                }

                /* Read the optional description */
                Substructure_GetDescription(InFile, i, Description);

                for (j = 0; j < Count_Type; j++) {
                    CNodes->Sub[i + j].SimStruct = (void*) calloc((size_t) 1, sizeof(ExactSimESP_t));
                    if (CNodes->Sub[i + j].SimStruct == NULL) {
                        exit(EXIT_FAILURE);
                    }

                    Substructure_ExactSolutionESP_Init(ftemp[0], ftemp[1], ftemp[2], DeltaTSub, Description, (ExactSimESP_t*) CNodes->Sub[i + j].SimStruct);

                    Print_Header( INFO);
                    printf("Simulating the substructure in the coupling node %d as an exact integration method (ESP).\n", CNodes->Array[i + j]);
                }
                free(ftemp);
            }
            break;
        case SIM_NEWMARK:
            /* Ignore coma */
            fscanf(InFile, "%*[,] %i", &itemp);
            if (itemp != NEWMARK_NUMPARAM_INIT) {
                Print_Header( ERROR);
                fprintf( stderr, "Wrong number of parameters for the substructue number %i of type Exact_ESP.\n", i);
                fprintf( stderr, "The number of init parameters should be %i\n", NEWMARK_NUMPARAM_INIT);
                exit( EXIT_FAILURE);
            } else {
                ftemp = NULL;
                ftemp = (hysl_float_t*) calloc((size_t) NEWMARK_NUMPARAM_INIT, sizeof( hysl_float_t ));

                /* Read the input parameters */
                for (j = 0; j < NEWMARK_NUMPARAM_INIT; j++) {
#if _FLOAT_
			 fscanf( InFile, "%f", &ftemp[j] );
#else
                    fscanf(InFile, "%lf", &ftemp[j]);
#endif
                }

                /* Read the optional description */
                Substructure_GetDescription(InFile, i, Description);

                for (j = 0; j < Count_Type; j++) {
                    CNodes->Sub[i + j].SimStruct = (void*) calloc((size_t) 1, sizeof(NewmarkSim_t));
                    if (CNodes->Sub[i + j].SimStruct == NULL) {
                        exit(EXIT_FAILURE);
                    }

                    Substructure_Newmark_Init(ftemp[0], ftemp[1], ftemp[2], DeltaTSub, InitCnt->Delta_t, InitCnt->TIntConst.Beta, InitCnt->TIntConst.Gamma, Description,
                            (NewmarkSim_t*) CNodes->Sub[i + j].SimStruct);

                    Print_Header( INFO);
                    printf("Simulating the substructure in the coupling node %d as an exact integration method (ESP).\n", CNodes->Array[i + j]);
                }
                free(ftemp);
            }
            break;
        case SIM_BOUCWEN:
            /* Ignore coma */
            fscanf(InFile, "%*[,] %i", &itemp);
            if ((itemp != BOUCWEN_NUMPARAM_INIT) && (itemp != BOUCWENDEG_NUMPARAM_INIT) && (itemp != BOUCWENBABERNOORI_NUMPARAM_INIT)) {
                Print_Header( ERROR);
                fprintf( stderr, "Wrong number of parameters for the substructue number %i of type Exact_ESP.\n", i);
                fprintf( stderr, "The number of init parameters should be %i for classic Bouc-Wen, %i for Bouc-Wen with material degradation and %i for the Bouc-Wen-Baber-Noori model\n",
                        BOUCWEN_NUMPARAM_INIT, BOUCWENDEG_NUMPARAM_INIT, BOUCWENBABERNOORI_NUMPARAM_INIT);
                exit( EXIT_FAILURE);
            } else {
                ftemp = NULL;
                ftemp = (hysl_float_t*) calloc((size_t) BOUCWENBABERNOORI_NUMPARAM_INIT, sizeof( hysl_float_t ));

                /* Read the input parameters */
                for (j = 0; j < itemp; j++) {
#if _FLOAT_
			 fscanf( InFile, "%f", &ftemp[j] );
#else
                    fscanf(InFile, "%lf", &ftemp[j]);
#endif
                }

                /* Read the optional description */
                Substructure_GetDescription(InFile, i, Description);

                for (j = 0; j < Count_Type; j++) {
                    CNodes->Sub[i + j].SimStruct = (void*) calloc((size_t) 1, sizeof(BoucWen_t));
                    if (CNodes->Sub[i + j].SimStruct == NULL) {
                        exit(EXIT_FAILURE);
                    }

                    if ( itemp == BOUCWEN_NUMPARAM_INIT ) {
                        Substructure_BoucWen_Init( ftemp[0], ftemp[1], ftemp[2], ftemp[3], ftemp[4], ftemp[5],       // Regular Bouc-Wen
                                1.0, 0.0, 1.0, 0.0, 1.0, 0.0,// Material degradation
                                0.0, 0.0, 0.0, 0.0, 0.0, 0.0,// Pitching
                                BOUC_WEN, Description, (BoucWen_t *) CNodes->Sub[i + j].SimStruct );
                    } else if ( itemp == BOUCWENDEG_NUMPARAM_INIT ) {
                        Substructure_BoucWen_Init( ftemp[0], ftemp[1], ftemp[2], ftemp[3], ftemp[4], ftemp[5],       // Regular Bouc-Wen
                                ftemp[6], ftemp[7], ftemp[8], ftemp[9], ftemp[10], ftemp[11],// Material degradation
                                0.0, 0.0, 0.0, 0.0, 0.0, 0.0,// Pitching
                                BOUC_WEN_DEG, Description, (BoucWen_t *) CNodes->Sub[i + j].SimStruct );
                    } else if ( itemp == BOUCWENBABERNOORI_NUMPARAM_INIT ) {
                        Substructure_BoucWen_Init( ftemp[0], ftemp[1], ftemp[2], ftemp[3], ftemp[4], ftemp[5],       // Regular Bouc-Wen
                                ftemp[6], ftemp[7], ftemp[8], ftemp[9], ftemp[10], ftemp[11],// Material degradation
                                ftemp[12], ftemp[13], ftemp[14], ftemp[15], ftemp[16], ftemp[17],// Pitching
                                BOUC_WEN_BABER_NOORI, Description, (BoucWen_t *) CNodes->Sub[i + j].SimStruct );
                    } else assert((itemp == BOUCWEN_NUMPARAM_INIT) && (itemp == BOUCWENDEG_NUMPARAM_INIT) && (itemp == BOUCWENBABERNOORI_NUMPARAM_INIT));

                    Print_Header( INFO);
                    printf("Simulating the substructure in the coupling node %d as Bouc-Wen.\n", CNodes->Array[i + j]);
                }
                free(ftemp);
            }
            break;
        case SIM_BOUCWEN_SURFACE:
            /* Ignore coma */
            fscanf(InFile, "%*[,] %i", &itemp);
            if (itemp != BOUCWENDEG_NUMPARAM_INIT) {
                Print_Header( ERROR);
                fprintf( stderr, "Wrong number of parameters for the substructue number %i of type Exact_ESP.\n", i);
                fprintf( stderr, "The number of init parameters should be %i for classic Bouc-Wen, %i for Bouc-Wen with material degradation and %i for the Bouc-Wen-Baber-Noori model\n",
                        BOUCWENDEG_NUMPARAM_INIT);
                exit( EXIT_FAILURE);
            } else {
                ftemp = NULL;
                ftemp = (hysl_float_t*) calloc((size_t) itemp, sizeof( hysl_float_t ));

                /* Read the input parameters */
                for (j = 0; j < itemp; j++) {
#if _FLOAT_
			 fscanf( InFile, "%f", &ftemp[j] );
#else
                    fscanf(InFile, "%lf", &ftemp[j]);
#endif
                }

                /* Read the optional description */
                Substructure_GetDescription(InFile, i, Description);

                for (j = 0; j < Count_Type; j++) {
                    CNodes->Sub[i + j].SimStruct = (void*) calloc((size_t) 1, sizeof(BoucWenSurface_t));
                    if (CNodes->Sub[i + j].SimStruct == NULL) {
                        exit(EXIT_FAILURE);
                    }

                    Substructure_BoucWenSurface_Init(ftemp[0], ftemp[1], ftemp[2], ftemp[3], ftemp[4], ftemp[5], ftemp[6], ftemp[7], ftemp[8], ftemp[9], ftemp[10], ftemp[11], BOUC_WEN_DEG,
                            Description, (BoucWenSurface_t*) CNodes->Sub[i + j].SimStruct);

                    Print_Header( INFO);
                    printf("Simulating the substructure in the coupling node %d as Bouc-Wen (Surface).\n", CNodes->Array[i + j]);
                }
                free(ftemp);
            }
            break;
        case SIM_UHYDE:
            /* Ignore coma */
            fscanf(InFile, "%*[,] %i", &itemp);
            if (itemp != UHYDE_NUMPARAM_INIT) {
                Print_Header( ERROR);
                fprintf( stderr, "Wrong number of parameters for the substructue number %i of type UHYDE.\n", i);
                fprintf( stderr, "The number of init parameters should be %i\n", UHYDE_NUMPARAM_INIT);
                exit( EXIT_FAILURE);
            } else {
                ftemp = NULL;
                ftemp = (hysl_float_t*) calloc((size_t) UHYDE_NUMPARAM_INIT, sizeof( hysl_float_t ));

                for (j = 0; j < UHYDE_NUMPARAM_INIT; j++) {
#if _FLOAT_
			 fscanf( InFile, "%f", &ftemp[j] );
#else
                    fscanf(InFile, "%lf", &ftemp[j]);
#endif
                }

                /* Read the optional description */
                Substructure_GetDescription(InFile, i, Description);

                for (j = 0; j < Count_Type; j++) {
                    CNodes->Sub[i + j].SimStruct = (void*) malloc(sizeof(UHYDEfbrSim_t));
                    Substructure_SimUHYDE_1D_Init(ftemp[0], ftemp[1], ftemp[2], Description, (UHYDEfbrSim_t*) CNodes->Sub[i + j].SimStruct);
                    Print_Header( INFO);
                    printf("Simulating the substructure in the coupling node %d as a UHYDE-fbr device.\n", CNodes->Array[i + j]);
                }
                free(ftemp);
            }
            break;
        case SIM_MEASURED:
            fscanf(InFile, "%*[,] %[^,]", FileMeas);

            /* Read the optional description */
            Substructure_GetDescription(InFile, i, Description);

            for (j = 0; j < Count_Type; j++) {
                CNodes->Sub[i + j].SimStruct = (void*) malloc(sizeof(MeasuredSim_t));
                Substructure_SimMeasured_Init(FileMeas, InitCnt->NStep, InitCnt->NSubstep, Description, (MeasuredSim_t*) CNodes->Sub[i + j].SimStruct);
                Print_Header( INFO);
                printf("Simulating the substructure in the coupling node %d using time history measured forces.\n", CNodes->Array[i + j]);
            }
            break;
        case SIM_STONEDRUMS:
            fscanf(InFile, "%*[,] %i", &itemp);
            if (itemp != STONEDRUMS_NUMPARAM_INIT) {
                Print_Header( ERROR);
                fprintf( stderr, "Wrong number of parameters for the substructue number %i of type UHYDE.\n", i);
                fprintf( stderr, "The number of init parameters should be %i\n", STONEDRUMS_NUMPARAM_INIT);
                exit( EXIT_FAILURE);
            } else {
                ftemp = NULL;
                ftemp = (hysl_float_t*) calloc((size_t) STONEDRUMS_NUMPARAM_INIT, sizeof( hysl_float_t ));

                intArray[0] = BOUC_WEN_DEG;
                intArray[1] = BOUC_WEN_DEG;
                intArray[2] = BOUC_WEN_DEG;
                for (j = 0; j < STONEDRUMS_NUMPARAM_INIT; j++) {
#if _FLOAT_
			 fscanf( InFile, "%f", &ftemp[j] );
#else
                    fscanf(InFile, "%lf", &ftemp[j]);
#endif
                }

                /* Read the optional description */
                Substructure_GetDescription(InFile, i, Description);

                for (j = 0; j < Count_Type; j++) {
                    CNodes->Sub[i + j].SimStruct = (void*) malloc(sizeof(StoneDrums_t));

                    Substructure_StoneDrums_Init((int) ftemp[0], (int) ftemp[1], (int) ftemp[2], (int) ftemp[3], ftemp + 4, intArray, Description, (StoneDrums_t*) CNodes->Sub[i + j].SimStruct);
                    Print_Header( INFO);
                    printf("Simulating the substructure in the coupling node %d as a Stone Drum.\n", CNodes->Array[i + j]);
                }
                free(ftemp);
            }
            break;
        case EXP_ADWIN:
            CNodes->OrderADwin = CNodes->OrderADwin + Count_Type;
            /* Read the optional description */
            Substructure_GetDescription(InFile, i, Description);

            for (j = 0; j < Count_Type; j++) {
                CNodes->Sub[i + j].SimStruct = (void*) malloc(sizeof(ExpSub_t));
                Substructure_Experimental_Init(Description, (ExpSub_t*) CNodes->Sub[i + j].SimStruct);
                Print_Header( INFO);
                printf("The substructure in the coupling node %d is computed in ADwin.\n", CNodes->Array[i + j]);
            }
            break;
        case REMOTE:
            /* Read IP Address and Port */
            fscanf(InFile, "%*[,] %s %s %[^,]", RemoteType, IPAddress, Port);

            if (strcmp(RemoteType, "NSEP") == 0) {
                fscanf(InFile, "%*[,] %s %[^,]", Account_Name, Account_Password);
            }

            /* Read the optional description */
            Substructure_GetDescription(InFile, i, Description);

            for (j = 0; j < Count_Type; j++) {
                CNodes->Sub[i + j].SimStruct = (void*) malloc(sizeof(Remote_t));
                Substructure_Remote_Init(RemoteType, IPAddress, Port, Account_Name, Account_Password, (unsigned int) Count_Type, &CNodes->Array[i], Description,
                        (Remote_t*) CNodes->Sub[i + j].SimStruct);
                Print_Header( INFO);
                printf("The substructure in the coupling node %d is computed is computed at %s:%s using %s protocol.\n", CNodes->Array[i + j], IPAddress, Port, RemoteType);
            }
            break;
        }
        i = i + Count_Type;
    }

    /* Close the file */
    fclose(InFile);

    /* Sort the coupling nodes in ascending order */
    Substructure_SortCouplingNodes(CNodes);
}

void Substructure_SortCouplingNodes (CouplingNode_t *const CNodes) {

    int i, j, Tmp;
    Substructure_t TmpSub;

    for (i = 0; i < CNodes->Order; i++) {
        for (j = 0; j < CNodes->Order - 1; j++) {
            if (CNodes->Array[j] > CNodes->Array[j + 1]) {
                Tmp = CNodes->Array[j + 1];
                TmpSub = CNodes->Sub[j + 1];

                CNodes->Array[j + 1] = CNodes->Array[j];
                CNodes->Sub[j + 1] = CNodes->Sub[j];

                CNodes->Array[j] = Tmp;
                CNodes->Sub[j] = TmpSub;
            }
        }
    }
}

void Substructure_Identify (char *const Type, int *const Identity_Num) {
    int ID; /* A counter */
    bool Found = false;

    ID = 0;
    /* Identify with substructure is in Type. Exit when found */
    while (ID < NUM_TYPE_SUB && !Found) {
        if (strcmp(Substructure_Type[ID], Type) == 0) {
            Found = true;
        } else {
            ID = ID + 1;
        }
    }

    /* The substructure in Type is not supported.*/
    if (!Found) {
        Print_Header( ERROR);
        fprintf( stderr, "Substructure_Identify: The substructure type '%s' is not supported. Valid substructures are:\n", Type);
        for (ID = 0; ID < NUM_TYPE_SUB; ID++) {
            fprintf( stderr, "[......] %d) %s.\n", ID + 1, Substructure_Type[ID]);
        }
        exit( EXIT_FAILURE);
    } else {
        /* Assign the identity */
        (*Identity_Num) = ID;
    }
}

void Substructure_GetDescription (FILE *const InFile, const int LineNum, char *const Description) {

    fscanf(InFile, "%*[,]");
    fgets(Description, MAX_DESCRIPTION, InFile);

    if (Description[strlen(Description) - 1] != '\n' && !feof(InFile)) {
        Print_Header( ERROR);
        fprintf( stderr, "Substructure_ReadCouplingNodes: Maximum description length (%d) exceeded in line %d.\n", MAX_DESCRIPTION, LineNum + 1);
        exit( EXIT_FAILURE);
    }

    if (Description[strlen(Description) - 2] != ';' && Description[strlen(Description) - 1] != ';') {
        Print_Header( ERROR);
        fprintf( stderr, "Substructure_ReadCouplingNodes: Line number %d should terminate with ';'.\n", LineNum + 1);
        exit( EXIT_FAILURE);
    }
}

void Substructure_DeleteCouplingNodes (CouplingNode_t *CNodes) {
    int i;

    for (i = 0; i < CNodes->Order; i++) {
        switch (CNodes->Sub[i].Type) {
        case SIM_EXACT_MDOF:
            Substructure_ExactSolutionMDOF_Destroy((ExactSim_t*) CNodes->Sub[i].SimStruct);
            break;
        case SIM_EXACT_SDOF:
            Substructure_ExactSolutionSDOF_Destroy((ExactSim_t*) CNodes->Sub[i].SimStruct);
            break;
        case SIM_EXACT_ESP:
            Substructure_ExactSolutionESP_Destroy((ExactSimESP_t*) CNodes->Sub[i].SimStruct);
            break;
        case SIM_NEWMARK:
            Substructure_Newmark_Destroy((NewmarkSim_t*) CNodes->Sub[i].SimStruct);
            break;
        case SIM_BOUCWEN:
            Substructure_BoucWen_Destroy((BoucWen_t*) CNodes->Sub[i].SimStruct);
            break;
        case SIM_BOUCWEN_SURFACE:
            Substructure_BoucWenSurface_Destroy((BoucWenSurface_t*) CNodes->Sub[i].SimStruct);
        case SIM_UHYDE:
            Substructure_SimUHYDE_Destroy((UHYDEfbrSim_t*) CNodes->Sub[i].SimStruct);
            break;
        case SIM_MEASURED:
            Substructure_SimMeasured_Destroy((MeasuredSim_t*) CNodes->Sub[i].SimStruct);
            break;
        case SIM_STONEDRUMS:
            Substructure_StoneDrums_Destroy((StoneDrums_t*) CNodes->Sub[i].SimStruct);
            break;
        case EXP_ADWIN:
            Substructure_Experimental_Destroy((ExpSub_t*) CNodes->Sub[i].SimStruct);
            break;
        case REMOTE:
            Substructure_Remote_Destroy((Remote_t*) CNodes->Sub[i].SimStruct);
            break;
        }
        free(CNodes->Sub[i].SimStruct);
    }

    CNodes->Order = 0;
    free(CNodes->Sub);
    free(CNodes->Array);
    free(CNodes->VecTdT0_c0);
}
