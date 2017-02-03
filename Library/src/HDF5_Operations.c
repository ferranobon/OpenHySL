#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "HDF5_Operations.h"
#include "MatrixVector.h"
#include "Algorithm_Aux.h"

#include "Print_Messages.h"

#include "Definitions.h"
#include "Substructure_CouplingNodes.h"
#include "Substructure_Exact.h"
#include "Substructure_BoucWen.h"
#include "Substructure_Newmark.h"
#include "Substructure_UHYDEfbr.h"
#include "Substructure_StoneDrums.h"
#include "Substructure_SimMeasured.h"
#include "Substructure_Experimental.h"

#if _ADWIN_
#include <libadwin.h>        /* ADwin routines Boot(), Set_DeviceNo(), ... */
#include <libadwin/errno.h>  /* ADwin error handling */
#include "ADwin_Routines.h"
#endif

/* hdf5 header files */
#include "hdf5.h"
#include "hdf5_hl.h"

hid_t HDF5_CreateFile( const char *Filename )
{

     hid_t    file_id;

     /* Create a new file using default properties. */
     file_id = H5Fcreate( Filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
     
     return file_id;
}

void HDF5_CreateGroup_Parameters( const hid_t hdf5_file, const AlgConst_t *const InitCnt,			  
				  const CouplingNode_t *const CNodes, const HYSL_FLOAT *const Acc1,
				  const HYSL_FLOAT *const Vel1, const HYSL_FLOAT *const Disp1,
				  const HYSL_FLOAT *const Acc2, const HYSL_FLOAT *const Vel2,
				  const HYSL_FLOAT *const Disp2, const HYSL_FLOAT *const Acc3,
				  const HYSL_FLOAT *const Vel3, const HYSL_FLOAT *const Disp3 )
{
     hid_t    group_id, dataset_id, dataspace_id;
     hsize_t  dims[2];
     herr_t   status;
     HYSL_FLOAT   *dArray;
     int      *iArray, i;
     char     **Entry_Names;
     
     Entry_Names = (char **) malloc( (size_t) 6*sizeof(char*) );
     dArray = (HYSL_FLOAT *) malloc( (size_t) 4*sizeof(HYSL_FLOAT) );

     /* Create a group named "/Time History" in the file. */
     group_id = H5Gcreate( hdf5_file, "/Test Parameters", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );     
     status = H5Gclose( group_id );
     if( status < 0 ){
	  Print_Header( ERROR );
	  fprintf( stderr, "Error creating group Test Parameters in HDF5 file.\n" );
	  exit( EXIT_FAILURE );
     }
     
     /* Add input load */
     group_id = H5Gcreate( hdf5_file, "/Test Parameters/Input X", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
     dims[0] = InitCnt->NStep;
     dims[1] = 1;
     dataspace_id = H5Screate_simple( 2, dims, NULL );
          if( Acc1 != NULL ){
	  dataset_id = H5Dcreate ( group_id, "/Test Parameters/Input X/Acceleration", H5T_NATIVE_HYSL_FLOAT,
				   dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	  status = H5LTset_attribute_string( hdf5_file, "Test Parameters/Input X/Acceleration", "Units", "m/s^2" );
	  status = H5Dwrite( dataset_id, H5T_NATIVE_HYSL_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, Acc1 );
	  status = H5Dclose( dataset_id );
     }

     if( Vel1 != NULL ){
	  dataset_id = H5Dcreate ( group_id, "/Test Parameters/Input X/Velocity", H5T_NATIVE_HYSL_FLOAT,
				   dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	  status = H5LTset_attribute_string( hdf5_file, "Test Parameters/Input X/Velocity", "Units", "m/s" );
	  status = H5Dwrite( dataset_id, H5T_NATIVE_HYSL_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, Vel1 );
	  status = H5Dclose( dataset_id );
     }

     if( Disp1 != NULL ){
	  dataset_id = H5Dcreate ( group_id, "/Test Parameters/Input X/Displacement", H5T_NATIVE_HYSL_FLOAT,
				   dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	  status = H5LTset_attribute_string( hdf5_file, "Test Parameters/Input X/Displacement", "Units", "m" );
	  status = H5Dwrite( dataset_id, H5T_NATIVE_HYSL_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, Disp1 );
	  status = H5Dclose( dataset_id );
     }
     status = H5Gclose( group_id );

     group_id = H5Gcreate( hdf5_file, "/Test Parameters/Input Y", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
     dims[0] = InitCnt->NStep;
     dims[1] = 1;
     dataspace_id = H5Screate_simple( 2, dims, NULL );
     if( Acc2 != NULL ){
	  dataset_id = H5Dcreate ( group_id, "/Test Parameters/Input Y/Acceleration", H5T_NATIVE_HYSL_FLOAT,
				   dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	  status = H5LTset_attribute_string( hdf5_file, "Test Parameters/Input Y/Acceleration", "Units", "m/s^2" );
	  status = H5Dwrite( dataset_id, H5T_NATIVE_HYSL_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, Acc2 );
	  status = H5Dclose( dataset_id );
     }

     if( Vel2 != NULL ){
	  dataset_id = H5Dcreate ( group_id, "/Test Parameters/Input Y/Velocity", H5T_NATIVE_HYSL_FLOAT,
				   dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	  status = H5LTset_attribute_string( hdf5_file, "Test Parameters/Input Y/Velocity", "Units", "m/s" );
	  status = H5Dwrite( dataset_id, H5T_NATIVE_HYSL_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, Vel2 );
	  status = H5Dclose( dataset_id );
     }

     if( Disp2 != NULL ){
	  dataset_id = H5Dcreate ( group_id, "/Test Parameters/Input Y/Displacement", H5T_NATIVE_HYSL_FLOAT,
				   dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	  status = H5LTset_attribute_string( hdf5_file, "Test Parameters/Input Y/Displacement", "Units", "m" );
	  status = H5Dwrite( dataset_id, H5T_NATIVE_HYSL_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, Disp2 );
	  status = H5Dclose( dataset_id );
     }
     status = H5Gclose( group_id );

     group_id = H5Gcreate( hdf5_file, "/Test Parameters/Input Z", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
     dims[0] = InitCnt->NStep;
     dims[1] = 1;
     dataspace_id = H5Screate_simple( 2, dims, NULL );
     if( Acc3 != NULL ){
	  dataset_id = H5Dcreate ( group_id, "/Test Parameters/Input Z/Acceleration", H5T_NATIVE_HYSL_FLOAT,
				   dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	  status = H5LTset_attribute_string( hdf5_file, "Test Parameters/Input Z/Acceleration", "Units", "m/s^2" );
	  status = H5Dwrite( dataset_id, H5T_NATIVE_HYSL_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, Acc2 );
	  status = H5Dclose( dataset_id );
     }

     if( Vel3 != NULL ){
	  dataset_id = H5Dcreate ( group_id, "/Test Parameters/Input Z/Velocity", H5T_NATIVE_HYSL_FLOAT,
				   dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	  status = H5LTset_attribute_string( hdf5_file, "Test Parameters/Input Z/Velocity", "Units", "m/s" );
	  status = H5Dwrite( dataset_id, H5T_NATIVE_HYSL_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, Vel2 );
	  status = H5Dclose( dataset_id );
     }

     if( Disp3 != NULL ){
	  dataset_id = H5Dcreate ( group_id, "/Test Parameters/Input Z/Displacement", H5T_NATIVE_HYSL_FLOAT,
				   dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	  status = H5LTset_attribute_string( hdf5_file, "Test Parameters/Input Z/Displacement", "Units", "m" );
	  status = H5Dwrite( dataset_id, H5T_NATIVE_HYSL_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, Disp2 );
	  status = H5Dclose( dataset_id );
     }
     status = H5Gclose( group_id );

     /* Add Integration parameters */
     Entry_Names[0] = strdup( "Beta" );         dArray[0] = InitCnt->TIntConst.Beta;
     Entry_Names[1] = strdup( "Gamma" );        dArray[1] = InitCnt->TIntConst.Gamma;
     Entry_Names[2] = strdup( "Hilber Alpha" ); dArray[2] = InitCnt->TIntConst.HilberAlpha;
     Entry_Names[3] = strdup( "Time step" );    dArray[3] = InitCnt->Delta_t;
     
     HDF5_AddFloatArray_AsTable( hdf5_file, "/Test Parameters/Time Integration", Entry_Names, dArray, 4, 1 );
     free( Entry_Names[0] );
     free( Entry_Names[1] );
     free( Entry_Names[2] );
     free( Entry_Names[3] );

     /* Add Rayleigh damping values */
     Entry_Names[0] = strdup( "Alpha" ); dArray[0] = InitCnt->Rayleigh.Alpha;
     Entry_Names[1] = strdup( "Beta" );  dArray[1] = InitCnt->Rayleigh.Beta;
     HDF5_AddFloatArray_AsTable( hdf5_file, "/Test Parameters/Rayleigh", Entry_Names, dArray, 2, 1 );
     free( Entry_Names[0] ); free( Entry_Names[1] );

     /* Add PID parameters */
     Entry_Names[0] = strdup( "P" ); dArray[0] = InitCnt->PID.P;
     Entry_Names[1] = strdup( "I" ); dArray[1] = InitCnt->PID.I;
     Entry_Names[2] = strdup( "D" ); dArray[2] = InitCnt->PID.D;
     HDF5_AddFloatArray_AsTable( hdf5_file, "/Test Parameters/PID", Entry_Names, dArray, 3, 1 );
     for( i = 0; i < 3; i++ ){
	  free( Entry_Names[i] );
     }
     /* Free the HYSL_FLOAT array */
     free( dArray );

     iArray = (int *) malloc( (size_t) 6*sizeof(int) );

     Entry_Names[0] = strdup( "Num. DOF" );           iArray[0] = (int) InitCnt->Order;
     Entry_Names[1] = strdup( "Num. steps" );         iArray[1] = (int) InitCnt->NStep; 
     Entry_Names[2] = strdup( "Num. sub-steps" );     iArray[2] = (int) InitCnt->NSubstep;
     Entry_Names[3] = strdup( "Num. substructures" ); iArray[3] = (int) InitCnt->OrderSub;      
     HDF5_AddIntArray_AsTable( hdf5_file, "/Test Parameters/Misc", Entry_Names, iArray, 4 );
     free( iArray );
     
     /* Save the coupling nodes and information */
     if( InitCnt->OrderSub > 0 ){
	  Save_InformationCNodes( hdf5_file, "/Test Parameters/Substructures", InitCnt, CNodes );
     }

     for( i = 0; i < 4; i++ ){
	  free( Entry_Names[i] );
     }

     free( Entry_Names );
}


void Save_InformationCNodes( const hid_t file_id, const char *Name_path, const AlgConst_t *const InitCnt, const CouplingNode_t *const CNodes )
{
     int i, j, k, count, is_adwin, is_exact, is_newmark, is_uhyde, is_boucwen, is_stonedrum, is_measured;
     hid_t    group_id, strtype, memtype, space, dset;
     hsize_t  dims[1];
     herr_t   status;
     HDF5_Exact_UHYDE_t *Nodes;
     HDF5_Newmark_t *Nodes_Newmark;
     HDF5_Exp_Meas_t *Nodes_Exp;
     HDF5_BoucWen_t *Nodes_BoucWen;
     HDF5_StoneDrum_t *Nodes_StoneDrum;
     ExactSimESP_t *ExSim;
     NewmarkSim_t *Newmark;
     UHYDEfbrSim_t *UHYDE;
     ExpSub_t *Experimental;
     BoucWen_t *BoucWen;
     StoneDrums_t *StoneDrum;

   
     group_id = H5Gcreate( file_id, Name_path, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );     

     strtype = H5Tcopy( H5T_C_S1 );
     status = H5Tset_size( strtype, MAX_DESCRIPTION );

     if( status < 0 ){
	  Print_Header( ERROR );
	  fprintf( stderr, "Error during H5Tset_size().\n" );
	  exit( EXIT_FAILURE );
     }

     for( i = 0; i <= NUM_TYPE_SUB; i++ ){
	  count = 0;
	  is_adwin = 0; is_exact = 0; is_newmark = 0; is_uhyde=0; is_measured = 0;
	  is_boucwen = 0; is_stonedrum  = 0;
	  
	  for( j = 0; j < CNodes->Order; j++ ){
	       if ( i == EXP_ADWIN && CNodes->Sub[j].Type == EXP_ADWIN ){
		    count = count + 1;
		    is_adwin = 1;	
	       } else if ( i == SIM_EXACT_ESP && CNodes->Sub[j].Type == SIM_EXACT_ESP ){
		    count = count + 1;
		    is_exact = 1;
	       } else if ( i == SIM_NEWMARK && CNodes->Sub[j].Type == SIM_NEWMARK ){
		    count = count + 1;
		    is_newmark = 1;
	       } else if ( i == SIM_UHYDE && CNodes->Sub[j].Type == SIM_UHYDE ){
		    count = count + 1;
		    is_uhyde = 1;
	       } else if ( i == SIM_MEASURED && CNodes->Sub[j].Type == SIM_MEASURED ){
		    count = count + 1;		    
		    is_measured = 1;
	       } else if ( i == SIM_BOUCWEN && CNodes->Sub[j].Type == SIM_BOUCWEN ){
		    count = count + 1;
		    is_boucwen = 1;
	       } else if ( i == SIM_STONEDRUMS && CNodes->Sub[j].Type == SIM_STONEDRUMS ){
		    count = count + 1;
		    is_stonedrum = 1;
	       }
	  }

	  /* Allocate space */
	  Nodes = NULL; Nodes_Newmark = NULL; Nodes_Exp = NULL;
	  if( is_uhyde || is_exact ){
	       Nodes = (HDF5_Exact_UHYDE_t *) calloc( (size_t) count, sizeof( HDF5_Exact_UHYDE_t ) );
	  } else if ( is_newmark ){
	       Nodes_Newmark = (HDF5_Newmark_t *) calloc( (size_t) count, sizeof( HDF5_Newmark_t ) );
	  } else if ( is_adwin || is_measured ){
	       Nodes_Exp = (HDF5_Exp_Meas_t *) calloc( (size_t) count, sizeof( HDF5_Exp_Meas_t ) );
	  } else if ( is_boucwen ){
	       Nodes_BoucWen = (HDF5_BoucWen_t *) calloc( (size_t) count, sizeof( HDF5_BoucWen_t ) );
	  } else if ( is_stonedrum ){
	       Nodes_StoneDrum = (HDF5_StoneDrum_t *) calloc( (size_t) count, sizeof( HDF5_StoneDrum_t ) );
	  }

	  /* Copy the matching entries into the newly created structure */
	  k = 0;
	  for( j = 0; j < CNodes->Order; j++ ){
	       if ( i == EXP_ADWIN && CNodes->Sub[j].Type == EXP_ADWIN ){

		    Nodes_Exp[k].Position = CNodes->Array[j] - 1;  /* 0-based index */
		    Experimental = (ExpSub_t *) CNodes->Sub[j].SimStruct;
		    strcpy( Nodes_Exp[k].Description, Experimental->Description );
		    k = k + 1;
	       } else if ( i == SIM_EXACT_ESP && CNodes->Sub[j].Type == SIM_EXACT_ESP ){
		    Nodes[k].Position = CNodes->Array[j] - 1;      /* 0-based index */

		    ExSim = (ExactSimESP_t *) CNodes->Sub[j].SimStruct;
		    Nodes[k].InitValues[0] = ExSim->Mass;
		    Nodes[k].InitValues[1] = ExSim->Damp;
		    Nodes[k].InitValues[2] = ExSim->Stiff;

		    strcpy( Nodes[k].Description, ExSim->Description );	
		    k = k + 1;
	       } else if ( i == SIM_NEWMARK && CNodes->Sub[j].Type == SIM_NEWMARK ){
		    Nodes_Newmark[k].Position = CNodes->Array[j] - 1;      /* 0-based index */

		    Newmark = (NewmarkSim_t *) CNodes->Sub[j].SimStruct;
		    Nodes_Newmark[k].InitValues[0] = Newmark->Mass;
		    Nodes_Newmark[k].InitValues[1] = Newmark->Damp;
		    Nodes_Newmark[k].InitValues[2] = Newmark->Stiff;
		    Nodes_Newmark[k].InitValues[3] = InitCnt->TIntConst.Beta;
		    Nodes_Newmark[k].InitValues[4] = InitCnt->TIntConst.Gamma;
		    Nodes_Newmark[k].InitValues[5] = InitCnt->Delta_t/(HYSL_FLOAT) InitCnt->NSubstep;

		    strcpy( Nodes_Newmark[k].Description, Newmark->Description );	
		    k = k + 1;
	       } else if ( i == SIM_UHYDE && CNodes->Sub[j].Type == SIM_UHYDE ){
		    Nodes[k].Position = CNodes->Array[j] - 1;      /* 0-based index */

		    UHYDE = (UHYDEfbrSim_t *) CNodes->Sub[j].SimStruct;
		    Nodes[k].InitValues[0] = UHYDE->qyield;
		    Nodes[k].InitValues[1] = UHYDE->qyield/UHYDE->qplastic;
		    Nodes[k].InitValues[2] = UHYDE->qplastic*UHYDE->k;

		    strcpy( Nodes[k].Description, UHYDE->Description );
		    k = k + 1;
	       } else if ( i == SIM_MEASURED && CNodes->Sub[j].Type == SIM_MEASURED ){
		    Nodes_Exp[k].Position = CNodes->Array[j] - 1; /* 0-based index */
		    Experimental = (ExpSub_t *) CNodes->Sub[j].SimStruct;
		    strcpy( Nodes_Exp[k].Description, Experimental->Description );
		    k = k + 1;
	       } else if ( i == SIM_BOUCWEN && CNodes->Sub[j].Type == SIM_BOUCWEN ){
		    Nodes_BoucWen[k].Position = CNodes->Array[j] - 1;      /* 0-based index */

		    BoucWen = (BoucWen_t*) CNodes->Sub[j].SimStruct;
		    Nodes_BoucWen[k].InitValues[0] = BoucWen->alpha;
		    Nodes_BoucWen[k].InitValues[1] = BoucWen->ko;
		    Nodes_BoucWen[k].InitValues[2] = BoucWen->beta;
		    Nodes_BoucWen[k].InitValues[3] = BoucWen->gamma;
		    Nodes_BoucWen[k].InitValues[4] = BoucWen->n;
		    Nodes_BoucWen[k].InitValues[5] = BoucWen->A0;
		    Nodes_BoucWen[k].InitValues[6] = BoucWen->deltaA;
		    Nodes_BoucWen[k].InitValues[7] = BoucWen->nu0;
		    Nodes_BoucWen[k].InitValues[8] = BoucWen->deltaNu;
		    Nodes_BoucWen[k].InitValues[9] = BoucWen->eta0;
		    Nodes_BoucWen[k].InitValues[10] = BoucWen->deltaEta;
		    Nodes_BoucWen[k].InitValues[11] = BoucWen->vs0;
		    Nodes_BoucWen[k].InitValues[12] = BoucWen->p;
		    Nodes_BoucWen[k].InitValues[13] = BoucWen->q;
		    Nodes_BoucWen[k].InitValues[14] = BoucWen->psi0;
		    Nodes_BoucWen[k].InitValues[15] = BoucWen->deltaPsi;
		    Nodes_BoucWen[k].InitValues[16] = BoucWen->lambda;
		    strcpy( Nodes_BoucWen[k].Description, BoucWen->Description );

		    k = k + 1;
	       } else if ( i == SIM_STONEDRUMS && CNodes->Sub[j].Type == SIM_STONEDRUMS ){
		    Nodes_StoneDrum[k].Position = CNodes->Array[j] - 1;      /* 0-based index */

		    StoneDrum = (StoneDrums_t*) CNodes->Sub[j].SimStruct;
		    strcpy( Nodes_StoneDrum[k].Description, StoneDrum->Description );

		    k = k + 1;
	       }
	  }

	  if ( (i == SIM_EXACT_ESP) && is_exact ){
	       memtype = H5Tcreate (H5T_COMPOUND, sizeof( HDF5_Exact_UHYDE_t) );
	       H5Tinsert( memtype, "Eq. column (0-based index)", HOFFSET( HDF5_Exact_UHYDE_t, Position), H5T_NATIVE_INT );
	       H5Tinsert( memtype, "Mass [kg]", HOFFSET( HDF5_Exact_UHYDE_t, InitValues[0]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "Damping", HOFFSET( HDF5_Exact_UHYDE_t, InitValues[1]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "Stiffness [N/m]", HOFFSET( HDF5_Exact_UHYDE_t, InitValues[2]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "Description", HOFFSET( HDF5_Exact_UHYDE_t, Description), strtype );

	       dims[0] = (hsize_t) count; 
	       space = H5Screate_simple (1, dims, NULL);
	       /* Create the dataset */
	       dset = H5Dcreate ( group_id, "/Test Parameters/Substructures/Exact", memtype, space, H5P_DEFAULT,
				  H5P_DEFAULT, H5P_DEFAULT);
	       status = H5Dwrite ( dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, Nodes );
	       if( status < 0 ){
		    Print_Header( ERROR );
		    fprintf( stderr, "Error writing to HDF5 file.\n" );
		    exit( EXIT_FAILURE );
	       }
	       
	       free( Nodes );

	       status = H5Dclose( dset );
	       if( status < 0 ){
		    Print_Header( ERROR );
		    fprintf( stderr, "Error closing dataset in HDF5 file.\n" );
		    exit( EXIT_FAILURE );
	       }
	       status = H5Sclose( space );
	       if( status < 0 ){
		    Print_Header( ERROR );
		    fprintf( stderr, "Error closing dataspace in HDF5 file.\n" );
		    exit( EXIT_FAILURE );
	       }
	       status = H5Tclose( memtype );
	       if( status < 0 ){
		    Print_Header( ERROR );
		    fprintf( stderr, "Error closing memtype in HDF5 file.\n" );
		    exit( EXIT_FAILURE );
	       }
	  } else if ( (i == SIM_NEWMARK) && is_newmark ){
	       memtype = H5Tcreate (H5T_COMPOUND, sizeof( HDF5_Newmark_t) );
	       H5Tinsert( memtype, "Eq. column (0-based index)", HOFFSET( HDF5_Newmark_t, Position), H5T_NATIVE_INT );
	       H5Tinsert( memtype, "Mass [kg]", HOFFSET( HDF5_Newmark_t, InitValues[0]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "Damping", HOFFSET( HDF5_Newmark_t, InitValues[1]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "Stiffness [N/m]", HOFFSET( HDF5_Newmark_t, InitValues[2]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "Beta", HOFFSET( HDF5_Newmark_t, InitValues[3]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "Gamma", HOFFSET( HDF5_Newmark_t, InitValues[4]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "Delta T [s]", HOFFSET( HDF5_Newmark_t, InitValues[5]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "Description", HOFFSET( HDF5_Newmark_t, Description), strtype );

	       dims[0] = (hsize_t) count; 
	       space = H5Screate_simple (1, dims, NULL);
	       /* Create the dataset */
	       dset = H5Dcreate ( group_id, "/Test Parameters/Substructures/Newmark", memtype, space, H5P_DEFAULT,
				  H5P_DEFAULT, H5P_DEFAULT);
	       status = H5Dwrite ( dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, Nodes_Newmark );
	       if( status < 0 ){
		    Print_Header( ERROR );
		    fprintf( stderr, "Error writing to HDF5 file.\n" );
		    exit( EXIT_FAILURE );
	       }
	       
	       free( Nodes_Newmark );

	       status = H5Dclose( dset );
	       if( status < 0 ){
		    Print_Header( ERROR );
		    fprintf( stderr, "Error closing dataset in HDF5 file.\n" );
		    exit( EXIT_FAILURE );
	       }
	       status = H5Sclose( space );
	       if( status < 0 ){
		    Print_Header( ERROR );
		    fprintf( stderr, "Error closing dataspace in HDF5 file.\n" );
		    exit( EXIT_FAILURE );
	       }
	       status = H5Tclose( memtype );
	       if( status < 0 ){
		    Print_Header( ERROR );
		    fprintf( stderr, "Error closing memtype in HDF5 file.\n" );
		    exit( EXIT_FAILURE );
	       }
	  } else if ( (i == SIM_UHYDE) && is_uhyde ){
	       memtype = H5Tcreate (H5T_COMPOUND, sizeof( HDF5_Exact_UHYDE_t) );
	       H5Tinsert( memtype, "Eq. column (0-based index)", HOFFSET( HDF5_Exact_UHYDE_t, Position), H5T_NATIVE_INT );
	       H5Tinsert( memtype, "qyield", HOFFSET( HDF5_Exact_UHYDE_t, InitValues[0]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "yield factor", HOFFSET( HDF5_Exact_UHYDE_t, InitValues[1]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "Friction force [N]", HOFFSET( HDF5_Exact_UHYDE_t, InitValues[2]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "Description", HOFFSET( HDF5_Exact_UHYDE_t, Description), strtype );

	       dims[0] = (hsize_t) count; 
	       space = H5Screate_simple (1, dims, NULL);
	       /* Create the dataset */
	       dset = H5Dcreate ( file_id, "/Test Parameters/Substructures/UHYDE-fbr", memtype, space, H5P_DEFAULT,
				  H5P_DEFAULT, H5P_DEFAULT);
	       status = H5Dwrite ( dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, Nodes );
	       
	       free( Nodes );

	       status = H5Dclose( dset );
	       status = H5Sclose( space );
	       status = H5Tclose( memtype );

	  } else if ( (i == SIM_MEASURED) && is_measured ){
	       memtype = H5Tcreate (H5T_COMPOUND, sizeof(HDF5_Exp_Meas_t) );
	       H5Tinsert( memtype, "Eq. column (0-based index)", HOFFSET( HDF5_Exp_Meas_t, Position), H5T_NATIVE_INT );
	       H5Tinsert( memtype, "Description", HOFFSET( HDF5_Exp_Meas_t, Description), strtype );

	       dims[0] = (hsize_t) count; 
	       space = H5Screate_simple (1, dims, NULL);
	       /* Create the dataset */

	       dset = H5Dcreate ( group_id, "/Test Parameters/Substructures/Measured", memtype, space, H5P_DEFAULT,
				  H5P_DEFAULT, H5P_DEFAULT);
	       status = H5Dwrite ( dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, Nodes_Exp );
	       
	       free( Nodes_Exp );

	       status = H5Dclose( dset );
	       status = H5Sclose( space );
	       status = H5Tclose( memtype );
	  } else if ((i == EXP_ADWIN) && is_adwin ){

	       memtype = H5Tcreate (H5T_COMPOUND, sizeof(HDF5_Exp_Meas_t) );
	       H5Tinsert( memtype, "Eq. column (0-based index)", HOFFSET( HDF5_Exp_Meas_t, Position), H5T_NATIVE_INT );
	       H5Tinsert( memtype, "Description", HOFFSET( HDF5_Exp_Meas_t, Description), strtype );

	       dims[0] = (hsize_t) count; 
	       space = H5Screate_simple (1, dims, NULL);
	       /* Create the dataset */

	       dset = H5Dcreate ( group_id, "/Test Parameters/Substructures/Experimental", memtype, space, H5P_DEFAULT,
				  H5P_DEFAULT, H5P_DEFAULT);
	       status = H5Dwrite ( dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, Nodes_Exp );
	       
	       free( Nodes_Exp );

	       status = H5Dclose( dset );
	       status = H5Sclose( space );
	       status = H5Tclose( memtype );

	  } else if ((i == SIM_BOUCWEN) && is_boucwen ){
	       memtype = H5Tcreate (H5T_COMPOUND, sizeof( HDF5_BoucWen_t) );
	       H5Tinsert( memtype, "Eq. column (0-based index)", HOFFSET( HDF5_BoucWen_t, Position), H5T_NATIVE_INT );
	       H5Tinsert( memtype, "Alpha", HOFFSET( HDF5_BoucWen_t, InitValues[0]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "Elastic stiffness ko [N/m]", HOFFSET( HDF5_BoucWen_t, InitValues[1]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "Beta", HOFFSET( HDF5_BoucWen_t, InitValues[2]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "Gamma", HOFFSET( HDF5_BoucWen_t, InitValues[3]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "n", HOFFSET( HDF5_BoucWen_t, InitValues[4]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "A0", HOFFSET( HDF5_BoucWen_t, InitValues[5]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "Delta A", HOFFSET( HDF5_BoucWen_t, InitValues[6]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "nu0", HOFFSET( HDF5_BoucWen_t, InitValues[7]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "Delta nu", HOFFSET( HDF5_BoucWen_t, InitValues[8]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "Eta0", HOFFSET( HDF5_BoucWen_t, InitValues[9]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "Delta Eta", HOFFSET( HDF5_BoucWen_t, InitValues[10]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "vs0", HOFFSET( HDF5_BoucWen_t, InitValues[11]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "p", HOFFSET( HDF5_BoucWen_t, InitValues[12]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "q", HOFFSET( HDF5_BoucWen_t, InitValues[13]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "Psi0", HOFFSET( HDF5_BoucWen_t, InitValues[14]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "DeltaPsi", HOFFSET( HDF5_BoucWen_t, InitValues[15]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "Lambda", HOFFSET( HDF5_BoucWen_t, InitValues[16]), H5T_NATIVE_HYSL_FLOAT );

	       H5Tinsert( memtype, "Description", HOFFSET( HDF5_BoucWen_t, Description), strtype );

	       dims[0] = (hsize_t) count; 
	       space = H5Screate_simple (1, dims, NULL);
	       /* Create the dataset */
	       dset = H5Dcreate ( group_id, "/Test Parameters/Substructures/Bouc Wen", memtype, space, H5P_DEFAULT,
				  H5P_DEFAULT, H5P_DEFAULT);
	       status = H5Dwrite ( dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, Nodes_BoucWen );
	       if( status < 0 ){
		    Print_Header( ERROR );
		    fprintf( stderr, "Error writing to HDF5 file.\n" );
		    exit( EXIT_FAILURE );
	       }
	       
	       free( Nodes_BoucWen );

	       status = H5Dclose( dset );
	       if( status < 0 ){
		    Print_Header( ERROR );
		    fprintf( stderr, "Error closing dataset in HDF5 file.\n" );
		    exit( EXIT_FAILURE );
	       }
	       status = H5Sclose( space );
	       if( status < 0 ){
		    Print_Header( ERROR );
		    fprintf( stderr, "Error closing dataspace in HDF5 file.\n" );
		    exit( EXIT_FAILURE );
	       }
	       status = H5Tclose( memtype );
	       if( status < 0 ){
		    Print_Header( ERROR );
		    fprintf( stderr, "Error closing memtype in HDF5 file.\n" );
		    exit( EXIT_FAILURE );
	       }
	  } else if ((i == SIM_STONEDRUMS) && is_stonedrum ){

	       memtype = H5Tcreate (H5T_COMPOUND, sizeof( HDF5_StoneDrum_t) );
	       H5Tinsert( memtype, "Eq. column (0-based index)", HOFFSET( HDF5_StoneDrum_t, Position), H5T_NATIVE_INT );
	       H5Tinsert( memtype, "Alpha", HOFFSET( HDF5_StoneDrum_t, InitValues[0]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "Elastic stiffness ko [N/m]", HOFFSET( HDF5_StoneDrum_t, InitValues[1]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "Beta", HOFFSET( HDF5_StoneDrum_t, InitValues[2]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "Gamma", HOFFSET( HDF5_StoneDrum_t, InitValues[3]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "n", HOFFSET( HDF5_StoneDrum_t, InitValues[4]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "A0", HOFFSET( HDF5_StoneDrum_t, InitValues[5]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "Delta A", HOFFSET( HDF5_StoneDrum_t, InitValues[6]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "nu0", HOFFSET( HDF5_StoneDrum_t, InitValues[7]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "Delta nu", HOFFSET( HDF5_StoneDrum_t, InitValues[8]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "Eta0", HOFFSET( HDF5_StoneDrum_t, InitValues[9]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "Delta Eta", HOFFSET( HDF5_StoneDrum_t, InitValues[10]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "vs0", HOFFSET( HDF5_StoneDrum_t, InitValues[11]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "p", HOFFSET( HDF5_StoneDrum_t, InitValues[12]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "q", HOFFSET( HDF5_StoneDrum_t, InitValues[13]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "Psi0", HOFFSET( HDF5_StoneDrum_t, InitValues[14]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "DeltaPsi", HOFFSET( HDF5_StoneDrum_t, InitValues[15]), H5T_NATIVE_HYSL_FLOAT );
	       H5Tinsert( memtype, "Lambda", HOFFSET( HDF5_StoneDrum_t, InitValues[16]), H5T_NATIVE_HYSL_FLOAT );

	       H5Tinsert( memtype, "Description", HOFFSET( HDF5_StoneDrum_t, Description), strtype );

	       dims[0] = (hsize_t) count; 
	       space = H5Screate_simple (1, dims, NULL);
	       /* Create the dataset */
	       dset = H5Dcreate ( group_id, "/Test Parameters/Substructures/Stone Drum", memtype, space, H5P_DEFAULT,
				  H5P_DEFAULT, H5P_DEFAULT);
	       status = H5Dwrite ( dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, Nodes_StoneDrum );
	       if( status < 0 ){
		    Print_Header( ERROR );
		    fprintf( stderr, "Error writing to HDF5 file.\n" );
		    exit( EXIT_FAILURE );
	       }
	       
	       free( Nodes_StoneDrum );

	       status = H5Dclose( dset );
	       if( status < 0 ){
		    Print_Header( ERROR );
		    fprintf( stderr, "Error closing dataset in HDF5 file.\n" );
		    exit( EXIT_FAILURE );
	       }
	       status = H5Sclose( space );
	       if( status < 0 ){
		    Print_Header( ERROR );
		    fprintf( stderr, "Error closing dataspace in HDF5 file.\n" );
		    exit( EXIT_FAILURE );
	       }
	       status = H5Tclose( memtype );
	       if( status < 0 ){
		    Print_Header( ERROR );
		    fprintf( stderr, "Error closing memtype in HDF5 file.\n" );
		    exit( EXIT_FAILURE );
	       }
	  }

	  
	  /* Close and release resources. */
     }    

     status = H5Tclose( strtype );
     status = H5Gclose( group_id );
}
     

void HDF5_CreateGroup_TimeIntegration( const hid_t hdf5_file, const AlgConst_t *const InitCnt )
{

     hid_t    group_id;
     herr_t   status;


     /* Create the Time integration group */
     group_id = H5Gcreate( hdf5_file, "/Time Integration", H5P_DEFAULT,
			  H5P_DEFAULT, H5P_DEFAULT);
     /* Acceleration */
     HDF5_Create_Dataset( hdf5_file, "/Time Integration/Acceleration", (int) InitCnt->NStep,
			  (int) InitCnt->Order );
     status = H5LTset_attribute_string( hdf5_file, "Time Integration/Acceleration", "Units", "m/s^2" );
     /* Velocity */
     HDF5_Create_Dataset( hdf5_file, "/Time Integration/Velocity", (int) InitCnt->NStep,
			  (int)InitCnt->Order );
     status = H5LTset_attribute_string( hdf5_file, "Time Integration/Velocity",	"Units", "m/s" );
     /* Displacement */
     HDF5_Create_Dataset( hdf5_file, "/Time Integration/Displacement", (int) InitCnt->NStep, 
			  (int) InitCnt->Order );
     status = H5LTset_attribute_string( hdf5_file, "Time Integration/Displacement", "Units", "m" );

     /* If the number of substructures is greater than zero, store the results */
     if( InitCnt->OrderSub > 0 ){
	  HDF5_Create_Dataset( hdf5_file, "/Time Integration/Coupling force", (int) InitCnt->NStep,
			       (int) InitCnt->Order );
	  status = H5LTset_attribute_string( hdf5_file, "Time Integration/Coupling force", "Units", "N" );

	  HDF5_Create_Dataset( hdf5_file, "/Time Integration/Bouc-Wen Forces", (int) InitCnt->NStep, (int) InitCnt->Order );
	  status = H5LTset_attribute_string( hdf5_file, "/Time Integration/Bouc-Wen Forces", "Units", "N" );
     }

     /* If the PID is used, then store the values */
     if( InitCnt->PID.P != 0.0 || InitCnt->PID.I != 0.0 || InitCnt->PID.D != 0.0 ){
	  HDF5_Create_Dataset( hdf5_file, "/Time Integration/Error force", (int) InitCnt->NStep,
			       (int) InitCnt->Order );
	  status = H5LTset_attribute_string( hdf5_file, "Time Integration/Error force", "Units", "N" );
     }

     status = H5Gclose( group_id );

}

void HDF5_Store_TimeHistoryData( const hid_t hdf5_file, const MatrixVector_t *const Acc,
				 const MatrixVector_t *const Vel, const MatrixVector_t *const Disp,
				 const MatrixVector_t *const fc, const MatrixVector_t *const fu,
				 const int istep, const AlgConst_t *const InitCnt )
{
     /* Create the Time integration group */

     HDF5_AddResults_to_Dataset( hdf5_file, "/Time Integration/Acceleration", Acc, istep );
     HDF5_AddResults_to_Dataset( hdf5_file, "/Time Integration/Velocity", Vel, istep );
     HDF5_AddResults_to_Dataset( hdf5_file, "/Time Integration/Displacement", Disp, istep );

     /* If the number of substructures is greater than zero, store the results */
     if( InitCnt->OrderSub > 0 ){
	  HDF5_AddResults_to_Dataset( hdf5_file, "/Time Integration/Coupling force", fc, istep );
     }

     /* If the PID is used, then store the values */
     if( InitCnt->PID.P != 0.0 || InitCnt->PID.I != 0.0 || InitCnt->PID.D != 0.0 ){
	  HDF5_AddResults_to_Dataset( hdf5_file, "/Time Integration/Error force", fu, istep );
     }
}

void HDF5_Store_TMD( const hid_t hdf5_file, const double *const Acc, const double *const Vel, const double *const Disp, const int istep )
{
     HYSL_FLOAT TMD[3];

     TMD[0] = *Acc; TMD[1] = *Vel; TMD[2] = *Disp;
     
     hid_t   dataset_id, filespace_id, memspace_id;
     herr_t  status;
     hsize_t size[2], offset[2], dims[2];

     size[0] = (hsize_t) istep; /* Num_Steps starts at 1 */
     size[1] = (hsize_t) 3;

     /* Create the data space for the dataset. */
     dataset_id = H5Dopen( hdf5_file, "/Time Integration/TMD", H5P_DEFAULT );
     H5Dextend( dataset_id, size );

     filespace_id = H5Dget_space( dataset_id );
     offset[0] = ((hsize_t) istep - 1);
     offset[1] = 0;
     dims[0] = 1; dims[1] = (hsize_t) 3;
     status = H5Sselect_hyperslab( filespace_id, H5S_SELECT_SET, offset, NULL, dims, NULL );
     
     memspace_id = H5Screate_simple( 2, dims, NULL );
     
     status = H5Dwrite( dataset_id, H5T_NATIVE_HYSL_FLOAT, memspace_id, filespace_id, H5P_DEFAULT, TMD );

     status = H5Dclose( dataset_id );
     status = H5Sclose( memspace_id );
     status = H5Sclose( filespace_id );
}

void HDF5_Store_BoucWen( const hid_t hdf5_file, const MatrixVector_t *const HistDisp, const MatrixVector_t *const HistLoop, const int istep)
{
     HDF5_AddResults_to_Dataset( hdf5_file, "/Time Integration/Bouc-Wen Displacements", HistDisp, istep );
     HDF5_AddResults_to_Dataset( hdf5_file, "/Time Integration/Bouc-Wen Forces", HistLoop, istep );
}

void HDF5_Store_Time( const hid_t hdf5_file, const SaveTime_t *const Time )
{

     hid_t    memtype, space, dset, strtype;
     hsize_t  dims[1] = {1};
     herr_t   status;


     strtype = H5Tcopy( H5T_C_S1 );
     status = H5Tset_size( strtype, H5T_VARIABLE );

     memtype = H5Tcreate (H5T_COMPOUND, sizeof( SaveTime_t ) );
     H5Tinsert( memtype, "Date", HOFFSET( SaveTime_t, Date_time), strtype );
     H5Tinsert( memtype, "Duration [ms]", HOFFSET( SaveTime_t, Elapsed_time), H5T_NATIVE_HYSL_FLOAT );
     
     space = H5Screate_simple (1, dims, NULL);
     /* Create the dataset */
     dset = H5Dcreate ( hdf5_file, "/Test Parameters/Date and Duration", memtype, space, H5P_DEFAULT,
			H5P_DEFAULT, H5P_DEFAULT);
     status = H5Dwrite ( dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, Time );
    
     status = H5Dclose( dset );
     status = H5Sclose( space );
     status = H5Tclose( memtype );
     status = H5Tclose( strtype );
}

#if _ADWIN_
void HDF5_StoreADwinData( const hid_t hdf5_file, const HYSL_FLOAT *Array, char **Entry_Names, const int Length )
{
     hid_t   group_id;
     herr_t  status;

     group_id = H5Gcreate( hdf5_file, "/ADwin measurements", H5P_DEFAULT,
			  H5P_DEFAULT, H5P_DEFAULT);     
     status = H5Gclose(group_id);

     HDF5_AddFloatArray_AsTable( hdf5_file, "/ADwin measurements/Recorded data", Entry_Names,
				  Array, NUM_CHANNELS, Length );
}

void ADwin_SaveData_HDF5( const hid_t hdf5_file, const int Num_Steps, const int Num_Sub,
			  const int Num_Channels, char **Chan_Names, const int DataIndex )
{
     int Length;
     HYSL_FLOAT *Data = NULL;

     Length = Num_Sub*Num_Steps*Num_Channels;
     Data = (HYSL_FLOAT *) calloc( (size_t) Length, sizeof( HYSL_FLOAT ) );
     if( Data == NULL ){
	  Print_Header( WARNING );
	  fprintf( stderr, "ADwin_SaveData_HDF5: Out of memory. Manual extraction of the data required.\n" );
     }

     /* Get the data from ADwin */
#if _FLOAT_
     GetData_Float( (int32_t) DataIndex, Data, 1, (int32_t) Length);
#else
     GetData_Double( (int32_t) DataIndex, Data, 1, (int32_t) Length);
#endif
  
     /* Save the data into an HDF5 file */
     HDF5_StoreADwinData( hdf5_file, Data, Chan_Names, Num_Sub*Num_Steps );

     /* Free allocated memory */
     free( Data );
}
#endif

void HDF5_AddFloatArray_AsTable( const hid_t file_id, const char *Name_path, char **Names, const HYSL_FLOAT *Array, const int Num_param, const int Length )
{
     int      i;
     hid_t    memtype, filetype, space, dset;
     hsize_t  dims[1] = { (hsize_t) Length};
     herr_t   status;


     /* Create the compound datatype for memory. */
     memtype = H5Tcreate( H5T_COMPOUND, (size_t) Num_param*H5Tget_size(H5T_NATIVE_HYSL_FLOAT) );
     for( i = 0; i < Num_param; i++ ){
	  status = H5Tinsert( memtype, Names[i], (size_t) i*H5Tget_size(H5T_NATIVE_HYSL_FLOAT),
			      H5T_NATIVE_HYSL_FLOAT );
     }
    
     /* Create the compound datatype for the file.*/
     filetype = H5Tcreate( H5T_COMPOUND, (size_t) Num_param*H5Tget_size(H5T_NATIVE_HYSL_FLOAT) );
     for( i = 0; i < Num_param; i++ ){
	  status = H5Tinsert( filetype, Names[i], (size_t) i*H5Tget_size(H5T_NATIVE_HYSL_FLOAT),
			      H5T_NATIVE_HYSL_FLOAT );
     }
 
     /* Create dataspace.  Setting maximum size to NULL sets the maximum size to be the current size. */    
     space = H5Screate_simple (1, dims, NULL);

     /* Create the dataset and write the compound data to it. */
     dset = H5Dcreate( file_id, Name_path, filetype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
     status = H5Dwrite( dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, Array );

     /* Close and release resources. */
     status = H5Tclose( filetype );
     status = H5Tclose( memtype);
     status = H5Dclose( dset );
     status = H5Sclose( space );
}

void HDF5_AddIntArray_AsTable( const hid_t file_id, const char *Name_path, char **Names, const int *Array,
			       const int Num_param )
{
     int      i;
     hid_t    memtype, filetype, space, dset;
     hsize_t  dims[1] = {1};
     herr_t   status;


     /* Create the compound datatype for memory. */
     memtype = H5Tcreate( H5T_COMPOUND, (size_t) Num_param*H5Tget_size(H5T_NATIVE_INT) );
     for( i = 0; i < Num_param; i++ ){
	  status = H5Tinsert( memtype, Names[i], (size_t) i*H5Tget_size(H5T_NATIVE_INT),
			      H5T_NATIVE_INT );
     }

    
     /* Create the compound datatype for the file.*/
     filetype = H5Tcreate( H5T_COMPOUND, (size_t) Num_param*H5Tget_size(H5T_NATIVE_INT) );
     for( i = 0; i < Num_param; i++ ){
	  status = H5Tinsert( filetype, Names[i], (size_t) i*H5Tget_size(H5T_NATIVE_INT),
			      H5T_NATIVE_INT);
     }
 
     /* Create dataspace. Setting maximum size to NULL sets the maximum size to be the current size. */    
     space = H5Screate_simple (1, dims, NULL);

     /* Create the dataset and write the compound data to it. */
     dset = H5Dcreate( file_id, Name_path, filetype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
     status = H5Dwrite( dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, Array );

     /* Close and release resources. */
     status = H5Dclose (dset);
     status = H5Sclose (space);
}

void HDF5_Create_Dataset( const hid_t file_id, const char *Path_name, const int NStep, const int Order )
{
     hsize_t dims[2] = {1, (hsize_t) Order };
     hsize_t max_dims[2] = {(hsize_t) NStep, (hsize_t) Order };

     hid_t cparms, dataspace_id, dataset_id;
     herr_t status;

     dataspace_id = H5Screate_simple( 2, dims, max_dims );

     cparms = H5Pcreate( H5P_DATASET_CREATE);
     status = H5Pset_chunk( cparms, 2, dims );

     dataset_id = H5Dcreate( file_id, Path_name, H5T_NATIVE_HYSL_FLOAT, dataspace_id, H5P_DEFAULT, cparms, H5P_DEFAULT);

     status = H5Pclose( cparms );
     status = H5Dclose( dataset_id );
     status = H5Sclose( dataspace_id );
}

void HDF5_AddResults_to_Dataset( const hid_t file_id, const char *Path_name, const MatrixVector_t *const Data,
				 const int Step_count )
{
     hid_t   dataset_id, filespace_id, memspace_id;
     herr_t  status;
     hsize_t size[2], offset[2], dims[2];

     size[0] = (hsize_t) Step_count; /* Num_Steps starts at 1 */
     size[1] = (hsize_t) Data->Rows;

     /* Create the data space for the dataset. */
     dataset_id = H5Dopen( file_id, Path_name, H5P_DEFAULT );
     H5Dextend( dataset_id, size );

     filespace_id = H5Dget_space( dataset_id );
     offset[0] = ((hsize_t) Step_count - 1);
     offset[1] = 0;
     dims[0] = 1; dims[1] = (hsize_t) Data->Rows;
     status = H5Sselect_hyperslab( filespace_id, H5S_SELECT_SET, offset, NULL, dims, NULL );
     
     memspace_id = H5Screate_simple( 2, dims, NULL );
     
     status = H5Dwrite( dataset_id, H5T_NATIVE_HYSL_FLOAT, memspace_id, filespace_id, H5P_DEFAULT, Data->Array );

     status = H5Dclose( dataset_id );
     status = H5Sclose( memspace_id );
     status = H5Sclose( filespace_id );
}

void HDF5_CloseFile( hid_t *const hdf5_file )
{    
     herr_t status;

     /* Terminate access to the file. */
     status = H5Fclose( *hdf5_file );
}
