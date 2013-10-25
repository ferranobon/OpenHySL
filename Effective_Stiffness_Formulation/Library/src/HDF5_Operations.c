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
#include "Substructure_UHYDEfbr.h"
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

int HDF5_CreateFile( const char *Filename )
{

     hid_t    file_id;

     /* Create a new file using default properties. */
     file_id = H5Fcreate( Filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
     
     return (int) file_id;
}

void HDF5_CreateGroup_Parameters( const int hdf5_file, const AlgConst_t *const InitCnt,
				  const CouplingNode_t *const CNodes, const HYSL_FLOAT *const Acc,
				  const HYSL_FLOAT *const Vel, const HYSL_FLOAT *const Disp )
{
     hid_t    group_id, dataset_id, dataspace_id;
     hsize_t  dims[2];
     herr_t   status;
     HYSL_FLOAT   *dArray;
     int      *iArray, i;
     char     **Entry_Names;
     
     Entry_Names = (char **) malloc( (size_t) 6*sizeof(char*) );
     dArray = (HYSL_FLOAT *) malloc( (size_t) 3*sizeof(HYSL_FLOAT) );

     /* Create a group named "/Time History" in the file. */
     group_id = H5Gcreate( hdf5_file, "/Test Parameters", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );     
     status = H5Gclose( group_id );
     if( status < 0 ){
	  Print_Header( ERROR );
	  fprintf( stderr, "Error creating group Test Parameters in HDF5 file.\n" );
	  exit( EXIT_FAILURE );
     }
     
     /* Add input load */
     group_id = H5Gcreate( hdf5_file, "/Test Parameters/Input", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
     dims[0] = InitCnt->NStep;
     dims[1] = 1;
     dataspace_id = H5Screate_simple( 2, dims, NULL );
     if( Acc != NULL ){
	  /* assssssssssfa asfjajs asjfnakjsdfajksh fafa safsjdfaskjh asjdf ajsf aksjdf aksj aksjdfnak sjdfn a
	   *  */
	  dataset_id = H5Dcreate ( group_id, "/Test Parameters/Input/Acceleration", H5T_NATIVE_HYSL_FLOAT,
				   dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	  status = H5LTset_attribute_string( hdf5_file, "Test Parameters/Input/Acceleration", "Units", "m/s^2" );
	  status = H5Dwrite( dataset_id, H5T_NATIVE_HYSL_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, Acc );
	  status = H5Dclose( dataset_id );
     }

     if( Vel != NULL ){
	  dataset_id = H5Dcreate ( group_id, "/Test Parameters/Input/Velocity", H5T_NATIVE_HYSL_FLOAT,
				   dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	  status = H5LTset_attribute_string( hdf5_file, "Test Parameters/Input/Velocity", "Units", "m/s" );
	  status = H5Dwrite( dataset_id, H5T_NATIVE_HYSL_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, Vel );
	  status = H5Dclose( dataset_id );
     }

     if( Disp != NULL ){
	  dataset_id = H5Dcreate ( group_id, "/Test Parameters/Input/Displacement", H5T_NATIVE_HYSL_FLOAT,
				   dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	  status = H5LTset_attribute_string( hdf5_file, "Test Parameters/Input/Displacement", "Units", "m" );
	  status = H5Dwrite( dataset_id, H5T_NATIVE_HYSL_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, Disp );
	  status = H5Dclose( dataset_id );
     }
     status = H5Gclose( group_id );

     /* Add Integration parameters */
     Entry_Names[0] = strdup( "Beta" );      dArray[0] = InitCnt->Newmark.Beta;
     Entry_Names[1] = strdup( "Gamma" );     dArray[1] = InitCnt->Newmark.Gamma;
     Entry_Names[2] = strdup( "Time step" ); dArray[2] = InitCnt->Delta_t;
     
     HDF5_AddFloatArray_AsTable( hdf5_file, "/Test Parameters/Time Integration", Entry_Names, dArray, 3, 1 );
     free( Entry_Names[0] );
     free( Entry_Names[1] );
     free( Entry_Names[2] );

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
	  Save_InformationCNodes( hdf5_file, "/Test Parameters/Substructures", CNodes );
     }

     for( i = 0; i < 4; i++ ){
	  free( Entry_Names[i] );
     }

     free( Entry_Names );
}


void Save_InformationCNodes( const hid_t file_id, const char *Name_path, const CouplingNode_t *const CNodes )
{
     int i, j, k, count, is_adwin, is_exact, is_uhyde, is_measured;
     hid_t    group_id, strtype, memtype, space, dset;
     hsize_t  dims[1];
     herr_t   status;
     HDF5_Exact_UHYDE_t *Nodes;
     HDF5_Exp_Meas_t *Nodes_Exp;
     ExactSim_t *TMD;
     ExactSimESP_t *ExSim;
     UHYDEfbrSim_t *UHYDE;
     ExpSub_t *Experimental;

   
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
	  is_adwin = 0; is_exact = 0; is_uhyde=0; is_measured = 0;

	  for( j = 0; j < CNodes->Order; j++ ){
	       if ( i == EXP_ADWIN && CNodes->Sub[j].Type == EXP_ADWIN ){
		    count = count + 1;
		    is_adwin = 1;	
	       } else if ( i == SIM_EXACT_ESP && CNodes->Sub[j].Type == SIM_EXACT_ESP ){
		    count = count + 1;
		    is_exact = 1;
	       } else if ( i == SIM_UHYDE && CNodes->Sub[j].Type == SIM_UHYDE ){
		    count = count + 1;
		    is_uhyde = 1;
	       } else if ( i == SIM_MEASURED && CNodes->Sub[j].Type == SIM_MEASURED ){
		    count = count + 1;		    
		    is_measured = 1;
	       }
	  }

	  /* Allocate space */
	  Nodes = NULL;
	  if( is_uhyde || is_exact ){
	       Nodes = (HDF5_Exact_UHYDE_t *) calloc( (size_t) count, sizeof( HDF5_Exact_UHYDE_t ) );
	  } else {
	       Nodes_Exp = (HDF5_Exp_Meas_t *) calloc( (size_t) count, sizeof( HDF5_Exp_Meas_t ) );
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
	       
	       for( j = 0; j < count; j++ ){
		    free( Nodes[j].Description );	     
	       }
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
	       
	       for( j = 0; j < count; j++ ){
		    free( Nodes_Exp[j].Description );	     
	       }
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
	       
	       for( j = 0; j < count; j++ ){
		    free( Nodes_Exp[j].Description );	     
	       }
	       free( Nodes_Exp );

	       status = H5Dclose( dset );
	       status = H5Sclose( space );
	       status = H5Tclose( memtype );

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
void HDF5_StoreADwinData( const int hdf5_file, const HYSL_FLOAT *Array, char **Entry_Names, const int Length )
{
     hid_t   group_id;
     herr_t  status;

     group_id = H5Gcreate( hdf5_file, "/ADwin measurements", H5P_DEFAULT,
			  H5P_DEFAULT, H5P_DEFAULT);     
     status = H5Gclose(group_id);

     HDF5_AddFloatArray_AsTable( hdf5_file, "/ADwin measurements/Recorded data", Entry_Names,
				  Array, 24, Length );
}

void ADwin_SaveData_HDF5( const int hdf5_file, const int Num_Steps, const int Num_Sub,
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
