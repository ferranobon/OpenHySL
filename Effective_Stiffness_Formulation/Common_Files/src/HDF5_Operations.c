#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Initiation.h"
#include "Substructure.h"
#include "MatrixVector.h"
#include "HDF5_Operations.h"
#include "RoutinesADwin.h"
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

void HDF5_CreateGroup_Parameters( int hdf5_file, AlgConst *const InitCnt, Coupling_Node *const CNodes )
{
     hid_t    file_id, group_id;
     herr_t   status;
     double   *dArray;
     int      *iArray, i;
     char **Entry_Names;
     
     file_id = (hid_t) hdf5_file;
     Entry_Names = (char **) malloc( (size_t) 6*sizeof(char*) );
     dArray = (double *) malloc( (size_t) 3*sizeof(double) );

     /* Create a group named "/Time History" in the file. */
     group_id = H5Gcreate(file_id, "/Test Parameters", H5P_DEFAULT,
			  H5P_DEFAULT, H5P_DEFAULT);     
     status = H5Gclose(group_id);

     /* Add Integration parameters */
     Entry_Names[0] = strdup( "Beta" ); Entry_Names[1] = strdup( "Gamma" );
     Entry_Names[2] = strdup( "Time step" );
     dArray[0] = InitCnt->Newmark.Beta; dArray[1] = InitCnt->Newmark.Gamma;
     dArray[2] = InitCnt->Delta_t;
     HDF5_AddDoubleArray_AsTable( file_id, "/Test Parameters/Time Integration", Entry_Names,
				    dArray, 3 );
     free( Entry_Names[0] ); free( Entry_Names[1] ); free( Entry_Names[2] );

     /* Add Rayleigh damping values */
     Entry_Names[0] = strdup( "Alpha" ); Entry_Names[1] = strdup( "Beta" );
     dArray[0] = InitCnt->Rayleigh.Alpha; dArray[1] = InitCnt->Rayleigh.Beta;
     HDF5_AddDoubleArray_AsTable( file_id, "/Test Parameters/Rayleigh", Entry_Names,
				    dArray, 2 );
     free( Entry_Names[0] ); free( Entry_Names[1] );

     /* Add PID parameters */
     Entry_Names[0] = strdup( "P" ); Entry_Names[1] = strdup( "I" );
     Entry_Names[2] = strdup( "D" );
     dArray[0] = InitCnt->PID.P; dArray[1] = InitCnt->PID.I;
     dArray[2] = InitCnt->PID.D;
     HDF5_AddDoubleArray_AsTable( file_id, "/Test Parameters/PID", Entry_Names,
				    dArray, 3 );

     /* Free the double array */
     free( dArray );

     iArray = (int *) malloc( (size_t) 6*sizeof(int) );

     Entry_Names[0] = strdup( "Num. DOF" ); Entry_Names[1] = strdup( "Num. steps" );
     Entry_Names[2] = strdup( "Num. sub-steps" ); Entry_Names[3] = strdup( "Num. substructures" );
     iArray[0] = (int) InitCnt->Order; iArray[1] = (int) InitCnt->Nstep; 
     iArray[2] = (int) InitCnt->NSubstep; iArray[3] = (int) InitCnt->OrderSub;
     HDF5_AddIntArray_AsTable( file_id, "/Test Parameters/Misc", Entry_Names,
				 iArray, 4 );

     free( iArray );
     
     /* Save the coupling nodes and information */
     if( InitCnt->OrderSub > 0 ){
	  Save_InformationCNodes( file_id, "/Test Parameters/Substructures", CNodes );
     }

     for( i = 0; i < 4; i++ ){
	  free( Entry_Names[i] );
     }

     free( Entry_Names );
}

typedef struct{
     int Position;
     double InitValues[3];
     char *Description;
} HDF5_Exact_UHYDE_t;

typedef struct{
     int Position;
     char *Description;
} HDF5_Exp_Meas_t;

void Save_InformationCNodes( hid_t file_id, const char *Name_path, Coupling_Node *const CNodes )
{
     int i, j, k, count, is_adwin, is_exact, is_uhyde, is_measured;
     hid_t    group_id, strtype, memtype, space, dset;
     hsize_t  dims[1];
     herr_t   status;
     HDF5_Exact_UHYDE_t *Nodes;
     HDF5_Exp_Meas_t *Nodes_Exp;
     TMD_Sim *TMD;
     UHYDE_Sim *UHYDE;
     EXP_Sub *Experimental;

   
     group_id = H5Gcreate(file_id, "/Test Parameters/Substructures", H5P_DEFAULT,
			  H5P_DEFAULT, H5P_DEFAULT);     

     strtype = H5Tcopy( H5T_C_S1 );
     status = H5Tset_size( strtype, H5T_VARIABLE );

     for( i = 0; i <= NUM_TYPE_SUB; i++ ){
	  count = 0;
	  is_adwin = 0; is_exact = 0; is_uhyde=0; is_measured = 0;

	  for( j = 0; j < CNodes->Order; j++ ){
	       if ( i == USE_ADWIN && CNodes->Sub[j].Type == USE_ADWIN ){
		    count = count + 1;
		    is_adwin = 1;	
	       } else if ( i == USE_EXACT && CNodes->Sub[j].Type == USE_EXACT ){
		    count = count + 1;
		    is_exact = 1;
	       } else if ( i == USE_UHYDE && CNodes->Sub[j].Type == USE_UHYDE ){
		    count = count + 1;
		    is_uhyde = 1;
	       } else if ( i == USE_MEASURED && CNodes->Sub[j].Type == USE_MEASURED ){
		    count = count + 1;		    
		    is_measured = 1;
	       }
	  }

	  /* Allocate space */
	  Nodes = NULL;
	  if( is_uhyde || is_exact ){
	       Nodes = (HDF5_Exact_UHYDE_t *) malloc( (size_t) count*sizeof( HDF5_Exact_UHYDE_t ) );
	  } else {
	       Nodes_Exp = (HDF5_Exp_Meas_t *) malloc( (size_t) count*sizeof( HDF5_Exp_Meas_t ) );
	  }

	  /* Copy the matching entries into the newly created structure */
	  k = 0;
	  for( j = 0; j < CNodes->Order; j++ ){
	       if ( i == USE_ADWIN && CNodes->Sub[j].Type == USE_ADWIN ){

		    Nodes_Exp[k].Position = CNodes->Array[j];
		    Experimental = (EXP_Sub *) CNodes->Sub[j].SimStruct;
		    Nodes_Exp[k].Description = strdup( Experimental->Description );
		    k = k + 1;
	       } else if ( i == USE_EXACT && CNodes->Sub[j].Type == USE_EXACT ){
		    Nodes[k].Position = CNodes->Array[j];

		    TMD = (TMD_Sim *) CNodes->Sub[j].SimStruct;
		    Nodes[k].InitValues[0] = TMD->Mass;
		    Nodes[k].InitValues[1] = TMD->Damp;
		    Nodes[k].InitValues[2] = TMD->Stiff;

		    Nodes[k].Description = strdup( TMD->Description );		
		    k = k + 1;
	       } else if ( i == USE_UHYDE && CNodes->Sub[j].Type == USE_UHYDE ){
		    Nodes[k].Position = CNodes->Array[j];

		    UHYDE = (UHYDE_Sim *) CNodes->Sub[j].SimStruct;
		    Nodes[k].InitValues[0] = UHYDE->qyield;
		    Nodes[k].InitValues[1] = UHYDE->qyield/UHYDE->qplastic;
		    Nodes[k].InitValues[2] = UHYDE->qplastic*UHYDE->k;

		    Nodes[k].Description = strdup( UHYDE->Description );
		    k = k + 1;
	       } else if ( i == USE_MEASURED && CNodes->Sub[j].Type == USE_MEASURED ){
		    Nodes_Exp[k].Position = CNodes->Array[j];
		    Experimental = (EXP_Sub *) CNodes->Sub[j].SimStruct;
		    Nodes_Exp[k].Description = strdup( Experimental->Description );
		    k = k + 1;
	       }
	  }

	  if ( (i == USE_ADWIN) && is_adwin ){

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

	  } else if ( (i == USE_EXACT) && is_exact ){
	       memtype = H5Tcreate (H5T_COMPOUND, sizeof( HDF5_Exact_UHYDE_t) );
	       H5Tinsert( memtype, "Eq. column (0-based index)", HOFFSET( HDF5_Exact_UHYDE_t, Position), H5T_NATIVE_INT );
	       H5Tinsert( memtype, "Mass [kg]", HOFFSET( HDF5_Exact_UHYDE_t, InitValues[0]), H5T_NATIVE_DOUBLE );
	       H5Tinsert( memtype, "Damping", HOFFSET( HDF5_Exact_UHYDE_t, InitValues[1]), H5T_NATIVE_DOUBLE );
	       H5Tinsert( memtype, "Stiffness [N/m]", HOFFSET( HDF5_Exact_UHYDE_t, InitValues[2]), H5T_NATIVE_DOUBLE );
	       H5Tinsert( memtype, "Description", HOFFSET( HDF5_Exact_UHYDE_t, Description), strtype );

	       dims[0] = (hsize_t) count; 
	       space = H5Screate_simple (1, dims, NULL);
	       /* Create the dataset */
	       dset = H5Dcreate ( group_id, "/Test Parameters/Substructures/Exact", memtype, space, H5P_DEFAULT,
				  H5P_DEFAULT, H5P_DEFAULT);
	       status = H5Dwrite ( dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, Nodes );
	       
	       for( j = 0; j < count; j++ ){
		    free( Nodes[j].Description );	     
	       }
	       free( Nodes );

	       status = H5Dclose( dset );
	       status = H5Sclose( space );
	       status = H5Tclose( memtype );

	  } else if ( (i == USE_UHYDE) && is_uhyde ){
	       memtype = H5Tcreate (H5T_COMPOUND, sizeof( HDF5_Exact_UHYDE_t) );
	       H5Tinsert( memtype, "Eq. column (0-based index)", HOFFSET( HDF5_Exact_UHYDE_t, Position), H5T_NATIVE_INT );
	       H5Tinsert( memtype, "qyield", HOFFSET( HDF5_Exact_UHYDE_t, InitValues[0]), H5T_NATIVE_DOUBLE );
	       H5Tinsert( memtype, "yield factor", HOFFSET( HDF5_Exact_UHYDE_t, InitValues[1]), H5T_NATIVE_DOUBLE );
	       H5Tinsert( memtype, "Friction force [N]", HOFFSET( HDF5_Exact_UHYDE_t, InitValues[2]), H5T_NATIVE_DOUBLE );
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

	  } else if ( (i == USE_MEASURED) && is_measured ){
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
	  }		  
	  /* Close and release resources. */
     }    

     status = H5Tclose( strtype );
     status = H5Gclose( group_id );
}
     

void HDF5_CreateGroup_TimeIntegration( int hdf5_file, AlgConst *const InitCnt )
{

     hid_t    file_id, group_id;
     herr_t   status;

     file_id = (hid_t) hdf5_file;
     /* Create the Time integration group */
     group_id = H5Gcreate(file_id, "/Time Integration", H5P_DEFAULT,
			  H5P_DEFAULT, H5P_DEFAULT);
     /* Acceleration */
     HDF5_Create_Dataset( file_id, "/Time Integration/Acceleration", (int) InitCnt->Nstep, (int) InitCnt->Order );
     status = H5LTset_attribute_string( file_id, "Time Integration/Acceleration",
					"Units", "m/s^2" );
     /* Velocity */
     HDF5_Create_Dataset( file_id, "/Time Integration/Velocity", (int) InitCnt->Nstep, (int)InitCnt->Order );
     status = H5LTset_attribute_string( file_id, "Time Integration/Velocity",
					"Units", "m/s" );
     /* Displacement */
     HDF5_Create_Dataset( file_id, "/Time Integration/Displacement", (int) InitCnt->Nstep, (int) InitCnt->Order );
     status = H5LTset_attribute_string( file_id, "Time Integration/Displacement",
					"Units", "m" );

     /* Input load */
     //HDF5_Create_Dataset( file_id, "/Time Integration/Input Load", (int) InitCnt->Nstep, (int) InitCnt->Order );
     //status = H5LTset_attribute_string( file_id, "Time Integration/Input Load",
	//				"Units", "N" );

     /* If the number of substructures is greater than zero, store the results */
     if( InitCnt->OrderSub > 0 ){
	  HDF5_Create_Dataset( file_id, "/Time Integration/Coupling force", (int) InitCnt->Nstep, (int) InitCnt->Order );
	  status = H5LTset_attribute_string( file_id, "Time Integration/Coupling force",
					     "Units", "N" );
     }

     /* If the PID is used, then store the values */
     if( InitCnt->PID.P != 0.0 || InitCnt->PID.I != 0.0 || InitCnt->PID.D != 0.0 ){
	  HDF5_Create_Dataset( file_id, "/Time Integration/Error force", (int) InitCnt->Nstep, (int) InitCnt->Order );
	  status = H5LTset_attribute_string( file_id, "Time Integration/Error force",
					     "Units", "N" );
     }

     status = H5Gclose( group_id );

}

void HDF5_Store_TimeHistoryData( int hdf5_file, MatrixVector *const Acc, MatrixVector *const Vel, MatrixVector *const Disp, MatrixVector *const InLoad,
				 MatrixVector *const fc, MatrixVector *const fu, int istep, AlgConst *InitCnt )
{

     hid_t    file_id, group_id;
     herr_t   status;

     file_id = (hid_t) hdf5_file;
     /* Create the Time integration group */

     HDF5_AddResults_to_Dataset( file_id, "/Time Integration/Acceleration", Acc, istep );
     HDF5_AddResults_to_Dataset( file_id, "/Time Integration/Velocity", Vel, istep );
     HDF5_AddResults_to_Dataset( file_id, "/Time Integration/Displacement", Disp, istep );

     //HDF5_AddResults_to_Dataset( file_id, "/Time Integration/Input Load", InLoad, istep );

     /* If the number of substructures is greater than zero, store the results */
     if( InitCnt->OrderSub > 0 ){
	  HDF5_AddResults_to_Dataset( file_id, "/Time Integration/Coupling force", fc, istep );
     }

     /* If the PID is used, then store the values */
     if( InitCnt->PID.P != 0.0 || InitCnt->PID.I != 0.0 || InitCnt->PID.D != 0.0 ){
	  HDF5_AddResults_to_Dataset( file_id, "/Time Integration/Error force", fu, istep );
     }
}

void HDF5_StoreTime( const int hdf5_file, const HDF5_time_t *Time )
{

     hid_t    file_id, memtype, filetype, space, dset, strtype;
     hsize_t  dims[1] = {1};
     herr_t   status;

     file_id = (hid_t) hdf5_file;

     strtype = H5Tcopy( H5T_C_S1 );
     status = H5Tset_size( strtype, H5T_VARIABLE );

     memtype = H5Tcreate (H5T_COMPOUND, sizeof( HDF5_time_t ) );
     H5Tinsert( memtype, "Date", HOFFSET( HDF5_time_t, Date_time), strtype );
     H5Tinsert( memtype, "Duration [ms]", HOFFSET( HDF5_time_t, Elapsed_time), H5T_NATIVE_DOUBLE );
     
     space = H5Screate_simple (1, dims, NULL);
     /* Create the dataset */
     dset = H5Dcreate ( file_id, "/Test Parameters/Date and Duration", memtype, space, H5P_DEFAULT,
			H5P_DEFAULT, H5P_DEFAULT);
     status = H5Dwrite ( dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, Time );
    
     status = H5Dclose( dset );
     status = H5Sclose( space );
     status = H5Tclose( memtype );
     status = H5Tclose( strtype );
}

void HDF5_StoreADwinData( const int hdf5_file, const double *Array, const int Length )
{
     int i;
     hid_t file_id;
     const char *Entry_Names[NUM_CHANNELS] = { "Sub-step",
				     "Control displacement actuator 1 [m]",
				     "Control displacement actuator 2 [m]",
				     "Measured displacement actuator 1 [m]",
				     "Measured displacement actuator 2 [m]",
				     "Control acceleration Actuator 1 [m/s^2]",
				     "Control acceleration Actuator 2 [m/s^2]",
				     "Measured acceleration actuator 1 [m/s^2]",
				     "Measured acceleration actuator 2 [m/s^2]",
				     "Acceleration TMD (y-direction) [m/s^2]",
				     "Acceleration TMD (x-direction) [m/s^2]",
				     "Acceleration Base-Frame (x-direction) [m/s^2]",
				     "Coupling Force 1 (y-direction) [N]",
				     "Coupling Force 2 (y-direction) [N]",
				     "Coupling Force 3 (x-direction) [N]",
				     "Displacement TMD (x-direction) [m]",
				     "Displacement TMD (y-direction) [m]",
				     "Control pressure [Pa]",
				     "Measured pressure [Pa]",
				     "Displacement at the coupling node [m]",
				     "Time spent doing the sub-step [ms]", 
				     "Sub-step time [ms]",
				     "Synchronisation time between PC and ADwin [ms]",
				     "Time between the first sub-step and arrival of the new update [ms]"
     };
     size_t *Offset;
     hid_t *Field_types;
     hsize_t Chunk_size;

     file_id = (hid_t) hdf5_file;
     Offset = (size_t *) malloc( (size_t) NUM_CHANNELS*sizeof(size_t) );
     Field_types = (hid_t *) malloc( (size_t) NUM_CHANNELS*sizeof(hid_t) );

     /* Chunk size will have the length of 100 sub-steps */
     Chunk_size = (hsize_t) 100*sizeof(double)*(hsize_t) NUM_CHANNELS;

     /* There are 24 channels defined in RoutinesADwin.h */
     
     for ( i = 0; i < NUM_CHANNELS; i++ ){
	  Offset[i] = (size_t) i*H5Tget_size(H5T_NATIVE_DOUBLE);
	  Field_types = H5T_NATIVE_DOUBLE;
     }

     H5TBmake_table( "/ADwin measurements", file_id, "Measurements", (hsize_t) NUM_CHANNELS, (hsize_t) Length, sizeof(double)*(size_t) NUM_CHANNELS, Entry_Names, Offset, Field_types, Chunk_size, NULL, 0, Array );

     free( Offset );
     free( Field_types );

}

void HDF5_AddDoubleArray_AsTable( hid_t file_id, const char *Name_path, char **Names, const double *Array, int Num_param )
{
     int      i;
     hid_t    memtype, filetype, space, dset;
     hsize_t  dims[1] = {1};
     herr_t   status;


     /* Create the compound datatype for memory. */
     memtype = H5Tcreate (H5T_COMPOUND,
			  Num_param*H5Tget_size(H5T_NATIVE_DOUBLE));
     for( i = 0; i < Num_param; i++ ){
	  status = H5Tinsert (memtype, Names[i],
			      (size_t) i*H5Tget_size(H5T_NATIVE_DOUBLE),
			      H5T_NATIVE_DOUBLE);
     }
    
     /* Create the compound datatype for the file.*/
     filetype = H5Tcreate (H5T_COMPOUND,
			   Num_param*H5Tget_size(H5T_NATIVE_DOUBLE));
     for( i = 0; i < Num_param; i++ ){
	  status = H5Tinsert (filetype, Names[i],
			      (size_t) i*H5Tget_size(H5T_NATIVE_DOUBLE),
			      H5T_NATIVE_DOUBLE);
     }
 
     /* 
      * Create dataspace.  Setting maximum size to NULL sets the maximum
      * size to be the current size.
      */    
     space = H5Screate_simple (1, dims, NULL);

     /* Create the dataset and write the compound data to it. */
     dset = H5Dcreate ( file_id, Name_path, filetype, space, H5P_DEFAULT,
			H5P_DEFAULT, H5P_DEFAULT);
     status = H5Dwrite (dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, Array);

     /* Close and release resources. */
     status = H5Tclose( filetype );
     status = H5Tclose( memtype);
     status = H5Dclose (dset);
     status = H5Sclose (space);
}

void HDF5_AddIntArray_AsTable( hid_t file_id, const char *Name_path, char **Names, const int *Array, int Num_param )
{
     int      i;
     hid_t    memtype, filetype, space, dset;
     hsize_t  dims[1] = {1};
     herr_t   status;


     /* Create the compound datatype for memory. */
     memtype = H5Tcreate (H5T_COMPOUND,
			  Num_param*H5Tget_size(H5T_NATIVE_INT));
     for( i = 0; i < Num_param; i++ ){
	  status = H5Tinsert (memtype, Names[i],
			      (size_t) i*H5Tget_size(H5T_NATIVE_INT),
			      H5T_NATIVE_INT);
     }

    
     /* Create the compound datatype for the file.*/
     filetype = H5Tcreate (H5T_COMPOUND,
			   Num_param*H5Tget_size(H5T_NATIVE_INT));
     for( i = 0; i < Num_param; i++ ){
	  status = H5Tinsert (filetype, Names[i],
			      (size_t) i*H5Tget_size(H5T_NATIVE_INT),
			      H5T_NATIVE_INT);
     }
 
     /* 
      * Create dataspace.  Setting maximum size to NULL sets the maximum
      * size to be the current size.
      */    
     space = H5Screate_simple (1, dims, NULL);

     /* Create the dataset and write the compound data to it. */
     dset = H5Dcreate ( file_id, Name_path, filetype, space, H5P_DEFAULT,
			H5P_DEFAULT, H5P_DEFAULT);
     status = H5Dwrite (dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, Array);

     /* Close and release resources. */
     status = H5Dclose (dset);
     status = H5Sclose (space);
}

void HDF5_Create_Dataset( hid_t file_id, const char *Path_name, int Nstep, int Order )
{
     hsize_t dims[2] = {1, (hsize_t) Order };
     hsize_t max_dims[2] = {(hsize_t) Nstep, (hsize_t) Order };

     hid_t cparms, dataspace_id, dataset_id;
     herr_t status;
	
     //  unsigned szip_options_mask;
     // unsigned szip_pixels_per_block;

     dataspace_id = H5Screate_simple( 2, dims, max_dims );

     cparms = H5Pcreate( H5P_DATASET_CREATE);
     status = H5Pset_chunk( cparms, 2, dims );

     //status = H5Pset_deflate (cparms, 6); 

     //szip_options_mask = H5_SZIP_NN_OPTION_MASK;
     //szip_pixels_per_block = 16;
     //status = H5Pset_szip (cparms, szip_options_mask, szip_pixels_per_block);

     dataset_id = H5Dcreate( file_id, Path_name, H5T_NATIVE_DOUBLE, dataspace_id,
			      H5P_DEFAULT, cparms, H5P_DEFAULT);

     status = H5Pclose( cparms );
     status = H5Dclose( dataset_id );
     status = H5Sclose( dataspace_id );
}

void HDF5_AddResults_to_Dataset( hid_t file_id, const char *Path_name, MatrixVector *const Data, int Step_count )
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
     offset[0] = (Step_count - 1);
     offset[1] = 0;
     dims[0] = 1; dims[1] = (hsize_t) Data->Rows;
     status = H5Sselect_hyperslab( filespace_id, H5S_SELECT_SET, offset, NULL, dims, NULL );
     
     memspace_id = H5Screate_simple( 2, dims, NULL );
     
     status = H5Dwrite( dataset_id, H5T_NATIVE_DOUBLE, memspace_id, filespace_id, H5P_DEFAULT, Data->Array );

     status = H5Dclose( dataset_id );
     status = H5Sclose( memspace_id );
     status = H5Sclose( filespace_id );
}

void HDF5_CloseFile( int hdf5_file )
{    
     herr_t status;
     hid_t file_id;

     file_id = (int) hdf5_file;
     /* Terminate access to the file. */
     status = H5Fclose( file_id );
}
