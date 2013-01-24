#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Initiation.h"

/* hdf5 header files */
#include "hdf5.h"
#include "hdf5_hl.h"

int main() {

     hid_t       file_id, group_id; 
     hsize_t     dims[1] = {1};
     
     double Array[5];
     char **Entry_Names;

     
 
          
     Entry_Names[0] = strdup( "Beta" );
     Entry_Names[1] = strdup( "Gamma" );
     Array[0] = 0.25; Array[1] = 0.5;
     Add_Array_as_Table( file_id, "/Test Parameters/Newmark", Entry_Names,
			 Array, 2 );
     free( Entry_Names[0] ); free( Entry_Names[1] );

     

     

     dims[0] = 4; dims[1] = 6; 

     /* Add atribute. */
     
     /* Close the group. */
     status = H5Gclose(group_id);

     

     free( Entry_Names );

     return 0;
}

hid_t HDF5_CreateFile( const char *Filename, AlgConst *InitCnt )
{
     int      i;
     hid_t    file_id, group_id;
     herr_t   status;
     hsize_t  dims[2];
     double   *dArray;
     int      *iArray;
     char     **Entry_Names;     

     Entry_Names = (char **) malloc( (size_t) 6*sizeof(char*) );
     dArray = (double *) malloc( (size_t) 3*sizeof(double) );

     /* Create a new file using default properties. */
     file_id = H5Fcreate( Filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
     
     /* Create a group named "/Time History" in the file. */
     group_id = H5Gcreate(file_id, "/Test Parameters", H5P_DEFAULT,
			  H5P_DEFAULT, H5P_DEFAULT);     
     status = H5Gclose(group_id);

     /* Add Integration parameters */
     Entry_Names[0] = strdup( "Beta" ); Entry_Names[1] = strdup( "Gamma" );
     Entry_Names[2] = strdup( "Time step" );
     dArray[0] = InitCnt->Newmark.Beta; dArray[1] = InitCnt->Newmark.Gamma;
     dArray[2] = InitCnt->Delta_t;
     HDF5_Add_DoubleArray_As_Table( file_id, "/Test Parameters/Integration parameters", Entry_Names,
				    dArray, 3 );
     free( Entry_Names[0] ); free( Entry_Names[1] ); free( Entry_Names[2] );

     /* Add Rayleigh damping values */
     Entry_Names[0] = strdup( "Alpha" ); Entry_Names[1] = strdup( "Beta" );
     dArray[0] = InitCnt->Rayleigh.Alpha; dArray[1] = InitCnt->Rayleigh.Beta;
     HDF5_Add_DoubleArray_As_Table( file_id, "/Test Parameters/Rayleigh", Entry_Names,
				    dArray, 2 );
     free( Entry_Names[0] ); free( Entry_Names[1] );

     /* Add PID parameters */
     Entry_Names[0] = strdup( "P" ); Entry_Names[1] = strdup( "I" );
     Entry_Names[2] = strdup( "D" );
     dArray[0] = InitCnt->PID.P; dArray[1] = InitCnt->PID.I;
     dArray[2] = InitCnt->PID.D;
     HDF5_Add_DoubleArray_As_Table( file_id, "/Test Parameters/PID", Entry_Names,
				    dArray, 3 );

     /* Free the double array */
     free( dArray );

     iArray = (int *) malloc( (size_t) 6*sizeof(int) );

     Entry_Names[0] = strdup( "Order" ); Entry_Names[1] = strdup( "Num. steps" );
     Entry_Names[2] = strdup( "Num. sub-steps" ); Entry_Names[3] = strdup( "Num. substructures" );
     iArray[0] = InitCnt->Order; iArray[1] = InitCnt->Nstep; 
     iArray[2] = InitCnt->NSubstep; iArray[3] = InitCnt->OrderSub;
     HDF5_Add_IntArray_As_Table( file_id, "/Test Parameters/Misc", Entry_Names,
				 iArray, 4 );
     for( i = 0; i < 4; i++ ){
	  free( Entry_Names[i] );
     }

     Entry_Names[0] = strdup( "X" ); Entry_Names[1] = strdup( "Y" );
     Entry_Names[2] = strdup( "Z" ); Entry_Names[3] = strdup( "Rot. X" );
     Entry_Names[2] = strdup( "Rot. Y" ); Entry_Names[3] = strdup( "Rot. Z" );
     iArray[0] = InitCnt->O 
     iArray[2] = InitCnt->NSubstep; iArray[3] = InitCnt->OrderSub;
     iArray[2] = InitCnt->NSubstep; iArray[3] = InitCnt->OrderSub;
     HDF5_Add_IntArray_As_Table( file_id, "/Test Parameters/Misc", Entry_Names,
				 iArray, 4 );
     for( i = 0; i < 4; i++ ){
	  free( Entry_Names[i] );
     }


     /* Create the Time integration group */
     group_id = H5Gcreate(file_id, "/Time Integration", H5P_DEFAULT,
			  H5P_DEFAULT, H5P_DEFAULT);

     dims[0] = (hsize_t) InitCnt->NStep; dims[1] = (hsize_t) InitCnt->Order;
     HDF5_AddDataAs2DDataset( file_id, "Time Integration/Acceleration", Array,
			      dims );
     status = H5LTset_attribute_string( file_id, "Time Integration/Acceleration",
					"Units", "m/s^2" );

     HDF5_AddDataAs2DDataset( file_id, "Time Integration/Velocity", Array,
			      dims );
     status = H5LTset_attribute_string( file_id, "Time Integration/Velocity",
					"Units", "m/s" );

     HDF5_AddDataAs2DDataset( file_id, "Time Integration/Displacement", Array,
			      dims );
     status = H5LTset_attribute_string( file_id, "Time Integration/Displacement",
					"Units", "m" );
     status = H5Gclose( group_id );

  
     free( Entry_Names );

     return file_id;
}

void HDF5_AddArrayAsTable( hid_t file_id, const char *Name_path, char **Names, const double *Array, int Num_param )
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
     status = H5Dclose (dset);
     status = H5Sclose (space);
}

void HDF5_AddDataAs2DDataset( hid_t file_id, const char *Path_name, const double *Array, hsize_t *const dims )
{
     hid_t   dataspace_id, dataset_id;
     herr_t  status;

     /* Create the data space for the dataset. */
     dataspace_id = H5Screate_simple(2, dims, NULL);
     
     /* Create the dataset. */
     dataset_id = H5Dcreate( file_id, Path_name, H5T_NATIVE_DOUBLE,
			     dataspace_id, H5P_DEFAULT, H5P_DEFAULT,
			     H5P_DEFAULT);

     /* End access to the dataset and release resources used by it. */
     status = H5Dclose(dataset_id);

     /* Terminate access to the data space. */ 
     status = H5Sclose(dataspace_id);
}

void HDF5_CloseFile( hid_t file_id )
{     
     /* Terminate access to the file. */
     status = H5Fclose( file_id );
}
