#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "HDF5_Operations.h"
#include "MatrixVector.h"
#include "Algorithm_Aux.h"

#include "Print_Messages.h"

#include "Substructure.h"
#include "Substructure_Exact.h"
#include "Substructure_UHYDEfbr.h"
#include "Substructure_SimMeasured.h"
#include "Substructure_Experimental.h"

#include "mpi.h"

#if _ADWIN_
#include <libadwin.h>        /* ADwin routines Boot(), Set_DeviceNo(), ... */
#include <libadwin/errno.h>  /* ADwin error handling */
#include "ADwin_Routines.h"
#endif

/* hdf5 header files */
#include "hdf5.h"
#include "hdf5_hl.h"

void HDF5_CreateFile_MPI( MPI_Comm Comm, const char *Filename, hid_t *const file_id )
{
     hid_t plist_id;

     MPI_Info info = MPI_INFO_NULL;

     plist_id = H5Pcreate ( H5P_FILE_ACCESS );
     H5Pset_fapl_mpio( plist_id, Comm, info );

     /* Create a new file using default properties. */
     *file_id = H5Fcreate( Filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id );


     /* Close property list */
     H5Pclose( plist_id );

}

void HDF5_Store_TimeHistoryData_MPI( const int hdf5_file, const PMatrixVector_t *const Acc,
				     const PMatrixVector_t *const Vel, const PMatrixVector_t *const Disp,
				     const PMatrixVector_t *const fc, const PMatrixVector_t *const fu,
				     int istep, const int nprow, const int myrow, const AlgConst_t *const InitCnt )
{

     hid_t    file_id;

     file_id = (hid_t) hdf5_file;
     /* Create the Time integration group */

     HDF5_AddResults_to_Dataset_MPI( file_id, "/Time Integration/Acceleration", Acc, istep, nprow, myrow );
     HDF5_AddResults_to_Dataset_MPI( file_id, "/Time Integration/Velocity", Vel, istep, nprow, myrow );
     HDF5_AddResults_to_Dataset_MPI( file_id, "/Time Integration/Displacement", Disp, istep, nprow, myrow );

     /* If the number of substructures is greater than zero, store the results */
     if( InitCnt->OrderSub > 0 ){
	  HDF5_AddResults_to_Dataset_MPI( file_id, "/Time Integration/Coupling force", fc, istep, nprow, myrow );
     }

     /* If the PID is used, then store the values */
     if( InitCnt->PID.P != 0.0 || InitCnt->PID.I != 0.0 || InitCnt->PID.D != 0.0 ){
	  HDF5_AddResults_to_Dataset_MPI( file_id, "/Time Integration/Error force", fu, istep, nprow, myrow );
     }
}

void HDF5_AddResults_to_Dataset_MPI( const hid_t file_id, const char *Path_name, const PMatrixVector_t *const Data,
				     const int Step_count, const int nprow, const int myrow )
{
     hid_t   plist_id, dataset_id, filespace_id, memspace_id;
     herr_t  status;
     hsize_t size[2], offset[2], dims[2], stride[2], count[2], block[2];

     size[0] = (hsize_t) Step_count; /* Num_Steps starts at 1 */
     size[1] = (hsize_t) Data->GlobalSize.Row;

     /* Create the data space for the dataset. */
     dataset_id = H5Dopen( file_id, Path_name, H5P_DEFAULT );
     H5Dextend( dataset_id, size );

     dims[0] = 1;
     dims[1] = (hsize_t) Data->LocalSize.Row*(hsize_t) Data->LocalSize.Col;

     offset[0] = (hsize_t) Step_count - 1;

     if( dims[1] != 0 ){
	  offset[1] = (hsize_t) Data->BlockSize.Row*(hsize_t) myrow;

	  count[0] = 1;
	  count[1] = (hsize_t) dims[1]/(hsize_t) Data->BlockSize.Row;

	  stride[0] = 1;
	  stride[1] = (hsize_t) Data->BlockSize.Row*(hsize_t) nprow;
	  
	  block[0] = 1;
	  block[1] = (hsize_t) Data->BlockSize.Row;
     } else {
	  offset[1] = 0;
	  count[0] = 0; count[1] = 0;
	  stride[0] = 1; stride[1] = 1;
	  block[0] = 0; block[1] = 0;
     }

     /*
      * Select hyperslab in the file.
      */
     filespace_id = H5Dget_space( dataset_id );
     memspace_id  = H5Screate_simple( 2, dims, NULL);

     if( dims[1] != 0 ){
	  status = H5Sselect_hyperslab( filespace_id, H5S_SELECT_SET, offset, stride, count, block );
     } else {
	  status = H5Sselect_none( memspace_id );
	  status = H5Sselect_none( filespace_id );
     }
    
     /*
     * Create property list for collective dataset write.
     */
     plist_id = H5Pcreate( H5P_DATASET_XFER );
     H5Pset_dxpl_mpio( plist_id, H5FD_MPIO_COLLECTIVE );

     status = H5Dwrite( dataset_id, H5T_NATIVE_DOUBLE, memspace_id, filespace_id, plist_id, Data->Array );

     status = H5Dclose( dataset_id );
     status = H5Sclose( memspace_id );
     status = H5Pclose( plist_id );
     status = H5Sclose( filespace_id );
}

void HDF5_Store_Time_MPI( const hid_t hdf5_file, const SaveTime_MPI_t *const Time )
{

     hid_t    memtype, space, dset, strtype;
     hsize_t  dims[1] = {1};
     herr_t   status;
  
     strtype = H5Tcopy( H5T_C_S1 );
     status = H5Tset_size( strtype, MPI_TIME_SLENGTH );

     memtype = H5Tcreate (H5T_COMPOUND, sizeof( SaveTime_MPI_t ) );
     H5Tinsert( memtype, "Date", HOFFSET( SaveTime_MPI_t, Date_time), strtype );
     H5Tinsert( memtype, "Duration [ms]", HOFFSET( SaveTime_MPI_t, Elapsed_time), H5T_NATIVE_DOUBLE );
   
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
