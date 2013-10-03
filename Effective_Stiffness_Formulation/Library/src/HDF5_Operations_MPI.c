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

void HDF5_CreateFile_MPI( MPI_Comm Comm, const char *Filename, hid_t *file_id, hid_t *plist_id )
{

     MPI_Info info = MPI_INFO_NULL;

     *plist_id = H5Pcreate ( H5P_FILE_ACCESS );
     H5Pset_fapl_mpio( *plist_id, Comm, info );

     /* Create a new file using default properties. */
     *file_id = H5Fcreate( Filename, H5F_ACC_TRUNC, H5P_DEFAULT, *plist_id );
}

void HDF5_Store_TimeHistoryData_MPI( int hdf5_file, PMatrixVector_t *const Acc, PMatrixVector_t *const Vel, PMatrixVector_t *const Disp, PMatrixVector_t *const InLoad,
				 PMatrixVector_t *const fc, PMatrixVector_t *const fu, int istep, AlgConst_t *InitCnt )
{

     hid_t    file_id;

     file_id = (hid_t) hdf5_file;
     /* Create the Time integration group */

     HDF5_AddResults_to_Dataset_MPI( file_id, "/Time Integration/Acceleration", Acc, istep );
     HDF5_AddResults_to_Dataset_MPI( file_id, "/Time Integration/Velocity", Vel, istep );
     HDF5_AddResults_to_Dataset_MPI( file_id, "/Time Integration/Displacement", Disp, istep );

     // HDF5_AddResults_to_Dataset( file_id, "/Time Integration/Input Load", InLoad, istep );

     /* If the number of substructures is greater than zero, store the results */
     if( InitCnt->OrderSub > 0 ){
	  HDF5_AddResults_to_Dataset_MPI( file_id, "/Time Integration/Coupling force", fc, istep );
     }

     /* If the PID is used, then store the values */
     if( InitCnt->PID.P != 0.0 || InitCnt->PID.I != 0.0 || InitCnt->PID.D != 0.0 ){
	  HDF5_AddResults_to_Dataset_MPI( file_id, "/Time Integration/Error force", fu, istep );
     }
}

void HDF5_AddResults_to_Dataset_MPI( hid_t file_id, const char *Path_name, PMatrixVector_t *const Data, const int Step_count )
{
     hid_t   dataset_id, filespace_id, memspace_id;
     herr_t  status;
     hsize_t size[2], offset[2], dims[2];

     size[0] = (hsize_t) Step_count; /* Num_Steps starts at 1 */
     size[1] = (hsize_t) Data->GlobalSize.Row;

     /* Create the data space for the dataset. */
     dataset_id = H5Dopen( file_id, Path_name, H5P_DEFAULT );
     H5Dextend( dataset_id, size );

     filespace_id = H5Dget_space( dataset_id );
     offset[0] = ((hsize_t) Step_count - 1);
     offset[1] = 0;
     dims[0] = 1; dims[1] = (hsize_t) Data->GlobalSize.Row;
     status = H5Sselect_hyperslab( filespace_id, H5S_SELECT_SET, offset, NULL, dims, NULL );
     
     memspace_id = H5Screate_simple( 2, dims, NULL );
     
     status = H5Dwrite( dataset_id, H5T_NATIVE_DOUBLE, memspace_id, filespace_id, H5P_DEFAULT, Data->Array );

     status = H5Dclose( dataset_id );
     status = H5Sclose( memspace_id );
     status = H5Sclose( filespace_id );
}

void HDF5_CloseFile_MPI( int *const hdf5_file, int *const hdf5_plist )
{

     /* Close property list */
     H5Pclose( *hdf5_plist );

     /* Close the file */
     H5Fclose( *hdf5_file );
}
