/**
 * \file HDF5_Operations.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 23rd of January 2013
 *
 * \brief HDF5 file operations.
 *
 * Operations used in creating and manipulating HDF5 files to save the data of the experiment.
 */

#ifndef HDF5_OPERATIONS_H
#define HDF5_OPERATIONS_H_

#include "hdf5.h"
#include "Algorithm_Aux.h"
#include "MatrixVector.h"
#include "Substructure.h"

#if _MPI_
#include <mpi.h>
#endif

typedef struct HDF5_Exact_UHYDE {
     int Position;
     HYSL_FLOAT InitValues[3];
     char *Description;
} HDF5_Exact_UHYDE_t;

typedef struct HDF5_Exp_Meas {
     int Position;
     char *Description;
} HDF5_Exp_Meas_t;

int HDF5_CreateFile( const char *Filename );

#if _MPI_
void HDF5_CreateFile_MPI( MPI_Comm Comm, const char *Filename, hid_t *const file_id );
void HDF5_Store_TimeHistoryData_MPI( const int hdf5_file, const PMatrixVector_t *const Acc, const PMatrixVector_t *const Vel, const PMatrixVector_t *const Disp,
				    const PMatrixVector_t *const fc, const PMatrixVector_t *const fu, int istep, const int nprow, const int myrow, const AlgConst_t *const InitCnt );
void HDF5_AddResults_to_Dataset_MPI( const hid_t file_id, const char *Path_name, const PMatrixVector_t *const Data, const int Step_count, const int nprow, const int myrow );
void HDF5_Store_Time_MPI( const hid_t hdf5_file, const SaveTime_MPI_t *const Time );
#endif

void HDF5_CreateGroup_Parameters( const hid_t hdf5_file, const AlgConst_t *const InitCnt, const CouplingNode_t *const CNode, const HYSL_FLOAT *const Acc, const HYSL_FLOAT *const Vel, const HYSL_FLOAT *const Disp );
void Save_InformationCNodes( const hid_t file_id, const char *Name_path, const CouplingNode_t *const CNodes );
void HDF5_CreateGroup_TimeIntegration( const hid_t hdf5_file, const AlgConst_t *const InitCnt );
void HDF5_Store_TimeHistoryData( const hid_t hdf5_file, const MatrixVector_t *const Acc, const MatrixVector_t *const Vel, const MatrixVector_t *const Disp, const MatrixVector_t *const fc, const MatrixVector_t *const fu, const int istep, const AlgConst_t *const InitCnt );
void HDF5_Store_Time( const hid_t hdf5_file, const SaveTime_t *const Time );
void HDF5_StoreADwinData( const int hdf5_file, const HYSL_FLOAT *Array, char **Entry_Names, const int Length );


/**
 * \brief Reads the data from ADwin and stores it in a HDF5 file.
 *
 * The data stored in the specified array ADwin is read and stored into a HDF5 file. The dataset will have as
 * many columns as the specified number of chanels and as many rows as \f$N_{steps}\cdot N_{substeps}\f$.
 *
 * \pre
 * - ADwin must be properly booted and the substructure test must be finished.
 * - Data must be a properly initialised array of length \f$L \geq N_{steps}\cdot N_{substeps}\cdot
 *   N_{channels}\f$.
 * - The hdf5_file identifier must point to a open hdf5 file.
 *
 * \param[in] hdf5_file HDF5 file identifier.
 * \param[in] Num_Steps Number of steps.
 * \param[in] Num_Sub Number of sub-steps.
 * \param[in] Num_Channels Number of data acquisition channels.
 * \param[in] Chan_Names List of channel names.
 * \param[in] DataIndex points out which data array is to be accessed in ADwin.
 *
 * \post The new dataset will have the the data stored in the array identified by \c Data_Index in ADwin with:
 * - A row containing a brief description of each channel.
 * - Number of columns equal to the number of channels.
 * - Number of rows equal to \f$N_{step}\cdot N_{substep}\f$.
 */
void ADwin_SaveData_HDF5( const int hdf5_file, const int Num_Steps, const int Num_Sub,
			  const int Num_Channels, char **Chan_Names, const int DataIndex );

void HDF5_AddFloatArray_AsTable( const hid_t file_id, const char *Name_path, char **Names, const HYSL_FLOAT *Array, const int Num_param, const int Length );
void HDF5_AddIntArray_AsTable( const hid_t file_id, const char *Name_path, char **Names, const int *Array, const int Num_param );

void HDF5_Create_Dataset( const hid_t file_id, const char *Path_name, const int Nstep, const int Order );
void HDF5_AddResults_to_Dataset( const hid_t file_id, const char *Path_name, const MatrixVector_t *const Data, const int Step_count );
void HDF5_CloseFile( hid_t *const hdf5_file );

#endif
