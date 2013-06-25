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

#include <sys/time.h> /* For struct timeval */
#include <time.h> /* For time_t */

typedef struct {
     time_t Date_start;
     char *Date_time;
     double Elapsed_time;
     struct timeval start;
     struct timeval end;
} HDF5time_t;

int HDF5_CreateFile( const char *Filename );
void HDF5_CreateGroup_Parameters( int hdf5_file, AlgConst_t *const InitCnt, CouplingNode_t *const CNode );
void Save_InformationCNodes( hid_t file_id, const char *Name_path, CouplingNode_t *const CNodes );
void HDF5_CreateGroup_TimeIntegration( int hdf5_file, AlgConst_t *const InitCnt );
void HDF5_Store_TimeHistoryData( int hdf5_file, MatrixVector_t *const Acc, MatrixVector_t *const Vel, MatrixVector_t *const Disp, MatrixVector_t *const InLoad, MatrixVector_t *const fc, MatrixVector_t *const fu, int istep, AlgConst_t *InitCnt );
void HDF5_StoreTime( int hdf5_file, const HDF5time_t *Time );
void HDF5_StoreADwinData( const int hdf5_file, const double *Array, const char **Entry_Names, const int Length );


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
void ADwin_SaveData_HDF5( const int hdf5_file, const unsigned int Num_Steps, const unsigned int Num_Sub,
			  const unsigned short int Num_Channels, const char **Chan_Names, const int DataIndex );

void HDF5_AddDoubleArray_AsTable( hid_t file_id, const char *Name_path, char **Names, const double *Array, int Num_param, int Length );
void HDF5_AddIntArray_AsTable( hid_t file_id, const char *Name_path, char **Names, const int *Array, int Num_param );

void HDF5_Create_Dataset( hid_t file_id, const char *Path_name, int Nstep, int Order );
void HDF5_AddResults_to_Dataset( hid_t file_id, const char *Path_name, MatrixVector_t *const Data, int Step_count );
void HDF5_CloseFile( int hdf5_file );
#endif
