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
void HDF5_StoreADwinData( const int hdf5_file, const double *Array, const int Length );
void HDF5_AddDoubleArray_AsTable( hid_t file_id, const char *Name_path, char **Names, const double *Array, int Num_param );
void HDF5_AddIntArray_AsTable( hid_t file_id, const char *Name_path, char **Names, const int *Array, int Num_param );

void HDF5_Create_Dataset( hid_t file_id, const char *Path_name, int Nstep, int Order );
void HDF5_AddResults_to_Dataset( hid_t file_id, const char *Path_name, MatrixVector_t *const Data, int Step_count );
void HDF5_CloseFile( int hdf5_file );
#endif
