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
#include "Initiation.h"
#include "MatrixVector.h"
#include <sys/time.h> /* For struct timeval */
#include <time.h> /* For time_t */

typedef struct {
     time_t Date_start;
     char *Date_time;
     double Elapsed_time;
     struct timeval start;
     struct timeval end;
} HDF5_time_t;

int HDF5_CreateFile( const char *Filename );
void HDF5_CreateGroup_Parameters( int hdf5_file, AlgConst *const InitCnt, Coupling_Node *const CNode );
void Save_InformationCNodes( hid_t file_id, const char *Name_path, Coupling_Node *const CNodes );
void HDF5_CreateGroup_TimeIntegration( int hdf5_file, AlgConst *const InitCnt );
void HDF5_Store_TimeHistoryData( int hdf5_file, MatrixVector *const Acc, MatrixVector *const Vel, MatrixVector *const Disp, MatrixVector *const InLoad, MatrixVector *const fc, MatrixVector *const fu, int istep, AlgConst *InitCnt );
void HDF5_StoreTime( int hdf5_file, const HDF5_time_t *Time );
void HDF5_StoreADwinData( const int hdf5_file, const double *Array, const int Length );
void HDF5_AddDoubleArray_AsTable( hid_t file_id, const char *Name_path, char **Names, const double *Array, int Num_param );
void HDF5_AddIntArray_AsTable( hid_t file_id, const char *Name_path, char **Names, const int *Array, int Num_param );

void HDF5_Create_Dataset( hid_t file_id, const char *Path_name, int Nstep, int Order );
void HDF5_AddResults_to_Dataset( hid_t file_id, const char *Path_name, MatrixVector *const Data, int Step_count );
void HDF5_CloseFile( int hdf5_file );
#endif
