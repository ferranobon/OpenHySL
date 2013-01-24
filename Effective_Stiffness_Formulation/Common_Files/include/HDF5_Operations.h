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

void HDF5_Add_DoubleArray_As_Table( hid_t file_id, const char *Name_path, char **Names, const double *Array, int Num_param );
void HDF5_Add_IntArray_As_Table( hid_t file_id, const char *Name_path, char **Names, const double *Array, int Num_param );
void HDF5_AddDataAs2DDataset( hid_t file_id, const char *Path_name, const double *Array, hsize_t *const dims );
void HDF5_CloseFile( hid_t file_id );
#endif
