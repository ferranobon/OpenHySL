/**
 * \file Cblacs.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 2nd of April 2013
 *
 * \brief Prototypes of the function for the Cblacs routines.
 *
 * This file contains the prototypes of the functions for the Cblacs routines. They should only be used with
 * the MPI version of the library.
 */

#ifndef CBLACS_H_
#define CBLACS_H_

#include <mpi.h>

void Cblacs_get( int icontxt, int what, int *val );
void Cblacs_gridinit(int *icontxt, char *order, int nprow, int npcol);
void Cblacs_gridinfo( int icontxt, int *nprow, int *npcol, int *myrow, int *mycol );
void Cblacs_gridmap( int *icontxt, int *usermap, int ldumap, int nprow, int npcol );
int Cblacs_pnum( int icontxt, int prow, int pcol );

int Csys2blacs_handle(MPI_Comm SysCtxt);
int Cfree_blacs_system_handle(int BlacsHandle);
void Cblacs_barrier(int icontxt, char *scope);
int Cblacs_gridexit( int icontxt );

#endif /* CBLACS_H_ */
