/**
 * \file Scalapack_Aux.h
 * \author Ferran Obón Santacana
 * \version 1.0
 * \date 2nd of April 2013
 *
 * \brief Prototypes of the function for the auxiliary ScaLapack routines.
 *
 * This file contains the prototypes of the functions for the auxiliary ScaLapack routines. They should only
 * be used with the MPI version of the library.
 */

#ifndef SCALAPACK_AUX_H_
#define SCALAPACK_AUX_H_

void descinit_(int *desc, int *m, int *n, int *mb, int *nb, int *irsrc, int *icsrc, int *ictxt, int *lld,
	       int *info);
int numroc_(int *n, int *nb, int *iproc, int *isrcproc, int *nprocs);
void infog2l_( int *grindx, int *gcindx, int *desc, int *nprow, int *npcol, int *myrow, int *mycol, int *lrindx,
	       int *lcindx, int *rsrc, int *csrc );

#endif /* CBLACS_H_ */
