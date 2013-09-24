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

/**
 * \brief Initialises the descriptor vector (1D array) necessary for distributed matrices/vectors (ScaLAPACK
 * routine) with the 8 input arguments: <em>m, n, mb, nb, irsrc, icsrc, ictxt</em> and \em lld.
 *
 * \param[out] desc  The array descriptor of a distributed matrix/vector.
 * \param[in]  m     (global) Number of rows in the distributed matrix/vector. \f$m \geq 0\f$.
 * \param[in]  n     (global) Number of columns in the distributed matrix/vector. \f$m \geq 0\f$.
 * \param[in]  mb    (global) Blocking factor used to distribute the rows among processes. \f$mb \geq 1\f$.
 * \param[in]  nb    (global) Blocking factor used to distribute the columns among processes. \f$nb \geq 1\f$.
 * \param[in]  irsrc (global) The process row over which the first row of the matrix is distributed. \f$0 \leq
 *                   irsrc < nprow\f$, where \f$nprow\f$ is the number of process rows in the grid.
 * \param[in]  icsrc (global) The process column over which the first column of the matrix is
 *                   distributed. \f$0 \leq icsrc < npcol\f$, where \f$npcol\f$ is the number of process
 *                   columns in the grid.
 * \param[in]  ictxt (global) The BLACS context handle, indicating the global context of the operation on the
 *                   matrix. The context itself is glogal.
 * \param[in]  lld   (local) The leading dimension of the local array storing the local blocks of the
 *                   distributed matrix. \f$lld \geq max(1,LOCr(m))\f$, where LOCr() denotes the number of
 *                   elements of the matrix that a process would perceive if the matrix were distributed over
 *                   \f$p\f$ processes of its process column. It can be determined using numroc_(): \f$ LOCr(
 *                   m ) = numroc_( m, mb, myrow, irsrc, nprow)\f$ Where \f$myrow\f$ is contains the local
 *                   process (local value).
 * \param[out] info  On output if:
 *                   - \f$info = 0\f$: successful exit;
                     - \f$info < 0\f$: if \f$info = -i\f$ the \em i-th argument had an illegal value.
 *
 * For a complete documentation the reader should refer to \cite ScaLAPACK_webpage \cite SLUG.
 */ 
void descinit_(int *desc, int *m, int *n, int *mb, int *nb, int *irsrc, int *icsrc, int *ictxt, int *lld,
	       int *info);

/**
 * \brief Computes the number of rows or columns of a distributed matrix that is owned by a particular process
 * (ScaLAPACK routine).
 *
 * \param[in] n        (global) The number of rows/columns in the distributed matrix.
 * \param[in] nb       (global) Block size, size of the blocks the distributed matrix is split into.
 * \param[in] iproc    (local) The coordinate of the process whose local array row or column is to be determined.
 * \param[in] isrcproc (global) The coordinate of the process that possesses the first row or column of the
 *                     distributed matrix.
 * \param[in] nprocs   (global) The total number of processes over which the matrix is distributed.
 * 
 * \returns Number of rows or columns of a distributed matrix that is owned by a particular process.
 *
 * For a complete documentation the reader should refer to \cite ScaLAPACK_webpage \cite SLUG.
 */
int numroc_(int *n, int *nb, int *iproc, int *isrcproc, int *nprocs);

/**
 * \brief Returns the local indexes corresponding the distributed matrix given their global coordinates. It
 * also returns globally the coordinates in the grid of the process owning the given matrix entry.coordinates
 * in the grid (local coordinates) of the process owning the given matrix entry (ScaLAPACK routine).

 *
 * \param[in]  grindx (global) The global row starting index of the submatrix.
 * \param[in]  gcindx (global) The global column starting index of the submatrix.
 * \param[in]  desc   (global) The array descriptor for the underlying distributed matrix.
 * \param[in]  nprow  (global) The total number of process rows over which the distributed matrix is
 *                    distributed.
 * \param[in]  npcol  (global) The total number of process columns over which the distributed matrix is
 *                    distributed.
 * \param[in]  myrow  (local) The row coordinate of the process calling this routine.
 * \param[in]  mycol  (local) The column coordinate of the process calling this routine.
 * \param[out] lrindx (local) The local rows starting index of the submatrix.
 * \param[out] lcindx (local) The local columns starting index of the submatrix.
 * \param[out] rsrc   (global) The row coordinate of the process that possesses the first row and column of
 *                    the submatrix.
 * \param[out] csrc   (global) The column coordinate of the process that possesses the first row and column of
 *                    the submatrix.
 *
 * For a complete documentation the reader should refer to \cite ScaLAPACK_webpage \cite SLUG.
 */
void infog2l_( int *grindx, int *gcindx, int *desc, int *nprow, int *npcol, int *myrow, int *mycol, int *lrindx,
	       int *lcindx, int *rsrc, int *csrc );

#endif /* CBLACS_H_ */
