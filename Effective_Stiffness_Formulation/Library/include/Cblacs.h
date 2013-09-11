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

/**
 * \brief This routine holds up execution of all processess within the indicated scope until they all have
 * called the routine (BLACS routine).
 *
 * \param[in] icontxt Integer handle indicating the BLACS context.
 * \param[in] scope   Indicates whether a process row (\c scope = 'R'), column (\c scope = 'C'), or entire
 *                    grid (\c scope = 'A') will participate in the barrier.
 * 
 * For a complete documentation the reader should refer to \cite{BLACS_webpage}.
 */
void Cblacs_barrier(int icontxt, char *scope);

/**
 * \brief Gets the values the BLACS are using for internal defaults (BLACS routine).
 *
 * Gets the values the BLACS are using for internal defaults (BLACS routine). The most common use is in
 * retrieving a default system context for input into CBlacs_gridinit() or CBlabs_gridmap(). The contents of
 * the output variable \c val will depend on \c what:
 *
 * - \c what = 0: Handle indicating the system context;
 * - \c what = 1: The BLACS message ID range;
 * - \c what = 2: The BLACS debug level the library was compiled with;
 * - \c what = 10: Handle indicating the system context used to define the BLACS context whose handle is \c
 *      icontxt.
 * - \c what = 11: Number of rings multiring topology is presently using;
 * - \c what = 12: Number of branches general tree topology is presently using
 *
 * \param[in]  icontxt Integer handle indicating the BLACS context if a value of \c what is tied to a
 *                     particular context. It is ignored otherwise. 
 * \param[in]  what    The BLACS internal that should be returned in \c val.
 * \param[out] val     The value of the BLACS internal.
 * 
 * For a complete documentation the reader should refer to \cite{BLACS_webpage}.
 */
void Cblacs_get( int icontxt, int what, int *val );

/**
 * \brief Frees a context and liberates the used resources. After beeing freed it can be re-used if new
 * contexts are defined (BLACS routine).
 *
 * \param[in] icontxt Integer handle indicating the BLACS context. 
 *
 * For a complete documentation the reader should refer to \cite{BLACS_webpage}.
 */
int Cblacs_gridexit( int icontxt );

/**
 * \brief This routine takes the available processes and assign, or map, them into a BLACS process grid (BLACS
 * routine).
 *
 * This routine takes the available processes and assign, or map, them into a BLACS process grid (BLACS
 * routine). In other words, they establish how the BLACS coordinate system will map into the native machine's
 * process numbering system. Each BLACS grid is contained in a context (its own message passing universe), so
 * that it does not interfere with distributed operations which occur within other grids/contexts. These grid
 * creation routines may be called repeatedly in order to define additional contexts/grids. All BLACS codes
 * must call this routine, or its sister routine CBlacs_gridmap().
 *
 * This routine creates a simple NPROW x NPCOL process grid. This process grid will use the first NPROW
 * NPCOL processes, and assign them to the grid in a row- or column-major natural ordering. If these
 * process-to-grid mappings are unacceptable, BLACS_GRIDINIT's more complex sister routine Cblacs_gridmap()
 * must be called instead.
 *
 * Depending on the value of \c order, the processes are mapped to the BLACS grid using following ordering:
 * - \c order = 'R': row-major natural ordering;
 * - \c order = 'C': column-major natural ordering;
 * - Any other value will result in row-major natural ordering.
 *
 * \param[in,out] icontxt On input an integer handle indicating the system context to be used in creating the
 *                        BLACS context. On output, the integer handle to the created BLACS context. 
 * \param[in]     order   Indicates how to map processes to BLACS grid.
 * \param[in]     nprow   Indicates how many process rows the process grid should contain.
 * \param[out]    npcol   Indicates how many process columns the processs grid should contain.
 *
 * For a complete documentation the reader should refer to \cite{BLACS_webpage}.
 */
void Cblacs_gridinit(int *icontxt, char *order, int nprow, int npcol);

/**
 * \brief Returns information on the current grid (BLACS routine). All values are returned as -1 if the
 * given context is not valid.
 *
 * \param[in]  icontxt Integer handle indicating the BLACS context. 
 * \param[out] nprow   It contains the number of process rows in the current process grid.
 * \param[out] npcol   It contains the number of process columns in the current proces grid.
 * \param[out] myrow   It contains the calling process's row coordinate in the process grid.
 * \param[out] mycol   It contains the calling process's column coordinate in the process grid.
 *
 * For a complete documentation the reader should refer to \cite{BLACS_webpage}.
 */
void Cblacs_gridinfo( int icontxt, int *nprow, int *npcol, int *myrow, int *mycol );

/**
 * \brief This routine takes the available processes and assign, or map, them into a BLACS process grid in an
 * arbitrary manner (BLACS routine).
 *
 * This routine takes the available processes, and assign, or map, them into a BLACS process grid in an
 * arbitrary manner. In other words, they establish how the BLACS coordinate system will map into the native
 * machine's process numbering system. Each BLACS grid is contained in a context (its own message passing
 * universe), so that it does not interfere with distributed operations which occur within other
 * grids/contexts. These grid creation routines may be called repeatedly in order to define additional
 * contexts/grids. All BLACS codes must call this routine, or its sister routine Cblacs_gridinit().
 *
 * \param[in,out] icontxt On input an integer handle indicating the system context to be used in creating the
 *                        BLACS context. On output, the integer handle to the created BLACS context. 
 * \param[in]     usermap Array of dimension (\c ldumap, \c npcol) indicating the process-to-grid mapping. \c
 *                        usermap(\em i,\em j) holds the process number of the process to be placed in \f${i,
 *                        j}\f$ of the process grid.
 * \param[in]     ldumap  The leading dimension of the 2D array \c usermap.
 * \param[in]     nprow   Indicates how many process rows the process grid should contain.
 * \param[out]    npcol   Indicates how many process columns the processs grid should contain.
 *
 * For a complete documentation the reader should refer to \cite{BLACS_webpage}.
 */
void Cblacs_gridmap( int *icontxt, int *usermap, int ldumap, int nprow, int npcol );

/**
 * \brief Returns the system process number of the process at the given coordinates in the process grid (BLACS
 * routine).
 *
 * \param[in] icontxt Integer handle indicating the BLACS context.
 * \param[in] prow    The row coordinate of the process who's system process number is to be determined.
 * \param[in] pcol    The column coordinate of the process who's sytem process number is to be determined.
 *
 * \return The sustem process number of the process at \f$(prow, pcol)\f$ in the process grid.
 *
 * For a complete documentation the reader should refer to \cite{BLACS_webpage}.
 */
int Cblacs_pnum( int icontxt, int prow, int pcol );

/**
 * \brief Translates a MPI communicator into a BLACS handle.
 * 
 * This function returns an integer handle which can be passed into the grid creation routines to
 * indicate the MPI communicator (system context) (BLACS routine).
 *
 * \param[in] SysCtxt The system context to be mapped to an integer BLACS handle. 
 *
 * For a complete documentation the reader should refer to \cite{MPI_BLACS_Issues}.
 */
int Csys2blacs_handle(MPI_Comm SysCtxt);

/**
 * \brief Free a BLACS handle associated with a MPI communicator.
 * 
 * This routine may be called to dissassociate a BlacsHandle with a system context. There is no need to
 * do this, other than to keep memory usage down (BLACS routine).
 *
 * \param[in] BlacsHandle Integer handle who's mapping is no longer needed. 
 *
 * For a complete documentation the reader should refer to \cite{BLACS_webpage}.
 */
int Cfree_blacs_system_handle(int BlacsHandle);

#endif /* CBLACS_H_ */
