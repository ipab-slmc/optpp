
#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#ifdef WITH_MPI
#include "mpi.h"
#endif

#include <string.h>

#include "pds.h"
#include "common.h"

extern struct pdscon pdscon;

#ifdef WITH_MPI
MPI_Op pdswapOpNum;
#endif

int pdscom(char *emesg)
{
  /*******************************************************************
   *
   * Service subroutine to fire up PVM or MPI, if either is in use,
   * and to initialize the global variables needed to handle
   * communications for the distributed memory versions.
   *
   * Written by Virginia Torczon (with help from R. M. Lewis)
   * MPI version designed and implemented by David Serafini
   *
   * Last modification: February 1995. (MPI version)
   *
   * Parameters
   *
   *    Output
   *
   *       ERROR  Error flag.  This is effectively a no-op for all but
   *              PVM  and MPI, where it is used to flag difficulties
   *              in firing up the parallel virtual machine.
   *
   *******************************************************************/

#ifdef WITH_MPI

  int error, flag, resultlen;
  extern MPI_User_function pdswap;
  char buffer[MPI_MAX_ERROR_STRING];

  /* Setup code for MPI */

  error = MPI_Initialized(&flag);
  if (error != MPI_SUCCESS) goto L99;

  if (flag == 0) {
      printf("\npdscom: MPI has not been initialized.\n");
      strcpy(emesg, "pdscom: MPI has not been initialized");
      return 15;
    }

  error = MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
  if (error != MPI_SUCCESS) goto L99;

  error = MPI_Comm_size(MPI_COMM_WORLD, &pdscon.nproc);
  if (error != MPI_SUCCESS) goto L99;

  error = MPI_Comm_rank(MPI_COMM_WORLD, &pdscon.me);
  if (error != MPI_SUCCESS) goto L99;

  error = MPI_Op_create(&pdswap, 1, &pdswapOpNum);
  if (error != MPI_SUCCESS) goto L99;

  return 0;

L99:
  MPI_Error_string(error, buffer, &resultlen);
  printf("\npdscom: MPI Error - %s\n", buffer);
  strcpy(emesg, "pdscom: MPI initialization failed");
  return 15;

#else

  /* For scalar version, a big no-op.  Nothing can possibly go
   * wrong. */

  pdscon.me = 0;
  pdscon.nproc = 1;

  strcpy(emesg, "");

  return 0;

#endif
}

