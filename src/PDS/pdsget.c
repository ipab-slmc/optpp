
#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>

#ifdef WITH_MPI
#include "mpi.h"
#endif

#include "pds.h"
#include "common.h"

#ifndef READ_TYPE
#define READ_TYPE (void *)
#endif

extern struct pdscon pdscon;

int pdsget(int ndim, FILE *fpin, int *sss, double *factor, int *beta,
	   char *emesg)
{
  /*******************************************************************
   *
   * This is the subroutine used to read in the search scheme for the
   * parallel direct search methods.
   *
   * Written by Virginia Torczon. 
   * MPI version written by David Serafini. 
   *
   * Last modification:  February, 1995. 
   *
   * Parameters
   *
   *    Input 
   *
   *       N             dimension of the problem to be solved 
   *
   *       IN            unit number from which the scheme is to be
   *                     read
   *
   *       SSS           size of the search scheme (the number of
   *                     points to be considered at each iteration)
   *
   *    Output
   *
   *       SSS           on a parallel computer, this contains the
   *                     number of points "my" node is to compute
   *                     before engaging in the global
   *                     communication/synchronization conducted at
   *                     each iteration
   *
   *       SCHEME        the int tuples used to define each point in
   *                     the search scheme
   *
   *       FACTOR        the scaling factor used to generate the
   *                     search scheme; note that this scaling factor
   *                     is necessary to reconstruct the real (or
   *                     double precision) counterparts to the int
   *                     tuples
   *
   *       RESIZE        the size of the smallest (complete) shrink
   *                     step seen during a single iteration when
   *                     using a search scheme of this size (this
   *                     factor is used to accelerate the shrinking
   *                     process if no improvement is seen during the
   *                     course of a single iteration)
   *
   *       ERROR         flag to signal whether or not the scheme file
   *                     was incompatible with the problem specified
   *                     by the user.  Return codes are: 
   *                        PDSAOK   no incompatibility detected 
   *                        PDSNNN   dimension mismatch 
   *                        PDSSGT   the size of the search scheme
   *                                  requested is greater than the
   *                                  total number of points in the
   *                                  file SCHEME
   *
   *******************************************************************/

  /* Variables */

  int error;
  int header[4];

  error = 0;

  /* Read in 4 bytes */

  fread(READ_TYPE &header[0], INT_SIZE, 4, fpin);

  /* HEADER should hold:
   *    N, 
   *    TOTAL (number of points in the search file),
   *    FACTOR (in integer form), and
   *    BETA. */
    
  /* Ensure the search scheme is for a problem of the same dimension. */

  if (header[0] == ndim) {

    /* Make sure the size of the search scheme requested does not
     * exceed the total number of points generated.  */

    if (*sss <= header[1]) {

      /* Convert the factor used to construct the search scheme from
       * its integer representation to a floating point number. */

      *factor = (double) header[2];

      /* Assign the contraction factor used to construct the search
       * scheme.  */

      *beta = header[3];

      /* In the shared-memory world, this is just a simple read. */

      *sss = ceil((double) *sss/pdscon.nproc);

    }
    else {
      error = 10;            /* sss is too big */
      strcpy(emesg, "pdsget: not enough points in search scheme");
    }
  }
  else {
    error = 11;              /* ndim is wrong */
    strcpy(emesg, "pdsget: search scheme of wrong dimension");
  }

  return error;
}

