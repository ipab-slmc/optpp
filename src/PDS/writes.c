/*
//--------------------------------------------------------------------
// Copyright (C) 1993,1994: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//--------------------------------------------------------------------
*/

#include <unistd.h>

#include "pds.h"

#define A 	-1
#define BEST 	0

int writes(int out, int ndim, int unique, int factor, int beta, 
	   int *scheme, int *list)
{
  /*******************************************************************
   *
   * Write the points in the search scheme out to an unformatted
   * fortran file for later use by the parallel direct search methods.
   * note that before the points in the search strategy are written
   * out to a file, we write out four pieces of "header" information:
   * n, unique, factor, and beta.  the first two pieces are for error
   * checking when the file is later used by the parallel direct
   * search methods; the second two pieces of information are used to
   * accelerate the contraction step when the best vertex is not
   * replaced.
   *  
   * Arguments:
   *
   *    out           file descriptor to which the scheme is to be
   *                  written
   *
   *    n             dimension used to generate the search scheme
   *
   *    total         the total number of points created (based on the
   *                  amount of workspace allocated in the main
   *                  program)
   * 
   *    unique        the total number of unique points created after
   *                  checking for duplicate n-tuples
   *
   *    factor        the scaling factor used to generate the search
   *                  scheme; note that this scaling factor is
   *                  necessary to reconstruct the real (or double
   *                  precision) counterparts to the integer tuples
   *
   *    beta          the contraction factor used to generate the
   *                  points in the search scheme (information
   *                  required by the optimization to accelerate
   *                  contractions when the best vertex is not
   *                  replaced)
   *
   *    scheme        the integer tuples used to define each point in
   *                  the search scheme
   *
   *    list          index array used to point to the unique n-tuples
   *                  in the search scheme; the ones actually to be
   *                  written out to the file
   *
   * Note that on the intel distributed memory machines it is *much*
   * faster to aggregate data on the node and write it out to the file
   * as a single block than to write out a piece at a time when using
   * the intel library routines for i/0.
   *
   *******************************************************************/

  /* System generated locals. */

    int scheme_dim1, scheme_offset;

    /* Local variables. */

    int i, j;
    int error = 0;

    /* Parameter adjustments. */

    scheme_dim1   = ndim + 2;
    scheme_offset = scheme_dim1 * (-ndim) - 1;
    scheme       -= scheme_offset;
    --list;

    /* Write out header information. */

    error = write (out, WRITE_TYPE &ndim,   INT_SIZE);
    error = write (out, WRITE_TYPE &unique, INT_SIZE);
    error = write (out, WRITE_TYPE &factor, INT_SIZE);
    error = write (out, WRITE_TYPE &beta,   INT_SIZE);

    /* Write out the SCHEME data. */
    
    for (i = 1; i <= unique; ++i) {
      error = write(out, WRITE_TYPE &scheme[list[i] * scheme_dim1 - 1],4);
      error = write(out, WRITE_TYPE &scheme[list[i] * scheme_dim1],4);

      for (j = 1; j <= ndim; ++j) {
	error = write(out,WRITE_TYPE &scheme[j + list[i] * scheme_dim1],4);
      }
    }

    if (error == -1) 
      return error;
    else
      return 0;
}

