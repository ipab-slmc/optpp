
#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef WITH_MPI
#include "mpi.h"
#endif

int pdshrk(int ndim, int beta, int *resize, int *num)
{
  /*******************************************************************
   *
   * If the best vertex is not replaced at the current step, the
   * convergence theorem requires us to shrink the simplex towards the
   * best vertex.  Since the choice of search schemes may allow us to
   * to "look-ahead" several iterations, we would like to know the
   * size of the smallest simplex seen during the step (which, by the
   * way the scheme is generated, is guaranteed to be a contraction
   * towards the best vertex) and use this smallest simplex to start
   * the search at the next iteration so that we avoid as much as
   * possible needless repetition of points already investigated
   * during the current iteration.
   *
   * This is easily accomplished by simply searching the list of
   * scalars necessary to reconstruct the simplex associated with a
   * given point I (contained in the -1 row of SCHEME) for the
   * smallest value.  We also count the number of times this smallest
   * value has been seen to insure that the entire shrink was
   * computed.
   *
   * Written by Virginia Torczon.
   *
   * Last modification:  January 6, 1994.
   *
   * Parameters
   *
   *    Input
   *
   *       N             dimension of the problem to be solved
   *
   *       SSS           size of the search scheme
   *
   *       SCHEME        matrix used to hold the processor's piece of
   *                     the template to be used for the search scheme
   *
   *    Output
   *
   *       RESIZE        the size of the smallest shrink factor found
   *                     in the given search scheme (at least on this
   *                     processor)
   *
   *       NUM           the number of times this shrink factor was
   *                     seen
   *******************************************************************/

#ifdef WITH_MPI
  int my_resize, my_num, resultlen, error=0;
  char buffer[MPI_MAX_ERROR_STRING], emesg[256];

  my_resize = *resize;
  my_num    = *num;

  error = MPI_Allreduce(&my_resize, resize, 1, MPI_INT, MPI_MIN,
			MPI_COMM_WORLD);
  if (error != MPI_SUCCESS) {
    MPI_Error_string(error, buffer, &resultlen);
    printf("\npdshrk:  MPI Error - %s\n", buffer);
    strcpy(emesg, "pdshrk:  error returned by MPI_Allreduce\n");
    return 15;
  }

  if (my_resize != *resize) my_num = 0;

  error = MPI_Allreduce(&my_num, num, 1, MPI_INT, MPI_SUM,
			MPI_COMM_WORLD);
  if (error != MPI_SUCCESS) {
    MPI_Error_string(error, buffer, &resultlen);
    printf("\npdshrk:  MPI Error - %s\n", buffer);
    strcpy(emesg, "pdshrk:  error returned by MPI_Allreduce\n");
    return 15;
  }
#endif

  if (*num < ndim) *resize *= beta;

  return 0;
}

