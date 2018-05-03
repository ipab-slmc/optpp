
#include <string.h>

#include "mpi.h"

#include "pds.h"
#include "cblas.h"
#include "common.h"

extern struct pdscon pdscon;
extern MPI_Op pdswapOpNum;

int pdsgop(double *mybest, int sizee, double *yourbest, char *emesg)
{
  /*******************************************************************
   *
   * This subroutine is used for doing global communication operations
   * with a user-defined function (in this case PDSWAP) to handle the
   * global decision-making.
   *
   * In the MPI version, PDSWAP must be a MPI_OP, which is a value
   * returned by the MPI_CREATE_OP() routine.
   *
   * Intel makes this process easy by providing the communication
   * library routine GOPF.
   *
   * In the case of PVM, we can make no assumptions about the
   * interconnection between processors, so we resort to the simplest
   * possible exchange of information: "my" processor broadcasts the
   * best point it has seen to all other processors and then sits in a
   * loop and waits to receive the same information from all other
   * processors.  The call to PDSWAP ensures that if a point with a
   * lower function value (or a point with the same function value but
   * from a processor with a lower number, i.e., earlier in the search
   * scheme generation strategy) is received, it is substituted for
   * "my" best point.
   *
   * Written by Virginia Torczon.
   * MPI version by David Serafini.
   *
   * Last modification:  May, 1995.
   *
   * Parameters
   *
   *    Input
   *
   *       MYBEST        Vector containing the point with the lowest
   *                     function value my node has yet seen.  Note
   *                     that the vertex is contained in positions 1
   *                     through N; its function value is contained
   *                     in position N+1; positions -3 through 0 are
   *                     used to pass additional information:
   *                       (-3) = N;            (-2) = processor#;
   *                       (-1) = scale factor;  (0) = BEST value
   *
   *       SIZEE         Number of entries in MYBEST.
   *
   *       YOURBEST      A work vector, which must be of the same
   *                     length as MYBEST, that is used to receive
   *                     the contributions from other nodes.
   *
   *       PDSWAP        Function to evaluate the contributions from
   *                     other processors, it effects the replacement
   *                     if another processor has found a point with
   *                     a lower function value.  In the MPI version,
   *                     this is an integer handle for a function
   *                     returned by MPI_OP_CREATE().
   *    Output
   *
   *       MYBEST        Vector containing the point with the lowest
   *                     function value my node has yet seen.
   *
   *******************************************************************/

  int resultlen;
  int error;
  int incx = 1, incy = 1;
  char buffer[MPI_MAX_ERROR_STRING];

  /* Execute the global exchange of information by using the MPI
   * global reduce with the user-defined reduction operation PDSWAP,
   * which allows for comparison on the function value only, but swaps
   * the vertex, the function value, the scalar and pointer necessary
   * to update the simplex, and a flag for the processor which
   * produced the point.  The message size for MPI functions is
   * defined by the number of elements.  The .TRUE. indicates PDSWAP
   * is commutative, and MPIOPNUM is a flag value that MPI_OP_CREATE()
   * assigned to the PDSWAP function. */

  error =  MPI_Allreduce(mybest, yourbest, sizee,
			 MPI_DOUBLE, pdswapOpNum,
			 MPI_COMM_WORLD);

  if (error != MPI_SUCCESS)
    {
      MPI_Error_string(error, buffer, &resultlen);
      printf("\npdsgop: MPI Error - %s\n", buffer);
      strcpy(emesg, "pdsgop: error returned by MPI_Allreduce\n");
      return 15;
    }

  /* The result is left in YOURBEST, so copy it to MYBEST */

  dcopy(&sizee, yourbest, &incx, mybest, &incy);

  return 0;
}
