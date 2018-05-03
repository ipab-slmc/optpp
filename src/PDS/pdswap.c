
#include <math.h>

#include "mpi.h"

#include "pds.h"
#include "cblas.h"
#include "common.h"

extern struct pdscon pdscon;

void pdswap(double *yourbest, double *mybest, int *len,
	    MPI_Datatype *type)
{
  /*******************************************************************
   *
   * Unique to distributed-memory machines.
   *
   * This is the function used to allow us to compare the two vectors
   * MYBEST and YOURBEST by looking only at the function value stored
   * in position N+1, but it is also used to allow us to exchange the
   * entire vector from positions -2 to N+1 if we have received a
   * vertex with a lower function value.
   *
   * Written by Virginia Torczon.
   * MPI version written by David Serafini.
   *
   * Last modification:  February, 1995. (MPI version)
   *
   * Input
   *
   *    MYBEST          Vector containing the point with the lowest
   *                    function value my node has yet seen.  Note
   *                    that the vertex is contained in positions 1
   *                    through N; its function value is contained in
   *                    position N+1; positions -3 through 0 are used
   *                    to pass additional information.
   *
   *    YOURBEST        A work vector, which must be of the same
   *                    length as MYBEST, that is used to receive the
   *                    contributions from other nodes.
   *
   *    LEN (MPI only)  Length of argument vectors
   *
   *    TYPE (MPI only) MPI datatype flag of argument vectors.  (both
   *                    of these can be ignored because they are
   *                    known)
   * Output
   *
   *    MYBEST          Vector containing the point with the lowest
   *                    function value my node has yet seen.
   *
   *******************************************************************/

  int error, n;
  int incx = 1, incy = 1;

  n = floor(mybest[0]);

  /* MPI_ALLREDUCE() may break the vector into pieces, which won`t
   * work with this routine, so make sure it isn`t happening. */

  if (*len != n+5)
    error = MPI_Abort(MPI_COMM_WORLD, 11);

  /* Copy YOURBEST to MYBEST if it is strictly less than MYBEST or if
   * they`re equal and YOURBEST has a lower processor number */

  if ((yourbest[n+4] < mybest[n+4]) || ((yourbest[n+4] == mybest[n+4])
				    && (yourbest[1] < mybest[1]))) {
    dcopy(len, yourbest, &incx, mybest, &incy);
  }
}
