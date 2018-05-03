
#include "pds.h"
#include "common.h"

extern struct pdscon pdscon;

int pdsglb(int ndim, double *mybest, double *yourbest, char *emesg)
{
  /*******************************************************************
   *
   * Unique to distributed-memory machines.
   *
   * This is a (the!) routine used to initiate the global
   * communication needed to run the parallel direct search schemes on
   * distributed- memory machines.
   *
   * Written by Virginia Torczon.  C C
   * Last modification: January 20, 1995.
   *
   * Input
   *
   *    N             dimension of the problem to be solved
   *
   *    MYBEST        vector containing the vertex computed by my node
   *                  during this iteration with the lowest function
   *                  value.  Note that the vertex is contained in
   *                  positions 1 through n; it`s function value is
   *                  contained in position n+1; positions -3 through
   *                  0 are used to pass additional information.
   *
   *    YOURBEST      a work vector, which must be of the same length
   *                  as MYBEST, that is used to receive the 
   *                  contributions from other nodes.
   *
   *    DEBUG         flag to signal debugging output
   *
   *    LPR           unit number for debugging output
   *
   * Output
   *
   *    MYBEST        vector containing the vertex with the least
   *                  function value to be found across all processors.
   *                  Note that it also contains the scalar ALPHA and
   *                  the pointer BEST needed to reconstruct the 
   *                  associated simplex, as well as the actual
   *                  function value.
   *
   *******************************************************************/

  int error, sizee;

  /* Initialize the message buffer with the dimension of the problem
   * and the number of "my" node---information needed during the global
   * exchange. */

  mybest[0] = ndim;
  mybest[1] = pdscon.me;

  sizee = ndim + 5;

  error = pdsgop(mybest, sizee, yourbest, emesg);

  return(error);
}
