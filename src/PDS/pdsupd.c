
#include "pds.h"

int pdsupd(long int worked, int best, int n, int *index, double *s,
	   double *plus, double alpha)
{
  /*******************************************************************
   *
   * This is the routine for updating the current simplex to produce
   * the simplex for the next iteration.
   *
   * Written by Virginia Torczon.
   *
   * Last modification:  January 5, 1994.
   *
   * Input
   *
   *    WORKED        flag to signal whether or not the best vertex
   *                  was replaced
   *
   *    BEST          index tracking which vertex in the current
   *                  simplex produced the new best vertex
   *
   *    N             dimension of the problem to be solved 
   *
   *    INDEX         a permutation array used to track the best
   *                  vertex in the simplex
   *
   *    S             an N by N+1 matrix containing the N+1 vertices
   *                  that define the current simplex
   *
   *    PLUS          vector containing the new best iterate
   *
   *    ALPHA         the scaling factor needed to resize and reorient
   *                  the current simplex to produce the correct
   *                  simplex for the next iteration
   *
   *    DEBUG         flag to signal whether or not to dump debugging
   *                  information to a file
   *                     0   no debugging output 
   *                     1   log starting simplex and each iterate
   *                     2   include some intermediate debugging
   *                     3   include almost all intermediate debugging
   *                     4   include all debugging (including SCHEME)
   *
   *    LPR           unit number for the debugging file
   *
   *    FBEST         the function value at the best vertex produced
   *                  by the search
   * Output
   *
   *    S             an N by N+1 matrix containing the N+1 vertices
   *                  that define the simplex for the next iteration
   *
   *******************************************************************/

  /* Local variables */

  static int temp;
  static double a;
  static int i, j, vj;

  /* Function Body */

  if (worked) {

    /* Update the simplex.  Note that we use the edges in the original
     * simplex to perform the update from the new best vertex. */

    for (j = 0; j < best; j++) {
      vj = index[j];

      for (i = 0; i < n; i++) {
	s[i + vj * n] = plus[i+4] + alpha * (s[i + vj * n] -
					     s[i + index[best] *  n]);
      }
    }

    for (j = best + 1; j <= n; j++) {
      vj = index[j];

      for (i = 0; i < n; i++) {
	s[i + vj * n] = plus[i+4] + alpha * (s[i + vj * n] - 
					     s[i + index[best] * n]);
      }
    }

    /* Now, copy the new best vertex into the simplex. */

    for (i = 0; i < n; i++) {
      s[i + index[best] * n] = plus[i+4];
    }

    /* Swap pointers to indicate the new best vertex. */

    temp = index[0];
    index[0] = index[best];
    index[best] = temp;
  }
  else {
    a = 1. - alpha;

    for (j = 0; j < n; j++) {
      vj = index[j+1];

      for (i = 0; i < n; i++) {
	s[i + vj * n] = a * s[i + index[0] * n] + 
	  alpha * s[i + vj * n];
      }
    }
  }

  return 0;
}

