
#include <math.h>

#include "pds.h"

int dqrdc(double *x, int ldx, int ndim, int p, double *qraux,
	  int *jpvt, double *work, int job);

int pdsdgn(int ndim, double *s, double *a, double *work1,
	   double *qraux, int *jpvt, double *rcond)
{
  /*******************************************************************
   *
   * Parameters
   *
   *    Input
   *
   *       N             dimension of the problem to be solved
   *
   *       S             two-dimensional array containing the vertex
   *                     specified by the user to begin the search
   *                     (as well as the remaining N vertices of the
   *                     initial simplex if they have already been
   *                     entered by the user) 
   *
   *       A             workspace vector of length N*N to hold the 
   *                     adjacency matrix, if needed
   *
   *       WORK1         workspace vector of length N used as
   *                     auxiliary storage during the QR
   *                     decomposition
   *
   *       QRAUX         workspace vector of length N used as
   *                     auxiliary storage during the QR
   *                     decomposition
   *
   *       JPVT          workspace vector of length N to control
   *                     pivoting during the QR decomposition
   *    Output
   *
   *       RCOND         an estimate of the inverse of the condition
   *                     number of the adjacency matrix.  This should
   *                     give some warning as to whether or not the
   *                     initial simplex is degenerate. 
   *
   *******************************************************************/

  /* Local variables */

  static int i, j, k;

  /* Form the matrix of edges adjacent to the initial vertex. */

  static int one = 1;
    
  for (j = 0; j < ndim; j++) {
    k = j + 1;

    for (i = 0; i < ndim; i++) {
      a[i + j * ndim] = s[i + k * ndim] - s[i];
    }

    jpvt[j] = 0;
  }

  /* Call QR with pivoting. */

  dqrdc(a, ndim, ndim, ndim, qraux, jpvt, work1, one);
  *rcond = fabs(a[ndim - 1 + (ndim-1)*ndim]) / fabs(a[0]);

  return 0;
}

