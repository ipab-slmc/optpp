
#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#include <math.h>
#include <string.h>

#ifdef WITH_MPI
#include "mpi.h"
#endif

#include "pds.h"
#include "cblas.h"

int pdsdone(int maxitr, int count, int n, double stol, double *length,
	    double *v0, double *conv, double finit, double fprev,
	    double fbest, double ftol, int max_fevals, int fevals,
	    char *emesg, int trpds)
{
  /*******************************************************************
   *
   * This is a service function used to test for convergence of the
   * parallel direct search methods.
   *
   * This function first tests for convergence by checking whether the
   * number of iterations has reached the maximum set by the user.  If
   * not, it then checks to see whether or not the following step
   * tolerance criterion has been met:
   *
   *      (1 / DELTA) * {    MAX     ||V - V || },
   *                         I,J        I   J  2
   * where
   *      DELTA = MAX(1, ||V || ),
   *                        0  2
   *
   * is less than some user-specified tolerance.
   *
   * The return values are:
   *
   *    TRUE    if either of the convergence criteria has been met
   *    FALSE   otherwise
   *
   * Written by Virginia Torczon.
   *
   * Last modification:  January 4, 1994.
   *
   *    Input
   *
   *       MAXITR        maximum number of iterations allowed (as
   *                     specified by the user during initialization)
   *
   *       COUNT         number of iterations completed
   *
   *       N             dimension of the problem to be solved
   *
   *       TOL           step tolerance criterion for completing the
   *                     search using the test described above (also
   *                     specified by the user during initialization)
   *
   *       LENGTH        length of the longest edge in the current
   *                     simplex
   * 
   *       V0            current best vertex
   *
   *    Output
   *
   *       CONV          computed convergence value (-TOL if MAXITR
   *                     exceeded)
   *
   *******************************************************************/

  /* System generated locals */

  int ret_val, total_fevals;
  int incx = 1;

  /* Local variables */

  static double norm;
  static double delta, deltaf, rftol;

#ifdef WITH_MPI

  char buffer[MPI_MAX_ERROR_STRING];
  int error, resultlen;

#endif

  /* Function Body */

  ret_val = 0;
  strcpy(emesg, "");

#ifdef WITH_MPI

  error = MPI_Allreduce(&fevals, &total_fevals, 1, MPI_INT, MPI_SUM,
			MPI_COMM_WORLD);

  if (error != MPI_SUCCESS)
    {
      MPI_Error_string(error, buffer, &resultlen);
      printf("\npdsdon: MPI Error - %s\n", buffer);
      strcpy(emesg, "pdsdone: error returned by MPI_Allreduce\n");
      return 15;
    }

#else

  total_fevals = fevals;

#endif

  /* check for maximum number of iterations, maximum number of
   * function evaluations, and minimum step size */

  if (count >= maxitr) {
    ret_val = 3;
    *conv = -stol;
    strcpy(emesg, "pdsdone: Maximum number of iterations exceeded");
  }
  else if (total_fevals >= max_fevals) {
    ret_val = 4;
    *conv = -stol;
    strcpy(emesg, "pdsdone: Maximum number of function evaluations exeeded");
  }
  else {
    norm = dnrm2(&n, v0, &incx);
    delta = max(1.,norm);
    *conv = *length / delta;
    ret_val = *conv <= stol;
    strcpy(emesg, "pdsdone: Step tolerance passed");
  }

  /* check for function decrease - If TRPDS, compare to initial
   * function value.  Otherwise, compare to previous function
   * value. */

  if (trpds) {

    if ((fbest < 0.0) && (finit < 0.0))
      ftol = 2-ftol;

    if (fbest <= (ftol*finit)) {
      ret_val = 2;
      strcpy(emesg, "pdsdone: Function tolerance passed");
    }
  }
  else {
    deltaf = fprev - fbest;
    rftol = ftol*max(1.0,fabs(fbest));

    if (deltaf <= rftol) {
      strcpy(emesg,"pdsdone: Function tolerance test passed");
      ret_val = 2;
    }
  }

  return ret_val;
}

