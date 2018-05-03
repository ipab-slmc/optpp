
#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#if (defined(__sgi) || defined(__xlc__) || defined(__xlC__))
#define WANT_MATH
#else
#define WANT_STREAM
#define WANT_MATH
#endif

#ifdef HAVE_STD
#include <cmath>
#include <cstdio>
#else
#include <math.h>
#include <stdio.h>
#endif

#include "OptPDS.h"
#include "pds.h"
#include "common.h"
#include "ioformat.h"

using namespace std;
using NEWMAT::ColumnVector;

extern struct pdscon pdscon;
extern struct conbcmni conbcmni;

namespace OPTPP {

void pdsquit(int, ofstream *, int *, double, int *, int, double);

int pdswork(NLP0* nlp, ostream *fout, ofstream *fplpr, int debug,
	    double tol, int maxitr, int sss, int *flag, double factor,
	    int beta, double *simplex, double *vscales,
	    int *pds_index, double *fbest, double *length, int *count,
	    int type, double scale, double *rcond, char *emesg,
	    double fcn_tol, double tr_size, int max_fevals,
	    int first, int trpds, double feas_tol, FILE *scheme_in)
{
  /*******************************************************************
   *
   * This is the controlling routine for the actual optimization.  For
   * a further description of the parallel direct search methods see
   * J. E. Dennis, Jr. and Virginia Torczon, "Direct Search Methods on
   * Parallel Machines," SIAM J. Optimization, Vol. 1, No. 4,
   * pp. 448--474, November 1991.
   *
   * Written by Virginia Torczon. 
   * MPI and PDSCONB/I versions written by David Serafini
   *
   * Last modification: 9Sept96, Added VSCALES array 
   *                    31May96, Added FLAG<>0 handling, termination
   *                             msg.
   *                    19Apr96, Added PDSLOG() to DEBUG=-1 handling.
   *
   * Input
   *
   *    N             dimension of the problem to be solved 
   *
   *    LPR           unit number for the debugging file(s) assumed
   *                  open for FORMATTED output if DEBUG.NE.-2 assumed
   *     	      open for UNFORMATTED output if DEBUG.EQ.-2
   *
   *    DEBUG         flag to signal whether or not to dump debugging
   *                  information to a file 
   *                     -3   dump timing info (in .01 seconds) [IBM
   *                          only]
   *                     -2   echo all function values and vertices to
   *                          debug unit
   *                     -1   interactive debugging output by calling
   *                          PDSLOG() (default PDSLOG() prints
   *                          (iter,objval,col#,resid,etc)
   *                      0   no debugging output 
   *                      1   log starting simplex and each iterate
   *                          (V0)
   *                      2   include some intermediate debugging 
   *                      3   include almost all intermediate
   *                          debugging 
   *                      4   include all debugging (including SCHEME)
   *
   *    TOL           tolerance used to check for convergence (this is
   *                  essentially a step tolerance test)
   *
   *    MAXITR        maximum number of iterations allowed
   *
   *    SSS           size of the search strategy (the number of
   *                  points to be considered at each iteration)
   *
   *    FCN           external subroutine to evaluate the objective
   *                  function
   *
   *    CON           external subroutine to evaluate the constraints
   *
   *    FLAG          used to pass back information about possible
   *                  failure in FCN to compute a function value for a
   *                  given point.  FLAG<>0 causes PDS to terminate.
   *
   *    FACTOR        the scaling factor used to reconstruct the real
   *                  (or double precision) counterparts to the int
   *                  tuples
   *
   *    SCHEME        the int tuples used to define each point in the
   *                  search strategy
   *
   *    RESIZE        the size of the smallest (complete) shrink step
   *                  seen during a single iteration when using a
   *                  search strategy of this size (this factor is
   *                  used to accelerate the shrinking process if no
   *                  improvement is seen during the course of a
   *                  single iteration)
   *
   *    S             an N by N+1 matrix containing the N+1 vertices
   *                  that define the simplex to be used to start the
   *                  search
   *
   *    VSCALES       vector of scale factors to multiply to the
   *                  vertices
   *
   *    PDS_INDEX     a permutation array used to track the best
   *                  vertex in the simplex
   *                  (i.e., V0 = S(*,PDS_INDEX(1)))
   *
   *    FBEST         the function value at the best vertex in the
   *                  initial simplex
   *
   *    LENGTH        length of the longest edge in the initial
   *                  simplex
   *
   *    UPARMDIM      sizes and dimensions of U*PARM() vectors
   *
   * Workspace
   *
   *    EDGE          workspace used to compute the length of any edge
   *                  in the simplex
   *
   *    C             workspace used both to compute the current point
   *                  in the list on each processor and, on the
   *                  distributed-memory machines, to provide
   *                  workspace for the global communication call
   *
   *    CS            workspace used for the current point * scale
   *                  factors
   *
   *    PLUS          workspace used to contain the best point seen
   *                  during the current iteration and, on the
   *                  distributed-memory machines, used as the message
   *                  buffer for the global communication call
   *
   *    UIPARM        for any int parameters required by the user to
   *                  evaluate the function
   *
   *    URPARM        for any real parameters required by the user to
   *                  evaluate the function
   * Output
   *
   *    S             the best vertex is contained in
   *                  S(*,PDS_INDEX(1))
   *
   *    PDS_INDEX     permutation array pointing to the current best
   *                  vertex
   *
   *    FBEST         the function value at the best vertex
   *
   *    LENGTH        length of the longest edge in the final simplex
   *
   *    COUNT(3)      number of iterations, function evals completed,   *                  constraint evals completed
   *
   *******************************************************************/

  /* Variables */

  double alpha, r, dist, fprev;
  int ndim = nlp->getDim();
  int point_count, feasible, file_end;
  ColumnVector x_curr = nlp->getXc();
  double finit = nlp->getF();
  int best, error, i, j, k, scheme_dim1, v0;
  bool worked, converged;
  int  done_code, num, resize = 999;
  ColumnVector x(ndim);

  scheme_dim1 = ndim + 2;

  /* Allocate work vectors */

  int *search_dir = new int[scheme_dim1];

  ColumnVector edge(ndim*ndim);
  ColumnVector c(ndim+5);
  ColumnVector cs(ndim+5);
  ColumnVector plus(ndim+5);

  /* Finish initializing simplex (if necessary), find the length of
   * the longest edge in the simplex, and determine the best vertex
   * and its function value.  */

  error = pdsinit(nlp, fout, debug, type, flag, count, scale,
		  simplex, vscales, length, pds_index, fbest, rcond,
		  edge.Store(), c.Store(), plus.Store(), emesg,
		  tr_size, first, trpds, feas_tol);

  if (*flag != 0)
      return(error);
	
  if (error != 0)
      return(error);
	
  v0 = pds_index[0];
  plus(3) = factor;
  plus(4) = 0.;
  for (i = 4; i <= ndim+3; i++)
    plus(i+1) = simplex[i - 4 + v0 * ndim];

  /* Print out iteration summary. */

  (*fout) << "\n\t\t PDS Iteration Summary\n";
  (*fout) << "\n  Iter    F(x)       ||edge||    step   fevals    jbest\n\n";
  (*fout) << d(count[0],5) << e(*fbest,13,4) << e(*length,13,4)
	  << "     " << d(count[1],10) << d(v0,8)<< endl;

  /* Main Loop */

  converged = false;
  while (! converged ) {

    /* Print status of previous (completed) iteration */

    if (debug) {
      (*fplpr) << " PDSWORK      ITERATION, FUNCTION, CONSTRAINT";
      (*fplpr) << "COUNT = " << d(count[0],4) << d(count[1],4)
	       << d(count[2],4) << "\n";
      (*fplpr) << " PDSWORK \n";
    }

    count[0]++;

    /* Calculate "my" vertices. */

    v0 = pds_index[0];
    best = 0;
    plus(ndim + 5) = *fbest;
    point_count = 0;
    file_end = 0;

    for (k=0; ((k<sss) && (!file_end)); k++) {

      feasible = 0;
      while ((!feasible) && (!file_end)) {
	if (debug) {
	  (*fplpr) << " PDSWORK \n";
	  (*fplpr) << " PDSWORK      VERTEX " << d(k,4) << "\n";
	}

	for (i = 0; i < ndim; i++)
	  c(i+5) = simplex[i + v0 * ndim];

	fseek(scheme_in, (4 + (pdscon.me +
			  point_count*pdscon.nproc)*scheme_dim1)*INT_SIZE,
	      SEEK_SET);
	fread(READ_TYPE search_dir, INT_SIZE, scheme_dim1, scheme_in);
	file_end = feof(scheme_in);

	if (!file_end) {
	  point_count++;

	  if (abs(search_dir[0]) < resize) {
	    resize = abs(search_dir[0]);
	    num = 1;
	  }
	  else if (abs(search_dir[0]) == resize)
	    num++;

	  for (j = 0; j < ndim; j++) {

	    for (i = 0; i < ndim; i++) {
	      edge(i+1) = simplex[i + pds_index[j + 1] * ndim] -
		simplex[i + v0 * ndim];
	    }

	    /* JCM: Add 2 to j pds_index to account for funny pds_indexing
	     * in Fortran code */

	    //	alpha = (double) scheme[j + 2 + k * scheme_dim1] / factor;
	    alpha = (double) search_dir[j+2]/factor;

	    for (i = 0; i < ndim; i++)
	      c(i+5) += alpha * edge(i+1);
	  }

	  /* Scale the current vertex to the values needed by FCN/CON */

	  for (i = 0; i < ndim; i++)
	    cs(i+1) = c(i+5) * vscales[i];

	  if (debug) {
	    for (i = 0; i < ndim; i++)
	      (*fplpr) << " PDSWORK      " << e(cs(i+1),30,14) << "\n";
	  }

	  /* Calculate the current vertex`s function value.  Compute
	   * the function only if the point is feasible */

	  feasible = pdschk(nlp, ndim, x_curr.Store(), cs.Store(), tr_size, &dist, trpds,
			    feas_tol);

	  /* A constraint eval is done except when the variable
	   * bounds are violated.*/

	  count[2]++;
	}
      }

      if (!file_end) {
	if (*flag != 0)
	  pdsquit(debug, fplpr, count, r, flag, maxitr, tol);

	for (i = 0; i < ndim; i++) {
	  x(i+1) = cs(i+1);
	}

	c(ndim+5) = nlp->evalF(x);
	count[1]++;

	if (*flag != 0)
	  pdsquit(debug, fplpr, count, r, flag, maxitr, tol);

	if (debug) {
	  (*fplpr) << " PDSWORK      WITH FUNCTION VALUE:"
		   << e(c(ndim+5),30,14) << "\n";
	}

	/* Synchronization/Communication point.
	 * Determine who has the new best vertex. [NOTE: this
	 * should be double-buffered, swapping a pointer instead
	 * of the whole vector. -dbs] */

	if (c(ndim + 5) < plus(ndim + 5)) {
	  best = k;
	  //	  plus[2] = (double)scheme[best*scheme_dim1];
	  //	  plus[3] = (double)scheme[best*scheme_dim1+1];
	  plus(3) = (double) search_dir[0];
	  plus(4) = (double) search_dir[1];

	  for (i = 4; i <= ndim+4; i++)
	    plus(i+1) = c(i+1);
	}
      }
      /*    }
	    else { */

      /* If the constraint violated is not bounds, then an
       * ineq constraint evaluation must have been done. */

    /*      if (debug) {
	(*fout) << " PDSWORK      VERTEX " << k
		<< " violates trust region constraint \n"
		<< " tr_size = " << tr_size
		<< " dist = " << dist << endl;

	for (i = 0; i < ndim; i++)
	  (*fout) << " PDSWORK      " << cs[i] << "\n";
	  }
	  }*/
    }

#ifdef WITH_MPI

    error = pdsglb(ndim, plus.Store(), c.Store(), emesg);
    if (error != 0) return error;

#endif

    /* Synchronization/Communication point.
     * Make sure that the function value of my best (new) vertex
     * is better than the function value at the current best
     * vertex. */

    if (plus(ndim + 5) < *fbest) {
      worked = true;

      if (debug) {
	(*fplpr) << " PDSWORK \n";
	(*fplpr) << " PDSWORK      REPLACED THE BEST VERTEX.\n";
	(*fplpr) << " PDSWORK \n";
      }

      /* Reinitialize FBEST, ALPHA, and BEST. */

      fprev = *fbest;
      *fbest = plus(ndim + 5);
      alpha = plus(3) / factor;
      best = (int)(plus(4));
    }
    else {
      worked = false;

      /* The best vertex has not been replaced.  Shrink the
       * simplex. */

      pdshrk(ndim, beta, &resize, &num);

      if (debug) {
	(*fplpr) << " PDSWORK \n";
	(*fplpr) << " PDSWORK      DID NOT REPLACE V0.\n";
	(*fplpr) << " PDSWORK      RESIZE = " << resize << endl;
	(*fplpr) << " PDSWORK \n";
      }

      /* Reinitialize ALPHA.  FBEST and BEST are unchanged.
       * Note that the "shrink" simplex is constructed using
       * RESIZE. */

      alpha = (double) (resize) / factor;
    }

    /* Update the simplex. */

    if (debug)
      (*fplpr) << " PDSWORK      ALPHA  = " << alpha  << endl;

    pdsupd(worked, best, ndim, pds_index, simplex, plus.Store(), alpha);

    /* Update the length of the longest edge for the convergence
     * test. */

    *length = fabs(alpha) * *length;

    /* Print iteration summary */

    (*fout) << d(count[0],5) 
	    << e(*fbest,13,4) << e(*length,13,4)
	    << d(worked,5) << d(count[1],10) << d(pds_index[0],8) << endl;

    done_code = pdsdone(maxitr, count[0], ndim, tol, length,
			&simplex[pds_index[0]*ndim], &r,
			finit, fprev, *fbest, fcn_tol, max_fevals,
			count[1], emesg, trpds);
   
    if (done_code > 0) {
      converged = true;
      (*fout) << emesg << "\n" << endl;
    }
  }

  (*fout).flush();

  pdsquit(debug, fplpr, count, r, flag, maxitr, tol);

  if (search_dir != NULL)
    delete[] search_dir;

  return(0);
}

void pdsquit(int debug, ofstream *fplpr, int *count, double r,
	     int *flag, int maxitr, double tol)
{
  /* Print last iteration info */

  if (debug) {
    (*fplpr) << " PDSWORK      ITERATION, FUNCTION, CONSTRAINT";
    (*fplpr) << " COUNT = " << d(count[0],4) << d(count[1],4)
	     << d(count[2],4) << "\n";
    (*fplpr) << " PDS \n";
  }

  /* Print the reason PDS terminated */

  if (debug) {
    (*fplpr) << " PDSWORK \n";

      if (*flag != 0) {
	(*fplpr) << " PDSWORK      FCN() OR CON() FAILED.  FLAG = "
		 << d(*flag,4) << "\n";
      }
      else if (count[0] >= maxitr) {
	(*fplpr) << " PDSWORK      MAXIMUM ITERATIONS REACHED.\n";
      }
      else {
	(*fplpr) << " PDSWORK      RESIDUAL < CONVERGENCE TOL. "
		 << e(r,12,4) << e(tol,12,4) << "\n";
      }
      (*fplpr) << " PDSWORK \n";
  }

  return;
}

} // namespace OPTPP

