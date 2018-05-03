
#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#ifdef HAVE_STD
#include <cstring>
#else
#include <string.h>
#endif

#ifdef WITH_MPI
#include "mpi.h"
#endif

#include "OptPDS.h"
#include "NLP0.h"
#include "pds.h"
#include "common.h"

using NEWMAT::ColumnVector;

// Structures for constraints and parallel configuration.

extern struct pdscon pdscon;
extern struct conbcmni conbcmni;

namespace OPTPP {

int pdsinit(NLP0* nlp, ostream *fplpr, int debug, int type,
	    int *flag, int *count, double scale, double *simplex,
	    double *vscales, double *length, int *pds_index,
	    double *fbest, double *rcond, double *work, double *work1,
	    double *work2, char *emesg, double tr_size, int first,
	    int trpds, double feas_tol) 
{
  /*******************************************************************
   *
   * This is a service subroutine used to finish constructing the
   * initial simplex (if one has not been entered by the user), to
   * determine the length of the longest edge in the initial simplex
   * (a value maintained throughout the search and used to test for
   * convergence), to compute the function value at the N+1 vertices
   * in the initial simplex, and to determine a "best" vertex (i.e.,
   * one with the least function value).  For large N on many
   * processors, this can be a significant cost, so the MPI version
   * parallelizes the function evaluations.  The logical function
   * PDSCHK() is used to identify infeasible vertices before they are
   * evaluated.
   *
   * Written by Virginia Torczon. 
   * MPI and Constrained versions written by David Serafini. 
   *
   * Last modifications: Sept96 v2.4epsilon, add variable scaling
   *                     VSCALES
   *                     July96 v2.4delta change TYPE=0 to TYPE=4
   *                     1Jun96, FLAG<>0 handling. 
   *                     1Jun95, PDSCONB/I mods. 
   * Parameters
   *
   *    Input
   *
   *       N             dimension of the problem to be solved
   * 
   *       FPLPR         File pointer of DEBUG file 
   *
   *       DEBUG         flag to signal whether or not to dump
   *                     debugging information to a file (see PDS()
   *                     for values)
   *
   *       TYPE          type of simplex chosen by the user: 
   *                         1      right-angled 
   *                         2      regular (edges of equal length) 
   *                         3      scaled right-angled 
   *                         4      user-defined (entered directly by
   *                                user)
   *
   *       FCN           subroutine to evaluate the objective function
   *
   *       CON           subroutine to evaluate the constraints
   *
   *       FLAG          used to pass back information about possible
   *                     failure in FCN to compute a function value
   *                     for a given point.  If FLAG<>0, PDS
   *                     terminates.
   *
   *       COUNT         Number of objective, constraint evaluations
   *                     performed; COUNT(2) is not used in the
   *                     unconstrained code.
   *
   *       SCALE         scale factor used if the simplex is to be
   *                     constructed (as specified during
   *                     initialization)
   *
   *       S             two-dimensional array containing the vertex
   *                     specified by the user to begin the search (as
   *                     well as the remaining N vertices of the
   *                     initial simplex if they have already been
   *                     entered by the user)
   *
   *       VSCALES       vector of scaling factors to be multiplied
   *                     into vertex coordinates before calling
   *                     FCN/CON.
   *
   *       UPARMDIM      Dimensions and actual sizes of U*PARM arrays.
   *                     (1) = # values in IPARMS,  (2) = # in RPARMS
   *                     (3) = dimension of IPARMS, (4) = dim. of
   *                                                      RPARMS
   *    Workspace
   *
   *       WORK          workspace vector of length Nndim (for DQRDC)
   * 
   *       WORK1         workspace vector of length N   (for DQRDC)
   *
   *       WORK2         workspace vector of length N   (for DQRDC)
   *
   *       UIPARM        for any int parameters required by the user 
   *                     to evaluate the function
   *
   *       URPARM        for any real parameters required by the user
   *                     to evaluate the function
   *
   *    Output
   *
   *       S             two-dimensional array containing the N+1
   *                     vertices of the simplex to be used to start
   *                     the search
   *
   *       LENGTH        length of the longest edge in the initial
   *                     simplex
   *
   *       PDS_INDEX     permutation array used to keep track of the 
   *                     current best vertex in the simplex
   *
   *       FBEST         the function value at the best vertex 
   *
   *******************************************************************/

  /* Variables */

  double temp, dist;
  int ndim  = nlp->getDim();
  ColumnVector x_curr = nlp->getXc();
  int i, i2, j, k, jbest, error;
  ColumnVector x(ndim);

  if(debug) (*fplpr) << "pdsinit: Entering\n";

#ifdef WITH_MPI

  char buffer[MPI_MAX_ERROR_STRING];
  int resultlen;
  struct {
    double value;
    int    loc;
  } localmin, globalmin;

#endif

  *flag = 0;

  /* Construct the simplex, if necessary. */

  if (type == 1)
    pdsrgt(ndim, scale, simplex);
  else if (type == 2)
    pdseql(ndim, scale, simplex);
  else if (type == 3)
    pdscld(ndim, scale, simplex);

  /* Determine the conditioning of the simplex.  Check, at the very
   * least, for numerical degeneracy in the simplex.  The simplex is
   * numerically degenerate.  To say that all is A-OK may be sanguine
   * as the simplex may be numerically ill-conditioned.  User beware! */

  pdsdgn(ndim, simplex, work, work1, work2, pds_index, rcond);

  if (*rcond + 1. == 1.) {
    strcpy(emesg, "pdsinit: initial simplex is degenerate");
    error = 12;
  }
  else
    error = 0;
    
  /* Determine the length of the longest edge in the simplex. */

  *length = pdslen(ndim, type, simplex, scale, work);
    
  /* get the current best function value.  If this is for TRPDS, this
   * value will be the function value at the current point.  JBEST
   * will be set to its position in the initial simplex.  If this is
   * stand-alone PDS, the best value will be +inf. */

  *fbest = nlp->getF();

  if (trpds) {

    if (first) {
      jbest = 1;
    }
    else {
      jbest = 2;
    }
  }
  else
    jbest = -1;

  error = 0;


  /* The first feasible vertex is JBEST.  Compute the rest and find
   * the min.  Distribute the rest of the vertices (cyclic) amongst
   * the processors to balance the load, then find the global min
   * value across all processors.  [NOTE: Processor numbers start at
   * 0, so shift them up by 1.]  [NOTE: The scalar code has ME=0,
   * NPROC=1, so this will work without change.] */
    
  i2 = pdscon.nproc;

  for (j = pdscon.me; i2 <= 0 ? j >= ndim : j <= ndim; j += i2) {

    for (i = 0; i < ndim; i++)
      work1[i] = simplex[i + j * ndim] * vscales[i];

    if (pdschk(nlp, ndim, x_curr.Store(), work1, tr_size, &dist, trpds, feas_tol)) {
      count[2]++;

      if (*flag != 0) {
	return 0;
      }

      for (k = 0; k < ndim; k++)
	x(k+1) = work1[k];

      temp = nlp->evalF(x);
      count[1]++;

      if (*flag != 0) {
	return 0;
      }

      if (temp < *fbest) {
	*fbest = temp;
	jbest = j;
      }
    }
    else {

      /* If the constraint violated is not bounds, then an ineq
       * constraint evaluation must have been done. */

      if (i == 0 || i > ndim)
	count[2]++;

      if (*flag != 0) {
	return 0;
      }

// PJW  12-01-99
//      cout << "pdsinit: P" << d(pdscon.me,2) << "-> returning early\n";
//      cout << "pdsinit: P" << d(pdscon.me,2) << "-> constraint violation\n";
// //     strcpy(emesg, "pdsinit: bound violation\n");
//      return 69;
// PJW  12-01-99
  }
  }

#ifdef WITH_MPI

  localmin.value = *fbest;
  localmin.loc   = jbest;

  /* The MPI_MINLOC operator takes a 2 element vector and finds
   * the minimum across the processors of the first element,
   * carrying the 2nd element along.  The ALLREDUCE function
   * makes sure every processor gets the answer. */

  error = MPI_Allreduce(&localmin, &globalmin, 1,
			MPI_DOUBLE_INT, MPI_MINLOC,
			MPI_COMM_WORLD);
  if (error != MPI_SUCCESS) {
    MPI_Error_string(error, buffer, &resultlen);
    printf("\npdsinit: MPI Error - %s\n", buffer);
    strcpy(emesg, "pdsinit: error returned by MPI_Allreduce\n");
    return 15;
  }

  *fbest = globalmin.value;
  pds_index[0] = globalmin.loc;

#else
    
  /* Scalar version: set the pointer to the best vertex in the
   * simplex. */
    
  pds_index[0] = jbest;

#endif

  /* In stand-alone PDS, return error if no vertex has a function
   * value less than +inf. */

  if (pds_index[0] == -1) {
    strcpy(emesg, "pdsinit: no vertex with f < +inf");
    return 13;
  }

  /* Initialize the rest of the PDS_INDEX array so every element
   * points to itself except the element which is the minimum,
   * which points */

  for (j = 1; j <= ndim; j++)
    pds_index[j] = j;

  pds_index[pds_index[0]] = 0;

  return error;
}

} // namespace OPTPP
