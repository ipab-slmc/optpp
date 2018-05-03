
#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#ifdef HAVE_STD
#include <cmath>
#include <cerrno>
#include <cstring>
#include <cstdio>
#else
#include <math.h>
#include <errno.h>
#include <string.h>
#include <stdio.h>
#endif

#ifdef WITH_MPI
#include "mpi.h"
#endif

#include "OptPDS.h"
#include "pds.h"
#include "common.h"
#include "ioformat.h"

using namespace std;

using NEWMAT::ColumnVector;

/* Structures for constraints and parallel configurations. */

extern struct pdscon pdscon;
extern struct conbcmni conbcmni;

namespace OPTPP {

int create_scheme(ostream *, int, int, char *, int *, int);

void pdslogerr(int, int, double *, int, ostream *, double, int,
	       double, double *, int, int, ColumnVector&, ColumnVector&);

/* Globals for this file only */

double rcond;
int flag, upper;
ofstream fpdebug;
FILE *fpscheme;

int pdsopt(NLP0* nlp, ostream *fout, double *simplex, int *pds_index,
	   int cflag, char *scheme_name, int debug, int reset_param,
	   double alpha, int maxitr, int sss, double scale,
	   double *vscales, double tol, int type, double *fbest,
	   int *iter, char *emesg, double fcn_tol, double tr_size,
	   double *length, int max_fevals, int first, int trpds,
	   double feas_tol)
{
  /*******************************************************************
   *
   * Main program to solve nonlinear optimization problems using the
   * parallel direct search methods.
   *      PDSOPT solves unconstrained problems.
   *      PDSCONB solves bounds constrained problems.
   *      PDSCONI solves nonlinear inequality constrained problems,
   *              including bounds constraints as a special case.
   *
   * All versions require a subroutine called FCN to evaluate the
   * objective function.  PDSCONI (inequality constrainted) requires a
   * subroutine called CON to evaluate the inequality constraints.
   * The program requires 2 input files to run.  One is a text file
   * called "INPUT" which contains the initial point and various PDS
   * control variables (see the documentation for details) and the
   * other is an unformatted data file written by the "crtpat"
   * program, which contains the pattern data (aka "scheme" in this
   * program).  If a solution 'restart' is requested in the input
   * (simplex-type <0), the restart file is an unformatted data file
   * containing the final simplex and various control variable values.
   * The file is named "RESTART.#", where # is the number of variables
   * in the problem (i.e. the number of dimensions).
   *
   * PDS works best if the values it works with are uniformly scaled.
   * To achieve this, the INPUT file must specify scales factors to
   * apply to the variable values that will (ideally) make them all
   * have nearly the same magnitude.  These scale factors are divided
   * into the actual values to compute the simplex points.  E.g. if
   * the
   *      initial variable values are:  1 10 100 1000 
   *      then scale factors might be:  1 10 100 1000 
   *      making the initial vertex:    1  1   1    1. 
   * It may not always be appropriate to scale the variables so the
   * initial vertex is (1,1,1,...).  It is the range of values that a
   * particular variable takes that determines what its scaling should
   * be.  If the range for the variables in the previous example was
   * such that each was expected to vary by +-10, then it would be
   * more appropriate to have scale factors 1,1,1,1.
   *
   * If interactive output is requested (via setting DEBUG=-1) a
   * subroutine called PDSLOG() is called at the end of each
   * iteration.  The default version of this routine can be used, or
   * an application specific version can be linked with the program to
   * perform whatever output is required.  The arguments to PDSLOG()
   * include the current best point and the objective function value,
   * the iteration number, the count of function evaluations, etc.
   * See PDS() for details.  PDSLOG() is also called after the simplex
   * has been initialized, with the best point in the initial simplex.
   * The iteration number for this call is 0.  This allows the initial
   * state to be logged as well as allowing PDSLOG() to initialize
   * whatever internal state it may have, since iteration 0 is always
   * the first call to PDSLOG(). PDSLOG() is also called at the end of
   * the program, with iteration number -1, to allow for any
   * termination actions to take place.  The argument values for this
   * call are not valid (except for Proc#,DbgUnit,N,F).  The default
   * version of PDSLOG() does nothing at termination.
   *
   * Written by Virginia Torczon. 
   * MPI and Constraint versions written by David Serafini.
   *
   * Last modifications: Sept96 v2.4z (Added caching of FCN,CON
   *                                   values) 
   *                     Sept96 v2.4e (Added variable scaling VSCALES
   *                                   array)
   *                     July96 v2.4d (Added restart capability)
   *                     June96 (Added FLAG<>0 handling) 
   *                     April96 (Added PDSLOG() for DEBUG=-1) 
   *
   * We include here a description of the workspace allocated in the
   * values used, and the declared dimension, for the int and real
   * PARMS arrays, respectively.  The PARM arrays are set by the user
   * in the input file.  The DEBUG flag signals whether or not to dump
   * debugging information
   *      -2   echo all function values and vertices to debug unit
   *      -1   interactive debugging output (invokes the PDSLOG()
   *           routine, which can be replaced by the application.  See
   *           PDS() for details.  Default PDSLOG() prints
   *           (iter,obj,col#,residual)
   *       0   no debugging output 
   *       1   log starting simplex and each iterate (V0) 
   *       2   include some intermediate debugging 
   *       3   include almost all intermediate debugging 
   *       4   include all debugging (including SCHEME) 
   *
   *******************************************************************/

  /* Variables */

  double factor, fprev, temp;
  int ndim = nlp->getDim();
  int i, j, ierr, beta, vbest;
  int count[3];

#ifdef WITH_MPI

  char buffer[MPI_MAX_ERROR_STRING], debug_name[100];
  int countsum[3], resultlen, status;

#endif

  ierr = 0;

  /* Verify that the information needed to define the problem is
   * valid.  UPPER is initialized to the size of the search scheme
   * array, and returned as the actual size used. */

  int limit = (ndim+2)*50*sss;
  upper = limit;

  if (ndim < 1) {
//    cout << "pdsopt: P" << d(pdscon.me,2) << "-> returning early\n";
//    cout << "pdsopt: P" << d(pdscon.me,2) << "-> ndim =" << d(ndim,2) << "\n";
    strcpy(emesg, "pdsopt: dimension < 1");
    return 1;
  }

  if (maxitr < 1) {
//    cout << "pdsopt: P" << d(pdscon.me,2) << "-> returning early\n";
//    cout << "pdsopt: P" << d(pdscon.me,2) << "-> maxitr =" << d(maxitr,2) << "\n";
    strcpy(emesg, "pdsopt: maximum # iterations < 1");
    return 2;
  }

  if (type < 1 || type > 4) {
//    cout << "pdsopt: P" << d(pdscon.me,2) << "-> returning early\n";
//    cout << "pdsopt: P" << d(pdscon.me,2) << "-> type =" << d(type,2) << "\n";
    strcpy(emesg, "pdsopt: simplex type < 1 or type > 4");
    return 3;
  }

  if (reset_param < 0) {

    if (pdscon.me == 0) {
      (*fout) <<"\npdsopt: WARNING --- reset_param < 0\n"
	      <<"pdsopt: reset_param will be set = 0\n\n";
    }

    reset_param = 0;
  }

  if (sss < 2*ndim) {

    if (pdscon.me == 0) {
      (*fout) <<"\npdsopt: WARNING --- sss < 2n.\n"
	      <<"pdsopt: PDS is not guaranteed to converge.\n\n";
    }
  }

// PJW  
  ColumnVector ltmp(ndim), utmp(ndim);
  bool hasConstraints = nlp->hasConstraints();
  if (hasConstraints){
    CompoundConstraint* constraints;
    constraints = nlp->getConstraints();
    conbcmni.ncb  = constraints->getNumOfVars(); 
    // Currently, OPT++ only supports bound-constrained PDS.
    conbcmni.ncni = 0; 
    ltmp = constraints->getLower();
    utmp = constraints->getUpper();
  }
  else{
    conbcmni.ncb  = 0;
    conbcmni.ncni = 0;
    ltmp = 0;
    utmp = 0;
  }
// PJW  

  /* Check that lower bounds are less than upper bounds. */

  for (i = 1; i <= conbcmni.ncb + conbcmni.ncni; i++) {

    if (ltmp(i) >= utmp(i)) {
//      cout << "pdsopt: P" << d(pdscon.me,2) << "-> returning early\n";
//      cout << "pdsopt: P" << d(pdscon.me,2) << "-> lowerbnd =" << e(ltmp(i),14,6) << "\n";
      strcpy(emesg, "pdsopt: lower bound exceeds corresponding upper bound");
      return 6;
    }
  }

  /* Calculate the upper limit on the size of the search scheme given
   * the amount of memory allocated in the main program.  Keep in
   * mind that each column (point) in the search scheme will be of
   * length n + 2. */

  upper /= ndim + 2;

  /* Open a file for debugging information if the user requested.
   * Initialize the function value output if that is what is
   * requested, i.e. (debug = -2) Header record is: words per
   * record (variables + function value + column pds_index +
   * iteration) #records (-1 means unknown) processor number for
   * this file */

  if (debug) {

#ifdef WITH_MPI

    sprintf (debug_name, "DEBUG.%d", pdscon.me);
    fpdebug.open(debug_name);

#else

    fpdebug.open("DEBUG");

#endif

//    fpdebug_fd = fpdebug.rdbuf()->fd();

    if (debug)
      fpdebug << d(ndim+3,4) << ", -1, " << d(pdscon.me,4) << "\n";
  }

  /* Read in the search scheme and determine the "shrink" factor
   * (which depends on the size of the search scheme that has been
   * specified).  Open a file for the search scheme. */

  int *scheme = new int[(ndim+2)*(ndim+1+50*sss)];

  if (cflag == true) {

#ifdef SHARED

    if (pdscon.me == 0) {
      ierr = create_scheme(fout, ndim, limit, scheme_name, scheme, debug);

      if (ierr != 0) {
// 	cout << "pdsopt: P" << d(pdscon.me,2) << "-> returning early\n";
// 	cout << "pdsopt: P" << d(pdscon.me,2) << "-> create_scheme failed\n";
	strcpy(emesg, "pdsopt: can't create scheme file");
	return 8;
      }
    }

#else

    ierr = create_scheme(fout, ndim, limit, scheme_name, scheme, debug);

    if (ierr != 0) {
//      cout << "pdsopt: P" << d(pdscon.me,2) << "-> returning early\n";
//      cout << "pdsopt: P" << d(pdscon.me,2) << "-> create_scheme failed\n";
      strcpy(emesg, "pdsopt: can't create scheme file");
      return 8;
    }

#endif

  }

#ifdef SHARED

  MPI_Bcast(scheme_name, 256, MPI_CHAR, 0, MPI_COMM_WORLD);

#endif

  if (scheme != NULL)
    delete[] scheme;

  if ((fpscheme = fopen (scheme_name, "r")) == NULL) {
//    cout << "pdsopt: P" << d(pdscon.me,2) << "-> returning early\n";
//    cout << "pdsopt: P" << d(pdscon.me,2) << "-> open scheme failed\n";
    strcpy(emesg, "pdsopt: can't open scheme file");
    return 9;
  }

  ierr = pdsget(ndim, fpscheme, &sss, &factor, &beta, emesg);
  (*fout) << "pdsopt: factor  = " << e(factor,12,4) << "\n";
  (*fout) << "        beta    = " << d(beta,11) << "\n";

  if (ierr != 0) {
//    cout << "pdsopt: P" << d(pdscon.me,2) << "-> returning early\n";
//    cout << "pdsopt: P" << d(pdscon.me,2) << "-> ierr =" << d(ierr,2) << "\n";
    return(ierr);
  }

  *iter = 0; fprev = 1.0e30;

  for (j = 0; j <= reset_param; j++) {

    /* Initialize counters for iterations, function and constraint
     * evals. */
      
    count[0] = 0;
    count[1] = 0;
    count[2] = 0;

    /* Optimize. */

    ierr = pdswork(nlp, fout, &fpdebug, debug, tol, maxitr, sss, &flag,
		   factor, beta, simplex, vscales,
		   pds_index, fbest, length, count, type, scale,
		   &rcond, emesg, fcn_tol, tr_size, max_fevals,
		   first, trpds, feas_tol, fpscheme);
    fclose(fpscheme);

    if (flag != 0) {
      ierr = -1;
      pdslogerr(ierr, ndim, simplex, type, fout, tol, maxitr,
		scale, vscales, debug, sss, ltmp, utmp);
      return(ierr);
    }

    if (ierr != 13) {

      /* Place best vertex in the first column of s for a possible
       * restart */ 

      vbest = pds_index[0];

      for (i = 0; i < ndim; i++) {
	temp = simplex[i];
	simplex[i] = simplex[i + ndim*vbest];
	simplex[i + ndim*vbest] = temp;
      }

      pds_index[0] = 0;
      *iter += count[0];

      if (*iter >= maxitr) {
	strcpy(emesg, "pdsopt: maximum # iterations reached");
	ierr = 14;
	break;
      }
	    
      if (j > 0)

	if (fabs(*fbest) > alpha*fabs(fprev)) {

	  if (pdscon.me == 0)
	    (*fout) <<"pdsopt: insufficient decrease in restart.\n";

	  break;
	}
	    
      fprev = *fbest;
    }
  }

  /* Done.  Normal termination.  Shut down gracefully. */

#ifdef WITH_MPI

  /* Parallel Version: first sum the function evaluations done by each
   * process.  Skip this for processors with errors.
   * Nothing to do if only one node */

  if (pdscon.nproc > 1) {

    /* Sum COUNT(2:3) on all processes into COUNTSUM on process 0
     * Copy back into COUNT(2:3) on processor 0. */

    status = MPI_Reduce(&count[1], countsum, 2, MPI_INT, MPI_SUM,
			0, MPI_COMM_WORLD);

    if (status != MPI_SUCCESS) {
      MPI_Error_string(status, buffer, &resultlen);
      printf("\npdsopt MPI Error - %s\n", buffer);
      strcpy(emesg, "pdsopt: error returned by MPI_Reduce\n");
      return 15;
    }

    if (pdscon.me == 0) {
      count[1] = countsum[0];
      count[2] = countsum[1];
    }
  }

#endif

  pdslogerr(ierr, ndim, simplex, type, fout, tol, maxitr, scale,
	    vscales, debug, sss, ltmp, utmp);

  return(ierr);
}
    
void pdslogerr(int error, int ndim, double *simplex, int simplex_type,
	       ostream *fout, double tol, int maxitr, double scale,
	       double *vscales, int debug, int sss,
	       ColumnVector& ltmp, ColumnVector& utmp)
{
  /* Variables */

  int i, j;

  if (pdscon.me != 0) {
    goto L99999;
  }

  if (error == -1) {
    (*fout) << "\nList of Parameters...\n\n";
    (*fout) << "     dimension                = " << d(ndim,11) << "\n";
    (*fout) << "     # bound constraints      = "
	    << e(conbcmni.ncb,30,14) << "\n";
    (*fout) << "     # inequality constraints = "
	    << e(conbcmni.ncni,30,14) << "\n";
    (*fout) << "     convergence tolerance    = " << e(tol,30,14)
	    << "\n"; 
    (*fout) << "     maximum # iterations     = " << d(maxitr,11)
	    << "\n"; 
    (*fout) << "     initial vertex           = "
	    << e(simplex[0],30,14) << "\n";
    for (i = 1; i < ndim; i++) {
      (*fout) << "                                "
	      << e(simplex[i],30,14) << "\n";
    }

    if (simplex_type == 4) {

      for (i = 1; i <= ndim; i++) {

	for (j = 0; j < ndim; j++) {
	  (*fout) << "                                "
		  << e(simplex[j + i*ndim],30,14) << "\n";
	}
      }
    }

    (*fout) << "     vertex scales            = " << e(vscales[0],30,14)
	    << "\n";

    for (i = 1; i < ndim; i++) {
      (*fout) << "                                "
	      << e(vscales[i],30,14) << "\n"; 
    }

    (*fout) << "     lower bounds             = \n";

    for (i = 1; i <= conbcmni.ncb; i++) {
      (*fout) << "                                "
	      << e(ltmp(i),30,14) << "\n";
    }

    (*fout) << "     upper bounds             = \n";

      for (i = 1; i <= conbcmni.ncb; i++) {
	(*fout) << "                                "
		<< e(utmp(i),30,14) << "\n";
      }

      (*fout) << "     simplex type             = " << d(simplex_type,11)
	      << "\n";
      (*fout) << "     simplex scale            = " << e(scale,30,14)
	      << "\n"; 
      (*fout) << "     debug flag               = " << d(debug,11) << "\n";
      (*fout) << "     # pattern points         = " << d(sss,11) << endl;
  }

  /* Write a termination message to the result file. */

  else if (error == 0) {

    /* Nothing left to do except terminate. */
  } 
  else if (error == 1) {
    (*fout) << "('                                             ')\n";
    (*fout) << "('                                             ')\n";
    (*fout) << "(' EITHER TRYING TO RUN PVM ON ONE PROCESSOR OR')\n";
    (*fout) << "(' NUMBER OF PROCESSORS EXCEEDS INTERNAL LIMIT.')\n";
    (*fout) << "(' EXITED WITHOUT CALLING PDS.                 ')\n";
    (*fout) << "(' CHECK PVM_NPROC_LIM IN FILE pdscon.h.       ')\n";
    (*fout) << "('                                             ')\n" << endl;
  }
  else if (error == 5) {
    (*fout) << "('                                             ')\n";
    (*fout) << "('                                             ')\n";
    (*fout) << "(' THE SIMPLEX IS NUMERICALLY DEGENERATE.      ')\n";
    (*fout) << "(' EXITED WITHOUT CALLING PDS.                 ')\n";
    (*fout) << "(' SEE DISCUSSION IN EXTERNAL DOCUMENTATION.   ')\n";
    (*fout) << "('                                             ')\n";
    (*fout) << "('                                             ')\n" << endl;
  }
  else if (error == 6) {
    (*fout) << "('                                             ')\n";
    (*fout) << "('                                             ')\n";
    (*fout) << "(' SEARCH SCHEME WAS OF THE WRONG DIMENSION.   ')\n";
    (*fout) << "(' EXITED WITHOUT CALLING PDS.                 ')\n";
    (*fout) << "(' RERUN PROGRAM TO CREATE SEARCH STRATEGY.    ')\n";
    (*fout) << "('                                             ')\n";
    (*fout) << "('                                             ')\n" << endl;
  }
  else if (error == 7) {
    (*fout) << "('                                             ')\n";
    (*fout) << "('                                             ')\n";
    (*fout) << "(' INSUFFICIENT NUMBER OF POINTS IN SCHEME.    ')\n";
    (*fout) << "(' EXITED WITHOUT CALLING PDS.                 ')\n";
    (*fout) << "(' SEE DISCUSSION IN EXTERNAL DOCUMENTATION.   ')\n";
    (*fout) << "('                                             ')\n";
    (*fout) << "('                                             ')\n" << endl;
  }
  else if (error == 9) {
    (*fout) << "('                                             ')\n";
    (*fout) << "('                                             ')\n";
    (*fout) << "(' EVERY VERTEX IN THE INITIAL SIMPLEX IS INFEASIBLE.')\n";
    (*fout) << "(' EXITED WITHOUT CALLING PDS.                 ')\n";
    (*fout) << "('                                             ')\n";
    (*fout) << "('                                             ')\n" << endl;
  }
  else if (error == 10) {
    (*fout) << "('                                             ')\n";
    (*fout) << "('                                             ')\n";
    (*fout) << "(' A LOWER BOUND EXCEEDS THE CORRESPONDING     ')\n";
    (*fout) << "UPPER BOUND AT PDS_INDEX " << d(upper,8) << "\n";
    (*fout) << "('                                             ')\n";
    (*fout) << "('                                             ')\n" << endl;
  }
  else if (error == 11) {
    (*fout) << "('                                             ')\n";
    (*fout) << "('                                             ')\n";
    (*fout) << "(' MPI ERROR IN PDS.  PDS TERMINATED.          ')\n";
    (*fout) << "('                                             ')\n";
    (*fout) << "('                                             ')\n" << endl;
  }

L99999:

  if (debug) {
    fpdebug.close();

    if (pdscon.me == 0)
      (*fout) << "pdsopt: exit\n";

  }

  return;
}

} // namespace OPTPP
