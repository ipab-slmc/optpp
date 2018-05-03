
//--------------------------------------------------------------------
// Copyright (C) 1993,1994: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//--------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#include <string>
#ifdef HAVE_STD
#include <cmath>
#include <cstdlib>
#include <cstring>
#else
#include <math.h>
#include <stdlib.h>
#include <string.h>
#endif

#include "OptPDS.h"
#include "precisio.h"
#include "pds.h"
#include "common.h"
#include "ioformat.h"

using namespace std;
using NEWMAT::ColumnVector;
using NEWMAT::Real;

// Structures for constraints and parallel configuration.

extern "C" {
  struct conbcmni conbcmni;
  struct pdscon pdscon;
}

namespace OPTPP {

void OptPDS::initOpt()
{
  ret_code = 0;

  if (!trpds) {
    nlp->initFcn();
    ColumnVector x_init = nlp->getXc();
    double perturb;
    for (int i=1; i <= nlp->getDim(); i++) {
      for (int j=1; j <= nlp->getDim()+1; j++) {
	simplex(i,j) = x_init(i);
      }
    }
    for (int i=1; i<= nlp->getDim(); i++) {
      perturb = x_init(i)*.01;
      simplex(i,i+1) = x_init(i) + perturb;
    }
  }

  readOptInput();

  if(nlp->hasConstraints()){
    CompoundConstraint* constraints = nlp->getConstraints();
    ColumnVector xstart = nlp->getXc();
    double feas_tol = tol.getCTol();
    bool feasible = constraints->amIFeasible(xstart, feas_tol);
    if (!feasible)
      *optout << "OptPDS WARNING:  Initial guess not feasible.\n"
	      << "PDS may be unable to make progress." << endl;
  }

  // If ret_code = 15, then MPI has not been initialized

  ret_code = pdscom(mesg);
}

void OptPDS::optimize()

  // Parallel Direct Search
{
  int i, j;
  
  // Allocate local vectors 

  int ndim = dim;
  ColumnVector x(ndim), xold(ndim);

  int count, ierr;
  double fbest, length;
  double alpha = 0.99;

  double scale     = 1.0;
  int restart        = 0;
  int *pds_index = new int[ndim+1];
  char scheme_name[256];    /* SCHEME default file name */
  char *tmpdir;
  int type, sss, cflag;

  int pds_debug;
  double pds_tol, pds_fcn_tol, feas_tol;
  int maxiter, max_fevals;
  int loc_first, loc_trpds;

  SpecOption SpecTmp = nlp->getSpecOption();

  // Initialize algorithm and set parameters.

  initOpt();

  nlp->setSpecOption(NoSpec);

  pds_tol     = tol.getStepTol();
  pds_fcn_tol = tol.getFTol();
  feas_tol    = tol.getCTol();
  maxiter     = tol.getMaxIter();
  max_fevals  = tol.getMaxFeval();

  pds_debug    = debug_;
  loc_first    = (int) first;
  loc_trpds    = (int) trpds;
  type         = getSimplexType();
  sss          = getSSS();
  cflag        = getCreateFlag();

  if (!trpds) {
    fbest = 1.e50;
    nlp->setF(fbest);
  }

  if (ret_code == 0) {

    ColumnVector pds_simplex(ndim*(ndim+1));

    tmpdir = getenv("PWD");
    if (tmpdir == NULL) {
      *optout << "pds WARNING: TMP environment variable not set./n"
	      << "Using /tmp..." << endl;
      strcpy(scheme_name, "/tmp");
    }
    else
      strcpy(scheme_name, tmpdir);
    strcat(scheme_name, "/");
    strcat(scheme_name, getSchemeFileName());
	
    xold = nlp->getXc();
	
    for (j = 0; j < ndim+1; j++) {
      for (i = 1; i <= ndim; i++) {
	pds_simplex(i+j*ndim) = simplex(i,j+1)/vscales(i);
      }
    }

    iter_taken = fcn_evals = 0;

    // Call main PDS routine.

    ierr = pdsopt(nlp, optout, pds_simplex.Store(), pds_index, cflag,
		  scheme_name, pds_debug, restart, alpha, maxiter,
		  sss, scale, vscales.Store(), pds_tol, type, &fbest,
		  &count, mesg, pds_fcn_tol, tr_size, &length,
		  max_fevals, loc_first, loc_trpds, feas_tol);

    // set output information.

    ret_code = ierr;
    setReturnCode(ret_code);

    if (ret_code != 13) {
      for (i = 1; i <= ndim; i++)
	x(i) = pds_simplex(i)*vscales(i);

      nlp->setX(x);
      nlp->setF(fbest);
      setSimplexSize(length);
      iter_taken = count;
      fcn_evals  = nlp->getFevals();
      nlp->setSpecOption(SpecTmp);
    }
  }
  if (pds_index != NULL)
    delete[] pds_index;
}

void OptPDS::printStatus(char *s)

  // set Message
{
  if (pdscon.me == 0) {
      *optout << "\n\n=========  " << s << "  ===========\n\n";
      *optout << "Optimization method       = " << method << "\n";
      *optout << "Dimension of the problem  = " << nlp->getDim()
	     <<	"\n";
      *optout << "Search Scheme Size        = " << search_scheme_size
	     << "\n";
      *optout << "Simplex type              = " << simplex_type
	     <<	"\n";
      *optout << "Return code               = " << ret_code << " ("
	     << mesg << ")\n";
      *optout << "No. iterations taken      = " << iter_taken  << "\n";
      *optout << "No. function evaluations  = " << fcn_evals << "\n";
      
      nlp->fPrintState(optout, s);
      tol.printTol(optout);
  }
}

int OptPDS::checkConvg()

  // Check convergence
{
  int    n;
  double stol, ftol, rftol;
  double xnorm, snorm;
  ColumnVector xc;

  double step_tol, fvalue;

  n  = nlp->getDim();
  xc = nlp->getXc();
  fvalue = nlp->getF();

  xnorm =  Norm2(xc);
  
  // Test 1. step tolerance 

  ColumnVector step(n);
  step = xc - xprev;
  step_tol = tol.getStepTol();
  snorm = Norm2(step);
  stol  = step_tol*max(1.0,xnorm);

  if (snorm  <= stol) {
    strcpy(mesg,"CheckConvg: Step tolerance test passed");
    *optout << "CheckConvg: snorm = " << e(snorm,12,4) 
	   << "  stol = " << e(stol,12,4) << "\n";
    return 1;
  }
  
  // Test 2. function tolerance

  Real deltaf = fprev - fvalue;
  ftol = tol.getFTol();
  rftol = ftol*max(1.0,fabs(fvalue));

  if (deltaf <= rftol) {
    strcpy(mesg,"Function tolerance test passed");
    *optout << "CheckConvg: deltaf = " << e(deltaf,12,4) 
	   << "  ftol = " << e(ftol,12,4) << "\n";
    return 2;
  }

  // Nothing to report 

  strcpy(mesg," ");

  return 0;
}

void OptPDS::readOptInput()

  // Read opt.input file if it exists
{
  /*******************************************************************
   *
   * A VERY simple routine for reading the optimization parameters *
   * We should really make this more general, but as a first pass this
   * * will have to do.
   * 
   * The input file should be of the form keyword = value
   * where keyword is one of the following
   * 
   *    search      = trustregion
   *    diff_option = forward
   *    max_iter    = 100
   *    maxfeval    = 1000
   *    grad_tol    = 1.e-6
   *    fcn_tol     = 1.e-9
   *    max_step    = 100.0
   *    fcn_accrcy  = 1.e-9
   *
   *******************************************************************/

  int  index, max_iter, max_feval, pds_max_feval;
  int  sss, s_type;
  real grad_tol,  max_step, fcn_accrcy;
  real pds_fcn_tol, feas_tol;

  char token[80], ignore[80], equals[1];

  //
  // Keywords allowed
  //

  string keyword;
  string cdebug("debug");
  string cdiff_option("diff_option");
  string cfcn_accrcy("fcn_accrcy");
  string cfcn_tol("pds.fcn_tol");
  string cgrad_tol("grad_tol");
  string cmaxfeval("maxfeval");
  string cpdsmaxfeval("pds.maxfeval");
  string cmaxiter("pds.maxiter");
  string cmax_step("max_step");
  string csss("sss");
  string csimplex("simplex");
  string cfeas_tol("feas_tol");

  string diff_option, debug_flag;

  int keyword_count = 0;

  // 
  // Default name of input file
  //

  char *opt_input  = {"opt.input"};

  //
  // Open opt.input file and check to see if we succeeded
  //

  ifstream optin(opt_input);

  if (!optin.rdbuf()->is_open()) {

    if (debug_) {
      *optout << "ReadOptInput: No opt.input file found\n";
      *optout << "ReadOptInput: Default values will be used\n";
    }

    return;
  }

  // Read opt.input file.

  if (debug_) *optout << "OptPDS::ReadOptInput: Reading opt.input file"
		    << "\n"; 

  optin >> token;

  while (!optin.eof()) {

    keyword = token;
    keyword_count++;

    if (keyword == cdiff_option) {
      optin >> equals >> token;
      diff_option = token;
    }
    else if (keyword == cdebug) {
      optin >> equals >> token;
      debug_flag = token;

      if ( debug_flag == "true") {
	setDebug();
	nlp->setDebug();
      }
    }    
    else if (keyword == cfcn_accrcy) {
      //optin >> equals >> fcn_accrcy;
      //nlp->setFcnAccrcy(fcn_accrcy);
      optin >> equals >> index >> fcn_accrcy;
      nlp->setFcnAccrcy(index, fcn_accrcy);
    }    
    else if (keyword == cfcn_tol) {
      optin >> equals >> pds_fcn_tol;
      setFcnTol(pds_fcn_tol);
    }    
    else if (keyword == cgrad_tol) {
      optin >> equals >> grad_tol;
      setGradTol(grad_tol);
    }    
    else if (keyword == cmaxfeval) {
      optin >> equals >> max_feval;
      setMaxFeval(max_feval);
    }    
    else if (keyword == cpdsmaxfeval) {
      optin >> equals >> pds_max_feval;
      setMaxFeval(pds_max_feval);
    }    
    else if (keyword == cmaxiter) {
      optin >> equals >> max_iter;
      setMaxIter(max_iter);
    }
    else if (keyword == cmax_step) {
      optin >> equals >> max_step;
      setMaxStep(max_step);
    }
    else if (keyword == csss) {
      optin >> equals >> sss;
      setSSS(sss);
    }
    else if (keyword == csimplex) {
      optin >> equals >> s_type;
      setSimplexType(s_type);
    }
    else if (keyword == cfeas_tol) {
      optin >> equals >> feas_tol;
      setConTol(feas_tol);
    }
    else {
      *optout << "pds: Unrecognized keyword '" << keyword << "'. "
	     << "Skipping the rest of this line\n";
      optin.getline(ignore, sizeof(ignore));
    }

    optin >> token;
  }

  if (debug_) {
    *optout << "\n\n======  Summary of input file  ======\n\n";

    *optout << cmaxiter     << " = " << max_iter << "\n";
    *optout << cmaxfeval    << " = " << max_feval << "\n";
    *optout << cpdsmaxfeval << " = " << pds_max_feval << "\n";
    *optout << cgrad_tol    << " = " << grad_tol << "\n";
    *optout << cfcn_tol     << " = " << pds_fcn_tol << "\n";
    *optout << cmax_step    << " = " << max_step << "\n";
    ColumnVector fcnacc = nlp->getFcnAccrcy();
    for(int i = 1; i <= fcnacc.Nrows(); i++)
        *optout << cfcn_accrcy  << " = " << fcnacc(i) << "\n";
    *optout << csss         << " = " << sss << "\n";
    *optout << csimplex     << " = " << s_type << "\n";
    
    tol.printTol(optout);
  }
}

} // namespace OPTPP
