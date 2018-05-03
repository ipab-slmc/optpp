//------------------------------------------------------------------------
// Copyright (C) 1993: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//------------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#ifdef HAVE_STD
#include <cstdio>
#include <cstring>
#include <ctime>
#else
#include <stdio.h>
#include <string.h>
#include <time.h>
#endif

#include <string>

using namespace std;

#include "OptBCNewtonLike.h"
#include "precisio.h"
#include "cblas.h"
#include "ioformat.h"

using NEWMAT::Real;
using NEWMAT::FloatingPointPrecision;
using NEWMAT::ColumnVector;
using NEWMAT::Matrix;
using NEWMAT::DiagonalMatrix;
using NEWMAT::LowerTriangularMatrix;
using NEWMAT::SymmetricMatrix;

#ifdef WITH_MPI
#include "mpi.h"
#endif

namespace OPTPP {

static char* class_name = "OptBCNewtonLike";

//------------------------------------------------------------------------
//
//   Constrained Newton Base Class functions
//   Notes:
//   These functions are first declared in OptBCNewtonLike.h as
//   virtual functions for the abstract base class OptBCNewtonLike
//   Therefore we need to define them so that the derived classes
//   can be instantiated by the end user.  Of course the user
//   can provide his/her own version of these.
//------------------------------------------------------------------------

//------------------------------------------------------------------------
//
// First the default functions
// defaultacceptStep
// defaultcomputeSearch
//
//------------------------------------------------------------------------
void OptBCNewtonLike::defaultAcceptStep(int iter, int step_type)
{
// Successful step
// Print out iteration summary and anything else
// you might want to do before the next iteration

  if (trace) 
    *optout << "\n***** OptBCNewtonLike:defaultacceptStep\n";

  NLP1* nlp = nlprob();
  int n     = nlp->getDim();

  static char *steps[] = {"C", "D", "N", "B"};
  ColumnVector xc(n), grad(n);
  double fvalue, gnorm;

  xc     = nlp->getXc();
  mem_step   = xc - xprev;
  step_length = Norm2(mem_step);

  fvalue = nlp->getF();

  grad   = nlp->getGrad();
  gnorm  = Norm2(grad);
  
  if (debug_) {
    *optout << "\n\t xc \t\t\t   grad \t\t   step\n";
    for(int i=1; i<=n; i++)
      *optout << i <<  e(xc(i),24,16) << e(grad(i),24,16) 
	   << e(mem_step(i),24,16) << "\n";
    *optout << "\nHessian";
    Print(Hessian);

//  Compute eigenvalues of Hessian

    DiagonalMatrix D(n);
    EigenValues(Hessian, D);
    *optout << "\nEigenvalues of Hessian";
    Print(D);

    *optout << "\n***************************************";
    *optout << "***************************************\n";


  }
//  Iteration summary
// 
  if(step_type >= 0){
  *optout 
    << d(iter,5)  << " " << e(fvalue,12,4) << " " << e(gnorm,12,4) << " "
    << e(step_length,12,4) << "  " << steps[step_type] << " " 
    << d(fcn_evals,5) << " " << d(grad_evals,5) << endl;
  }
  else{
  *optout 
    << d(iter,5)  << " " << e(fvalue,12,4) << " " << e(gnorm,12,4) << " "
    << e(step_length,12,4) << "  " << " "  << " " 
    << d(fcn_evals,5) << " " << d(grad_evals,5) << endl;
  }
}

ColumnVector OptBCNewtonLike::defaultComputeSearch(SymmetricMatrix& H)
{  
  NLP1* nlp = nlprob();
  int n     = nlp->getDim();

  ColumnVector sk(n);
  LowerTriangularMatrix L(n);

  L = MCholesky(H);
  sk = -(L.t().i()*(L.i()*gprev));
  return sk;

}

//------------------------------------------------------------------------
//
// Now all the other functions that can be generalized
// to all constrained Newton cases
//
// checkConvg
// CheckDeriv
// computeStep
// InitOpt
// initTrustRegionSize
// Optimize
//------------------------------------------------------------------------

int OptBCNewtonLike::checkConvg() // Check convergence
{
  NLP1* nlp = nlprob();
  ColumnVector xc(nlp->getXc());

// Test 1. step tolerance 

  double step_tol = tol.getStepTol();
  double snorm = stepTolNorm();
  double xnorm =  Norm2(xc);
  double stol  = step_tol*max(1.0,xnorm);
  if (snorm  <= stol) {
    strcpy(mesg,"Step tolerance test passed");
    *optout << "checkConvg: snorm = " << e(snorm,12,4) 
      << "  stol = " << e(stol,12,4) << "\n";
    return 1;
  }
  
// Test 2. function tolerance
  double ftol   = tol.getFTol();
  double fvalue = nlp->getF();
  double rftol  = ftol*max(1.0,fabs(fvalue));
  Real deltaf   = fprev - fvalue;
  if (deltaf <= rftol) {
    strcpy(mesg,"Function tolerance test passed");
    *optout << "checkConvg: deltaf = " << e(deltaf,12,4) 
         << "  ftol = " << e(ftol,12,4) << "\n";
    return 2;
  }
  

// Test 3. gradient tolerance 

  ColumnVector grad(nlp->getGrad());
  double gtol  = tol.getGTol();
  double rgtol = gtol*max(1.0,fabs(fvalue));
  double gnorm = Norm2(grad);
  if (gnorm <= rgtol) {
    strcpy(mesg,"Gradient tolerance test passed");
    *optout << "checkConvg: gnorm = " << e(gnorm,12,4) 
      << "  gtol = " << e(rgtol, 12,4) << "\n";
    return 3;
  }
  

// Test 4. absolute gradient tolerance 

  if (gnorm <= gtol) {
    strcpy(mesg,"Gradient tolerance test passed");
    *optout << "checkConvg: gnorm = " << e(gnorm,12,4) 
      << "  gtol = " << e(gtol, 12,4) << "\n";
    return 4;
  }
  
  // Nothing to report 
  return 0;

}

//---------------------------------------------------------------------------- 

int OptBCNewtonLike::checkAnalyticFDGrad() 
{
  int i, n = dim, retcode = GOOD; 
  double eta, gnorm, maxerr, third;
  ColumnVector error(n), fd_grad(n), grad(n);

  Real mcheps = FloatingPointPrecision::Epsilon();

  NLP1* nlp = nlprob();
  ColumnVector xc = nlp->getXc();
  double fx = nlp->getF();
  SpecOption tmpSpec = nlp->getSpecOption();

  nlp->setSpecOption(NoSpec);
  fd_grad   = nlp->FDGrad(sx, xc, fx, fd_grad); // Evaluate gradient using finite differences
  nlp->setSpecOption(tmpSpec);
  grad      = nlp->getGrad();  	// Now use the analytical functions
  third     = 0.33333;
  gnorm     = grad.NormInfinity();
  eta       = pow(mcheps,third)*max(1.0,gnorm);

  if(debug_){
     *optout << "Check_Deriv: Checking gradients versus finite-differences\n";
     *optout << "    i    gradient     fd grad       error\n";
     for (i=1; i<=n; i++) {
        error(i) = fabs(grad(i)-fd_grad(i));
        *optout << d(i,5) << e(grad(i),12,4) 
               << e(fd_grad(i),12,4) << e(error(i),12,4) << "\n";
      }
  }
  maxerr = error.NormInfinity();           
  if(debug_){
     *optout << "maxerror = " << e(maxerr, 12,4) 
            << "tolerance =  " << e(eta, 12,4) << "\n";
  }
  if (maxerr > eta) retcode = BAD;
  return retcode;
}


int OptBCNewtonLike::checkDeriv() // Check the analytic gradient with FD gradient
{return GOOD;}

//---------------------------------------------------------------------------- 
// 
// Compute a step along the direction sk using either a line search
// or model trust-region approach
//
//---------------------------------------------------------------------------- 
int OptBCNewtonLike::computeStep(ColumnVector sk)
{
  NLP1* nlp = nlprob();
  real stp_length = 1.0;
  real stptmp;
  real lstol  = tol.getLSTol();
//  real xtol   = tol.getStepTol();
//  real gtol   = tol.getGTol();
  real stpmax = tol.getMaxStep();
  real stpmin = tol.getMinStep();
  int  step_type;
  int  itnmax = tol.getMaxBacktrackIter();

  if (trace) *optout << class_name << ": computeStep\n";

//
// Compute the maximum step allowed
//
  stptmp = computeMaxStep(sk);
  stpmax = min(stpmax, stptmp);

  if (strategy == TrustRegion) {
    SymmetricMatrix H = Hessian;
    step_type = trustregion(nlp, optout, H, sk, sx, TR_size, stp_length, 
			    stpmax, stpmin);
  }
  else if (strategy == LineSearch) {
    step_type = linesearch(nlp, optout, sk, sx, &stp_length, stpmax, stpmin,
			   itnmax, lstol);
  }
  else if (strategy == TrustPDS) {
    SymmetricMatrix H = Hessian;
    step_type = trustpds(nlp, optout, H, sk, sx, TR_size, stp_length, 
			    stpmax, stpmin, searchSize);
  }
  else {
    return(-1);
  }
  
  if (step_type < 0) {
    setMesg("OptBCNewtonLike: Step does not satisfy sufficient decrease condition");
    ret_code = -1;
    setReturnCode(ret_code);
    return(ret_code);
  }
  fcn_evals   = nlp->getFevals();
  grad_evals  = nlp->getGevals();
  step_length = stp_length;
  return(step_type);
}

void OptBCNewtonLike::initHessian()
{ 
  int i;
  NLP1* nlp = nlprob();
  int ndim = nlp->getDim();

  if (WarmStart) {
    *optout << "OptBCNewton::initHessian: Warm Start specified\n";
  }
  else {
    Real typx, xmax, gnorm;
    ColumnVector grad(ndim), xc(ndim);
    xc     = nlp->getXc();
    grad   = nlp->getGrad();
    gnorm = Norm2(grad);
    DiagonalMatrix D(ndim);

    // Initialize xmax, typx and D to default values
    xmax   = -1.e30; typx   =  1.0; D      =  1.0;

    for (i=1; i <= ndim; i++) xmax = max(xmax,xc(i));
    if(xmax != 0.0) typx = xmax;
    if(gnorm!= 0.0) D    = gnorm/typx;
    if (debug_) {
      *optout << "OptBCNewton::initHessian: gnorm0 = " << gnorm
	<< "  typx = " << typx << "\n";
    }
    Hessian = 0.0;
    for (i=1; i <= ndim; i++) Hessian(i,i) = D(i);
   }
}

void OptBCNewtonLike::initOpt()
{
  double gnorm;
  NLP1* nlp = nlprob();
  int n = nlp->getDim();

  time_t t;
  char *c;

// get date and print out header

  t = time(NULL);
  c = asctime(localtime(&t));
  *optout << "**********************************************************\n";
  *optout << "OPT++ version " << OPT_GLOBALS::OPT_VERSION << "\n";
  *optout << "Job run at " << c << "\n";
  copyright();
  *optout << "**********************************************************\n";

//
// Read in OPT++ input file if it exists
// Be aware that anything in the input file will
// override any variable set so far
//

  nlp->initFcn();
  readOptInput();

  ret_code = 0;

  if(nlp->hasConstraints()){
    CompoundConstraint* constraints = nlp->getConstraints();
    ColumnVector xstart = nlp->getXc();
    double feas_tol = tol.getCTol();
    bool feasible = constraints->amIFeasible(xstart, feas_tol);
    if (!feasible) {
      *optout << "OptBCNewtonLike WARNING:  Initial guess not feasible.\n"
	      << "BCNewton may be unable to make progress." << endl;
    }
  }

  if (ret_code == 0) {
    // evaluate Function, gradient and compute initial Hessian

    nlp->eval();

    xprev = nlp->getXc();
    fprev = nlp->getF();
    gprev = nlp->getGrad();
    gnorm = Norm2(gprev);
  
    //  SymmetricMatrix Hk(n);
    //  Hessian = updateH(Hk,0);

    initHessian();
    setFcnScale(fprev);

    nlp->fPrintState(optout, "Initial state");

    if(strategy == TrustRegion) {
      *optout << "\n\t\t" << method << " Method with Trust Regions\n";
      TR_size = getTRSize();
      if (TR_size == 0.0) TR_size = getGradMult()*gnorm;
      *optout << "\t\t Initial Trust Region = " << e(TR_size,12,4) << "\n";
    }
    else if(strategy == TrustPDS) {
      *optout << "\n\t\t" << method << " Method with Trust Region / PDS\n";
      TR_size = getTRSize();
      if (TR_size == 0.0) TR_size = getGradMult()*gnorm;
      *optout << "\t\t Initial Trust Region = " << e(TR_size,12,4) << "\n";
    }
    else  
      *optout << "\n\t\t" << method << " Method with Line Search\n";

    *optout << "\n  Iter      F(x)       ||grad||     "
	    << "||step||      f/g\n\n"
	    << d(0,5) << " " << e(fprev,12,4) << " " << e(gnorm,12,4) << endl;
    if (debug_) {
      nlp->fPrintState(optout, "BCNewtonLike: Initial Guess");
      *optout << "xc, grad, step\n";
      for(int i=1; i<=n; i++)
	*optout << i << e(xprev(i),24,16) << e(gprev(i),24,16) << "\n";
      Print(Hessian);
    }
  }
}

double OptBCNewtonLike::initTrustRegionSize() const
{ 
  double init_tr;
//
// return minimum of 100||x||, tolerance default, or Maxstep
//
  init_tr = 100.0*Norm2(xprev);
  init_tr = min(init_tr, tol.getTRSize());    
  init_tr = min(init_tr, tol.getMaxStep());    

  return init_tr;
}

void OptBCNewtonLike::optimize()
//---------------------------------------------------------------------------- 
//
// Given a nonlinear operator nlp find the minimizer using a
// Newton-like method
//
//---------------------------------------------------------------------------- 
{
  int k;
  int convgd = 0;
  int maxiter, maxfev, myfevals, fevals, step_type;

// Allocate local vectors 

  int n = dim;
  ColumnVector sk(n);
  SymmetricMatrix Hk(n);
  NLP1* nlp = nlprob();

// Initialize iteration
// Evaluate Function, Gradient, and Hessian

  initOpt();

  if (ret_code == 0) {
    Hk = Hessian;
    maxiter = tol.getMaxIter();
    maxfev  = tol.getMaxFeval();

    // Check for convergence. Need to take into account that this is the
    // zeroth iteration
    //  convgd = objfcn.Check_Convg(tol,stat);
    //  if (convgd > 0) {
    //    stat.ret_code = convgd;
    //    return;
    //  }
  
    for (k=1; k <= maxiter; k++) {

      iter_taken = k;
      if (debug_)
	*optout << " **** OptBCNewtonLike : iteration count = " << k << "\n";

      //  Compute search direction

      sk = computeSearch(Hk);

      //  attempt to take a step in the direction sk from the current point. 
      //  The default method is to use a trust region

      if ((step_type = computeStep(sk)) >= 0) {
	acceptStep(k, step_type);
	convgd    = checkConvg();
        m_nconvgd = convgd;
      }

      //  Update Constraints

      ret_code = updateConstraints(step_type);

      //  Error checking 

      if (ret_code <= 0) { // constraints have not been modified
	if (step_type<0 && convgd==0) {//not converged and cannot take a step 
	  ret_code = step_type;
          setReturnCode(ret_code);
	  *optout << "OptBCNewtonLike : cannot take a step \n";
	  return;
	} else if (convgd > 0) { // converged
          setReturnCode(convgd);
	  *optout << "OptBCNewtonLike : convergence achieved. \n";
	  return;
	}
      }

      myfevals = nlp->getFevals();

#ifdef WITH_MPI

      char buffer[MPI_MAX_ERROR_STRING];
      int error, resultlen, flag;

      error = MPI_Initialized(&flag);
      if (error != MPI_SUCCESS)
	{
	  MPI_Error_string(error, buffer, &resultlen);
	  printf("\nOptNewtonLike: MPI Error - %s\n", buffer);
	  strcpy(mesg, "OptNewtonLike: error returned by MPI_Initialized\n");
	  ret_code = -14;
          setReturnCode(ret_code);
	}

      if (flag == 0) {
	fevals = myfevals;
      }
      else{
	error = MPI_Allreduce(&myfevals, &fevals, 1, MPI_INT, MPI_MAX,
			      MPI_COMM_WORLD);
	if (error != MPI_SUCCESS)
	  {
	    MPI_Error_string(error, buffer, &resultlen);
	    printf("\nOptNewtonLike: MPI Error - %s\n", buffer);
	    strcpy(mesg, "OptNewtonLike: error returned by MPI_Allreduce\n");
	    ret_code = -15;
            setReturnCode(ret_code);
	  }
      }

#else

      fevals = myfevals;

#endif

      if (fevals > maxfev) break;

      //  if not converged, update the Hessian 

      if (convgd <= 0 || ret_code > 0) {
	Hessian = updateH(Hk,k);
	Hk = Hessian;
	xprev = nlp->getXc();
	fprev = nlp->getF();
	gprev = nlp->getGrad();
      }
    }

    setMesg("OptBCNewtonLike: Maximum number of iterations or fevals");
    ret_code = -4;
    setReturnCode(ret_code);
  }
}

void OptBCNewtonLike::printStatus(char *s) // set Message
{
  NLP1* nlp = nlprob();

  *optout << "\n\n=========  " << s << "  ===========\n\n";
  *optout << "Optimization method       = " << method << "\n";
  *optout << "Dimension of the problem  = " << nlp->getDim()  << "\n";
  *optout << "No. of bound constraints  = " << nlp->getDim()  << "\n";
  *optout << "Return code               = " << ret_code << " ("
       << mesg << ")\n";
  *optout << "No. iterations taken      = " << iter_taken  << "\n";
  *optout << "No. function evaluations  = " << nlp->getFevals() << "\n";
  *optout << "No. gradient evaluations  = " << nlp->getGevals() << "\n";

  if (debug_) {
    Print(Hessian);
//  Compute eigenvalues of Hessian

    DiagonalMatrix D;
    SVD(Hessian, D);
    *optout << "\nEigenvalues of Hessian";
    Print(D);
  }

  nlp->fPrintState(optout, s);
  tol.printTol(optout);

}
void OptBCNewtonLike::readOptInput() // Read opt.input file if it exists
{
  NLP1* nlp = nlprob();

/* A VERY simple routine for reading the optimization parameters
 * We should really make this more general, but as a first pass this
 * will have to do.
 * 
 * The input file should be of the form keyword = value
 * where keyword is one of the following
 * 
 * search      = trustregion
 * diff_option = forward
 * max_iter    = 100
 * max_feval   = 1000
 * grad_tol    = 1.e-6
 * fcn_tol     = 1.e-9
 * max_step    = 100.0
 * fcn_accrcy  = 1.e-9
 *
 */
  
  int  index, max_iter, max_feval;
  real grad_tol,  fcn_tol, max_step, fcn_accrcy;

  char token[80], ignore[80], equals[1];
//
// Keywords allowed
//
  string keyword;
  string cdiff_option("diff_option");
  string cfcn_accrcy("fcn_accrcy");
  string cfcn_tol("fcn_tol");
  string cgrad_tol("grad_tol");
  string cmaxfeval("maxfeval");
  string cmaxiter("max_iter");
  string cmax_step("max_step");
  string csearch("search");

  string diff_option;
  string search;
  SearchStrategy s = TrustRegion;

  int keyword_count = 0;

// 
// default name of input file
//
  char *opt_input  = {"opt.input"};

//
// Open opt.input file and check to see if we succeeded
//

  ifstream optin(opt_input);
  if (!optin.rdbuf()->is_open()) {
    *optout << "readOptInput: No opt.input file found\n";
    *optout << "readOptInput: default values will be used\n";
    return;
  }

  *optout << "readOptInput: Reading opt.input file\n";

  optin >> token;

  while (!optin.eof()) {

    keyword = token;
    keyword_count++;

//debug    *optout << "keyword = " << keyword << "\n";

    if (keyword == cdiff_option) {

      optin >> equals >> token;
      diff_option = token;

      if ( diff_option == "forward")
	nlp->setDerivOption(ForwardDiff);
      else if ( diff_option == "backward")
	nlp->setDerivOption(BackwardDiff);
      else if ( diff_option == "central")
	nlp->setDerivOption(CentralDiff);
    }    
    else if (keyword == cfcn_accrcy) {
      //optin >> equals >> fcn_accrcy;
      //nlp->setFcnAccrcy(fcn_accrcy);
      optin >> equals >> index >> fcn_accrcy;
      nlp->setFcnAccrcy(index, fcn_accrcy);
    }    
    else if (keyword == cfcn_tol) {
      optin >> equals >> fcn_tol;
      setFcnTol(fcn_tol);
    }    
    else if (keyword == cgrad_tol) {
      optin >> equals >> grad_tol;
      setGradTol(grad_tol);
    }    
    else if (keyword == cmaxfeval) {
      optin >> equals >> max_feval;
      setMaxFeval(max_feval);
    }    
    else if (keyword == cmaxiter) {
      optin >> equals >> max_iter;
      setMaxIter(max_iter);
    }
    else if (keyword == cmax_step) {
      optin >> equals >> max_step;
      setMaxStep(max_step);
    }
    else if (keyword == csearch) {
      optin >> equals >> token;
      search = token;
      if ( search == "trustregion")
	s = TrustRegion;
      else if ( search == "linesearch")
	s = LineSearch;
      setSearchStrategy(s);
    }
    else {
      *optout << "Unrecognized keyword '" << keyword << "'. "
	<< "Skipping the rest of this line\n";
      optin.getline(ignore, sizeof(ignore));
    }
  optin >> token;
  }

  *optout << "\n\n======  Summary of input file  ======\n\n";

  *optout << csearch      << " = " << search << "\n";
  *optout << cdiff_option << " = " << diff_option << "\n";
  *optout << cmaxiter     << " = " << max_iter << "\n";
  *optout << cmaxfeval    << " = " << max_feval << "\n";
  *optout << cgrad_tol    << " = " << grad_tol << "\n";
  *optout << cfcn_tol     << " = " << fcn_tol << "\n";
  *optout << cmax_step    << " = " << max_step << "\n";
  ColumnVector fcnacc = nlp->getFcnAccrcy();
  for(int i = 1; i<= fcnacc.Nrows(); i++)
     *optout << cfcn_accrcy  << " = " << fcnacc(i) << "\n";

  tol.printTol(optout);

}

} // namespace OPTPP
