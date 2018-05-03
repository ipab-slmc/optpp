//------------------------------------------------------------------------
// Copyright (C) 1996: 
// J.C. Meza and Charles Tong
// Sandia National Laboratories
// meza@california.sandia.gov
//------------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#ifdef HAVE_STD
#include <cstring>
#include <ctime>
#else
#include <string.h>
#include <time.h>
#endif

#include <string>

using namespace std;

#include "OptBCEllipsoid.h"
#include "precisio.h"
#include "cblas.h"
#include "ioformat.h"

using NEWMAT::ColumnVector;
using NEWMAT::Matrix;
using NEWMAT::SymmetricMatrix;

namespace OPTPP {

//------------------------------------------------------------------------
// Bound Constrained Ellipsoid Class functions
//------------------------------------------------------------------------
// acceptStep
// checkConvg
// initOpt
// reset
// optimize
// printStatus
// readOptInput
// computeFeasibility
// computeConstraintSubgradient
// infeasibilityStep
// halfSpacStep
// computeGamma
//------------------------------------------------------------------------

//------------------------------------------------------------------------
// acceptStep - Print out iteration summary and anything else
//              you might want to do before the next iteration
//------------------------------------------------------------------------
void OptBCEllipsoid::acceptStep(int iter, int step_type)
{
  NLP1        *nlp = nlprob();
  int          i, n = nlp->getDim(), grad_evals;
  ColumnVector xc(n), grad(n);
  double       fvalue;

  xc          = nlp->getXc();
  mem_step    = xc - xprev;
  step_length = Norm2(mem_step);
  xprev       = xc;
  grad        = nlp->evalG();
  fvalue      = nlp->getF();

  if (debug_) {
    *optout << "\n\t xc \t\t\t grad \t\t\t step\n";
    for(i=1; i<=n; i++)
      *optout << i <<  e(xc(i),24,16) << e(grad(i),24,16)
	   << e(mem_step(i),24,16) << "\n";
  }

  fcn_evals  = nlp->getFevals();  
  grad_evals = nlp->getGevals();

  *optout << d(iter,5)      <<  e(fvalue,12,4)  << e(step_length,12,4) 
         << d(fcn_evals,5) << d(grad_evals,5) << "\n" << flush;
}

//------------------------------------------------------------------------
// checkConvg - check whether the distance between upper and lower bounds
//              is less than a prescribed threshold.
//------------------------------------------------------------------------
int OptBCEllipsoid::checkConvg() 
{
  NLP1         *nlp = nlprob();
  ColumnVector xc(nlp->getXc());
  double       fvalue = nlp->getF();
  double       ftol = tol.getFTol();
  double       delta;

  fval_upbound  = min(fval_upbound,fvalue);
  delta = fabs(fval_upbound - fval_lowbound);
  if (delta <= ftol) {
    strcpy(mesg,"Function tolerance test passed");
    ret_code = 2;
    setReturnCode(ret_code);
    return 1;
  } else 
    return 0;
}

//------------------------------------------------------------------------
// Initialization subroutine
//------------------------------------------------------------------------
void OptBCEllipsoid::initOpt()
{
  NLP1         *nlp = nlprob();
  int          i,n = nlp->getDim();
  double       dtmp=0.0;
  time_t       t;
  char         *c;

  // Get date and print out header
  t = time(NULL);
  c = asctime(localtime(&t));
  *optout << "**********************************************************\n";
  *optout << "OPT++ version " << OPT_GLOBALS::OPT_VERSION << "\n";
  *optout << "Job run at " << c << "\n";
  copyright();
  *optout << "**********************************************************\n";

  // Read in OPT++ input file if it exists. Be aware that anything in 
  // the input file will override any variable set so far

  nlp->initFcn();
  ColumnVector xc = nlp->getXc();
  readOptInput();

  ret_code = 0;

  if(nlp->hasConstraints()){
    CompoundConstraint* constraints = nlp->getConstraints();
    ColumnVector xstart = nlp->getXc();
    double feas_tol = tol.getCTol();
    bool feasible = constraints->amIFeasible(xstart, feas_tol);
    if (!feasible) {
      *optout << "OptBCEllipsoid WARNING:  Initial guess not feasible.\n"
	      << "Ellipsoid may be unable to make progress." << endl;
    }
  }

  if (ret_code == 0) {
    nlp->evalF();

    // if initial radius of ellipsoid is not specified, set it to something
    if (initial_radius < 0.0e0) {
      for (i=1; i<=n; i++) dtmp = max(dtmp, fabs(xc(i))); 
      initial_radius = 1.0e1 * dtmp + 1.0e5; 
    }
  
    *optout << "\n  Iter      F(x)   Steplength   "
	    << "fevals    gevals\n\n";
    
    if(debug_)
      *optout << "Radius of initial ellipsoid = " << initial_radius << "\n";  
  }
}

//---------------------------------------------------------------------------- 
// Reset optimization parameters 
//---------------------------------------------------------------------------- 
void OptBCEllipsoid::reset()
{
   NLP1* nlp = nlprob();
   int   n   = nlp->getDim();
   nlp->reset();
   OptimizeClass::defaultReset(n);
   initial_radius = -1.0e0;
   xscal_flag = deepcutflag = 0;
}

//---------------------------------------------------------------------------- 
// Given a nonlinear operator nlp find the minimizer using a
//---------------------------------------------------------------------------- 
void OptBCEllipsoid::optimize()
{
  NLP1*  nlp = nlprob();
  int convgd = 0;
  int       i,n = nlp->getDim(), step_type;
  ColumnVector xc = nlp->getXc(), xscale = getXScale(), xs(n);
  double          psi, dtmp;

  // Read input file and output initial message to file 
  initOpt();

  if (ret_code == 0) {
    iter_taken = 0;

    // Initialize convergence test variables 
    fval_lowbound = -FLT_MAX;
    fval_upbound  = FLT_MAX;

    // Initialize the A matrix
    SymmetricMatrix A(n);
    if (xscal_flag != 1) {xscale.ReSize(n); xscale = 1.0;}
    dtmp = initial_radius * initial_radius;
    A = 0.0;
    for (i=1; i<=n; i++) A(i,i) = dtmp / (xscale(i) * xscale(i));

    // scale the initial guess (if scaling is desired)
    for (i=1; i<=n; i++) xc(i) = xc(i) / xscale(i);

    // declare other vectors used in the iterations
    ColumnVector ghk(n), aghk(n), aghkscal(n);

    // assuming that the function has been evaluated, get the value
    fprev = nlp->getF();

    // Move the initial guess into the feasible region, if needed
    for (i=1; i<=n; i++) xs(i) = xc(i) * xscale(i);
    psi = computeFeasibility(xs);
    if (psi > 0.0) infeasibilityStep(xc,A,psi);

    while (convgd == 0) { 

      iter_taken++;
      //*optout << " **** OptBCEllipsoid : iteration count = " 
      //	 << iter_taken << "\n";

      // put away the last solution to prepare for current iteration 
      xprev = nlp->getXc();

      // perform one ellipsoid iteration (xc,A changed upon return)
      step_type = halfSpaceStep(xc,A,psi);

      // if the next solution is infeasible, do deep cut
      if (step_type == -1) infeasibilityStep(xc,A,psi);

      // update solution and update function value
      for (i=1; i<=n; i++) xs(i) = xc(i) * xscale(i);
      nlp->setX(xs);
      fprev = nlp->evalF();

      // check convergence
      acceptStep(iter_taken, 0);
      convgd = checkConvg();

      // debug  information - volume of ellipsoid
      //logdeterminant = A.LogDeterminant();
      //dtmp = 0.5 * n + 1.0;
      //determinant = sqrt(logdeterminant.Value()) * pow(pi,dtmp-1.0) 
      //					       / ComputeGamma(dtmp);
      //*optout << "Volume of current ellipsoid = " << determinant << "\n";
    }
  }
}

//---------------------------------------------------------------------------- 
// Print message to the output file
//---------------------------------------------------------------------------- 
void OptBCEllipsoid::printStatus(char *s) 
{
  NLP1         *nlp = nlprob();

  if (deepcutflag == 1)
    strcpy(method,"The Ellipsoid method with deep cut");
  else
    strcpy(method,"The Ellipsoid method ");

  *optout << "\n\n=========  " << s << "  ===========\n\n";
  *optout << "Optimization method       = " << method << "\n";
  *optout << "Dimension of the problem  = " << nlp->getDim()  << "\n";
  *optout << "Return code               = " << ret_code << " ("
       << mesg << ")\n";
  *optout << "No. iterations taken      = " << iter_taken  << "\n";
  *optout << "No. function evaluations  = " << nlp->getFevals() << "\n";
  *optout << "No. gradient evaluations  = " << nlp->getGevals() << "\n";

  tol.printTol(optout);
  nlp->fPrintState(optout, s);
}

//---------------------------------------------------------------------------- 
// read an input file called opt.input
//  - A VERY simple routine for reading the optimization parameters
//    We should really make this more general, but as a first pass this
//    will have to do.  The input file should be of the form keyword = value
//    where keyword is one of the following
//     max_iter    = 100
//     max_feval   = 1000
//     grad_tol    = 1.e-6
//     fcn_tol     = 1.e-9
//     fcn_accrcy  = 1.e-9
//---------------------------------------------------------------------------- 
void OptBCEllipsoid::readOptInput() 
{
  NLP1         *nlp = nlprob();
  int   index, max_iter, max_feval;
  real  grad_tol,  fcn_tol, fcn_accrcy;
  char  token[80], ignore[80], equals[1];

  // Keywords allowed
  string keyword;
  string cfcn_accrcy("fcn_accrcy");
  string cfcn_tol("fcn_tol");
  string cgrad_tol("grad_tol");
  string cmaxfeval("maxfeval");
  string cmaxiter("maxiter");

  // Default name of input file
  char *opt_input  = {"opt.input"};

  // Open opt.input file and check to see if we succeeded
  ifstream optin(opt_input);
  if (!optin.rdbuf()->is_open()) {
    *optout << "ReadOptInput: No opt.input file found\n";
    *optout << "ReadOptInput: Default values will be used\n";
    return;
  }

  *optout << "ReadOptInput: Reading opt.input file\n";

  fcn_tol    = tol.getFTol();
  grad_tol   = tol.getGTol();
  max_feval  = tol.getMaxFeval();
  max_iter   = tol.getMaxIter();
  while ((optin >> token)) {

    keyword = token;

    if (keyword == cfcn_accrcy) {
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
    else {
      *optout << "Unrecognized keyword '" << keyword << "'. "
	<< "Skipping the rest of this line\n";
      optin.getline(ignore, sizeof(ignore));
    }
  }

  *optout << "\n\n======  Summary of input file  ======\n\n";

  *optout << cmaxiter     << " = " << max_iter << "\n";
  *optout << cmaxfeval    << " = " << max_feval << "\n";
  *optout << cgrad_tol    << " = " << grad_tol << "\n";
  *optout << cfcn_tol     << " = " << fcn_tol << "\n";
  ColumnVector fcnacc     = nlp->getFcnAccrcy();
  for(int i = 1; i<= fcnacc.Nrows(); i++)
     *optout << cfcn_accrcy  << " = " << fcnacc(i) << "\n";

  tol.printTol(optout);
}

//------------------------------------------------------------------------
// ComputeFeasibility 
//------------------------------------------------------------------------
double OptBCEllipsoid::computeFeasibility(ColumnVector& vec)
{  
  NLP1       *nlp = nlprob();
  int          i, n; 
  ColumnVector upper, lower;
  double       feas = -FLT_MAX, fdiff;

  n     = nlp->getDim();
  upper = nlp->getConstraints()->getUpper();
  lower = nlp->getConstraints()->getLower();

  for (i=1; i<=n; i++) {
    if (lower(i) != -FLT_MAX) {
      fdiff = lower(i) - vec(i);
      feas = max(fdiff, feas);
    }
    /* 02/01/01 PJW changed the comparison from -FLT_MAX to FLT_MAX */
    if (upper(i) !=  FLT_MAX) {
      fdiff = vec(i) - upper(i);
      feas = max(fdiff, feas);
    }
  }
  return feas;
}

//------------------------------------------------------------------------
// ComputeConstraintSubgradient - pick the row corresponding to the most
//                                infeasible constraint
//------------------------------------------------------------------------
ColumnVector OptBCEllipsoid::computeConstraintSubgradient(ColumnVector& x)
{  
  NLP1        *nlp = nlprob();
  int           i, n, index;
  double        dtmp, dblemax = -FLT_MAX;
  ColumnVector  subgrad, upper, lower;  

  n     = nlp->getDim();
  upper = nlp->getConstraints()->getUpper();
  lower = nlp->getConstraints()->getLower();
  subgrad.ReSize(n);

  // find the constraint that gives the most infeasibility
  for (i=1; i<=n; i++) {
    dtmp = - x(i) + lower(i);
    if (dtmp > dblemax) {dblemax = dtmp; index = i;}
    dtmp =   x(i) - upper(i);
    if (dtmp > dblemax) {dblemax = dtmp; index = n + i;}
  }

  // set the subgradient 
  subgrad = 0.0;
  if (index <= n) subgrad(index) = -1.0;
  else            subgrad(index-n) = 1.0;
  return subgrad;
}

//------------------------------------------------------------------------
// A step to be taken in case the current x is infeasible
//------------------------------------------------------------------------
int OptBCEllipsoid::infeasibilityStep(ColumnVector& x,
	                              SymmetricMatrix &A, double &feas)
{  
  NLP1       *nlp = nlprob();
  int          i, n;
  ColumnVector ghk, aghk, aghkscal, xscale=getXScale(), xs;
  double       alpha, scale, fact, fact2, psi;
  Matrix       Atmp;

  n = nlp->getDim();
  ghk.ReSize(n);
  aghk.ReSize(n);
  aghkscal.ReSize(n);
  xs.ReSize(n);
  Atmp.ReSize(n,n);

  psi = feas;
  while (psi > 0.0) {
    //*optout << "===>InfeasibilityStep, psi = " << psi << "\n";

    // calculate a constraint subgradient g
    for (i=1; i<=n; i++) xs(i) = x(i) * xscale(i);
    ghk  = computeConstraintSubgradient(xs);
    for (i=1; i<=n; i++) ghk(i) = ghk(i) * xscale(i);

    // calculate sqrt(g'Ag)
    aghk = A * ghk;
    scale = Dot(ghk,aghk);
    if (scale <= 0.0) {
      *optout << "Error in OptBCEllipsoid : sqrt of negative number.\n";
      exit(-1);
    }
    scale = sqrt(scale);
    if (psi > scale) {
      *optout << "Error in OptBCEllipsoid : feasible set is empty.\n";
      exit(-1);
    }

    // scale the vector (A*g)
    aghkscal = aghk * (1.0/scale); 

    // depending on whether deep-cut is desired, set the alpha
    if (deepcutflag == 1) alpha = psi / scale;
    else                  alpha = 0.0;

    // update x
    fact2 = (1.0 + alpha * n) / (1.0 + n);
    x = x - aghkscal * fact2;

    // update A
    fact = (n * n) / (n * n - 1.0e0) * (1.0 - alpha * alpha);
    fact2 = 2.0 * fact2 / (1.0 + alpha); 
    Atmp = A;
    Atmp = (Atmp - (aghkscal * aghkscal.t()) * fact2) * fact;
    A << Atmp;

    // compute the feasibility measure of the updated x
    for (i=1; i<=n; i++) xs(i) = x(i) * xscale(i);
    psi = computeFeasibility(xs);
  }
  return 0;
}

//------------------------------------------------------------------------
// Deep cut step for upper bound 
//------------------------------------------------------------------------
int OptBCEllipsoid::halfSpaceStep(ColumnVector& x,
	                          SymmetricMatrix &A, double &psi)
{
  NLP1*      nlp = nlprob();
  int          i, n, numstep = 0;
  ColumnVector ghk, aghk, aghkscal, xs, xscale=getXScale();
  double       alpha, scale, fact, fact2, gnormA;
  Matrix       Atmp;

  n = nlp->getDim();
  ghk.ReSize(n);
  aghk.ReSize(n);
  aghkscal.ReSize(n);
  xs.ReSize(n);
  Atmp.ReSize(n,n);

  while (fprev > fval_upbound || numstep == 0) {
    numstep++;
    //*optout << "===>HalfSpaceStep, fprev, fval_up = " << fprev << " "
    //	   << fval_upbound << "\n";

    // compute the gradient (incorporate scaling) 
    for (i=1; i<=n; i++) xs(i) = x(i) * xscale(i);
    ghk = nlp->evalG(xs); 
    for (i=1; i<=n; i++) ghk(i) = ghk(i) * xscale(i);

    // compute A times the gradient and compute sqrt(g'Ag)
    aghk = A * ghk;
    scale = Dot(ghk,aghk);
    if (scale <= 0.0) {
      *optout << "Error in OptBCEllipsoid : sqrt of negative number.\n";
      exit(-1);
    }
    gnormA = scale = sqrt(scale);

    // update the function value lower bound (refer to formula in the text)
    fval_lowbound = max(fprev-gnormA,fval_lowbound);

    // compute a scaled version of A * tilde(g) where tilde(g)=g/scale
    aghkscal = aghk * (1.0/scale); 

    // if using deep cut and the current function value is above the upper 
    // bound, then do deep cut by setting alpha != 0
    if (fprev>fval_upbound && deepcutflag==1) alpha=(fprev-fval_upbound)/scale;
    else                                      alpha = 0.0;

    // update x
    fact2 = (1.0 + alpha * n) / (1.0 + n);
    x = x - aghkscal * fact2;

    // update A
    fact = (n * n) / (n * n - 1.0e0) * (1.0 - alpha * alpha);
    fact2 = 2.0 * fact2 / (1.0 + alpha); 
    Atmp = A;
    Atmp = (Atmp - (aghkscal * aghkscal.t()) * fact2) * fact;
    A << Atmp;

    // compute the feasibility measure of the updated x
    for (i=1; i<=n; i++) xs(i) = x(i) * xscale(i);
    psi = computeFeasibility(xs);

    // if the solution is not feasible, return with an error flag to
    // indicate that an infeasibility step is needed.
    if (psi > 0.0) return -1;

    // if there is no deep cut (for upper bound), do not loop 
    if (deepcutflag == 0) return 0;

    // otherwise, update the previous function value and repeat
    if (fprev > fval_upbound) {
      fprev = nlp->evalF(xs); 
    }
  }
  return 0;
}

//------------------------------------------------------------------------
// Calculation of gamma function given x
//    gamma(x) = int from 0 to infinity [t^x exp(-t)] dt
//------------------------------------------------------------------------
double OptBCEllipsoid::computeGamma(double x)
{
  int    i, n=1000;
  double t, a=0.0, b=100.0, h, gamma;

  if (x < 1.0) {
    *optout << "A gamma function of <1 is not supported.\n";
    exit(-1);
  }
  h   = (b - a) / n;
  gamma = 0.0;
  for (i=1; i<=n; i++) {
    t = i * h;
    gamma = gamma + pow(t,x) * exp(-t) * h;
  }
  return gamma;
}

} // namespace OPTPP
