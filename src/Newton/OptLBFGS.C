//------------------------------------------------------------------------
// Copyright (C) 2003: 
// R.A.Oliva, Lawrence Berkeley National Laboratory.
// raoliva@lbl.gov
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

#include "OptLBFGS.h"
#include "precisio.h"
#include "cblas.h"
#include "ioformat.h"

using namespace std;

using NEWMAT::Real;
using NEWMAT::ColumnVector;

namespace OPTPP {

int OptLBFGSLike::checkDeriv() // check the analytic gradient with FD gradient
{return GOOD;}

int OptLBFGSLike::checkConvg() // check convergence
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
  double ftol = tol.getFTol();
  double fvalue = nlp->getF();
  double rftol = ftol*max(1.0,fabs(fvalue));
  Real deltaf = fprev - fvalue;

  if (deltaf <= rftol) {
    strcpy(mesg,"Function tolerance test passed");
    *optout << "checkConvg: deltaf = " << e(deltaf,12,4) 
	    << "  ftol = " << e(ftol,12,4) << "\n";
    return 2;
  }
  

// Test 3. gradient tolerance 

  ColumnVector grad(nlp->getGrad());
  double gtol = tol.getGTol();
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
    strcpy(mesg,"Gradient absolute tolerance test passed");
    *optout << "checkConvg: gnorm = " << e(gnorm,12,4) 
      << "  gtol = " << e(gtol, 12,4) << "\n";
    return 4;
  }
  
  // Nothing to report 

  return 0;

}

void OptLBFGS::printStatus(char *s) // set Message
{

  *optout << "\n\n=========  " << s << "  ===========\n\n";
  *optout << "Optimization method       = " << method << "\n";
  *optout << "Dimension of the problem  = " << dim    << "\n";
  *optout << "Return code               = " << ret_code << " ("
       << mesg << ")\n";
  *optout << "No. iterations taken      = " << iter_taken  << "\n";
  *optout << "No. function evaluations  = " << fcn_evals << "\n";
  *optout << "No. gradient evaluations  = " << grad_evals << "\n";
  *optout << "Function Value            = " << nlp->getF() << "\n";
  *optout << "Norm of gradient          = " << Norm2(nlp->getGrad()) << "\n";

  tol.printTol(optout);

  if (printXs) nlp->fPrintState(optout, s);

}

void OptLBFGS::reset() // Reset parameters 
{
   NLP1* nlp = nlprob();
   int   n   = nlp->getDim();
   nlp->reset();
   OptimizeClass::defaultReset(n);
   grad_evals = 0;
   // Still to do. Reset memM.
}

int OptLBFGS::checkDeriv() // check the analytic gradient with FD gradient
{return GOOD;}

real OptLBFGS::stepTolNorm() const
{
  return Norm2(nlp->getXc()-xprev);
}

int OptLBFGS::computeStep(ColumnVector& sk, double stp)
//---------------------------------------------------------------------------- 
// 
// compute a step along the direction sk using either a backtrack line search
// or a More-Thuente search, starting with stp as suggested step length.
//
//---------------------------------------------------------------------------- 
{
  int  step_type;
  int  itnmax = tol.getMaxBacktrackIter();
  real stp_length = stp;
  real stpmax = tol.getMaxStep();
  real stpmin = tol.getMinStep();
  real ftol = 5.e-1;
  real xtol = tol.getStepTol();
  real gtol = 5.e-1;

  fprev   = nlp->getF();
  xprev   = nlp->getXc();
  gprev   = nlp->getGrad();  

  step_type = linesearch(nlp, optout, sk, sx, &stp_length, stpmax, stpmin,
			   itnmax, ftol, xtol, gtol);
  if (step_type < 0) {
    setMesg("OptLBFGS: Step does not satisfy sufficient decrease condition");
    ret_code = -1;
    setReturnCode(ret_code);
    return(-1);
  }
  fcn_evals   = nlp->getFevals();
  grad_evals  = nlp->getGevals();
  step_length = stp_length;
  return(step_type);
}

void OptLBFGS::initOpt()
{
  time_t t;
  char *c;

// get date and print out header

  t = time(NULL);
  c = asctime(localtime(&t));

  *optout << "************************************************************\n";
  *optout << "OPT++ version " << OPT_GLOBALS::OPT_VERSION << "\n";
  *optout << "Job run at " << c << "\n";
  copyright();
  *optout << "************************************************************\n";

  int n   = nlp->getDim();
  nlp->initFcn();
  nlp->eval();

  if(nlp->hasConstraints()){
//    cerr << "Error: OptLBFGS does not support bound, linear, or nonlinear "
//         << "constraints.\n       Please select a different method for "
//         << "constrained problems." << endl;
    abort_handler(-1);
  }


  fprev   = nlp->getF();
  xprev   = nlp->getXc();
  gprev   = nlp->getGrad();  

  *optout << "\n\t\tNonlinear LBFGS with m = " << memM
	  << "\n  Iter      F(x)      ||grad||    "
	  << "||step||       gtp      fevals  \n\n";

  if (debug_) {
    nlp->fPrintState(optout, "LBFGS: Initial Guess");
    *optout << "xc, grad, step\n";
    for(int i=1; i<=n; i++)
      *optout << d(i,6) << e(xprev(i),24,16) << e(gprev(i),24,16) << "\n";
  }
}

void OptLBFGS::optimize()
//------------------------------------------------------------------------
// Limited Memory BFGS Method for Large Scale Optimization
//
// Solves the unconstrained minimization problem
//
//          min F(x),    x = (x_1, x_2, ..., x_N)
//
// using the limited memory BFGS method, where N can be large.
//
// The inverse Hessian approximation Hk is computed via BFGS
// updates to a diagonal matrix H0.  The number of updates
// depend on the previous m steps kept in memory, as set by the
// user. H0 can be any symmetric positive definite matrix 
// specified by the user (else a default one is constructed).
//
// References:
//
// D. Liu and J. Nocedal, 
// "On the limited memory BFGS method for large scale optimization"
// Mathematical Programming B 45 (1989), 503-528
//
// Mathematical description of the algorithm
//
// f() is objective function, g its gradient at current x
//
// Input: 
// -----
//    x0    -- initial point
//    m     -- memory parameter (m<<N for N large)
//    beta' --     0 < beta' < 1/2
//    beta  -- beta' < beta  < 1
//    H0    -- A symmetric positive definite matrix (Hessian inv approx)
// 
// 0. set k=0, x[0] = x0; H[0]= H0.  
//              
// 1. compute
//         z[k] := H[k] g[k],  // using Strang formula 
//         d[0] := -z[0];      // the search direction
//
// 2. starting with alfa[k]=1, find alpha[k] that minimizes 
//        f(x[k] + alpha[k] * d[k]) subject to the Wolfe conditions
//      
//        f(x[k] + alpha[k] d[k]) <= f(x[k]) + beta' alfa[k] g[k]^T d[k]
//        g(x[k] + alpha[k] d[k])^T d[k] >= beta g[k]^T d[k]
//
// 3. test for convergence 
//
// 4. set m1 = min{k, m-1} and do m1 BFGS updates of H[k] using the pairs 
//     {y[j], s[j]}, j=k-m1,...,k; s[k]==x[k+1]-x[k]; y[k]==g[k+1]-g[k]:
//    
//          H[k+1] = V[k]^T  H[k] V[k]^T + p[k] s[k] s[k]^T
//
//     where p[k] = 1 / (y[k]^T s[k]),  V[k] = I - p[k] y[k] s[k]^T.
//
// 5. set k := k + 1 and go to 1.
//----------------------------------------------------------------------------
//
{ // BEGIN optimize() 

  //--
  // Initialize parent:
  //  this inits the nlp object and outputs the initial state
  //   but it is otherwise algorithm independent
  initOpt(); 

  //---------------------------------------------------------  
  // Init iteration:
  //  this should init all that's needed to compute a step
  //  and for a next step update
  //---------------------------------------------------------  
  int n = dim;
  ColumnVector xk(n), grad(n), W(n);
  double fvalue, gnorm, slope, step, stp1, stp, ginf;

  xk = nlp->getXc();
  grad = nlp->getGrad();
  gnorm = Norm2(grad); 
  ginf = grad.NormInfinity();

  // the initial step_length for the linesearch (mcsrch)
  stp1 = 1.0/gnorm;   

  int m = memM;

  // allocate the memory vectors
  ColumnVector * s = new ColumnVector[m];
  ColumnVector * y = new ColumnVector[m];
  for (int i=0; i<m; i++) { 
    y[i].ReSize(n);   // gradient differences
    s[i].ReSize(n);   // steps
  }
    
  // allocate other vectors used to compute Hg
  ColumnVector rho(m), alpha(m); 

  // check storage
  if (!s || !y || !rho.Storage() || !alpha.Storage() || !W.Storage() )  {
//    cerr << "memory error. " << endl;
    ret_code = -10;
    setReturnCode(ret_code);
    delete [] y; delete [] s;
    return;
  }

  // For now, we will use the default H0=I
  ColumnVector diag(n);
  diag = 1;

  // Init s[0]
  for (int i=1; i<=n; i++) {
    s[0](i) = -grad(i)*diag(i); // diag==1 for now, but not in gnrl
  }

  //--------------------------------------------------
  // Iteration loop:
  //--------------------------------------------------
  double ys, yy;
  int cp, npt = 0, point = 0, bound; // circular indices
  double truestep; // used for output
  int maxiter = tol.getMaxIter();

  printIter(0, nlp->getF(), gnorm, 0.0, 0.0, 0);  
  updateModel(0, n, nlp->getXc());
  for (int iter=1; iter < maxiter; iter++) {

    bound = iter - 1;
    
    if (iter == 1) goto L165;

    if (iter > m) bound = m;

    // update to diag(), npt indicates the previous iterate
    ys = Dot(y[npt], s[npt]); 
    yy = Dot(y[npt], y[npt]);
    
    for (int i=1; i<=n; i++)
      diag(i) = ys / yy;

    //--
    // computation of  -H*grad via Nocedal's (Strang's) formula
    //--
    double sq, yr, bet;

    // point is a circular index over the memory vectors
    // cp is a temporary circular index used in the loops
    cp = (point==0)? m : point ;  
    rho(cp) = 1 / ys;

    W = -grad;

    cp = point;
    for (int i=1; i<=bound; i++) {
      
      if (--cp == -1) cp = m-1;

      sq = Dot(s[cp], W);
      alpha(cp+1) = rho(cp+1) * sq; 
    
      W += -alpha(cp+1) * y[cp]; // daxpy
    }
    
    for (int i=1; i<=n; i++)
      W(i) *= diag(i); 

    for (int i=1; i<=bound; i++) {
   
      yr = Dot(y[cp], W);
      bet = rho(cp+1) * yr;
      bet = alpha(cp+1) - bet;
      
      W += bet * s[cp];

      if (++cp == m) cp = 0;
    }
    
    s[point] = W;
    
  L165:      
    
    stp = (iter==1)? stp1 : 1.0;
    W = grad;

    int step_rc = computeStep(s[point], stp);
      // computes step based on current data;
      // -- accepts step if good, 
    if (step_rc < 0) {
      setMesg("lbfgs: Step does not satisfy sufficient decrease condition");
      ret_code = step_rc;
      setReturnCode(ret_code);
      delete [] y; delete [] s; 
      return;
    }
    iter_taken = iter;
    step       = step_length;
    truestep       = Norm2(xprev - nlp->getXc()); // used for output
    fvalue     = nlp->getF();
    grad       = nlp->getGrad();
    gnorm      = sqrt(Dot(grad,grad));
    ginf       = grad.NormInfinity();
    slope      = Dot(grad,s[point]);

    //  Test for Convergence
    int convgd = checkConvg();
    if (convgd > 0) {
      ret_code = convgd;
      setReturnCode(ret_code);
      printIter(iter, fvalue, gnorm, truestep, slope, fcn_evals);
      updateModel(iter, n, nlp->getXc());
      delete [] y; delete [] s;
      return;
    }

    // --------------------------------------------
    // UPDATES
    // --------------------------------------------
    
    printIter(iter, fvalue, gnorm, truestep, slope, fcn_evals);
    updateModel(iter, n, nlp->getXc());

    // step and gradient changes
    s[point] *= step;
    y[point] = grad - W;

    npt = point;
    if (++point == m) point = 0;    
  }

  // too many iterations
  setMesg("Max numbers of iterations reached");
  ret_code = -4;
  setReturnCode(ret_code);

  // clean up;
  delete [] y; delete [] s;

} // END optimize() 

void OptLBFGS::printIter(int iter, double fvalue, double gnorm, 
			 double truestep, double slope, int nfev) 
{
    *optout 
      << d(iter,5) << " " << e(fvalue,12,4) << " "
      << e(gnorm,12,4) << " " << e(truestep,12,4) << " " 
      << e(slope,12,4) << " " << d(nfev,6)
      << endl;
}

} // namespace OPTPP
