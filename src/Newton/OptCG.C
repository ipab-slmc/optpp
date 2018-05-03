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
#include <cstring>
#include <ctime>
#else
#include <string.h>
#include <time.h>
#endif

#include "OptCG.h"
#include "cblas.h"
#include "ioformat.h"

using namespace std;

using NEWMAT::Real;
using NEWMAT::ColumnVector;

namespace OPTPP {

int OptCGLike::checkDeriv() // check the analytic gradient with FD gradient
{return GOOD;}

int OptCGLike::checkConvg() // check convergence
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
    strcpy(mesg,"Gradient tolerance test passed");
    *optout << "checkConvg: gnorm = " << e(gnorm,12,4) 
      << "  gtol = " << e(gtol, 12,4) << "\n";
    return 4;
  }
  
  // Nothing to report 

  return 0;

}

void OptCG::printStatus(char *s) // set Message
{

  *optout << "\n\n=========  " << s << "  ===========\n\n";
  *optout << "Optimization method       = " << method << "\n";
  *optout << "Dimension of the problem  = " << dim    << "\n";
  *optout << "Return code               = " << ret_code << " ("
       << mesg << ")\n";
  *optout << "No. iterations taken      = " << iter_taken  << "\n";
  *optout << "No. function evaluations  = " << fcn_evals << "\n";
  *optout << "No. gradient evaluations  = " << grad_evals << "\n";

  tol.printTol(optout);

  nlp->fPrintState(optout, s);
}

void OptCG::reset() // Reset parameters 
{
   NLP1* nlp = nlprob();
   int   n   = nlp->getDim();
   nlp->reset();
   OptimizeClass::defaultReset(n);
   grad_evals = 0;
}

int OptCG::checkDeriv() // check the analytic gradient with FD gradient
{return GOOD;}

real OptCG::stepTolNorm() const
{
  return Norm2(nlp->getXc()-xprev);
}

int OptCG::computeStep(ColumnVector sk)
//---------------------------------------------------------------------------- 
// 
// compute a step along the direction sk using either a backtrack line search
// or a More-Thuente search
//
//---------------------------------------------------------------------------- 
{
  int  step_type;
  int  itnmax = tol.getMaxBacktrackIter();
  real stp_length = 1.0;
  real stpmax = tol.getMaxStep();
  real stpmin = tol.getMinStep(); // 1.e-9;
  real ftol = 5.e-1;
  real xtol = 2.2e-16;
  real gtol = 5.e-1;

  step_type = linesearch(nlp, optout, sk, sx, &stp_length, stpmax, stpmin,
			   itnmax, ftol, xtol, gtol);
  if (step_type < 0) {
    setMesg("OptCG: Step does not satisfy sufficient decrease condition");
    ret_code = -1;
    setReturnCode(ret_code);
    return(-1);
  }
  fcn_evals   = nlp->getFevals();
  grad_evals  = nlp->getGevals();
  step_length = stp_length;
  return(step_type);
}

void OptCG::initOpt()
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

  nlp->initFcn();
  ret_code = 0;

  if(nlp->hasConstraints()){
    CompoundConstraint* constraints = nlp->getConstraints();
    ColumnVector xstart = nlp->getXc();
    double feas_tol = tol.getCTol();
    bool feasible = constraints->amIFeasible(xstart, feas_tol);
    if (!feasible) {
      *optout << "OptCG WARNING:  Initial guess not feasible.\n"
              << "CG may be unable to make progress." << endl;
    }
    //cerr << "Error: OptCG does not support bound, linear, or nonlinear "
    //     << "constraints.\n       Please select a different method for "
    //     << "constrained problems." << endl;
    //abort_handler(-1);
  }

  if (ret_code == 0) {
    double fvalue, gnorm;

    int n     = nlp->getDim();
    nlp->eval();

    fvalue  = nlp->getF();
    fprev   = fvalue;
    xprev   = nlp->getXc();
    gprev   = nlp->getGrad();  gnorm = Norm2(gprev);


    *optout << "\n\t\t\t\tNonlinear CG"
	    << "\n  Iter      F(x)       ||grad||    "
	    << "||step||     beta       gtp        fcn\n\n"
	    << d(0,5) << " " << e(fvalue,12,4) << " " << e(gnorm,12,4) << endl;

    if (debug_) {
      nlp->fPrintState(optout, "qnewton: Initial Guess");
      *optout << "xc, grad, step\n";
      for(int i=1; i<=n; i++)
	*optout << d(i,6) << e(xprev(i),24,16) << e(gprev(i),24,16) << "\n";
    }
  }
}

void OptCG::optimize()
//------------------------------------------------------------------------
// Nonlinear Preconditioned Conjugate Gradient
// 
// Given a nonlinear operator objfcn find the minimizer using a
// nonlinear conjugate gradient method
// This version uses the Polak-Ribiere formula.
// and a line search routine due to More and Thuente as implemented
// in the routines mcsrch and mcstep
//
// Notes: The parameters ftol and gtol should be set so that
//        0 < ftol < gtol < 0.5
//        Default values: ftol = 1.e-1, gtol = 5.e-1
//        This results in a fairly accurate line search
//
// Here is the mathematical description of the algorithm
// (g = grad f).
//                     -1
//        1.  set z = M  g, search = -z; 
//
//        2.  for i=0 until convergence
//
//                 find alpha that minimizes f(x + alpha*search)
//                 subject to the strong Wolfe conditions
//
//                 Test for convergence
//
//
//                 beta = -( g  ,  (z     - z ) ) / ( g  , z  )
//                            i+1    i + 1   i         i    i
//
//                 search     =  - z     +   beta * search
//                       i+1        i+1                   i
//
//----------------------------------------------------------------------------
     
{

  int i, nlcg_iter;
  int convgd = 0;
  int step_type;

  double beta;
  double delta, delta_old, delta_mid, delta_new;
  double slope, gnorm;

  double step;
  double zero = 0.;
  
// Allocate local vectors 

  int n = dim;
  int maxiter;
  double fvalue;
  ColumnVector search(n), grad(n), z(n), diag(n), xc(n);

// Initialize iteration

  maxiter = tol.getMaxIter();

  initOpt();

  if (ret_code == 0) {
    //  compute preconditioned gradient

    diag = getFcnScale();
    grad = nlp->getGrad();
    for (i=1; i<=n; i++) z(i) = grad(i)/diag(i);

    search    = -z;
    delta_old = delta_new = Dot(grad,z);
    gnorm     = sqrt(Dot(grad,grad));

    step    = 1.0/gnorm;

//---------------------------------------------------------------------------
//
//
//
//
//---------------------------------------------------------------------------

    for (nlcg_iter=1; nlcg_iter <= maxiter; nlcg_iter++) {

      iter_taken = nlcg_iter;

      //  compute a step along the direction search 

      if ((step_type = computeStep(search)) < 0) {
	setMesg("nlcg: Step does not satisfy sufficient decrease condition");
	ret_code = step_type;
        setReturnCode(ret_code);
	return;
      }
    
      //  Accept this step and update the nonlinear model

      acceptStep(nlcg_iter, step_type);
      updateModel(nlcg_iter, n, xprev);

      xc         = nlp->getXc();
      mem_step   = xc - xprev;
      step       = Norm2(mem_step);

      fvalue     = nlp->getF();
      grad       = nlp->getGrad();
      gnorm      = sqrt(Dot(grad,grad));
      slope      = Dot(grad,search);

      //  Test for Convergence

      convgd = checkConvg();
      if (convgd > 0) {
	ret_code = convgd;
        setReturnCode(ret_code);
        setMesg("OptCG: Algorithm converged");
	*optout  << d(nlcg_iter,5) << " " << e(fvalue,12,4)  << " "
		 << e(gnorm,12,4)  << e(step,12,4) << "\n";
	return;
      }

      //
      //  compute a new search direction
      //  1. compute preconditioned gradient,  z = grad;
      //  2. beta is computed using Polak-Ribiere Formula constrained 
      //     so that beta > 0
      //  3  Update search direction and norms 
  
      delta_old = delta_new; delta_mid = Dot(grad,z);

      for (i=1; i<=n; i++) z(i) = grad(i)/diag(i);

      delta_new = Dot(grad,z);
      delta     = delta_new - delta_mid;
      beta      = max(zero,delta/delta_old);

      search = -z + search*beta;

      xprev  = nlp->getXc();
      fprev  = fvalue;
      gprev  = grad;

      *optout 
	<< d(nlcg_iter,5) << " " << e(fvalue,12,4) << " " << e(gnorm,12,4) 
	<< e(step,12,4)   << " " << e(beta,12,4)   << " " << e(slope,12,4) 
	<< d(fcn_evals,4) << " " << d(grad_evals,4) << endl;
    }

    setMesg("Maximum number of iterations in nlcg");
    ret_code = 4;
    setReturnCode(ret_code);
  }
}

} // namespace OPTPP
