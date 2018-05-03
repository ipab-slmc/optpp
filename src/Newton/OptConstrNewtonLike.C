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

#include "OptConstrNewtonLike.h"
#include "precisio.h"
#include "cblas.h"
#include "ioformat.h"

using NEWMAT::Real;
using NEWMAT::FloatingPointPrecision;
using NEWMAT::ColumnVector;
using NEWMAT::Matrix;
using NEWMAT::DiagonalMatrix;
using NEWMAT::LowerTriangularMatrix;
using NEWMAT::UpperTriangularMatrix;
using NEWMAT::SymmetricMatrix;

#ifdef WITH_MPI
#include "mpi.h"
#endif

namespace OPTPP {

static char* class_name = {"OptConstrNewtonLike"};

//------------------------------------------------------------------------
//
//   Newton-Like Base Class functions
//   Notes:
//   These functions are first declared in Opt.h as
//   virtual functions for the abstract base class optimize
//   Therefore we need to define them so that the derived classes
//   can be instantiated by the end user.  Of course the user
//   can provide his/her own version of these.
//------------------------------------------------------------------------

//------------------------------------------------------------------------
//
// First the default functions
// defaultAcceptStep
// defaultComputeSearch
//
//------------------------------------------------------------------------

void OptConstrNewtonLike::defaultAcceptStep(int iter, int step_type)
{
// Successful step
// Print out iteration summary and anything else
// you might want to do before the next iteration

  double condh;

  if (trace) 
    *optout << "\n***** OptConstrNewtonLike:DefaultAcceptStep\n";

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
    FPrint(optout, Hessian);

//  Compute eigenvalues of Hessian

    DiagonalMatrix D(n);
    EigenValues(Hessian, D);
    *optout << "\nEigenvalues of Hessian";
    FPrint(optout, D);
    condh = D(n,n) / D(1,1);

    *optout << "Reciprocal Condition Number of H = " << condh << "\n";
    *optout << "\n***************************************";
    *optout << "***************************************\n";


  }
//  Iteration summary
// 

  if(step_type >= 0){
  *optout 
	<< d(iter,5)  << " " << e(fvalue,12,4) << " " << e(gnorm,12,4) << " "
	<< e(step_length,12,4) << "  " << steps[step_type] << " " 
	<< d(fcn_evals,5) << " " << d(grad_evals,5) << "\n" << flush;
  }
  else{
  *optout 
	<< d(iter,5)  << " " << e(fvalue,12,4) << " " << e(gnorm,12,4) << " "
	<< e(step_length,12,4) << "  " << "  " << " " 
	<< d(fcn_evals,5) << " " << d(grad_evals,5) << "\n" << flush;
  }


}
ColumnVector OptConstrNewtonLike::defaultComputeSearch(SymmetricMatrix& H)
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
// to all Newton cases
//
// checkAnalyticFDGrad
// checkConvg
// checkDeriv
// computeStep
// initOpt
// initHessian
// initTrustRegionSize
// optimize
// reset
// readOptInput
// printStatus
//------------------------------------------------------------------------

int OptConstrNewtonLike::checkAnalyticFDGrad()
{
  Real mcheps = FloatingPointPrecision::Epsilon();

  int i;
  int retcode = GOOD;
  int n = dim;

  Real third = 0.33333;
  ColumnVector error(n);

  NLP1* nlp = nlprob();
  ColumnVector xc = nlp->getXc();
  double fx = nlp->getF();
  SpecOption tmpSpec = nlp->getSpecOption();
  ColumnVector fd_grad(n);
  nlp->setSpecOption(NoSpec);
  fd_grad = nlp->FDGrad(sx, xc, fx, fd_grad);    // Evaluate gradient using finite differences
  nlp->setSpecOption(tmpSpec);
  ColumnVector grad(nlp->getGrad());        // Now use the analytical functions

  double gnorm = grad.NormInfinity();
  double eta   = pow(mcheps,third)*max(1.0,gnorm);

  *optout << "Check_Deriv: Checking gradients versus finite-differences\n";
  *optout << "    i    gradient     fd grad       error\n";
  for (i=1; i<=n; i++) {
    error(i) = fabs(grad(i)-fd_grad(i));
    *optout << d(i,5) << e(grad(i),12,4) 
	   << e(fd_grad(i),12,4) << e(error(i),12,4);
  }
  Real maxerr = error.NormInfinity();		
  *optout << "maxerror = " << e(maxerr, 12,4) 
     << "tolerance =  " << e(eta, 12,4) << "\n";
  if (maxerr > eta) retcode = BAD;

  return retcode;
}

int OptConstrNewtonLike::checkConvg() // Check convergence
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
    *optout << "CheckConvg: snorm = " << e(snorm,12,4) 
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
    *optout << "CheckConvg: deltaf = " << e(deltaf,12,4) 
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
    *optout << "CheckConvg: gnorm = " << e(gnorm,12,4) 
      << "  gtol = " << e(rgtol, 12,4) << "\n";
    return 3;
  }
  

// Test 4. absolute gradient tolerance 

  if (gnorm <= gtol) {
    strcpy(mesg,"Gradient tolerance test passed");
    *optout << "CheckConvg: gnorm = " << e(gnorm,12,4) 
      << "  gtol = " << e(gtol, 12,4) << "\n";
    return 4;
  }
  
  // Nothing to report 

  return 0;

}

// Check the analytic gradient with FD gradient
int OptConstrNewtonLike::checkDeriv()
{return GOOD;}

real OptConstrNewtonLike::computeMaxStep(ColumnVector sk)
{return FLT_MAX;}

int OptConstrNewtonLike::computeStep(ColumnVector sk)
//---------------------------------------------------------------------------- 
// 
// Compute a step along the direction sk using either a line search
// or model trust-region approach
//
//---------------------------------------------------------------------------- 
{
  NLP1*  nlp = nlprob();
  real stp_length = 1.0;
  real lstol  = tol.getLSTol();
//  real xtol   = tol.getStepTol();
//  real gtol   = tol.getGTol();
  real stpmax = tol.getMaxStep();
  real stpmin = tol.getMinStep();
  int  step_type;
  int  itnmax = tol.getMaxBacktrackIter();

  if (trace) *optout << class_name << ": ComputeStep\n";

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
    setMesg("OptConstrNewtonLike: Step does not satisfy sufficient decrease condition");
    ret_code = -1;
    setReturnCode(ret_code);
    return(ret_code);
  }
  fcn_evals   = nlp->getFevals();
  grad_evals  = nlp->getGevals();
  step_length = stp_length;
  return(step_type);
}

void OptConstrNewtonLike::initOpt()
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
      *optout << "OptConstrNewtonLike WARNING:  Initial guess not feasible.\n"
	      << "ConstrNewton may be unable to make progress." << endl;
    }
  }

  if (ret_code == 0) {
    // Evaluate Function, gradient and compute initial Hessian

    nlp->eval();

    xprev = nlp->getXc();
    fprev = nlp->getF();
    gprev = nlp->getGrad();
    gnorm = Norm2(gprev);
  
    //  SymmetricMatrix Hk(n);
    //  Hessian = updateH(Hk,0);

    initHessian();
    setFcnScale(fprev);

    // get Optimization parameters

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
      nlp->fPrintState(optout, "OptConstrNewtonLike: Initial Guess");
      *optout << "xc, grad, step\n";
      for(int i=1; i<=n; i++)
	*optout << i << e(xprev(i),24,16) << e(gprev(i),24,16) << "\n";
      FPrint(optout, Hessian);
    }
  }
}
void OptConstrNewtonLike::initHessian()
{ 
  int i;
  NLP1* nlp = nlprob();
  int ndim = nlp->getDim();

  if (WarmStart) {
    *optout << "OptNewtonlike::initHessian: Warm Start specified\n";
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

    for (i=1; i <= ndim; i++) xmax = max(xmax,fabs(xc(i)));
    if(xmax != 0.0) typx = xmax;
    if(gnorm!= 0.0) D    = gnorm/typx;
    if (debug_) {
      *optout << "OptNewtonlike::initHessian: gnorm0 = " << gnorm
	<< "  typx = " << typx << "\n";
    }
    Hessian = 0.0;
    for (i=1; i <= ndim; i++) Hessian(i,i) = D(i);
   }
}
double OptConstrNewtonLike::initTrustRegionSize() const
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

//---------------------------------------------------------------------------- 
// Reset optimization parameters 
//---------------------------------------------------------------------------- 
void OptConstrNewtonLike::reset()
{
   NLP1* nlp = nlprob();
   int   n   = nlp->getDim();
   nlp->reset();
   OptimizeClass::defaultReset(n);
   me = mi = grad_evals = 0;
   TR_size = cost       = 0.0;
   gradMult= 0.1;
   searchSize = 64;
   gprev      = 0;
   gradl      = 0;
   gradlprev  = 0;
   constraintResidual     = 0;
   constraintGradient     = 0;
   constraintGradientPrev = 0;
   
}

void OptConstrNewtonLike::optimize()
//---------------------------------------------------------------------------- 
//
// Given a nonlinear operator nlp find the minimizer using a
// Newton-like method
//
//---------------------------------------------------------------------------- 
{
  int k;
  int maxiter, maxfev, myfevals, fevals;
  int convgd = 0;
  int step_type;

// Allocate local vectors 

  int n = dim;
  ColumnVector sk(n);
  SymmetricMatrix Hk(n);

// Initialize iteration
// Evaluate Function, Gradient, and Hessian

  initOpt();

  if (ret_code == 0) {
    maxiter = tol.getMaxIter();
    maxfev  = tol.getMaxFeval();

    Hk = Hessian;

    // Check for convergence. Need to take into account that this is the
    // zeroth iteration
    //  convgd = objfcn.check_Convg(tol,stat);
    //  if (convgd > 0) {
    //    stat.ret_code = convgd;
    //    return;
    //  }
  
    for (k=1; k <= maxiter; k++) {

      iter_taken = k;

      //  Solve for the Newton direction
      //  H * step = -grad;

      sk = computeSearch(Hk);

      //  ComputeStep will attempt to take a step in the direction sk 
      //  from the current point. 
      //  The default method is to use a trust region

      if ((step_type = computeStep(sk)) < 0) {
	*optout << "step_type = " << step_type << "\n";
	setMesg("OptConstrNewtonLike: Step does not satisfy sufficient decrease condition");
	ret_code = step_type;
        setReturnCode(ret_code);
	return;
      }

      //  Accept this step and update the nonlinear model

      acceptStep(k, step_type);

      //  Test for Convergence

      convgd = checkConvg();
      if (convgd > 0) {
	ret_code = convgd;
        setReturnCode(ret_code);
	return;
      }

      NLP1* nlp = nlprob();
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

      // Update state
      Hessian = updateH(Hk,k);
      Hk = Hessian;

      xprev = nlp->getXc();
      fprev = nlp->getF();
      gprev = nlp->getGrad();

      updateModel(k, n, xprev);
    }

    setMesg("OptConstrNewtonLike: Maximum number of iterations or fevals");
    ret_code = -4;
    setReturnCode(ret_code);
  }
}

void OptConstrNewtonLike::printStatus(char *title) // set Message
{
  NLP1* nlp = nlprob();

  *optout << "\n\n=========  " << title << "  ===========\n\n";
  *optout << "Optimization method       = " << method << "\n";
  *optout << "Dimension of the problem  = " << nlp->getDim()  << "\n";
  *optout << "Return code               = " << ret_code << " ("
       << mesg << ")\n";
  *optout << "No. iterations taken      = " << iter_taken  << "\n";
  *optout << "No. function evaluations  = " << nlp->getFevals() << "\n";
  *optout << "No. gradient evaluations  = " << nlp->getGevals() << "\n";

  if (debug_) {
    *optout << "\nHessian";
    FPrint(optout, Hessian);
//  Compute eigenvalues of Hessian

    DiagonalMatrix D;
    SVD(Hessian, D);
    *optout << "\nEigenvalues of Hessian";
    FPrint(optout, D);
  }

  tol.printTol(optout);

  nlp->fPrintState(optout, title);
  fPrintMultipliers(optout, title);
}

void OptConstrNewtonLike::printMultipliers(char *cs) // set Message
{
  int i;
   // Print out current state: Multipliers 
  cout << "\n\n=========  " << cs << "  ===========\n\n";
  cout << "\n    i\t   y    \n\n";
  for (i=1; i < me + 1; i++)
     cout << d(i,5) << e(y(i),12,4) << "\n";
  cout <<"\n\n=====================================\n\n";
  cout << "\n    i\t    z \t      s\n\n";
  for (i = 1; i < mi+1; i++)
     cout << d(i,5) << e(z(i),12,4) <<  e(s(i),12,4) << "\n";
}

void OptConstrNewtonLike::fPrintMultipliers(ostream *nlpout,char *cs) 
{
  int i;
   // Print out current state: Multipliers 
  (*nlpout) << "\n\n=========  " << cs << "  ===========\n\n";
  (*nlpout) << "\n    i\t   y    \n\n";
  for (i=1; i < me + 1; i++)
     (*nlpout) << d(i,5) << e(y(i),12,4) << "\n";
  (*nlpout) <<"\n\n=====================================\n\n";
  (*nlpout) << "\n    i\t    z \t      s\n\n";
  for (i = 1; i < mi+1; i++)
     (*nlpout) << d(i,5) << e(z(i),12,4) <<  e(s(i),12,4) << "\n";
}

void OptConstrNewtonLike::fPrintSecSuff(ostream *nlpout, ColumnVector& info) 
{
  int i;
  int ncActive = (int) info(dim+mi+1);
  int rank     = (int) info(dim+mi+2);
  (*nlpout) << "\n\n=========  Second-Order Sufficiency Test   ===========\n\n";
  (*nlpout) << "Number of active constraints         =  " << d(ncActive,5) << "\n";
  (*nlpout) << "Approx rank of gradient set (active) =  " << d(rank,5) << "\n";
  (*nlpout) << "List of active/non-active constraints " << "\n";
  (*nlpout) << "      Active( 0 = N, 1= YES)          " << "\n";
  for (i=1; i < mi + 1; i++)
     (*nlpout) << d(i,5) << e(info(dim+i),3,1) << "\n";
  (*nlpout) << "Eigenvalues of the projected hessian " << "\n";
  for (i=1; i < dim - rank; i++)
     (*nlpout) << d(i,5) << e(info(i),3,1) << "\n";
  (*nlpout) <<"\n\n===================================================\n\n";
}

ColumnVector OptConstrNewtonLike::computeFFK1Ind(const ColumnVector& xc)
{
/*  
 * For details, see Facchinei, Fischer & Kanzow's "On the Identification 
 * Of Active Constraints", Siam J. Optim. Vol.9, No. 1, pp. 14-32.
 */ 

   // Local variables
    int i;
    real r, rhoxz, dotprod = 0.0e0;
    const real zero        = 0.0e0;
    const real point_nine  = 0.9e0;
    ColumnVector conresid(me+mi), g(mi), zplus(mi), activeSet(mi);

//
//   Evaluate the  constraints
//
    conresid   = getConstraintResidual();

    for(i = 1; i <= mi; i++){
       g(i)     = max(0.0, -conresid(me+ i));
       zplus(i) = max(0.0, -z(i));
       dotprod += conresid(me+ i)*z(i);
    }
    
/*  r(x,z) = || gradl || + |z'g(x)| + ||max(0,-z)|| + || max(0, -g(x)|| */ 
    r  = Norm2(getGradL()) + Norm2(g)+ Norm2(zplus) + fabs(dotprod); 
    
    if( r == zero )
       rhoxz = 0; 
    if( r < point_nine  &&  r > zero)
       rhoxz = -1/log(r); 
    if( r >= point_nine )
       rhoxz = -1/log(point_nine); 

    for(i = 1; i <= mi; i++){
       if(conresid(me+ i) <= rhoxz)
          activeSet(i)  = 1;
       else
          activeSet(i)  = 0;
    }
    return activeSet;
}


ColumnVector OptConstrNewtonLike::computeFFK2Ind(const ColumnVector& xc)
{
/*  
   For details, see Facchinei, Fischer & Kanzow's "On the Identification 
   Of Active Constraints", Siam J. Optim. Vol.9, No. 1, pp. 14-32.
 */ 

   // Local variables
    int i;
    real rhoxz;
    ColumnVector conresid(me+mi),  min_gz(mi), activeSet(mi);
//
//   Evaluate the constraints
//
    conresid   = getConstraintResidual();

    for(i = 1; i <= mi; i++)
       min_gz(i)  = min(conresid(me+ i),z(i));
    
//  r(x,z) = || gradl | min_gz|| 
    min_gz &= getGradL();
    
    rhoxz = sqrt(Norm2(min_gz)); 

    for(i = 1; i <= mi; i++){
       if(conresid(me+i) <= rhoxz)
          activeSet(i)  = 1;
       else
          activeSet(i)  = 0;
    }
    return activeSet;
}

ColumnVector OptConstrNewtonLike::computeTapiaInd(const ColumnVector& sk)
{
/*  
 * For details, see Tapia's "On the role of slack variables in quasi-Newton
 * methods for constrained optimization", in Numerical Optimization of
 * Dynamic Systems,  L.C.W. Dixon and G.P. Szego, eds., pp. 235-246.
 *
 * See also, El-Bakry, Tapia, and Zhang's "A Study of Indicators For 
 * Identifying Zero Varialbes in Interior-Point Methods", SIAM Review
 * Vol. 36, No. 1, pp 45-72, March 1994
 */ 

   // Local variables
    NLP1* nlp  = nlprob();
    int i, n   = nlp->getDim();
    ColumnVector Ts(mi), Tz(mi), activeSet(mi);

    for(i = 1; i <= mi; i++){
       Ts(i)  = ( s(i) + sk(n+me+mi+i) )/s(i);
       Tz(i)  = ( z(i) + sk(n+me+i) )/z(i);

       if( fabs(Ts(i)) + fabs(1-Tz(i)) <= .2e0 )
          activeSet(i)  = 1;
       else
          activeSet(i)  = 0;
    }

    return activeSet;
}

void OptConstrNewtonLike::readOptInput() // Read Opt.input file if it exists
{
  NLP1* nlp = nlprob();

/* A VERY simple routine for reading the Optimization parameters
 * We should really make this more general, but as a first pass this
 * will have to do.
 * 
 * The input file should be of the form keyword = value
 * where keyword is one of the following
 * 
 * search      = trustregion
 * diff_Option = forward
 * max_iter    = 100
 * maxfeval   = 1000
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
  string cdebug("debug");
  string cdiff_Option("diff_Option");
  string cfcn_accrcy("fcn_accrcy");
  string cfcn_tol("fcn_tol");
  string cgrad_tol("grad_tol");
  string cmaxfeval("maxfeval");
  string cmaxiter("max_iter");
  string cmax_step("max_step");
  string csearch("search");

  string diff_Option, debug_flag;
  string search;
  SearchStrategy search_strategy = TrustRegion;

  int keyword_count = 0;

// 
// Default name of input file
//
  char *Opt_input  = {"Opt.input"};

//
// Open Opt.input file and check to see if we succeeded
//

  ifstream Optin(Opt_input);
  if (!Optin.rdbuf()->is_open()) {
    if (debug_) {
      *optout << "OptConstrNewtonLike::ReadOptInput: No Opt.input file found\n";
      *optout << "OptConstrNewtonLike::ReadOptInput: Default values will be used\n";
    }
    return;
  }

  if (debug_) *optout << "OptConstrNewtonLike::ReadOptInput: Reading Opt.input file\n";

  Optin >> token;

  while (!Optin.eof()) {

    keyword = token;
    keyword_count++;

    if (keyword == cdiff_Option) {

      Optin >> equals >> token;
      diff_Option = token;

      if ( diff_Option == "forward")
	nlp->setDerivOption(ForwardDiff);
      else if ( diff_Option == "backward")
	nlp->setDerivOption(BackwardDiff);
      else if ( diff_Option == "central")
	nlp->setDerivOption(CentralDiff);
    }    
    else if (keyword == cdebug) {
      Optin >> equals >> token;
      debug_flag = token;
      if ( debug_flag == "true") {
	setDebug();
	nlp->setDebug();
      }
    }    
    else if (keyword == cfcn_tol) {
      Optin >> equals >> fcn_tol;
      setFcnTol(fcn_tol);
    }    
    else if (keyword == cgrad_tol) {
      Optin >> equals >> grad_tol;
      setGradTol(grad_tol);
    }    
    else if (keyword == cmax_step) {
      Optin >> equals >> max_step;
      setMaxStep(max_step);
    }
    else if (keyword == csearch) {
      Optin >> equals >> token;
      search = token;
      if ( search == "trustregion")
	search_strategy = TrustRegion;
      else if ( search == "linesearch")
	search_strategy = LineSearch;
      else if ( search == "trustpds")
	search_strategy = TrustPDS;
      setSearchStrategy(search_strategy);
    }
    else if (keyword == cfcn_accrcy) {
      Optin >> equals >> index >> fcn_accrcy;
      nlp->setFcnAccrcy(index, fcn_accrcy);
    }    
    else if (keyword == cmaxfeval) {
      Optin >> equals >> max_feval;
      setMaxFeval(max_feval);
    }    
    else if (keyword == cmaxiter) {
      Optin >> equals >> max_iter;
      setMaxIter(max_iter);
    }
    else {
      *optout << "Unrecognized keyword '" << keyword << "'. "
	<< "Skipping the rest of this line\n";
      Optin.getline(ignore, sizeof(ignore));
    }
  Optin >> token;
  }

  *optout << "\n\n======  Summary of input file  ======\n\n";

  *optout << csearch      << " = " << search << "\n";
  *optout << cdiff_Option << " = " << diff_Option << "\n";
  *optout << cmaxiter     << " = " << max_iter << "\n";
  *optout << cmaxfeval    << " = " << max_feval << "\n";
  *optout << cgrad_tol    << " = " << grad_tol << "\n";
  *optout << cfcn_tol     << " = " << fcn_tol << "\n";
  *optout << cmax_step    << " = " << max_step << "\n";
  ColumnVector fcnacc = nlp->getFcnAccrcy();
  for(int i=1; i <= fcnacc.Nrows(); i++)
     *optout << cfcn_accrcy  << " = " << fcnacc(i) << "\n";

  tol.printTol(optout);

}

} // namespace OPTPP
