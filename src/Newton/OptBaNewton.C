//------------------------------------------------------------------------
// Copyright (C) 1996:
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

#include "OptBaNewton.h"
#include "precisio.h"
#include "cblas.h"
#include "ioformat.h"

using NEWMAT::ColumnVector;
using NEWMAT::Matrix;
using NEWMAT::DiagonalMatrix;
using NEWMAT::SymmetricMatrix;
using NEWMAT::LowerTriangularMatrix;

namespace OPTPP {

//------------------------------------------------------------------------
// Notes: Some of these functions are first declared in opt.h as virtual
//        functions for the abstract base class Optimize.
//------------------------------------------------------------------------
// initHessian
// initOpt
// optimize
// checkInnerConvg
// checkConvg
// updateBarrierMultiplier
// compute_Barrier_Gradient
// compute_Barrier_Hessian
// computeSearch2
// computeStep
// computeMaxStep
// scalarNewton
//------------------------------------------------------------------------

//---------------------------------------------------------------------------- 
// Initialize the barrier Hessian
//---------------------------------------------------------------------------- 
void OptBaNewton::initHessian()
{
  NLP2*       nlp = nlprob2();
  int           n = nlp->getDim();
  ColumnVector xc = nlp->getXc();

  // get the Hessian and Compute the barrier Hessian 

  Hessian   = nlp->getHess();
  Hess_barrier.ReSize(n);
  //JCM: get rid of warning message
  //  Hess_barrier = compute_Barrier_Hessian(nlp->getHess(),xc);
  Hess_barrier = compute_Barrier_Hessian(Hessian,xc);
}

//---------------------------------------------------------------------------- 
// Initialize the other barrier parameters
//---------------------------------------------------------------------------- 
void OptBaNewton::initOpt()
{
  // initialize mu - the multiplier for the barrier term

  mu = 1.0e-2;

  // perform whatever that was done in OptBCNewtonLike::InitOpt

  OptBCNewtonLike::initOpt();

  if (ret_code == 0) {
    // get the nonlinear problem 

    NLP2* nlp = nlprob2();
    int   n = nlp->getDim();

    // Compute the barrier function value

    double   fvalue = nlp->getF();
    ColumnVector xc = nlp->getXc();
    fvalue_barrier  = compute_Barrier_Fvalue(fvalue,xc);

    // Compute the barrier gradient 

    ColumnVector local_grad = nlp->getGrad();
    grad_barrier.ReSize(n,1);
    grad_barrier = compute_Barrier_Gradient(local_grad,xc);
  }
}

//---------------------------------------------------------------------------- 
// Given a nonlinear operator nlp find the minimizer using the Newton's method
//---------------------------------------------------------------------------- 
void OptBaNewton::optimize()
{
  int     inner_convgd, outer_convgd, step_type;
  int     inner_iter_taken, outer_iter_taken;

  // get the nonlinear problem and its dimension

  NLP2* nlp = nlprob2();
  int     n = nlp->getDim();

  // declare a local search vector of length n 

  ColumnVector    search_vector(n); 

  // Initialize Function , Gradient, and Hessian

  initOpt();

  if (ret_code == 0) {
    // Initialize other iteration parameters 

    outer_convgd     = false;
    outer_iter_taken = 0;
    iter_taken       = 0;

    // The main outer loop

    while (!outer_convgd) {
    
      // initialize for inner iterations

      outer_iter_taken++; 
      inner_convgd = false;
      inner_iter_taken = 0;
      fprev_outer = nlp->getF();

      while (!inner_convgd) {
	inner_iter_taken++; 
	if(debug_)
	  *optout << "OptBaNewton::Optimize: iteration count = " << iter_taken << "\n";
	iter_taken++;

	// temporarily put variables aside to accommodate new data

	setAsideCurrentVariables();

	// calculate the search direction

	search_vector = computeSearch2(Hess_barrier,grad_barrier);

	// compute the step length using quadratic-logarithmic interpolation

	step_type = computeStep(search_vector);
	if( debug_) *optout << "step_type = " << step_type << "\n";

	// if successful, accept the step; otherwise terminate inner iterations

	if (step_type < 0) inner_convgd = true;
	else {
	  acceptStep(iter_taken, step_type);
	  inner_convgd = checkInnerConvg(outer_iter_taken);
	} 
      } // while - inner loop 

      // Compute the next mu and check for convergence

      updateBarrierMultiplier();
      outer_convgd = checkConvg();

    } // while - outer loop
  }
}

//---------------------------------------------------------------------------- 
// Check for convergence in the inner iterations
//---------------------------------------------------------------------------- 
int OptBaNewton::checkInnerConvg(int iter) 
{
  NLP2*     nlp = nlprob2();
  ColumnVector xc = nlp->getXc();
  double      dtmp, epik, xnorm, gnorm;

  epik  = pow(10.0e0,-(iter+1.0e0));
  epik  = max(1.0e-5,epik);
  xnorm = Norm2(xc);
  dtmp  = max(1.0e0,xnorm);
  gnorm = Norm2(grad_barrier);
  dtmp  = gnorm / dtmp;
  if( debug_) 
     *optout << "CheckInnerConvg : " << dtmp << " < " << epik << " ? \n";
  if (dtmp < epik) return true; else return false;
}

//---------------------------------------------------------------------------- 
// Check for convergence in the outer iterations
//---------------------------------------------------------------------------- 
int OptBaNewton::checkConvg() 
{
  NLP2*        nlp = nlprob2();
  ColumnVector xc = nlp->getXc();
  ColumnVector grad = nlp->getGrad();
  ColumnVector upper = nlp->getConstraints()->getUpper();
  ColumnVector lower = nlp->getConstraints()->getLower();
  double       fvalue, deltaf, rftol;
  double       gnorm, xnorm, q1, q2, qtmp;
  int          i, n = nlp->getDim();

  // Test 1. function tolerance

  if(mu < 1.0e-12){
    strcpy(mesg,"Mu is TOO SMALL to continue");
    return 3;
  }

  fvalue = nlp->getF();
  deltaf = fprev_outer - fvalue;
  if (deltaf == 0.0) return 0;

  rftol = 1.0e-6 * (1.0e0 + fabs(fprev));
  if (deltaf <= rftol) {
    *optout << "CheckConvg: deltaf = " << e(deltaf,12,4) 
         << " rftol = " << e(rftol,12,4) << "\n";
    return 1;
  }
  
  // Test 2. gradient tolerance 

  xnorm = Norm2(xc);
  for (i=1; i<=n; i++) {
    if (fabs(xc(i)-lower(i))<1.0e-4 || fabs(upper(i)-xc(i))<1.0e-4) grad(i) = 0.0;
  }
  gnorm = Norm2(grad_barrier);
  q1    = gnorm / (1.0e0 + xnorm);
  *optout << "CheckConvg: gnorm/(1+xnorm) = " << e(q1,12,4) << "\n"; 
  q2    = FLT_MAX;
  for (i=1; i<=n; i++) {
    qtmp = xc(i) - lower(i); q2 = (qtmp < q2) ? qtmp : q2;
    qtmp = upper(i) - xc(i); q2 = (qtmp < q2) ? qtmp : q2;
  }
  q2 = - q2;
  qtmp = max(q1, q2);
  if (qtmp < 1.0e-4) {
    strcpy(mesg,"Function and gradient tolerance test passed");
    return 2;
  } else
    return 0;
}

//---------------------------------------------------------------------------- 
// Update the Lagrange multipliers and mu
//---------------------------------------------------------------------------- 
void OptBaNewton::updateBarrierMultiplier()
{
  NLP2*        nlp = nlprob2();
  ColumnVector    xc = nlp->getXc();
  int           i, n = nlp->getDim();
  ColumnVector  upper = nlp->getConstraints()->getUpper();
  ColumnVector  lower = nlp->getConstraints()->getLower();
  double       maxmu, dtmp;

  maxmu = 10.0;
  for (i=1; i<=n; i++) {
    if (lower(i) != -FLT_MAX) {
      dtmp = (xc(i) - lower(i)) / mu;
      if (dtmp < 0.0) maxmu = min(maxmu, 1.0e0/dtmp);
    }
  }
  for (i=1; i<=n; i++) {
    if (upper(i) != FLT_MAX) {
      dtmp = (upper(i) - xc(i)) / mu;
      if (dtmp < 0.0) maxmu = min(maxmu, 1.0e0/dtmp);
    }
  }
  maxmu = mu / min(maxmu, 1.0e1);
  mu = maxmu;
  *optout << "UpdateBarrierMultiplier: new mu = " << mu << "\n";
}

//---------------------------------------------------------------------------- 
// Compute the barrier part of the function value
//---------------------------------------------------------------------------- 
double OptBaNewton::compute_Barrier_Fvalue(double fcurrent, ColumnVector &xc)
{
  NLP2*         nlp   = nlprob2();
  int           i, n  = nlp->getDim();
  ColumnVector  upper = nlp->getConstraints()->getUpper();
  ColumnVector  lower = nlp->getConstraints()->getLower();
  double        dtmp1, dtmp2, fval;

  fval = fcurrent;
  for (i=1; i<=n; i++) {
    if (lower(i) != -FLT_MAX) {
      dtmp1 = xc(i) - lower(i);
      dtmp1 = log(dtmp1);
    } else dtmp1 = 0.0;
    if (upper(i) !=  FLT_MAX) {
      dtmp2 = upper(i) - xc(i);
      dtmp2 = log(dtmp2);
    } else dtmp2 = 0.0;
    fval = fval - mu * (dtmp2 + dtmp1);
  }
  return fval;
}

//---------------------------------------------------------------------------- 
// Compute the barrier part of the gradient 
//---------------------------------------------------------------------------- 
ColumnVector OptBaNewton::compute_Barrier_Gradient(ColumnVector &ingrad,
						    ColumnVector &xc)
{
  NLP2*         nlp = nlprob2();
  int           i, n= nlp->getDim();
  ColumnVector  upper = nlp->getConstraints()->getUpper();
  ColumnVector  lower = nlp->getConstraints()->getLower();
  ColumnVector  gk(n);
  double        dtmp1, dtmp2;
  
  gk = ingrad;
  for (i=1; i<=n; i++) {
    if (lower(i) != -FLT_MAX) dtmp1 = 1.0 / (xc(i) - lower(i)); else dtmp1=0.0;
    if (upper(i) !=  FLT_MAX) dtmp2 = 1.0 / (upper(i) - xc(i)); else dtmp2=0.0;
    gk(i) = gk(i) + mu * (dtmp2 - dtmp1);
  }
  return gk;
}

//---------------------------------------------------------------------------- 
// Compute the barrier part of the Hessian 
//---------------------------------------------------------------------------- 
SymmetricMatrix OptBaNewton::compute_Barrier_Hessian(SymmetricMatrix &H,
						      ColumnVector &xc)
{
  NLP2*           nlp = nlprob2();
  int              i, n= nlp->getDim();
  ColumnVector    upper = nlp->getConstraints()->getUpper();
  ColumnVector    lower = nlp->getConstraints()->getLower();
  double          dtmp1, dtmp2;
  SymmetricMatrix H2(n);

  H2 = H;
  for (i=1; i<=n; i++) {
    if (lower(i) != -FLT_MAX) {
      dtmp1  = xc(i) - lower(i);
      dtmp1  = 1.0 / (dtmp1 * dtmp1);
    } else dtmp1 = 0.0;
    if (upper(i) != FLT_MAX) {
      dtmp2  = upper(i) - xc(i);
      dtmp2  = 1.0 / (dtmp2 * dtmp2);
    } else dtmp2 = 0.0;
    H2(i,i) = H2(i,i) + mu * (dtmp1 + dtmp2);
  }
  return H2;
}

//---------------------------------------------------------------------------- 
// Compute the Search direction 
//---------------------------------------------------------------------------- 
ColumnVector OptBaNewton::computeSearch2(SymmetricMatrix &H, ColumnVector &g)
{
  NLP2*                 nlp = nlprob2();
  int                   n   = nlp->getDim();
  ColumnVector          sk(n);
  LowerTriangularMatrix L(n);

  L = MCholesky(H);
  sk = -(L.t().i()*(L.i()*g));
  return sk;
}

//---------------------------------------------------------------------------- 
// Compute the step length along pk
// The algorithm goes as follow :
//  1. initialize the interval of uncertainty
//  2. Calculate alpha(0) based on feasibility and the barrier function
//        alpha(0) = max(ComputeMaxStep(pk)+mu/(g'*pk), 0.5*ComputeMaxStep(pk))
//
//  Reference : "Line search procedures for the logarithmic barrier function"
//              by Murray and Wright (SIAM J. Optimization, May 1994)
//---------------------------------------------------------------------------- 
int OptBaNewton::computeStep(ColumnVector pk)
{
  NLP2*       nlp = nlprob2();
  int          n = nlp->getDim();
  double       alpha_bar, inner_gp, alpha_b, alpha_bar_plus, alpha;
  double       alpha_u, fplus, initslope, ftol, fnext, b, d, y;
  double       inner_gpnew, a, c, dtmp1, dtmp2, dtmp3;
  ColumnVector gplus(n), gnext(n);
  ColumnVector xc = nlp->getXc(), xplus(n);

  // Initialization

  ftol     = tol.getFTol();
  alpha_u  = 1.0;

  // compute alpha_bar (max step that can be taken w/o violating constraints)

  alpha_bar = computeMaxStep(pk);
  if(debug_)
    *optout << "ComputeStep : max alpha that can be taken = " << alpha_bar << "\n";

  // choose a reasonable step on the barrier function based on alpha_bar

  inner_gp = Dot(grad_barrier,pk);
  alpha_bar_plus = alpha_bar + mu / inner_gp;
  if (alpha_bar < FLT_MAX && alpha_bar_plus < 0.0) 
     alpha_b = max(alpha_bar_plus,0.5*alpha_bar);
  else if (alpha_bar < FLT_MAX && alpha_bar_plus >= 0.0) 
     alpha_b = 0.95 * alpha_bar;
  else alpha_b = FLT_MAX;
  if(debug_)
    *optout << "ComputeStep : best alpha that can be taken = " << alpha_b << "\n";

  // set initial upper bound for step length 

  alpha = (alpha_b < alpha_u) ? alpha_b : alpha_u;
  if(debug_) *optout << "ComputeStep : initial alpha = " << alpha << "\n";

  // Solve quadratic-logarithmic function to find optimal alpha
  // Q(x) = a + bx + cx^2 - mu * log(d-x)

  xplus = xc + pk * alpha;
  fnext = nlp->evalF(xplus);
  fplus = compute_Barrier_Fvalue(fnext,xplus);

  initslope = -Dot(grad_barrier,grad_barrier);
  if (fplus < fvalue_barrier + initslope * ftol) {
    nlp->setX(xplus);
    nlp->setF(fnext);
    nlp->evalG();
    nlp->evalH();
    Hessian = nlp->getHess();
    fcn_evals   = nlp->getFevals();
    grad_evals  = nlp->getGevals();
    step_length = alpha;
    return 0 ;
  } 

  gnext = nlp->evalG(xplus);
  gplus = compute_Barrier_Gradient(gnext,xplus);
  inner_gpnew = Dot(gplus,pk);

  if(debug_){
    *optout << "ComputeStep : fval (old, new) = " << fvalue_barrier << " " << fplus << "\n";
    *optout << "ComputeStep : g'p  (old, new) = " << inner_gp << " " << inner_gpnew << "\n";
  }
  y = scalarNewton(fvalue_barrier, inner_gp, fplus, inner_gpnew, alpha);
  if( debug_) *optout << "ComputeStep : y = " << y << "\n";
  if (y == 1) return -1;

  d = alpha / (1.0e0 - y);
  c = (inner_gpnew - inner_gp + mu / d - mu / (d - alpha)) / (2.0e0 * alpha);
  b = inner_gp - mu / d;
  a = fvalue_barrier + mu * log(d);
  if( debug_) *optout << "ComputeStep : a,b,c,d = " << a << " " << b << " " << c << " " << d << "\n";
  dtmp1 = 2.0 * c * d - b;
  dtmp2 = dtmp1 * dtmp1 + 8.0 * c * (mu + b * d);
  dtmp2 = sqrt(dtmp2);
  dtmp3 = 4 * c;
  if (c == 0.0) {
    *optout << "ComputeStep: error - divide by 0. \n";
    return -1;
  }
  alpha = (dtmp1 - dtmp2) / dtmp3;
  if(debug_){
    *optout << "ComputeStep : alpha chosen    = " << alpha << "\n";
    *optout << "ComputeStep : the other alpha = " << (dtmp1+dtmp2)/dtmp3 << "\n";
  }

  // Check to see if step OK

  xplus = xc + pk * alpha;
  fnext = nlp->evalF(xplus);
  fplus = compute_Barrier_Fvalue(fnext,xplus);
  if (fplus < fvalue_barrier + initslope * ftol) {
    nlp->setX(xplus);
    nlp->setF(fnext);
    nlp->evalG();
    nlp->evalH();
    Hessian = nlp->getHess();
    fcn_evals   = nlp->getFevals();
    grad_evals  = nlp->getGevals();
    step_length = alpha;
    return 0 ;
  } else {
    setMesg("OptBaNewton: Step does not satisfy sufficient decrease condition");
    return -1;
  }
}

//---------------------------------------------------------------------------- 
// Compute the maximum step allowed along the search direction sk
// before we hit a constraint
//--------------------------------------------------------------------------- 
double OptBaNewton::computeMaxStep(ColumnVector &sk)
{
  NLP2* nlp = nlprob2();
  int i, n= nlp->getDim();
  double gamma=FLT_MAX, delta;
  ColumnVector lower = nlp->getConstraints()->getLower();
  ColumnVector upper = nlp->getConstraints()->getUpper();
  ColumnVector xc    = nlp->getXc();

  double feas_tol = 1.e-3;

  for (i=1; i<=n; i++) {
    if      (sk(i) > 0.0e0) {
      delta = (upper(i)-xc(i)) / sk(i);
      if (delta <= feas_tol) {
	if (debug_)
	  *optout << "OptBaNewton::ComputeMaxStep: variable " << i << " hits upper constraint.\n";
      }
    } else if (sk(i) < 0.0e0) {
      delta = (lower(i)-xc(i)) / sk(i);
      if (delta <= feas_tol) {
	if (debug_)
	  *optout << "OptBaNewton::ComputeMaxStep: variable " << i << " hits lower constraint.\n";
      }
    }
    if (delta < 0.0e0) delta = 0.0e0;
    gamma = min(gamma,delta);
  }
  if( debug_) 
    *optout << "OptBaNewton::ComputeMaxStep: maximum step allowed = " << gamma << "\n";
  return gamma;
}

//---------------------------------------------------------------------------- 
// Use Newton's method to find the root of the equation
// f(z) = ln(z) + 0.5 * (1/z - z) - kappa = 0 where
//    kappa = (0.5*alpha*(phi1'+phi2')-(phi2-phi1)) / mu
//--------------------------------------------------------------------------- 
double OptBaNewton::scalarNewton(double phi1, double phi1_prime, double phi2, 
			          double phi2_prime, double alpha)
{
  int    convgd = 0;
  double y=1.0e-6, fval, f_prime, kappa;

  if(debug_){
    *optout << "ScalarNewton: phi1       = " << phi1 << "\n";
    *optout << "ScalarNewton: phi1_prime = " << phi1_prime << "\n";
    *optout << "ScalarNewton: phi2       = " << phi2 << "\n";
    *optout << "ScalarNewton: phi2_prime = " << phi2_prime << "\n";
    *optout << "ScalarNewton: alpha      = " << alpha << "\n";
  }
  kappa = (0.5 * alpha * (phi1_prime + phi2_prime)- phi2 + phi1) / mu;
  if (debug_) *optout << "ScalarNewton: kappa = " << kappa << "\n";
  if (kappa <= 0.0) {
    if (debug_) *optout << "ScalarNewton: Error - interpolant inadequate. \n";
    return 1.0e0;
  }
  while (convgd == 0) {
    fval  = log(y) + 0.5 * (1/y - y) - kappa;
    if (fabs(fval) < 1.0e-4) convgd = 1;
    else {
      f_prime = 1/y - 1/(2*y*y) - 0.5;
      y       = y - fval / f_prime;
    } 
  } 
  if (debug_) *optout << "ScalarNewton: y, f       = " << y << " " << fval << "\n";
  return y;
} 

//---------------------------------------------------------------------------- 
// Accept the step and update the parameters
//--------------------------------------------------------------------------- 
void OptBaNewton::acceptStep(int iter, int steptype)
{
  OptBCNewtonLike::acceptStep(iter, steptype);

  NLP2*      nlp   = nlprob2();
  ColumnVector xc  = nlp->getXc();
  ColumnVector gg  = nlp->getGrad();
  double       fv  = nlp->getF();

  Hess_barrier   = compute_Barrier_Hessian(Hessian,xc);
  grad_barrier   = compute_Barrier_Gradient(gg,xc);
  fvalue_barrier = compute_Barrier_Fvalue(fv,xc);
}

//---------------------------------------------------------------------------- 
// update the previous variables
//--------------------------------------------------------------------------- 
void OptBaNewton::setAsideCurrentVariables()
{
  NLP2* nlp   = nlprob2();
  xprev       = nlp->getXc();
  fprev       = nlp->getF();
  gprev       = nlp->getGrad();
  fprev_barrier = fvalue_barrier;
  gprev_barrier = grad_barrier;
}

SymmetricMatrix OptBaNewton::updateH(SymmetricMatrix&, int) 
{// JCM: Not used by optbanewton
  //     Default to evaluating the Hessian
  return nlprob()->evalH();
}

} // namespace OPTPP
