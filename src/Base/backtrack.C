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
#include <cmath>
#else
#include <math.h>
#endif

#include "Opt.h"

using NEWMAT::ColumnVector;

const   int MAXITER = 20;
//------------------------------------------------------------------------
// external subroutines referenced by this module 
//------------------------------------------------------------------------

#if !(defined(__GNUC__) && __GNUC__ >= 3)
extern "C" {
  double copysign(double, double);
}
#else
#ifdef CYGWIN
extern "C" {
  extern double copysign _PARAMS((double, double));
}
#endif
#endif

namespace OPTPP {

int backtrack(NLP1* nlp, ostream *fout, 
	      ColumnVector& search_dir, ColumnVector& sx,
	      double *stp, int itnmax, double ftol, double stpmax, double stpmin)
{
/****************************************************************************
 * subroutine backtrack
 * 
 *    FIND the next ITERATE BY LINE SEARCH. 

 * PARAMETERS 
 * ---------- 
 * N            --> DIMENSION OF PROBLEM 
 * X(N)         --> OLD ITERATE:   X[K-1] 
 * F            --> FUNCTION VALUE AT OLD ITERATE, F(X) 
 * G(N)         --> GRADIENT AT OLD ITERATE, G(X), OR APPROXIMATE 
 * P(N)         --> NON-ZERO NEWTON STEP 
 * XPLUS(N)     <--  NEW ITERATE X[K] 
 * FPLUS        <--  FUNCTION VALUE AT NEW ITERATE, F(XPLUS) 
 * FCN          --> NAME OF SUBROUTINE TO EVALUATE FUNCTION 
 * MXTAKE      <--  BOOLEAN FLAG INDICATING STEP OF MAXIMUM LENGTH USED 

 * STEPMX       --> MAXIMUM ALLOWABLE STEP SIZE 
 * STEPTL       --> RELATIVE STEP SIZE AT WHICH SUCCESSIVE ITERATES 
 *                  CONSIDERED CLOSE ENOUGH TO TERMINATE ALGORITHM 
 * SX(N)        --> DIAGONAL SCALING MATRIX FOR X 
 * IPR          --> DEVICE TO WHICH TO SEND OUTPUT 

 * INTERNAL VARIABLES 
 * ------------------ 
 * SLN              NEWTON LENGTH 
 * RELLENGTH              RELATIVE LENGTH OF NEWTON STEP 
 ***************************************************************************/
  /* Local variables */
  real disc;
  real a, b;
  int i;
  real fplus, t1, t2, t3, tlmbda, minlambda;
  real scl, rellength, sln, initslope;

  real lambda = 1.;
  real plmbda = -1.0;
  real pfplus = -1.0;
  double one;

  int    n = nlp->getDim();
  ColumnVector work(n), xc(n), xplus(n), grad(n);
  ColumnVector p(n);
  real  fx, tmp1, tmp2;
  int   iter;
  bool   debug = nlp->getDebug();
  bool modeOverride = nlp->getModeOverride();

  xc     = nlp->getXc();
  fx     = nlp->getF();
  grad   = nlp->getGrad();


  p      = search_dir;
  work   = sx.AsDiagonal()*p;
  sln    = Norm2(work);
 
    
  if (sln >= stpmax && sln != 0.0) { // STEP LONGER THAN MAXIMUM ALLOWED
    scl = stpmax / sln;
    p   = p   * scl;
    sln = stpmax;
  }

  initslope = Dot(grad, p);
  if (initslope >= 0.0) {
    *fout <<"backtrack: Initial search direction not a descent direction\n";
    *fout <<"backtrack: Replacing search direction with negative gradient\n";
    search_dir = -grad;
    initslope  = -Dot(grad, grad);
  }

  
  if (debug) {
    *fout <<"initslope = " << initslope << "\n";
    *fout <<"sln       = " << sln       << "\n";
    *fout <<"search dir\n";
    for( i=1; i<=n; i++) *fout << i << "\t" << p(i) << "\n";
  }
  
  rellength = 0.;
  
  for (i = 1; i <= n; ++i) {
    tmp1 = fabs(p(i));
    tmp2 = max(fabs(xc(i)),1.0/sx(i));
    rellength = max(rellength,tmp1)/tmp2;
  }
  
  minlambda = stpmin / rellength;
  
  /* CHECK IF NEW ITERATE SATISFACTORY.  GENERATE NEW LAMBDA IF NECESSARY. */
  iter = 0;
  //const  int MAXITER = itnmax;
  while (iter < MAXITER) {
    iter++;
    xplus = xc + p*(lambda);
    if (modeOverride) {
      nlp->setX(xplus);
      nlp->eval();
      fplus = nlp->getF();
    }
    else
      fplus = nlp->evalF(xplus);
    if (debug) {
      *fout 
	<< "iter:" << iter 
	<< "fplus  = " << fplus
        << "lambda = " << lambda << "\n";
    }
    
    // Is this an acceptable step ?

    if (fplus <= fx + initslope * ftol * lambda) { 
      // Yes !
      if (debug) *fout <<"Accept\n";
      *stp = lambda;
      if (!modeOverride) {
	nlp->setX(xplus);
	nlp->setF(fplus);
	nlp->evalG();
      }
      if (iter == 1) return(Newton_Step); // Newton Step
      else return(Backtrack_Step); // Backtrack Step
    }

    // Insufficient decrease. Try to backtrack
      
    if (lambda < minlambda) { // Step size smaller than min allowed
      *stp = lambda;
      nlp->setX(xc);
      nlp->setF(fx);
      nlp->setGrad(grad);
      return (-1); // Error
    }  
    else { // Compute new lambda 

      if (iter == 1) {/* FIRST BACKTRACK: QUADRATIC FIT */
	if (debug) *fout <<"First Backtrack\n";
	tlmbda = -initslope / ((fplus - fx - initslope) * 2.);
      }
      else {/* ALL SUBSEQUENT BACKTRACKS: CUBIC FIT */
	if (debug) *fout << "More Backtrack\n";
	t1 = fplus - fx - lambda * initslope;
	t2 = pfplus - fx - plmbda * initslope;
	t3 = 1. / (lambda - plmbda);
	a  = t3 * (t1 / (lambda * lambda) - t2 / (plmbda * plmbda));
	b  = t3 * (t2 * lambda / (plmbda * plmbda) - 
		   t1 * plmbda / (lambda * lambda));
	disc = b * b - a * 3. * initslope;
	if (disc > b * b) {
	  // ONLY ONE POSITIVE CRITICAL POINT, MUST BE MINIMUM
	  one = copysign(1.0, a);
	  tlmbda = (-b + one * sqrt(disc)) / (a * 3.);
	}
	else {
	  // BOTH CRITICAL POINTS POSITIVE, FIRST IS MINIMUM 
	  one = copysign(1.0, a);
	  tlmbda = (-(real)b - one * sqrt(disc)) / (a * 3.);
	}
	if (tlmbda > lambda * .5) tlmbda = lambda * .5;
      }
      
      plmbda = lambda;
      pfplus = fplus;
      if (tlmbda < .1 * lambda) lambda *= .1;
      else lambda = tlmbda;
    }
  }

  nlp->setX(xc);
  nlp->setF(fx);
  nlp->setGrad(grad);
  return (-1); // Too many iterations
}

} // namespace OPTPP
