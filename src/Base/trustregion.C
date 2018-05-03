//------------------------------------------------------------------------
// Copyright (C) 1993: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//------------------------------------------------------------------------

#include "Opt.h"
#include "ioformat.h"

using NEWMAT::ColumnVector;
using NEWMAT::SymmetricMatrix;

namespace OPTPP {

int trustregion(NLP1* nlp, ostream *fout, 
		SymmetricMatrix& H, ColumnVector& search_dir, 
		ColumnVector& sx, real& TR_size, real& step_length, 
		real stpmax, real stpmin)
{
/****************************************************************************
 *   subroutine trustregion
 *
 *   Purpose
 *   find a step which satisfies the Goldstein-Armijo line search conditions
 *
 *        Compute the dogleg step
 *        Compute the predicted reduction, pred, of the quadratic model
 *        Compute the actual reduction, ared
 *        IF ared/pred > eta
 *          THEN x_vec = x_vec + d_vec
 *               TR_size >= TR_size
 *               Compute the gradient g_vec at the new point
 *          ELSE TR_size < ||d_vec||
 *
 *   Parameters
 *     nlp          -->  pointer to nonlinear problem object 
 *
 *     search_dir   -->  Vector of length n which specifies the 
 *                       newton direction on input. On output it will
 *                       contain the step 
 *
 *     step_length  <-- is a nonnegative variable. 
 *                      On output step_length contains the step size taken
 *
 *     ftol  -->  default Value = 1.e-4
 *                ftol should be smaller than 5.e-1
 *                suggested value = 1.e-4 for newton methods
 *                                = 1.e-1 for more exact line searches
 *     xtol  -->  default Value = 2.2e-16
 *     gtol  -->  default Value = 0.9
 *                gtol should be greater than 1.e-4 
 *
 *                termination occurs when the sufficient decrease 
 *                condition and the directional derivative condition are 
 *                satisfied. 
 *
 *     stpmin and TR_size are nonnegative input variables which 
 *       specify lower and upper bounds for the step.
 *       stpmin Default Value = 1.e-9
 *
 *   Initial version    Juan Meza November 1994
 *
 *
 *****************************************************************************/

  // Local variables

  int n = nlp->getDim();
  bool debug = nlp->getDebug();
  bool modeOverride = nlp->getModeOverride();

  ColumnVector tgrad(n), newton_dir(n), tvec(n), xc(n), xtrial(n);
  real fvalue, fplus, dnorm;
  real eta1 = .001;
  real eta2 = .1;
  real eta3 = .75;
  real rho_k;
  int iter = 0;
  int iter_max = 100;
  real dd1, dd2;
  real ared, pred;
  int dog_step;
  real TR_MAX = stpmax;
  static char *steps[] = {"C", "D", "N", "B"};
  static bool accept;

  //
  // Initialize variables
  //

  fvalue = nlp->getF();
  xc     = nlp->getXc();
  step_length = 1.0;
  tgrad      = nlp->getGrad();
  newton_dir = search_dir;

  if (debug) {
    *fout << "\n***************************************";
    *fout << "***************************************\n";
    *fout << "\nComputeStep using trustregion\n";
    *fout << "\tStep   ||step||       ared          pred        TR_size \n";
  }

  while (iter < iter_max) {
    iter++;
    //
    // Compute the dogleg step 
    //
    search_dir = newton_dir;
    dog_step = dogleg(nlp, fout, H, tgrad, search_dir, sx, dnorm, TR_size, stpmax);
    step_length = dnorm;
    //
    // Compute pred = -g'd - 1/2 d'Hd
    //
    dd1  = Dot(tgrad,search_dir);
    tvec = H * search_dir;
    dd2  = Dot(search_dir,tvec);
    pred = -dd1 - dd2/2.0;
    
    //
    // Compute objective function at trial point
    //

    xtrial = xc + search_dir;
    if (modeOverride) {
      nlp->setX(xtrial);
      nlp->eval();
      fplus = nlp->getF();
    }
    else
      fplus  = nlp->evalF(xtrial);
    ared   = fvalue - fplus;

    //
    // Should we take this step ?
    //

    rho_k = ared/pred;

    accept = false;
    if (rho_k >= eta1) accept = true;

    //
    // Update the trust region
    //

    if (accept) {

      //
      // Do the standard updating
      //
      
      if (rho_k <= eta2) {

	// Model is just sort of bad
	// New trust region will be TR_size/2
	
	if (debug) {
	  *fout << "trustregion: rho_k  = " << e(rho_k,14,4) 
	      << " eta2 = " << e(eta2,14,4) << "\n";
	}
	TR_size = step_length / 2.0;

	if (TR_size < stpmin) {
	  *fout << "***** Trust region too small to continue.\n";
	  nlp->setX(xc);
	  nlp->setF(fvalue);
	  nlp->setGrad(tgrad);
	  return(-1);
	}
      }

      else if ((eta3 <= rho_k) && (rho_k <= (2.0 - eta3))) {

	// Model is PRETTY good
	// Double trust region

	if (debug) {
	  *fout << "trustregion: rho_k = " << e(rho_k,14,4) 
	      << " eta3 = " << e(eta3,14,4) << "\n";
	}
	TR_size = min(2.0*TR_size, TR_MAX);
      }

      else {

	//
	// All other cases
	//

	TR_size = max(2.0*step_length,TR_size);
	TR_size = min(TR_size, TR_MAX);
	if (debug) {
	  *fout << "trustregion: rho_k = " << e(rho_k,14,4) << "\n";
	}
      }
    }
    
    else {

      // Model is REALLY bad
      //

      TR_size = step_length/10.0;

      if (debug) {
	*fout << "trustregion: rho_k = " << e(rho_k,14,4) << "n"
	    << " eta1 = " << e(eta1,14,4) << "\n";
      }
      if (TR_size < stpmin) {
	*fout << "***** Trust region too small to continue.\n";
	nlp->setX(xc);
	nlp->setF(fvalue);
	nlp->setGrad(tgrad);
	return(-1);
      }
    }

    //
    //  Accept/Reject Step
    //      

    if (accept) {
      //
      // Update x, f, and grad
      //

      if (!modeOverride) {
	nlp->setX(xtrial);
	nlp->setF(fplus);
	nlp->evalG();
      }

      if (debug) {
	*fout << "\t Step     ||step||       ared          pred        TR_size \n";
	*fout << "Accept  " << steps[dog_step] << e(step_length,14,4) 
	      << e(ared,14,4) << e(pred,14,4) << e(TR_size,14,4) << "\n";
      }
      return(dog_step);
    }

    else {
      //
      // Reject step
      //
	
      if (debug) {
	*fout << "\t Step     ||step||       ared          pred        TR_size \n";
	*fout << "Reject  " << steps[dog_step] << e(step_length,14,4) 
	      << e(ared,14,4)  << e(pred,14,4)  << e(TR_size,14,4) << "\n";
      }
    }
  }
  nlp->setX(xc);
  nlp->setF(fvalue);
  nlp->setGrad(tgrad);
  return(-1); // Too many iterations
}

} // namespace OPTPP
