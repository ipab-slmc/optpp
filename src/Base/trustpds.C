
//--------------------------------------------------------------------
// Copyright (C) 1993: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//--------------------------------------------------------------------

#include "Opt.h"
#include "ioformat.h"

using NEWMAT::ColumnVector;
using NEWMAT::SymmetricMatrix;

namespace OPTPP {

int trustpds(NLP1* nlp, ostream *fout, SymmetricMatrix& H,
	     ColumnVector& search_dir, ColumnVector& sx, real&
	     TR_size, real& step_length, real stpmax, real stpmin,
	     int searchSize)
{
  /*******************************************************************
   *
   * subroutine trustpds
   *
   * Purpose
   *
   *    find a step which satisfies the Goldstein-Armijo line search
   *    conditions
   *
   *    Compute the dogleg step
   *    Compute the predicted reduction, pred, of the quadratic model
   *    Compute the actual reduction, ared
   *    IF ared/pred > eta
   *       THEN x_vec = x_vec + d_vec
   *          TR_size >= TR_size
   *          Compute the gradient g_vec at the new point
   *    ELSE TR_size < ||d_vec||
   *
   * Parameters
   *
   *    nlp          -->  pointer to nonlinear problem object 
   *
   *    search_dir   -->  Vector of length n which specifies the
   *                      newton direction on input. On output it will
   *                      contain the step 
   *
   *    step_length  <-- is a nonnegative variable. 
   *                     On output step_length contains the step size
   *                     taken
   *
   *    ftol  -->  default Value = 1.e-4
   *               ftol should be smaller than 5.e-1
   *               suggested value = 1.e-4 for newton methods
   *                               = 1.e-1 for more exact line
   *                                 searches
   *
   *    xtol  -->  default Value = 2.2e-16
   *
   *    gtol  -->  default Value = 0.9
   *               gtol should be greater than 1.e-4 
   *
   *               termination occurs when the sufficient decrease 
   *               condition and the directional derivative condition
   *               are satisfied. 
   *
   *    stpmin and TR_size are nonnegative input variables which 
   *      specify lower and upper bounds for the step.
   *      stpmin Default Value = 1.e-9
   *
   * Initial version    Juan Meza November 1994
   *
   *******************************************************************/

  // Local variables
  int n = nlp->getDim();
  bool debug = nlp->getDebug();
  SpecOption SpecTmp = nlp->getSpecOption();

  ColumnVector tgrad(n), newton_dir(n), tvec(n), xtrial(n), xsave(n);
  real fsave, fplus, dnorm;
  real eta1 = .001;
  real eta2 = .1;
  real eta3 = .75;
  real rho_k;

  int iter = 0;
  int iter_max = 100;
  real dd1, dd2;
  real ared, pred;
  int pds_step;
  real TR_MAX = stpmax;
  static char *steps[] = {"C", "P", "N"};
  double pds_radius;
  static bool init_value = true;
  static bool accept;

  //
  // Initialize variables
  //

  step_length = 1.0;
  tgrad      = nlp->getGrad();
  newton_dir = search_dir;

  if (debug) {
    *fout << "\n***************************************";
    *fout << "***************************************\n";
    *fout << "\nComputeStep using trustpds\n";
  }

  xsave      = nlp->getXc();
  fsave      = nlp->getF();

  while (iter < iter_max) {
    iter++;

    //
    // Compute the PDS step 
    //

    search_dir = newton_dir;

//  PJW 12/13/99
    pds_step = pdsstep(nlp, fout, H, tgrad, search_dir, sx, dnorm,
		       TR_size, stpmax, pds_radius, init_value, searchSize);
//  PJW 12/13/99

    if (debug ) {
      (*fout) << "\ntrustpds: Returned from pds_step\n";
      (*fout) << "trustpds: step_length = " << e(dnorm, 14,4) << "\n";
    }

    step_length = dnorm;

    if (step_length == 0.0) {
      xtrial = xsave;
      fplus  = fsave;
      ared   = fsave - fplus;
      rho_k  = 0.0;
    }
    else {

      //
      // Compute Actual and Predicted reduction
      //

      xtrial = nlp->getXc();
      fplus  = nlp->getF();
      ared   = fsave - fplus;
      
      // Compute pred = -g'd - 1/2 d'Hd
      //

      dd1  = Dot(tgrad,search_dir);
      tvec = H * search_dir;
      dd2  = Dot(search_dir,tvec);
      pred = -dd1 - dd2/2.0;
    }
    
    if (debug ) {
      (*fout) << "trustpds:Trial point\n";
      FPrint(fout, xtrial);
      (*fout) << "\ntrustpds: fsave  = " << e(fsave,14,4) <<"\n";
      (*fout) << "trustpds: fplus  = " << e(fplus,14,4)  <<"\n";
    }

    //
    // Should we take this step ?
    // Note: Because of the step that PDS takes it is possible
    //       for pred to be negative. However, if we actually get
    //       a reduction, we should accept the step
    //       Don't think this won't happen, because it DID in our test cases

    accept = false;

    if (pred <= 0.0) {
      if (ared > 0.0) accept = true;
    }
    else {
      rho_k = ared / pred;
      if (rho_k >= eta1) accept = true;
    }

    //
    // Update the trust region
    //

    if (accept) {

      //
      // Need to check for special case of negative pred
      // For lack of anything better to do, just leave
      // the trust region size alone !

      if (pred <= 0.0) {
	TR_size = max(2.0*step_length,TR_size);
	TR_size = min(TR_size, TR_MAX);
	if (debug) {
	  *fout << "trustpds: pred <= 0.0, ared > 0. \n";
	}
      }

      //
      // Predicted decrease positive so do the standard updating
      //
      
      else if (rho_k <= eta2) {

	// Model is just sort of bad
	// New trust region will be the minimum of TR_size/2 and 
	// length of longest edge returned by PDS
	
	if (debug) {
	  *fout << "trustpds: rho_k  = " << e(rho_k,14,4) 
	  << " eta2 = " << e(eta2,14,4) << "\n";
	}
	
	TR_size = step_length / 2.0;

	if (TR_size < stpmin) {
	  *fout << "***** Trust region too small to continue.\n";
	  return(-1);
	}
      }
      else if ((eta3 <= rho_k) && (rho_k <= (2.0 - eta3))) {

	// Model is PRETTY good
	// Double trust region

	if (debug) {
	  *fout << "trustpds: rho_k = " << e(rho_k,14,4) 
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
	  *fout << "trustpds: rho_k = " << e(rho_k,14,4) << "\n";
	}
      }
    }
    else {

      // Model is REALLY bad or 
      // New trust region will be the minimum of TR_size/10 and 
      // length of longest edge returned by PDS
      //

      if (step_length == 0.0)
	TR_size = min(TR_size, pds_radius)/10.0;
      else
	TR_size = step_length/10.0;

      if (debug) {
	*fout << "trustpds: ared = " << e(ared,14,4) 
	    << " pred = " << e(pred,14,4) << "\n"
	    << " eta1 = " << e(eta1,14,4) << "\n";
      }
      if (TR_size < stpmin) {
	*fout << "***** Trust region too small to continue.\n";
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

      if (init_value) init_value = false;

      nlp->setX(xtrial);
      nlp->setF(fplus);
      nlp->setSpecOption(NoSpec);
      nlp->evalG();
      nlp->setSpecOption(SpecTmp);

      if (debug) {
	*fout << "\t Step     ||step||       ared          pred        TR_size \n";
	*fout << "Accept  " << steps[pds_step] << e(step_length,14,4) 
	  << e(ared,14,4) << e(pred,14,4) << e(TR_size,14,4) << "\n";
      }

      return(pds_step);
    }
    else {

      //
      // Reject step
      //

      nlp->setX(xsave);
      nlp->setF(fsave);
	
      if (debug) {
	*fout << "\t Step     ||step||       ared          pred        TR_size \n";
	*fout << "Reject  " << steps[pds_step] << e(step_length,14,4) 
	  << e(ared,14,4)  << e(pred,14,4)  << e(TR_size,14,4) << "\n";
      }
    }
  }

  return(-1); // Too many iterations
}

} // namespace OPTPP
