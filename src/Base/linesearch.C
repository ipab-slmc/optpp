//------------------------------------------------------------------------
// Copyright (C) 1993: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//------------------------------------------------------------------------
// Modified by Ricardo Oliva, raoliva@lbl.gov
// Last modification date: March 2004.
//------------------------------------------------------------------------
#define WANT_MATH

#include "Opt.h"
using NEWMAT::ColumnVector;

namespace OPTPP {

int linesearch(NLP1* nlp, ostream *optout, 
	       ColumnVector& search_dir, ColumnVector& sx,
	       double *step_length, double stpmax, double stpmin,
	       int itnmax, double ftol, double xtol, double gtol)
{
/****************************************************************************
 *   subroutine linesearch
 *
 *   Purpose
 *   find a step which satisfies the Goldstein-Armijo line search conditions
 *
 *   Parameters
 *     nlp  -->  pointer to nonlinear problem object 
 *
 *     s    -->  Vector of length n which specifies the 
 *               search direction
 *
 *     step_length  <--> is a nonnegative variable. on input stp contains an 
 *               initial estimate of a satisfactory step. 
 *               On output stp contains the final estimate. 
 *
 *     *nfev <--  variable set to the number of calls to fcn. 
 *
 *     maxfev --> positive input variable. termination 
 *                occurs when the number of calls to fcn is at least 
 *                maxfev by the end of an iteration
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
 *     stpmin and stpmax are nonnegative input variables which 
 *       specify lower and upper bounds for the step. (in this reverse 
 *       communication implementatin they are defined in a common 
 *       statement). 
 *       stpmin Default Value = 1.e-9
 *       stpmax Default Value = 1.e3
 *
 *     itnmax  -->  iteration limit for the linesearch algorithm
 *                default Value = 5
 *
 *
 *   Initial version    Juan Meza November 1994
 *
 *   Modified by Ricardo Oliva. March 2004.
 *   -- The step_length is no longer set to 1.0 here. Calling programs
 *      now pass the suggested initial step size that mcsrch expects.
 *   -- The scaled stpmax is adjusted to prevent stpmax < stpmin
 *
 *****************************************************************************/

  int expensive_function;

  /* local variables */
  
  int step_type;
  //  *step_length = 1.0; 
  expensive_function = nlp->getIsExpensive();

  
  if (expensive_function)
    step_type = backtrack(nlp, optout, search_dir, sx, step_length, 
			  itnmax, ftol, stpmax, stpmin);
  else {

    // Since stpmax is the absolute max step size allowed 
    // we need to scale stpmax because mcsrch needs stpmax to
    // be relative to the search_dir

    // But the scaling
    //       stpmax = stpmax / snorm
    // could make stpmax < stpmin, which is rejected by mcsrch().
    //. -RAO

    double snorm = Norm2(search_dir);
    stpmax = stpmax / snorm;
    if (stpmax<stpmin) {
//      cerr << "Warning: in linesearch(): stpmax/stpnorm < stpmin\n";
      stpmax = 10*stpmin;      
    }
    step_type = mcsrch(nlp, search_dir, optout, step_length, itnmax,
	               ftol, xtol, gtol, stpmax, stpmin);
  }
  return(step_type);
}

} // namespace OPTPP
