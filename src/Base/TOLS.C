//------------------------------------------------------------------------
// Copyright (C) 1993: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//------------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#if (defined(__sgi) || defined(__xlc__) || defined(__xlC__))
#define WANT_MATH
#else
#define WANT_STREAM
#define WANT_MATH
#endif

#ifdef HAVE_STD
#include <cfloat>
#include <cmath>
#include <cstring>
#else
#include <float.h>
#include <math.h>
#include <string.h>
#endif

#include "precisio.h"

#include "cblas.h"
#include "TOLS.h"

using namespace std;

namespace OPTPP {

//-------------------------------------------------------------------------
// Routines for Tolerance setting 
//-------------------------------------------------------------------------

void TOLS::setDefaultTol()  // set the default tolerances
{
  mcheps         = DBL_EPSILON;          // Machine epsilon
  max_step       = 1.e3;                 // Maximum step allowed
  min_step       = sqrt(mcheps);         // Minimum step allowed
  step_tol       = sqrt(mcheps);         // Convergence tolerance for step size
  fcn_tol        = sqrt(mcheps);         // Convergence tolerance for function
  con_tol        = sqrt(mcheps);         // Convergence tolerance for constraint
  grad_tol       = pow(mcheps,1.0/3.0);  // Convergence tolerance for gradient
  linesearch_tol = 1.e-4;                // Line search tolerance
  tr_size        = max_step;             // Initial trust region size;
  max_iter       = 100;                  // Maximum number of iterations
  max_backiter   = 5;                   // Maximum number of backtracks
  max_feval      = 1000;                 // Maximum number of functions 
}

void TOLS::setMaxStep(double x) {max_step = x;} // set the maximum step
void TOLS::setMinStep(double x) {min_step = x;} // set the maximum step
void TOLS::setStepTol(double x) {step_tol = x;} // set the step tolerance
void TOLS::setFTol(double x) {fcn_tol  = x;}    // set the Function tolerance
void TOLS::setCTol(double x) {con_tol  = x;}    // set the Constraint tolerance
void TOLS::setGTol(double x) {grad_tol = x;}    // set the gradient tolerance
void TOLS::setLSTol(double x) {linesearch_tol = x;} // set the step tolerance
void TOLS::setMaxIter(int k) {max_iter = k;}    // set the max number of iters
void TOLS::setMaxBacktrackIter(int k) {max_backiter = k;}    
           // set the max number of backtracks in the linesearch routine
void TOLS::setMaxFeval(int k) {max_feval = k;}  // set the max number of fevals

void TOLS::printTol()  // set the default tolerances
{
  cout << "\n\n==========  Tolerances  ===========\n\n";
  cout << "Machine Epsilon      = " <<  mcheps << "\n" ;
  cout << "Maximum Step         = " <<  max_step << "\n" ;
  cout << "Minimum Step         = " <<  min_step << "\n" ;
  cout << "Maximum Iter         = " <<  max_iter << "\n" ;
  cout << "Maximum Backtracks   = " <<  max_backiter << "\n" ;
  cout << "Maximum Fcn Eval     = " <<  max_feval << "\n" ;
  cout << "Step Tolerance       = " <<  step_tol << "\n" ;
  cout << "Function Tolerance   = " <<  fcn_tol << "\n" ;
  cout << "Constraint Tolerance = " <<  con_tol << "\n" ;
  cout << "Gradient Tolerance   = " <<  grad_tol << "\n" ;
  cout << "LineSearch Tolerance = " <<  linesearch_tol << "\n" ; 

}
void TOLS::printTol(ostream *tolout)  // set the default tolerances
{
  (*tolout) << "\n\n==========  Tolerances  ===========\n\n";
  (*tolout) << "Machine Epsilon      = " <<  mcheps << "\n" ;
  (*tolout) << "Maximum Step         = " <<  max_step << "\n" ;
  (*tolout) << "Minimum Step         = " <<  min_step << "\n" ;
  (*tolout) << "Maximum Iter         = " <<  max_iter << "\n" ;
  (*tolout) << "Maximum Backtracks   = " <<  max_backiter << "\n" ;
  (*tolout) << "Maximum Fcn Eval     = " <<  max_feval << "\n" ;
  (*tolout) << "Step Tolerance       = " <<  step_tol << "\n" ;
  (*tolout) << "Function Tolerance   = " <<  fcn_tol << "\n" ;
  (*tolout) << "Constraint Tolerance = " <<  con_tol << "\n" ;
  (*tolout) << "Gradient Tolerance   = " <<  grad_tol << "\n" ;
  (*tolout) << "LineSearch Tolerance = " <<  linesearch_tol << "\n" ; 

}

} // namespace OPTPP
