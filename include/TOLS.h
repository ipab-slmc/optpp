#ifndef tols_h
#define tols_h


/*----------------------------------------------------------------------
 Copyright (c) 2001, Sandia Corporation.   Under the terms of Contract 
 DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 work by or on behalf of the U.S. Government.
 
 J.C. Meza, Sandia National Laboratories, meza@ca.sandia.gov
 ----------------------------------------------------------------------*/

#include <iostream>
#include <fstream>

using std::ostream;

namespace OPTPP {

/**
 * TOLS is the Base Class for Tolerances which will be used 
 * in the optimization methods.
 *
 * @author J. C. Meza, Sandia National Laboratories, meza@ca.sandia.gov
 * @note Modified by P.J. Williams, pwillia@sandia.gov 
 * @date Last modified 02/2006
 */

class TOLS{
  double  mcheps;         ///< Machine epsilon
  double  max_step;       ///< Maximum step allowed in computestep
  double  min_step;       ///< Minimum step allowed in computestep
  double  step_tol;       ///< Step tolerance used for convergence test
  double  fcn_tol;        ///< Function tolerance used for convergence test
  double  con_tol;        ///< Constraint tolerance used for convergence test
  double  grad_tol;       ///< Gradient tolerance used for convergence test
  double  linesearch_tol; ///< Line search tolerance
  double  tr_size;        ///< Initial trust region size
  int     max_iter;       ///< Maximum number of iterations allowed
  int     max_backiter;   ///< Maximum number of backtracks allowed in lnsrch
  int     max_feval;      ///< Maximum number of function evaluations allowed

public:
  /// Default Constructor
  TOLS() {}           		
  /// Copy Constructor
  TOLS(TOLS& t) {*this = t;}    
  /// Destructor
  virtual ~TOLS() {}       	

  /// Assignment Operator 
  void operator= (TOLS& t) {
    mcheps         = t.mcheps;
    max_step       = t.max_step;
    min_step       = t.min_step;
    step_tol       = t.step_tol;
    fcn_tol        = t.fcn_tol;
    con_tol        = t.con_tol ;
    grad_tol       = t.grad_tol;
    linesearch_tol = t.linesearch_tol;
    max_iter       = t.max_iter;
    max_backiter   = t.max_backiter;
    max_feval      = t.max_feval;
  }
  /// Set all default tolerances
  void    setDefaultTol();	
  /// Set maximum allowable steplength
  void    setMaxStep(double x);	
  /// Set minimum allowable steplength
  void    setMinStep(double x);	
  /// Set tolerance used for step convergence test
  void    setStepTol(double x); 
  /// Set tolerance used in function convergence test
  void    setFTol(double x); 	
  /// Set tolerance used for constraint feasibility test
  void    setCTol(double x);	
  /// Set tolerance used in gradient convergence test
  void    setGTol(double x);	
  /// Set linesearch tolerance
  void    setLSTol(double x);	
  /// Set trust-region radius
  void    setTRSize(double x);	
  /// Set maximum number of iterations
  void    setMaxIter(int k);	
  /// Set maximum backtrack iterations
  void    setMaxBacktrackIter(int k);	
  /// Set maximum allowable function evaluations
  void    setMaxFeval(int k);	

  /**
   * @return Maximum allowable steplength
   */
  double  getMaxStep()  const {return max_step;}

  /**
   * @return Minimum allowable steplength
   */
  double  getMinStep()  const {return min_step;}

  /**
   * @return Tolerance used for step convergence test
   */
  double  getStepTol()  const {return step_tol;}

  /**
   * @return Tolerance used in function convergence test
   */
  double  getFTol()     const {return fcn_tol;}

  /**
   * @return Tolerance used for constraint feasibility convergence test
   */
  double  getCTol()     const {return con_tol;}

  /**
   * @return Tolerance used in gradient convergence test
   */
  double  getGTol()     const {return grad_tol;}

  /**
   * @return Linesearch tolerance
   */
  double  getLSTol()    const {return linesearch_tol;}

  /**
   * @return Trust-region radius 
   */
  double  getTRSize()   const {return tr_size;}

  /**
   * @return Maximum number of iterations allowed
   */
  int     getMaxIter()  const {return max_iter;}

  /**
   * @return Maximum number of backtrack iterations allowed
   */
  int     getMaxBacktrackIter()  const {return max_backiter;}

  /**
   * @return Maximum number of function evaluations allowed
   */
  int     getMaxFeval() const {return max_feval;}

  /**
   * Please use setOutputFile with version 1.6 and higher
   * See OptimizeClass in Opt.h
   * void setOutput(ofstream& f) {cerr << "setOutput no longer supported\n";}
   */


  /**
   * Print the current tolerances being used
   */
  void printTol();

  /**
   * Print the current tolerances to file ofstream 
   */
  void printTol(ostream *);
};

} // namespace OPTPP

#endif
