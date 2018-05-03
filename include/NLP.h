#ifndef NLP_h
#define NLP_h

/*----------------------------------------------------------------------
 Copyright (c) 2001, Sandia Corporation.   Under the terms of Contract 
 DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 work by or on behalf of the U.S. Government.

 P.J. Williams, Sandia National Laboratories, pwillia@sandia.gov
 ----------------------------------------------------------------------*/


#include "NLPBase.h"
#include "OptppSmartPtr.h"

using std::ostream;

namespace OPTPP {

/**
 * NLP is a handle class for NLPBase.
 * This class is an interface to NLP0-2 for evaluating nonlinear functions.
 *
 * @author P.J. Williams, Sandia National Laboratories, pwillia@sandia.gov
 * @date Last modified 03/2007
 */


class NLP{

private:
  SmartPtr<NLPBase> ptr_;  ///< Pointer to an NLPBase object

public:
  /**
   * Default Constructor
   */
  NLP();
  /**
   * @param base pointer to an NLPBase object 
   */
  NLP(NLPBase* base); 

  /// Set the ith component of the vector x 
  void setX(const int i, const real& x);

  /// Set the current point
  void setX(const NEWMAT::ColumnVector& x);

  /// Set the function value
  void setF(const real& fx);

  /**
   * e = 1, simple backtracking linesearch is used in opt. algorithm
   *
   * e = 0, More-Thuente linesearch is used in opt. algorithm
   */
  void setIsExpensive(const int e);

  /// Set the ith component of the function accuracy 
  void setFcnAccrcy(const int i, const real& accrcy);

  /// Set the function accuracy  
  void setFcnAccrcy(const NEWMAT::ColumnVector& accrcy);

  /**
   * @return Problem Dimension
   */
  int  getDim()         const;

  /**
   * @return Number of function evaluations taken 
   */
  int  getFevals()      const;

  /**
   * @return   1,  Function evaluation is expensive
   *
   * @return  = 0,  Function evaluation is inexpensive
   */
  int  getIsExpensive() const;

  /**
   * @return The current value of function 
   */
  real getF()           const;

  /**
   * @return User-specified function accuracy
   */
  NEWMAT::ColumnVector getFcnAccrcy()   const;

  /**
   * @return The current value of x 
   */
  NEWMAT::ColumnVector getXc()  const;

  /**
   * @return CPU time used 
   */
  real getFcnTime()     const;

  /**
   * @return Total number of constraints
   */
  int getNumOfCons() const; 

  /**
   * @return Total number of nonlinear constraints
   */
  int getNumOfNLCons() const; 

  /**
   * @return  = true,  Problem is constrained
   *
   * @return  = false, Problem is unconstrained
   */
  bool hasConstraints() ; 

  /// Print value of constraints to the screen 
  void printConstraints(); 

  /// Set debug parameter = true
  void setDebug();

  /**
   * @return  = true,  debug output statements printed 
   *
   * @return  = false, debug output statements are not printed 
   */
  bool getDebug() const;

  /// Reset parameter values 
  void reset();  			

  /// Initialize selected function 
  void initFcn();  			

  /// Evaluate the function 
  real evalF();  			

  /// Evaluate the function at x 
  real evalF(const NEWMAT::ColumnVector& x);  	

  /// Evaluate the gradient 
  NEWMAT::ColumnVector evalG();  		

  /// Evaluate the gradient at x 
  NEWMAT::ColumnVector evalG(const NEWMAT::ColumnVector& x);

  /// Evaluate Hessian 
  NEWMAT::SymmetricMatrix evalH();  		

  /// Evaluate Hessian at x 
  NEWMAT::SymmetricMatrix evalH(NEWMAT::ColumnVector& x);

  /// Evaluate the function, gradient, and Hessian 
  void eval();  			

  /// Evaluate the constraints at x 
  NEWMAT::ColumnVector evalCF(const NEWMAT::ColumnVector& x); 

  /// Evaluate the constraint gradient at x 
  NEWMAT::Matrix evalCG(const NEWMAT::ColumnVector& x);  

  /// Evaluate the constraint Hessian at x 
  NEWMAT::SymmetricMatrix evalCH(NEWMAT::ColumnVector& x);   

  /// Evaluate the constraint Hessian at x 
  OptppArray<NEWMAT::SymmetricMatrix> evalCH(NEWMAT::ColumnVector& x, int darg);   
  void evalC(const NEWMAT::ColumnVector& x); 

  /// Print status of the nonlinear function to the screen 
  void printState(char *); 
  /// Print status of the nonlinear function to file 
  void fPrintState(ostream *, char *); 

};

} // namespace OPTPP
#endif

