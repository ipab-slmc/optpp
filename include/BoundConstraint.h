#ifndef BoundConstraint_h
#define BoundConstraint_h

/*---------------------------------------------------------------------
  Copyright (c) 2001, Sandia Corporation.   Under the terms of Contract 
  DE-AC04-94AL85000, there is a non-exclusive license for use of this 
  work by or on behalf of the U.S. Government.
  -----------------------------------------------------------------------*/

#include "ConstraintBase.h"

/**
 * Simple Bounds Class
 * The standard form representation is x >= a.
 *
 * @author P.J. Williams, Sandia National Laboratories, pwillia@sandia.gov
 * @date Last modified 02/2007
 */

namespace OPTPP {

class BoundConstraint: public ConstraintBase {

protected:
  /// Number of constraints.
  int          	numOfCons_; 	
  /// Number of variables.
  int          	numOfVars_; 	
  /// Number of finite lower bounds
  int          	nnzl_;		
  /// Number of finite upper bounds
  int          	nnzu_;		

  /// Lower bounds on the variables 
  NEWMAT::ColumnVector 	lower_;   	
  /// Upper bounds on the variables 
  NEWMAT::ColumnVector 	upper_;   	
  /// Value of the variables 
  mutable NEWMAT::ColumnVector 	cvalue_;   	

  /// Indicator of a fixed variable 
  BoolVector   	fixedVar_;   	
  /// Indicator of a free variable 
  BoolVector   	freeVar_;   	

  /// Denotes whether a constraint is written in standard form or not
  const BoolVector   	stdForm_; 
  /// Type of constraint
  NEWMAT::ColumnVector		ctype_; 	
  
  /// Index vector of finite constraints
  OptppArray<int> 	constraintMappingIndices_;  

public:
/**
 * Default Constructor
 * @see BoundConstraint(int nc, const NEWMAT::ColumnVector& lower)
 * @see BoundConstraint(int nc, const NEWMAT::ColumnVector& bound, 
 *                            const BoolVector& bdFlag)
 * @see BoundConstraint(int nc, const NEWMAT::ColumnVector& lower, 
 *                            const NEWMAT::ColumnVector& upper)
 */
  BoundConstraint();

/**
 * @param nc an integer argument 
 * @param lower a NEWMAT::ColumnVector
 * @see BoundConstraint(int nc, const NEWMAT::ColumnVector& bound, 
 *                            const BoolVector& bdFlag)
 * @see BoundConstraint(int nc, const NEWMAT::ColumnVector& lower, 
 *                            const NEWMAT::ColumnVector& upper)
 * @note Assumes all constraints are in the form x >= c
 */
  BoundConstraint(int nc, const NEWMAT::ColumnVector& lower);
  
/**
 * @param nc an integer argument 
 * @param bound a NEWMAT::ColumnVector
 * @param bdFlag a BoolVector
 * @see BoundConstraint(int nc, const NEWMAT::ColumnVector& lower)
 * @see BoundConstraint(int nc, const NEWMAT::ColumnVector& lower, 
 *                            const NEWMAT::ColumnVector& upper)
 * @note User must specify whether the set of constraints is in standard form.  
 */
  BoundConstraint(int nc, const NEWMAT::ColumnVector& bound, const BoolVector& bdFlag);

/**
 * @param nc an integer argument 
 * @param lower a NEWMAT::ColumnVector
 * @param upper a NEWMAT::ColumnVector
 * @see BoundConstraint(int nc, const NEWMAT::ColumnVector& lower)
 * @see BoundConstraint(int nc, const NEWMAT::ColumnVector& bound, 
 *                            const BoolVector& bdFlag)
 * @note Includes lower and upper bounds
 */
  BoundConstraint(int nc, const NEWMAT::ColumnVector& lower, const NEWMAT::ColumnVector& upper);

/**
 * Destructor
 */
  ~BoundConstraint() {;}

/**
 * @return Number of constraints.
 */
  virtual int getNumOfCons() const {return numOfCons_;}

/**
 * @return Number of variables.
 */
  virtual int getNumOfVars() const {return numOfVars_;}

/**
 * @return Lower bounds on the variables.
 */
  virtual NEWMAT::ColumnVector  getLower()  const {return lower_;}

/**
 * Set lower bounds on the variables.
 */
  void setLower(NEWMAT::ColumnVector& x) {lower_ = x;}

/**
 * @return Upper bounds on the variables.
 */
  virtual NEWMAT::ColumnVector  getUpper()  const {return upper_;}

/**
 * Set upper bounds on the variables.
 */
  void  setUpper(NEWMAT::ColumnVector& x) {upper_ = x;}

/**
 * @return Type of constraint.
 */
  virtual NEWMAT::ColumnVector  getConstraintType()  const {return ctype_;}

/**
 * @return Value of constraint, in this case, the current iterate.
 */
  virtual NEWMAT::ColumnVector  getConstraintValue()  const {return cvalue_;}

/**
 * CPJW Placeholder!!
 * @return Constraint violation 
 */
  virtual NEWMAT::ColumnVector  getConstraintViolation()  const {return cvalue_;}

/**
 * @return Indices of constraints with finite bounds 
 */
  OptppArray<int> getConstraintMappingIndices() const 
	  { return constraintMappingIndices_; }  

/**
 * @return Fixed variable indicator. 
 */
  BoolVector getFixedVar() const {return fixedVar_;}

/**
 * @return Free variable indicator. 
 */
  BoolVector getFreeVar()  const {return freeVar_;}

/**
 * @return Standard form representation.
 */
  BoolVector getStdForm()  const {return stdForm_;}

/**
 * Takes two arguments and returns a bool.
 * @param xc a NEWMAT::ColumnVector
 * @param epsilon a real argument
 * @return true - constraints are feasible 
 * @return false - constraints are infeasible 
 */
  virtual bool amIFeasible(const NEWMAT::ColumnVector& xc, double epsilon) const;

/**
 * Takes no arguments and returns a bool.
 * @return true - lower  < upper 
 * @return false - lower > upper 
 */
  bool amIConsistent() const;

// Evaluation Methods

/**
 * Takes one argument and returns a NEWMAT::ColumnVector of reals.
 * @param xc a NEWMAT::ColumnVector
 * @return The residuals of the constraints.
 */
  virtual NEWMAT::ColumnVector evalResidual(const NEWMAT::ColumnVector& xc) const;
  virtual void evalCFGH(const NEWMAT::ColumnVector& xc) const;

private:
/**
 * Takes one argument and returns a real Matrix.
 * @param xc a ColumnVector
 * @return The gradient of the constraints.
 */
  virtual NEWMAT::Matrix evalGradient(const NEWMAT::ColumnVector& xc) const;
/**
 * Takes one argument and returns a SymmetricMatrix.
 * @param xc a ColumnVector
 * @return The Hessian of the constraints.
 */
  virtual NEWMAT::SymmetricMatrix evalHessian(NEWMAT::ColumnVector& xc) const;

/**
 * Takes one argument and returns an array of real SymmetricMatrices.
 * @param xc a ColumnVector
 * @param darg an integer argument
 * @return An array of constraint Hessians.
 */
  virtual OptppArray<NEWMAT::SymmetricMatrix> evalHessian(NEWMAT::ColumnVector& xc, int darg) const;
};

} // namespace OPTPP
#endif
