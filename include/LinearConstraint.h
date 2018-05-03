#ifndef LinearConstraint_h
#define LinearConstraint_h

/*----------------------------------------------------------------------
 Copyright (c) 2001, Sandia Corporation.   Under the terms of Contract 
 DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 work by or on behalf of the U.S. Government.

 P.J. Williams, Sandia National Laboratories, pwillia@sandia.gov
 ----------------------------------------------------------------------*/


#include "ConstraintBase.h"

/**
 * LinearConstraint is a derived class of ConstraintBase.
 * LinearConstraint is an abstract class, which
 * provides common data and functionality 
 * to LinearEquation and LinearInequality.
 *
 * @author P.J. Williams, Sandia National Laboratories, pwillia@sandia.gov
 * @date Last modified 02/06/2007
 */

namespace OPTPP {

class LinearConstraint: public ConstraintBase {
protected:
  /// Number of constraints	
  int          	numOfCons_;	
  /// Number of variables	
  int          	numOfVars_;	
  /// Number of finite lower bounds	
  int	  	nnzl_;		
  /// Number of finite upper bounds	
  int	  	nnzu_;		
  /// Matrix representation of the constraints
  NEWMAT::Matrix       A_;		

  /// Matrix-vector product 
  NEWMAT::ColumnVector Ax_;		
  /// Lower bounds on the variables
  NEWMAT::ColumnVector lower_;		
  /// Upper bounds on the variables
  NEWMAT::ColumnVector upper_;		
  /// The value of the linear constraints 
  mutable NEWMAT::ColumnVector cvalue_;		
  /// The constraint violation, zero if all constraints satisfied 
  mutable NEWMAT::ColumnVector cviolation_;		

  /// Indices of finite bounds
  OptppArray<int> constraintMappingIndices_;
  /// Standard form representation of constraints
  bool   	stdForm_;	

public:
// Constructors
/**
 * Default Constructor
 */
  LinearConstraint();
/**
 * @param A a real Matrix
 * @see LinearConstraint(const Matrix& A, const NEWMAT::ColumnVector& b)
 * @see LinearConstraint(const Matrix& A, const NEWMAT::ColumnVector& b, 
 *                                 const bool rowFlag)
 * @see LinearConstraint(const Matrix& A, const NEWMAT::ColumnVector& lower,
 *                            const NEWMAT::ColumnVector& upper)
 */
  LinearConstraint(const NEWMAT::Matrix& A);
/**
 * @param A a real Matrix
 * @param b a ColumnVector
 * @see LinearConstraint(const Matrix& A);
 * @see LinearConstraint(const Matrix& A, const NEWMAT::ColumnVector& b, 
 *                                 const bool rowFlag)
 * @see LinearConstraint(const Matrix& A, const NEWMAT::ColumnVector& lower,
 *                            const NEWMAT::ColumnVector& upper)
 */
  LinearConstraint(const NEWMAT::Matrix& A, const NEWMAT::ColumnVector& b); 
/**
 * @param A a real Matrix
 * @param b a ColumnVector
 * @param rowFlag a bool
 * @see LinearConstraint(const Matrix& A);
 * @see LinearConstraint(const Matrix& A, const NEWMAT::ColumnVector& b)
 * @see LinearConstraint(const Matrix& A, const NEWMAT::ColumnVector& lower,
 *                            const NEWMAT::ColumnVector& upper)
 */
  LinearConstraint(const NEWMAT::Matrix& A, const NEWMAT::ColumnVector& b, 
                   const bool rowFlag);
/**
 * @param A a real Matrix
 * @param lower a ColumnVector
 * @param upper a ColumnVector
 * @see LinearConstraint(const Matrix& A);
 * @see LinearConstraint(const Matrix& A, const NEWMAT::ColumnVector& b)
 * @see LinearConstraint(const Matrix& A, const NEWMAT::ColumnVector& b, 
 *                                 const bool rowFlag)
 */
  LinearConstraint(const NEWMAT::Matrix& A, const NEWMAT::ColumnVector& lower,
                   const NEWMAT::ColumnVector& upper);

/**
 *  Destructor
 */
  virtual ~LinearConstraint(){}


/**
 * @return The number of constraints.
 */
  virtual int getNumOfCons() const {return numOfCons_;}

/**
 * @return The number of variables.
 */
  virtual int getNumOfVars() const {return numOfVars_;}

/**
 * @return The lower bounds on the variables.
 */
  virtual NEWMAT::ColumnVector getLower() const {return lower_;}

/**
 * @return The upper bounds on the variables.
 */
  virtual NEWMAT::ColumnVector getUpper() const {return upper_;}

/**
 * @return Value of the linear constraints.
 */
  virtual NEWMAT::ColumnVector getConstraintValue() const {return cvalue_;}

/**
 * @return Constraint violation.
 */
  virtual NEWMAT::ColumnVector getConstraintViolation() const {return cviolation_;}

/**
 * @return Indices of constraints with finite bounds 
 */
  OptppArray<int> getConstraintMappingIndices() const 
  		{ return constraintMappingIndices_; }

  
/**
 * Assigns a value to the constraint matrix.
 */
  void setA(NEWMAT::Matrix & A);

/**
 * Pure Virtual Functions
 */

/**
 * @return Type of constraint - Leqn
 */
  virtual NEWMAT::ColumnVector getConstraintType() const = 0;

/**
 * Takes one argument and returns a ColumnVector.
 * @param xc a ColumnVector
 * @return Matrix-vector product of A and xc.
 */
  virtual NEWMAT::ColumnVector evalAx(const NEWMAT::ColumnVector& xc) const = 0;

/**
 * Takes one argument and returns a ColumnVector.
 * @param xc a ColumnVector
 * @return The residual of the linear equations
 * evaluated at xc.
 */
  virtual NEWMAT::ColumnVector evalResidual(const NEWMAT::ColumnVector& xc) const = 0;
  virtual void evalCFGH(const NEWMAT::ColumnVector& xc) const = 0;

/**
 * Takes one argument and returns a real Matrix.
 * @param xc a ColumnVector
 * @return The gradient of the linear equations
 * evaluated at xc.
 */
  virtual NEWMAT::Matrix evalGradient(const NEWMAT::ColumnVector& xc) const = 0;

/**
 * Takes two arguments and returns a bool.
 * @param xc a ColumnVector
 * @param epsilon a real argument
 * @return A bool
 */
  virtual bool amIFeasible(const NEWMAT::ColumnVector& xc, double epsilon) const = 0;

private:
/**
 * Takes one argument and returns a SymmetricMatrix
 * @param xc a ColumnVector
 * @return The constraint Hessain evaluated at xc 
 */
  virtual NEWMAT::SymmetricMatrix evalHessian(NEWMAT::ColumnVector& xc) const ;
/**
 * Takes two arguments and returns an array of real SymmetricMatrices.
 * @param xc a ColumnVector
 * @param darg an integer argument
 * @return An array of constraint Hessians.
 */
  virtual OptppArray<NEWMAT::SymmetricMatrix> evalHessian(NEWMAT::ColumnVector& xc, int darg) const;

/**
 * Takes one arguments and returns a bool.
 * @param A a Matrix
 * @return A bool
 */
  bool dimMatch(NEWMAT::Matrix& A);

};

} // namespace OPTPP
#endif
