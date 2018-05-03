#ifndef LinearEquation_h
#define LinearEquation_h


#include "LinearConstraint.h"

/**
 * LinearEquation is a derived class of LinearConstraint.
 *
 * @author P.J. Williams, Sandia National Laboratories, pwillia@sandia.gov
 * @date  modified 02/2006
 *
 */

//----------------------------------------------------------------------
// Linear Equations
//----------------------------------------------------------------------

namespace OPTPP {

class LinearEquation: public LinearConstraint {

protected:
  /// Right-hand side of equation
  NEWMAT::ColumnVector b_;		
  /// Type of constraint - Leqn
  NEWMAT::ColumnVector ctype_;		

public:

/**
 * Default Constructor
 * @see LinearEquation(const NEWMAT::Matrix& A, const NEWMAT::ColumnVector& rhs);
 */
  LinearEquation();
/**
 * @param A a real NEWMAT::Matrix
 * @param rhs NEWMAT::ColumnVector
 * @see LinearEquation()
 */
  LinearEquation(const NEWMAT::Matrix& A, const NEWMAT::ColumnVector& rhs);

/**
 * Destructor
 */
  virtual ~LinearEquation(){}

/**
 * @return Type of constraint - Leqn
 */
  virtual NEWMAT::ColumnVector getConstraintType() const { return ctype_;};

/**
 * @return The right-hand side of the equation.
 */
  NEWMAT::ColumnVector getB() const {return b_;}

/**
 * Takes one argument and returns a ColumnVector.
 * @param xc a ColumnVector
 * @return Matrix-vector product of A and xc.
 */
  virtual NEWMAT::ColumnVector evalAx(const NEWMAT::ColumnVector& xc) const;

  /**
   * Takes one argument and returns a ColumnVector.
   * @param xc a ColumnVector
   * @return The residual of the linear equations
   * evaluated at xc.
   */
  virtual NEWMAT::ColumnVector evalResidual(const NEWMAT::ColumnVector& xc) const;
  virtual void evalCFGH(const NEWMAT::ColumnVector& xc) const;

  /**
   * Takes one argument and returns a real Matrix.
   * @param xc a ColumnVector
   * @return The gradient of the linear equations
   * evaluated at xc.
   */
  virtual NEWMAT::Matrix evalGradient(const NEWMAT::ColumnVector& xc) const;

  /**
   * Takes two arguments and returns a bool.
   * @param xc a ColumnVector
   * @param epsilon a real argument.
   * @return The feasibility of the linear equations at xc.
   */
  virtual bool amIFeasible(const NEWMAT::ColumnVector& xc, double epsilon) const;
};

} // namespace OPTPP
#endif
