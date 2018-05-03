#ifndef CompoundConstraint_h
#define CompoundConstraint_h

/*----------------------------------------------------------------------
 Copyright (c) 2001, Sandia Corporation.   Under the terms of Contract 
 DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 work by or on behalf of the U.S. Government.
 ----------------------------------------------------------------------*/


#include "Constraint.h"
#include "OptppArray.h"

/**
 * CompoundConstraint is an array of constraints. 
 *
 * @author P.J. Williams, Sandia National Laboratories, pwillia@sandia.gov
 * @date   Last modified 02/2006
 */

namespace OPTPP {

class CompoundConstraint: public ConstraintBase{

private:
  OptppArray<Constraint> constraints_; 	///< Array of constraints
  int numOfSets_;			///< Number of constraint sets
  NEWMAT::ColumnVector lower_;		///< Lower bound on the constraints 
  NEWMAT::ColumnVector upper_;		///< Upper bound on the constraints

public:
 /**
  * Default Constructor
  * @see CompoundConstraint(const Constraint& c1)
  * @see CompoundConstraint(const Constraint& c1, const Constraint& c2)
  * @see CompoundConstraint(const OptppArray<Constraint>& constraints)
  * @see CompoundConstraint(const CompoundConstraint& cc)
  */
  CompoundConstraint(); 
 /**
  * @param c1 a Constraint
  * @see CompoundConstraint(const Constraint& c1, const Constraint& c2)
  * @see CompoundConstraint(const OptppArray<Constraint>& constraints)
  * @see CompoundConstraint(const CompoundConstraint& cc)
  */
  CompoundConstraint(const Constraint& c1); 
 /**
  * @param c1 a Constraint
  * @param c2 a Constraint
  * @see CompoundConstraint(const Constraint& c1)
  * @see CompoundConstraint(const OptppArray<Constraint>& constraints)
  * @see CompoundConstraint(const CompoundConstraint& cc)
  */
  CompoundConstraint(const Constraint& c1, const Constraint& c2); 
 /**
  * @param constraints an array of Constraints 
  * @see CompoundConstraint(const Constraint& c1)
  * @see CompoundConstraint(const Constraint& c1, const Constraint& c2)
  * @see CompoundConstraint(const CompoundConstraint& cc)
  */
  CompoundConstraint(const OptppArray<Constraint>& constraints); 
 /**
  * Copy constructor
  * @param cc a CompoundConstraint
  * @see CompoundConstraint(const Constraint& c1)
  * @see CompoundConstraint(const Constraint& c1, const Constraint& c2)
  * @see CompoundConstraint(const OptppArray<Constraint>& constraints)
  */
  CompoundConstraint(const CompoundConstraint& cc);

 /**
  * Destructor
  */
  ~CompoundConstraint() {;}

  /// Assignment Operator
  CompoundConstraint& operator=(const CompoundConstraint& cc);

  /**
   * @return Integer 
   */
  int compare(const Constraint& c1, const Constraint& c2); 

  /**
   * @return Number of different constraint sets.
   */
  int getNumOfSets() const {return numOfSets_;}

  /**
   * @return Number of nonlinear constraints.
   */
  int getNumOfNLCons() const;

  /**
   * @return Number of total constraints.
   */
  virtual int getNumOfCons() const;

  /**
   * @return Number of variables.
   */
  virtual int getNumOfVars() const;

  /**
   * @return Lower bounds on the constraints 
   */
  virtual NEWMAT::ColumnVector getLower() const;

  /**
   * @return Upper bounds on the constraints 
   */
  virtual NEWMAT::ColumnVector getUpper() const;

  /**
   * @return Constraint type 
   */
  virtual NEWMAT::ColumnVector getConstraintType() const;

  /**
   * @return Value of the entire constraint set 
   */
  virtual NEWMAT::ColumnVector getConstraintValue() const;

  /**
   * @return Violation of the entire constraint set 
   */
  virtual NEWMAT::ColumnVector getConstraintViolation() const;

  /**
   * @return Value of the nonlinear constraints only 
   */
  NEWMAT::ColumnVector getNLConstraintValue() const;

  /**
   * @return Indices of the constraints with finite bounds 
   */
  virtual OptppArray<int> getConstraintMappingIndices() const;

  /**
   * @return Reference to constraint i
   */
  Constraint& operator[]( int i ) { return constraints_[i];}

  /**
   * @return Const reference to constraint i
   */
  const Constraint& operator[](int i) const {return constraints_[i];}

  /**
   * @return Feasible vector with respect to bounds 
   */
  void computeFeasibleBounds(NEWMAT::ColumnVector& xcurrent, double epsilon);

  /**
   * @return Feasible vector with respect to linear and nonlinear inequalities 
   */
  void computeFeasibleInequalities(NEWMAT::ColumnVector& xcurrent, double ftol);

  /**
   * @return Sorted constraints - equations followed by inequalities 
   */
  void insertSort();

  /**
   * @return Output constraint values to the screen 
   */
  void printConstraints();

  /**
   * Takes one argument and returns a ColumnVector
   * @param xcurrent a ColumnVector
   * @return The residual of the constraints.
   */
  virtual NEWMAT::ColumnVector evalResidual(const NEWMAT::ColumnVector& xcurrent ) const ;
  virtual void evalCFGH(const NEWMAT::ColumnVector& xcurrent ) const ;

  /**
   * Takes one argument and returns a real Matrix.
   * @param xcurrent a ColumnVector
   * @return The gradient of the constraints.
   */
  virtual NEWMAT::Matrix evalGradient(const NEWMAT::ColumnVector& xcurrent ) const ;

  /**
   * Takes two arguments and returns a real SymmetricMatrix
   * @param xcurrent a ColumnVector
   * @param LagMultiplier a ColumnVector
   * @return The Hessian of the constraints multiplied by its assoc. multiplier.
   */
  NEWMAT::SymmetricMatrix evalHessian(NEWMAT::ColumnVector& xcurrent, 
                            const NEWMAT::ColumnVector& LagMultiplier) const ;
  /**
   * Takes one arguments and returns a NEWMAT::SymmetricMatrix
   * @param xcurrent a ColumnVector
   * @return The Hessian of the constraints.
   */
  virtual NEWMAT::SymmetricMatrix evalHessian(NEWMAT::ColumnVector& xcurrent ) const ;

  /**
   * Takes two arguments and returns an array of real SymmetricMatrices
   * @param xcurrent a ColumnVector
   * @param darg an integer argument
   * @return  An array of Hessians.
   */
  virtual OptppArray<NEWMAT::SymmetricMatrix> evalHessian(NEWMAT::ColumnVector& xcurrent, int darg) const ;

  /**
   * Takes two arguments and returns a bool.
   * @param xcurrent a ColumnVector
   * @param epsilon a real argument
   * @return true - constraints are feasible 
   *
   * @return false - constraints are infeasible
   */
  virtual bool amIFeasible(const NEWMAT::ColumnVector& xcurrent, double epsilon) const;

private:
  /**
   * Sorts an array of constraints.  
   * Equations are followed by inequalities.
   */
  void insertSort(const OptppArray<Constraint>& constraints);
};

} // namespace OPTPP
#endif
