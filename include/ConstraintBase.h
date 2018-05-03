#ifndef ConstraintBase_h
#define ConstraintBase_h

/*----------------------------------------------------------------------
 Copyright (c) 2001, Sandia Corporation.   Under the terms of Contract 
 DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 work by or on behalf of the U.S. Government.
 ----------------------------------------------------------------------*/

#include <iostream>

#include "globals.h"

#include "include.h"
#include "newmatap.h"
#include "precisio.h"

#include "BoolVector.h"
#include "OptppArray.h"
#include "OptppExceptions.h"
#include "OptppFatalError.h"

double const MIN_BND = -FLT_MAX;
double const MAX_BND =  FLT_MAX;
double const BIG_BND =  1.0e10;

/**
 * ConstraintBase is an abstract class.
 * All other constraint classes are derived from this one
 *
 * @author P.J. Williams, Sandia National Laboratories, pwillia@sandia.gov
 * @date Last modified 02/2006
 *
 */

namespace OPTPP {

class ConstraintBase {

public:
 /**
  * Destructor
  */
  virtual ~ConstraintBase() {}

  // Virtual Functions 
  
  /**
   * @return Number of constraints.
   */
   virtual int getNumOfCons() const = 0;

  /**
   * @return Number of variables.
   */
   virtual int getNumOfVars() const = 0;

  /**
   * @return Lower bounds on the constraints
   */
   virtual NEWMAT::ColumnVector getLower() const = 0;

  /**
   * @return Upper bounds on the constraints
   */
   virtual NEWMAT::ColumnVector getUpper() const = 0;

  /**
   * @return Type of constraints.
   */
   virtual NEWMAT::ColumnVector getConstraintType() const = 0;

  /**
   * @return Value of constraints.
   */
   virtual NEWMAT::ColumnVector getConstraintValue() const = 0;

  /**
   * @return The constraint violation.
   */
   virtual NEWMAT::ColumnVector getConstraintViolation() const = 0;

  /**
   * @return Indices of constraints with finite bounds.
   */
   virtual OptppArray<int> getConstraintMappingIndices() const = 0;

  /**
    * Takes one argument and returns a ColumnVector.
    * @param xcurrent a ColumnVector
    * @return The gradient of the constraints evaluated at xcurrent.
    */
   virtual NEWMAT::ColumnVector evalResidual(const NEWMAT::ColumnVector& xcurrent) const = 0;
   virtual void evalCFGH(const NEWMAT::ColumnVector& xcurrent) const = 0;

  /**
    * Takes one argument and returns a Matrix.
    * @param xcurrent a ColumnVector
    * @return The gradient of the constraints evaluated at xcurrent.
    */
   virtual NEWMAT::Matrix evalGradient(const NEWMAT::ColumnVector& xcurrent) const = 0;

   /**
    * Takes one arguments and returns a SymmetricMatrix
    * @param xcurrent a ColumnVector
    * @return The constraint Hessian evaluated at xcurrent.
    */
   virtual NEWMAT::SymmetricMatrix evalHessian(NEWMAT::ColumnVector& xcurrent) const = 0;

   /**
    * Takes two arguments and returns an array of real SymmetricMatrices.
    * @param xcurrent a ColumnVector
    * @param darg an integer argument
    * @return An array of constraint Hessians.
    */
   virtual OptppArray<NEWMAT::SymmetricMatrix> evalHessian(NEWMAT::ColumnVector& xcurrent, int darg) const = 0;

  /**
   * Takes two arguments and returns a bool.
   * @param xcurrent a ColumnVector
   * @param epsilon a real argument.
   * @return The feasibility of the nonlinear equations at xcurrent.
   */
   virtual bool amIFeasible(const NEWMAT::ColumnVector& xcurrent, double epsilon) const = 0;
};

} // namespace OPTPP
#endif
