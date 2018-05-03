#include "LinearEquation.h"

using NEWMAT::ColumnVector;
using NEWMAT::Matrix;
using NEWMAT::SymmetricMatrix;

//------------------------------------------------------------------------
// P.J. Williams
// Sandia National Laboratories
// pwillia@sandia.gov
// Last modified 03/20/2003 
//------------------------------------------------------------------------

namespace OPTPP {

// Constructors 
LinearEquation::LinearEquation():
      LinearConstraint(), b_(0), ctype_(0){;}

LinearEquation::LinearEquation(const Matrix& A, const ColumnVector& rhs):
      LinearConstraint(A, rhs), b_(rhs), ctype_(A.Nrows())
      {ctype_.ReSize(numOfCons_); ctype_ = Leqn;}

// Functions for computing various quantities 
ColumnVector LinearEquation::evalAx(const ColumnVector& xc) const 
{ 
      int i, index;
      ColumnVector Ax(numOfCons_);
      Matrix tmp(numOfCons_, numOfVars_);
      for( i = 1; i <= numOfCons_; i++){
         index = constraintMappingIndices_[i-1];
	 tmp.Row(i) = A_.Row(index);
      }
      Ax = tmp*xc;
      return Ax;
}

void LinearEquation::evalCFGH(const ColumnVector & xc) const
{
  return;
}

ColumnVector LinearEquation::evalResidual(const ColumnVector& xc) const 
{ 
      int i, index;
      cvalue_         = A_*xc;
      ColumnVector Ax = evalAx(xc);
      ColumnVector residual(numOfCons_);

      for( i = 1; i <= numOfCons_; i++){
         index = constraintMappingIndices_[i-1];
         residual(i) = Ax(i) - b_(index); 
      }
      return residual;
}

Matrix LinearEquation::evalGradient(const ColumnVector& xc) const 
{ 
      int i, index;
      Matrix tmp(numOfCons_, numOfVars_);

      for( i = 1; i <= numOfCons_; i++){
         index = constraintMappingIndices_[i-1];
         tmp.Row(i) = A_.Row(index);
      }
      return tmp.t();
}


bool LinearEquation::amIFeasible(const ColumnVector & xc, double epsilon) const
{
     int i;
     bool feasible = true;
     ColumnVector residual = evalResidual(xc);
     for(i = 1; i <= numOfCons_; i++){
        if( (residual(i) > epsilon) || (residual(i) < -epsilon) ){
           feasible = false;
           break;
        }
     }
     return feasible;
}

} // namespace OPTPP
