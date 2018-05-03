//------------------------------------------------------------------------
// P.J. Williams
// Sandia National Laboratories
// pwillia@sandia.gov
// Last modified 03/20/2003 
//------------------------------------------------------------------------

#include "LinearInequality.h"

using NEWMAT::ColumnVector;
using NEWMAT::Matrix;
using NEWMAT::SymmetricMatrix;

namespace OPTPP {

// Constructors
LinearInequality::LinearInequality():
      LinearConstraint( ), ctype_(0) {}

LinearInequality::LinearInequality(const Matrix& A, const ColumnVector& rhs):
      LinearConstraint(A,rhs,true), ctype_(A.Nrows())
      {ctype_.ReSize(numOfCons_); ctype_ = Lineq;}

LinearInequality::LinearInequality(const Matrix& A, const ColumnVector& rhs, 
                                   const bool rowFlag):
      LinearConstraint(A,rhs,rowFlag), ctype_(A.Nrows())
      {ctype_.ReSize(numOfCons_); ctype_ = Lineq;}

LinearInequality::LinearInequality(const Matrix& A, const ColumnVector& lower, 
                                   const ColumnVector& upper):
      LinearConstraint(A,lower,upper), ctype_(2*A.Nrows())
      {ctype_.ReSize(numOfCons_); ctype_ = Lineq;}

// Evaluation Methods
ColumnVector LinearInequality::evalAx(const ColumnVector& xc) const 
{

      int i, index, nnz = nnzl_ + nnzu_;
      ColumnVector Ax(numOfCons_);
      Matrix tmp(numOfCons_, numOfVars_);
      for( i = 1; i <= nnzl_; i++){
         index = constraintMappingIndices_[i-1];
	 tmp.Row(i) = A_.Row(index);
      }
      for( i = nnzl_+1; i <= nnz; i++){
         index = constraintMappingIndices_[i-1];
	 tmp.Row(i) = -A_.Row(index);
      }
      Ax = tmp*xc;
      return Ax;
}

void LinearInequality::evalCFGH(const ColumnVector & xc) const
{
  return;
}

ColumnVector LinearInequality::evalResidual(const ColumnVector& xc) const 
{
      int i, index, nnz = nnzl_ + nnzu_;
      cvalue_               = A_*xc;
      ColumnVector residual = evalAx(xc);
      for( i = 1; i <= nnzl_; i++){
         index = constraintMappingIndices_[i-1];
         residual(i)-= lower_(index); 
      }
      for( i = nnzl_+1; i <= nnz; i++){
         index = constraintMappingIndices_[i-1];
         residual(i)+= upper_(index); 
      }
      return residual;
}

Matrix LinearInequality::evalGradient(const ColumnVector& xc) const 
{
      int i, index, nnz = nnzl_ + nnzu_;
      Matrix tmp(numOfCons_, numOfVars_); 
      for( i = 1; i <= nnzl_; i++){
         index = constraintMappingIndices_[i-1];
         tmp.Row(i) = A_.Row(index); 
      }
      for( i = nnzl_+1; i <= nnz; i++){
         index = constraintMappingIndices_[i-1];
         tmp.Row(i) = -A_.Row(index); 
      }
      return tmp.t();
}

bool LinearInequality::amIFeasible(const ColumnVector& xc, double epsilon) const
{
     bool feasible = true;
     int i, index;
     ColumnVector residual = evalResidual(xc);

     for(i = 1; i <= numOfCons_; i++){
       index = constraintMappingIndices_[i-1];
       if( residual(i) < -epsilon ){
          cviolation_(index) = residual(i);
          feasible = false;
       }
     }
     return feasible;
}

} // namespace OPTPP
