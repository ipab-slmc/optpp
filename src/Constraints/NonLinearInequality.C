#include "NonLinearInequality.h"

using NEWMAT::ColumnVector;
using NEWMAT::Matrix;
using NEWMAT::SymmetricMatrix;

//------------------------------------------------------------------------
// P.J. Williams
// Sandia National Laboratories
// pwillia@sandia.gov
// Last modified 03/20/2003
// 
// Standard form g(x) >= 0
//------------------------------------------------------------------------

namespace OPTPP {

// Constructors
// PJW: The +/- infinity bound checked is done in the NonLinearConstraint
// constructor.  If any bound is infinite, then the vector ctype_
// must be resized accordingly.
NonLinearInequality::NonLinearInequality():
     NonLinearConstraint(), ctype_(0), oneSided_(true){}

NonLinearInequality::NonLinearInequality(NLP* nlprob, int numconstraints):
     NonLinearConstraint(nlprob,true,numconstraints), 
     ctype_( numconstraints ), oneSided_(true)
     {ctype_.ReSize(numOfCons_); ctype_ = NLineq;}

NonLinearInequality::NonLinearInequality(NLP* nlprob, const bool flag, 
                         int numconstraints):
     NonLinearConstraint(nlprob,flag,numconstraints), 
     ctype_( numconstraints ), oneSided_(true)
     {ctype_.ReSize(numOfCons_); ctype_ = NLineq;}

NonLinearInequality::NonLinearInequality(NLP* nlprob, const ColumnVector& b, 
                         int numconstraints):
     NonLinearConstraint(nlprob,b,true,numconstraints), 
     ctype_( numconstraints ), oneSided_(true)
     {ctype_.ReSize(numOfCons_); ctype_ = NLineq;}

NonLinearInequality::NonLinearInequality(NLP* nlprob, const ColumnVector& b, 
                                     const bool flag, int numconstraints):
     NonLinearConstraint(nlprob,b,flag,numconstraints), 
     ctype_( numconstraints ), oneSided_(true)
     {ctype_.ReSize(numOfCons_); ctype_ = NLineq;}

NonLinearInequality::NonLinearInequality(NLP* nlprob, 
                                     const ColumnVector& lower, 
                                     const ColumnVector& upper,
                                     int numconstraints):
     NonLinearConstraint(nlprob, lower, upper,numconstraints), 
     ctype_(2*numconstraints ), oneSided_(false)
     {ctype_.ReSize(numOfCons_); ctype_ = NLineq;}


// Methods
void NonLinearInequality::evalCFGH(const ColumnVector & xc) const
{
  nlp_->evalC(xc);
}

ColumnVector NonLinearInequality::evalResidual(const ColumnVector & xc) const 
{
      int i, index;
      ColumnVector resid( numOfCons_);
      cvalue_ = nlp_->evalCF(xc);
     
      for( i = 1; i <= nnzl_; i++){
          index = constraintMappingIndices_[i-1];
	  resid(i) = cvalue_(index) - lower_(index); 
      }
      for( i = nnzl_+1; i <= numOfCons_; i++){
          index = constraintMappingIndices_[i-1];
	  resid(i) = upper_(index) - cvalue_(index); 
      }
      return resid;
}

Matrix NonLinearInequality::evalGradient(const ColumnVector & xc) const 
{
      int j, index;
      Matrix grad(numOfVars_, numOfCons_);
      Matrix constraint_grad = nlp_->evalCG(xc);
     
      for( j = 1; j <= nnzl_; j++){
          index = constraintMappingIndices_[j-1];
	  grad.Column(j) = constraint_grad.Column(index);
      }
      for( j = nnzl_+1; j <= numOfCons_; j++){
          index = constraintMappingIndices_[j-1];
	  grad.Column(j) = -constraint_grad.Column(index);
      }
      return grad;
}

SymmetricMatrix NonLinearInequality::evalHessian(ColumnVector & xc) const 
{
      SymmetricMatrix hess, constraint_hess, nconstraint_hess;
      constraint_hess = nlp_->evalCH(xc);
     
      if(oneSided_){
        if(stdForm_) 
          return constraint_hess;
        else
          return -constraint_hess;
      }
      else{
         nconstraint_hess = -constraint_hess;
	 hess = constraint_hess & nconstraint_hess;
	 return hess;
      }
}

OptppArray<SymmetricMatrix> NonLinearInequality::evalHessian(ColumnVector& xc,
                                                        int darg)const 
{
      int i, index;
      OptppArray<SymmetricMatrix> hess(numOfCons_);
      OptppArray<SymmetricMatrix> constraint_hess = nlp_->evalCH(xc,darg);
       
      for( i = 1; i <= nnzl_; i++){
          index = constraintMappingIndices_[i-1];
	  hess[i-1] = constraint_hess[index-1];
      }
      for( i = nnzl_+1; i <= numOfCons_; i++){
          index = constraintMappingIndices_[i-1];
	  hess[i-1] = -constraint_hess[index-1];
      }
      return hess;
}


bool NonLinearInequality::amIFeasible(const ColumnVector & xc, double epsilon) const
{
      int i, index;
      bool feasible = true;
      ColumnVector residual = evalResidual(xc);
      for(i = 1; i <= numOfCons_; i++){
         index = constraintMappingIndices_[i-1];
         if(residual(i) < -epsilon ){
            cviolation_(index) = residual(i);
            feasible = false;
         }
      }
      return feasible;
}

} // namespace OPTPP
