#include "NonLinearEquation.h"

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
NonLinearEquation::NonLinearEquation():
   NonLinearConstraint(), b_(0), ctype_(0){} 

NonLinearEquation::NonLinearEquation(NLP* nlprob, int numconstraints ):
   NonLinearConstraint(nlprob,numconstraints), b_(numconstraints), 
   ctype_(numconstraints) 
   {b_ = 0.0; ctype_.ReSize(numOfCons_); ctype_ = NLeqn; } 

NonLinearEquation::NonLinearEquation(NLP* nlprob, const ColumnVector& b,
                                 int numconstraints ):
   NonLinearConstraint(nlprob, b, numconstraints), b_(b),
   ctype_(nlprob->getDim()) 
   {ctype_.ReSize(numOfCons_); ctype_ = NLeqn;} 

// Functions for computing various quantities 
void NonLinearEquation::evalCFGH(const ColumnVector & xc) const
{
  nlp_->evalC(xc);
}

ColumnVector NonLinearEquation::evalResidual(const ColumnVector& xc) const 
{ 
      int i, index;
      ColumnVector resid( numOfCons_);
      cvalue_ = nlp_->evalCF(xc);
     
      // 08/06/2001 PJW - Okay to only loop over nnzl.
      // We assume nnzl = total constraints and nnzu = 0
      // for the equality constrained case.
      for( i = 1; i <= nnzl_; i++){
          index = constraintMappingIndices_[i-1];
	  resid(i) = cvalue_(index) - b_(index); 
      }
      return resid;
}

Matrix NonLinearEquation::evalGradient(const ColumnVector& xc) const 
{ 
      int j, index;
      Matrix grad(numOfVars_, numOfCons_);
      Matrix constraint_grad = nlp_->evalCG(xc);
     
      for( j = 1; j <= nnzl_; j++){
          index = constraintMappingIndices_[j-1];
	  grad.Column(j) = constraint_grad.Column(index);
      }
      return grad;
}

SymmetricMatrix NonLinearEquation::evalHessian(ColumnVector& xc) const 
{ 
     SymmetricMatrix hess, constraint_hess;
     constraint_hess = nlp_->evalCH(xc);
     
     hess = constraint_hess;
     return hess;
}

OptppArray<SymmetricMatrix> NonLinearEquation::evalHessian(ColumnVector& xc,
							  int darg) const 
{ 
     int i, index;
     OptppArray<SymmetricMatrix> hess(numOfCons_);
     OptppArray<SymmetricMatrix> constraint_hess = nlp_->evalCH(xc,darg);
     
     for( i = 1; i <= nnzl_; i++){
         index = constraintMappingIndices_[i-1];
	 hess[i-1] = constraint_hess[index-1];
     }
     return hess;
}

bool NonLinearEquation::amIFeasible(const ColumnVector& xc, double epsilon) const
{
     int i;
     bool feasible = true;
     ColumnVector residual = evalResidual(xc);
     for(i = 1; i <= numOfCons_; i++){
        if( (residual(i) < -epsilon) || (residual(i) > epsilon) ){
           feasible = false;
           break;
        }
     }
     return feasible;
}

} // namespace OPTPP
