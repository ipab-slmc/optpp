
#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#include "NonLinearConstraint.h"

using NEWMAT::ColumnVector;
using NEWMAT::Matrix;
using NEWMAT::SymmetricMatrix;

//------------------------------------------------------------------------
// P.J. Williams
// Sandia National Laboratories
// pwillia@sandia.gov
// Last modified 10/19/1999
//------------------------------------------------------------------------

namespace OPTPP {

// Constructors
NonLinearConstraint::NonLinearConstraint():
  nlp_(0), lower_(0), upper_(0), cvalue_(0), cviolation_(0),
  numOfCons_(0), numOfVars_(0), nnzl_(0), nnzu_(0), 
  constraintMappingIndices_(0), stdForm_(true){;}

NonLinearConstraint::NonLinearConstraint(NLP* nlprob, int numconstraints):
   nlp_(nlprob), lower_(numconstraints), upper_(numconstraints), 
   cvalue_(numconstraints), cviolation_(numconstraints),
   numOfCons_(numconstraints), numOfVars_(nlprob->getDim()), 
   nnzl_(0), nnzu_(0), constraintMappingIndices_(0), stdForm_(true)
   { 
     cvalue_ = 1.0e30; cviolation_ = 0.0;
     lower_ = 0.0; upper_ = MAX_BND;
     nnzl_ = numconstraints;
     for(int i = 1; i <= numconstraints; i++)
        constraintMappingIndices_.append(i);
   }

NonLinearConstraint::NonLinearConstraint(NLP* nlprob, const bool conFlag, int numconstraints):
   nlp_(nlprob), lower_(numconstraints), upper_(numconstraints),
   cvalue_(numconstraints), cviolation_(numconstraints),
   numOfCons_(numconstraints), numOfVars_(nlprob->getDim()), 
   nnzl_(0), nnzu_(0), constraintMappingIndices_(0), stdForm_(conFlag)
   {
      cvalue_ = 1.0e30; cviolation_ = 0.0;
      if (stdForm_) { 
         lower_ = 0.0; upper_ = MAX_BND; 
         nnzl_ = numconstraints;
         for(int i = 1; i <= numconstraints; i++)
            constraintMappingIndices_.append(i);
      }
      else{ 
	 lower_ = MIN_BND; upper_ = 0.0; 
         nnzu_ = numconstraints;
         for(int i = 1; i <= numconstraints; i++)
            constraintMappingIndices_.append(i);
      }
   }

NonLinearConstraint::NonLinearConstraint(NLP* nlprob, const ColumnVector& rhs, 
                                                 int numconstraints):
   nlp_(nlprob), lower_(rhs), upper_(rhs), 
   cvalue_(numconstraints), cviolation_(numconstraints), 
   numOfCons_(numconstraints), numOfVars_(nlprob->getDim()), 
   nnzl_(0), nnzu_(0), constraintMappingIndices_(0), stdForm_(true) 
   {
      cvalue_ = 1.0e30; cviolation_ = 0.0;
      for(int i = 1; i <= numconstraints; i++){
         if(lower_(i) > -BIG_BND && upper_(i) < BIG_BND ){
            constraintMappingIndices_.append(i);
	    nnzl_++;
         }
      }
      numOfCons_ = nnzl_;
   }


NonLinearConstraint::NonLinearConstraint(NLP* nlprob, const ColumnVector& rhs,
                                           const bool conFlag, int numconstraints):
   nlp_(nlprob), lower_(numconstraints), upper_(numconstraints), 
   cvalue_(numconstraints), cviolation_(numconstraints),
   numOfCons_(numconstraints), numOfVars_(nlprob->getDim()), 
   nnzl_(0), nnzu_(0), constraintMappingIndices_(0), stdForm_(conFlag)
   {
      cvalue_ = 1.0e30; cviolation_ = 0.0;
      if (stdForm_) { 
	 lower_ = rhs;
         upper_ = MAX_BND; 
         for(int i = 1; i <= numconstraints; i++){
	    if( lower_(i) > -BIG_BND ){
               constraintMappingIndices_.append(i);
               nnzl_++;
            }
         }
      }
      else{ 
	 lower_ = MIN_BND; 
	 upper_ = rhs; 
         for(int i = 1; i <= numconstraints; i++){
	    if( upper_(i) < BIG_BND ){
              constraintMappingIndices_.append(i);
              nnzu_++;
            }
         }
      }
     numOfCons_ = nnzl_ + nnzu_;
   }


NonLinearConstraint::NonLinearConstraint(NLP* nlprob, const ColumnVector& lower,
                                       const ColumnVector& upper, int numconstraints):
   nlp_(nlprob), lower_(lower), upper_(upper), 
   cvalue_(numconstraints), cviolation_(numconstraints),
   numOfCons_(2*numconstraints), numOfVars_(nlprob->getDim()), 
   nnzl_(0), nnzu_(0), constraintMappingIndices_(0), stdForm_(true)
   { 
      cvalue_ = 1.0e30; cviolation_ = 0.0;
      for(int i = 1; i <= numconstraints; i++){
	 if( lower_(i) > -BIG_BND ){
           constraintMappingIndices_.append(i);
           nnzl_++;
         }
      }
      for(int i = 1; i <= numconstraints; i++){
	  if( upper_(i) < BIG_BND ){
            constraintMappingIndices_.append(i);
            nnzu_++;
          }
      }
      numOfCons_ = nnzl_ + nnzu_;
   }

#ifdef DAKOTA_OPTPP
NonLinearConstraint::NonLinearConstraint(NLP* nlprob, const ColumnVector& lower,
                                       const ColumnVector& upper, 
				       int ne, int ni):
   nlp_(nlprob), lower_(lower), upper_(upper), cvalue_(ne+ni), cviolation_(ne+ni),
   numOfCons_(ne+ni), numOfVars_(nlprob->getDim()) , 
   nnzl_(0), nnzu_(0), constraintMappingIndices_(0), 
   stdForm_(true), ctype_(ne +ni) 
   { 
      OptppArray<int> temp;
      cvalue_ = 1.0e30; cviolation_ = 0.0;
      for(int i = 1; i <= ne; i++){
	 if( lower_(i) > -BIG_BND ){
           constraintMappingIndices_.append(i);
	   temp.append(NLeqn);
           nnzl_++;
         }
      }
      for(int i = ne+1; i <= ne + ni; i++){
	 if( lower_(i) > -BIG_BND ){
           constraintMappingIndices_.append(i);
	   temp.append(NLineq);
           nnzl_++;
         }
      }
      for(int i = ne+1; i <= ne + ni; i++){
	 if( upper_(i) < BIG_BND ){
           constraintMappingIndices_.append(i);
	   temp.append(NLineq);
           nnzu_++;
         }
      }
      numOfCons_ = nnzl_ + nnzu_;
      ctype_.ReSize(numOfCons_);
      for(int i = 1; i <= numOfCons_; i++)
         ctype_(i) = temp[i-1];
   }

void NonLinearConstraint::evalCFGH(const ColumnVector & xc) const
{
  nlp_->evalC(xc);
}

ColumnVector NonLinearConstraint::evalResidual(const ColumnVector & xc) const 
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

Matrix NonLinearConstraint::evalGradient(const ColumnVector & xc) const 
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

SymmetricMatrix NonLinearConstraint::evalHessian(ColumnVector & xc) const 
{
   // 09/05/01 PJW Dummy routine 
      SymmetricMatrix hess, constraint_hess, nconstraint_hess;
      constraint_hess = nlp_->evalCH(xc);
     
      nconstraint_hess = -constraint_hess;
      hess = constraint_hess & nconstraint_hess;
      return hess;
}

OptppArray<SymmetricMatrix> NonLinearConstraint::evalHessian(ColumnVector& xc,
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


bool NonLinearConstraint::amIFeasible(const ColumnVector & xc, double epsilon) const
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
#endif // DAKOTA_OPTPP 

} // namespace OPTPP
