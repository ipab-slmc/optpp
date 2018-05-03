//------------------------------------------------------------------------
// P.J. Williams
// Sandia National Laboratories
// pwillia@sandia.gov
// Last modified 03/20/2003
//------------------------------------------------------------------------

#include "BoundConstraint.h"

using NEWMAT::ColumnVector;
using NEWMAT::Matrix;
using NEWMAT::SymmetricMatrix;

namespace OPTPP {

// Constructors
BoundConstraint::BoundConstraint(): numOfCons_(0), numOfVars_(0), 
   nnzl_(0), nnzu_(0), lower_(0), upper_(0), cvalue_(0), 
   fixedVar_(0), freeVar_(0), stdForm_(0), 
   ctype_(0), constraintMappingIndices_(0){}

BoundConstraint::BoundConstraint(int nc, const ColumnVector& lower):
   numOfCons_(0), numOfVars_(nc), nnzl_(0), nnzu_(0), 
   lower_(nc), upper_(nc), cvalue_(nc), fixedVar_(nc,false), 
   freeVar_(nc,true), stdForm_(nc,true), ctype_(nc), constraintMappingIndices_(0)
   { 
     cvalue_ = 1.0e30;
     lower_  = lower; 
     upper_  = MAX_BND;     
     for(int i = 1; i <= nc; i++){
        if (lower_(i) > -BIG_BND){
	    constraintMappingIndices_.append(i);
	    nnzl_++;
        }
     }
     numOfCons_ = nnzl_;
     ctype_.ReSize(nnzl_);
     ctype_ = Bound;
   }

BoundConstraint::BoundConstraint(int nc, const  ColumnVector& bound, 
                                 const BoolVector& bdFlag):
   numOfCons_(0), numOfVars_(nc), nnzl_(0), nnzu_(0), 
   lower_(nc), upper_(nc), cvalue_(nc), fixedVar_(nc,false), 
   freeVar_(nc,true), stdForm_(nc,bdFlag), ctype_(nc), constraintMappingIndices_(0)
   { 
     cvalue_ = 1.0e30;
     for(int i = 1; i <= nc; i++){
        if ( stdForm_(i) ){
            lower_(i)  = bound(i); 
	    upper_(i)  = MAX_BND; 
            if (lower_(i) > -BIG_BND){
	        nnzl_++;
	        constraintMappingIndices_.append(i);
            }
	}
     }
     for(int i = 1; i <= nc; i++){
        if ( !stdForm_(i) ){
	    lower_(i)  = MIN_BND; 
            upper_(i)  = bound(i); 
            if (upper_(i) < BIG_BND){
	        nnzu_++;
	        constraintMappingIndices_.append(i);
            }
	}
     }
     numOfCons_ = nnzl_ + nnzu_;
     ctype_.ReSize(numOfCons_);
     ctype_ = Bound;
   }

BoundConstraint::BoundConstraint(int nc, const ColumnVector& lower,
                                 const ColumnVector& upper):
   numOfCons_(0), numOfVars_(nc), nnzl_(0), nnzu_(0), 
   lower_(nc), upper_(nc), cvalue_(nc),
   fixedVar_(nc,false), freeVar_(nc,true), stdForm_(nc,true), 
   ctype_(2*nc), constraintMappingIndices_(0)
   { 
     cvalue_ = 1.0e30;
     lower_  = lower; 
     for(int i = 1; i <= nc; i++){
        if (lower_(i) > -BIG_BND){
	    nnzl_++;
	    constraintMappingIndices_.append(i);
        }
     }
     upper_  = upper; 
     for(int i = 1; i <= nc; i++){
        if (upper_(i) <  BIG_BND){
	    nnzu_++;
	    constraintMappingIndices_.append(i);
        }
     }
     numOfCons_ = nnzl_ + nnzu_;
     ctype_.ReSize(numOfCons_);
     ctype_ = Bound;
     if(!amIConsistent() ) 
       OptppmathError("Error in Constructor - Lower bound exceeds upper bound");  
   }

void BoundConstraint::evalCFGH(const ColumnVector & xc) const
{
  return;
}

ColumnVector BoundConstraint::evalResidual(const ColumnVector& xc) const 
{ 
    int i, index, nnz = nnzl_ + nnzu_;
    ColumnVector resid(nnz);
    for(i = 1; i <= nnzl_; i++){
        index = constraintMappingIndices_[i-1];
	resid(i) = xc(index) - lower_(index);
    }
    for(i = nnzl_+1; i <= nnz; i++){
        index = constraintMappingIndices_[i-1];
	resid(i) = upper_(index) - xc(index);
    }
    cvalue_ = xc;
    return resid;
}

Matrix BoundConstraint::evalGradient(const ColumnVector& xc) const 
{ 
    int i, j, nnz = nnzl_+ nnzu_;
    Matrix D(numOfVars_, nnz);
    D = 0.0;
    
    for(j = 1; j <= nnzl_; j++){
      i = constraintMappingIndices_[j-1];
      D(i,j) = 1.0;
    }
    for(j = nnzl_+1; j <= nnz; j++){
      i = constraintMappingIndices_[j-1];
      D(i,j) = -1.0;
    }
    return D;
}

SymmetricMatrix BoundConstraint::evalHessian(ColumnVector& xc) const 
{ 
    SymmetricMatrix H(numOfCons_);
    H = 0;
    return H;
}

OptppArray<SymmetricMatrix> BoundConstraint::evalHessian(ColumnVector& xc, 
                                           int darg) const 
{ 
    OptppArray<SymmetricMatrix> Hessian(1);
    SymmetricMatrix H(numOfCons_);
    H = 0;
    Hessian[0] = H;
    return Hessian;
}

bool BoundConstraint::amIFeasible(const ColumnVector& xc, double epsilon) const 
{
    int i;
    bool feasible = true;

    for(i = 1; i <= numOfVars_; i++)
      if(xc(i) <  lower_(i) || xc(i) > upper_(i)){
        feasible = false;  
        break;
      }
    return feasible;
}

bool BoundConstraint::amIConsistent() const
{
    int i;
    bool consistent = true;
    for(i = 1; i <= numOfVars_; i++)
       if(lower_(i) > upper_(i)){
	  consistent = false;
	  break;
      }
    return consistent;
}

} // namespace OPTPP
