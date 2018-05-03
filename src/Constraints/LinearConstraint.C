//------------------------------------------------------------------------
// P.J. Williams
// Sandia National Laboratories
// pwillia@sandia.gov
// Last modified 03/20/2003 
//------------------------------------------------------------------------

#include "LinearConstraint.h"

using NEWMAT::ColumnVector;
using NEWMAT::Matrix;
using NEWMAT::SymmetricMatrix;

namespace OPTPP {

// Constructors
LinearConstraint::LinearConstraint():
    numOfCons_(0), numOfVars_(0), nnzl_(0), nnzu_(0),
    A_(0,0), Ax_(0), lower_(0), upper_(0), cvalue_(0), 
    cviolation_(0), constraintMappingIndices_(0), stdForm_(true) {;}

LinearConstraint::LinearConstraint(const Matrix& A):
    numOfCons_( A.Nrows() ), numOfVars_( A.Ncols() ), nnzl_(0), nnzu_(0),
    A_(A), Ax_( A.Nrows() ), lower_( A.Nrows() ), upper_( A.Nrows() ),
    cvalue_( A.Nrows() ), cviolation_(0), 
    constraintMappingIndices_(0), stdForm_(true)
    { 
      cvalue_ = 1.0e30; cviolation_ = 0.0;
      lower_ = 0.0; upper_ = MAX_BND;
      for(int i = 1; i <= numOfCons_; i++){
 	 constraintMappingIndices_.append(i);
         nnzl_++;
      }
      numOfCons_ = nnzl_;
    }

LinearConstraint::LinearConstraint(const Matrix& A, const ColumnVector& b):
    numOfCons_( A.Nrows() ), numOfVars_( A.Ncols() ), nnzl_(0), nnzu_(0),
    A_(A), Ax_( A.Nrows() ), lower_( b ), upper_( b ), 
    cvalue_( A.Nrows() ), cviolation_( A.Nrows() ),
    constraintMappingIndices_(0), stdForm_(true)
    {
       cvalue_ = 1.0e30; cviolation_ = 0.0;
       for(int i = 1; i <= numOfCons_; i++){
          if(lower_(i) > -BIG_BND ){
             constraintMappingIndices_.append(i);
             nnzl_++;
         }
       }
       numOfCons_ = nnzl_;
     }

LinearConstraint::LinearConstraint(const Matrix& A, const ColumnVector& b,
                                   const bool rowFlag):
    numOfCons_( A.Nrows() ), numOfVars_( A.Ncols() ), nnzl_(0), nnzu_(0),
    A_(A), Ax_( A.Nrows() ), lower_( A.Nrows() ), upper_( A.Nrows() ), 
    cvalue_( A.Nrows()), cviolation_( A.Nrows() ),
    constraintMappingIndices_(0), stdForm_(rowFlag)
    {
      int i;

      cvalue_  = 1.0e30; cviolation_ = 0.0;
      if( stdForm_ ){
        lower_ = b;
        upper_ = MAX_BND;
        for(i = 1; i <= numOfCons_; i++){
            if (lower_(i) > -BIG_BND){
	        constraintMappingIndices_.append(i);
	        nnzl_++;
            }
	}
      }
      else{
        upper_ = b;
	lower_ = MIN_BND; 
        for(i = 1; i <= numOfCons_; i++){
            if (upper_(i) < BIG_BND){
	        constraintMappingIndices_.append(i);
	        nnzu_++;
            }
	}
      }

      numOfCons_ = nnzl_ + nnzu_;
   }


LinearConstraint::LinearConstraint(const Matrix& A, const ColumnVector& lower,
                                   const ColumnVector& upper):
    numOfCons_( 2*A.Nrows() ), numOfVars_( A.Ncols() ), nnzl_(0), nnzu_(0),
    A_(A), Ax_( A.Nrows() ), lower_( lower ), upper_( upper ),
    cvalue_( A.Nrows()), cviolation_( A.Nrows()) ,
    constraintMappingIndices_(0), stdForm_(true)
    {
       int i, numconstraints = A.Nrows();

       cvalue_  = 1.0e30; cviolation_ = 0.0;
       for(i = 1; i <= numconstraints; i++){
           if(lower_(i) > -BIG_BND){
	      constraintMappingIndices_.append(i);
	      nnzl_++;
           }
	}
        for(i = 1; i <= numconstraints; i++){
            if(upper_(i) < BIG_BND){
	       constraintMappingIndices_.append(i);
	       nnzu_++;
            }
	}
        numOfCons_ = nnzl_ + nnzu_;
   }

void LinearConstraint::setA(Matrix& A)
{
    if( dimMatch(A) ) 
      A_ = A;
    else 
      OptppmathError("Check matrix dimensions.  Error in the setA method. ");
}

bool LinearConstraint::dimMatch(Matrix& A)
{
    bool match = true;
    if (numOfCons_ != A.Nrows() || numOfVars_ !=  A.Ncols() )
		   match = false;
    return match;
}

SymmetricMatrix LinearConstraint::evalHessian(ColumnVector& xc) const
{
    SymmetricMatrix H(numOfVars_);
    H = 0;
    return H;
}

OptppArray<SymmetricMatrix> LinearConstraint::evalHessian(ColumnVector& xc, int darg) const
{
    OptppArray<SymmetricMatrix>  H(1);
    SymmetricMatrix Htmp(numOfVars_);
    Htmp = 0;
    H[0] = Htmp;
    return H;
}

} // namespace OPTPP
