//------------------------------------------------------------------------
// P.J. Williams
// Sandia National Laboratories
// pwillia@sandia.gov
// Last modified 03/20/2003
//------------------------------------------------------------------------

#include "CompoundConstraint.h"
#include "ioformat.h"

using NEWMAT::ColumnVector;
using NEWMAT::Matrix;
using NEWMAT::SymmetricMatrix;
using std::cout;


namespace OPTPP {

// Constructors
CompoundConstraint::CompoundConstraint(): constraints_(0), 
   numOfSets_(0), lower_(0) , upper_(0){}

CompoundConstraint::CompoundConstraint(const Constraint& c1): 
   constraints_(0), numOfSets_(1)
{ 
   constraints_.append(c1); 
   lower_ = getLower(); 
   upper_ = getUpper(); 
}

CompoundConstraint::CompoundConstraint(const Constraint& c1, 
                                       const Constraint& c2): 
   constraints_(0), numOfSets_(2)
{ 
   constraints_.append(c1);
   constraints_.append(c2);
   insertSort();
   lower_       = getLower();
   upper_       = getUpper();
}

CompoundConstraint::CompoundConstraint(const OptppArray<Constraint>& constraints): 
   constraints_(constraints), numOfSets_(constraints.length())
{  
   insertSort();
   lower_  = getLower();  
   upper_   = getUpper();  
}

CompoundConstraint::CompoundConstraint(const CompoundConstraint& cc):
   constraints_(0), numOfSets_(cc.numOfSets_), lower_(cc.lower_), 
   upper_(cc.upper_) 
{ 
   if(numOfSets_ > 0){
     for(int i = 0; i < numOfSets_; i++)
        constraints_.append(cc[i]);
   }
}

CompoundConstraint& CompoundConstraint::operator=(const CompoundConstraint& cc)
{ 
   if(this != &cc){
      numOfSets_ = cc.numOfSets_;
      lower_     = cc.lower_;
      upper_     = cc.upper_;
      for(int i = 0; i < numOfSets_; i++)
         constraints_.append(cc[i]);
   }                                
   return *this;                       
}

// Accessor Methods 
int CompoundConstraint::getNumOfNLCons() const
{
   int Mi, i, k = 0;
   Constraint test;

   for(i = 0; i < numOfSets_; i++){
     test = constraints_[i];
     ColumnVector temp =  test.getConstraintType();
     if (temp(1) == NLeqn || temp(1) == NLineq){
        Mi   = test.getNumOfCons();
        k   += Mi; 
     }
   }
   return k;
}

int CompoundConstraint::getNumOfCons() const
{
   int Mi, i, k = 0;
   Constraint test;

   for(i = 0; i < numOfSets_; i++){
     test = constraints_[i];
     Mi   = test.getNumOfCons();
     k   += Mi; 
   }
   return k;
}

int CompoundConstraint::getNumOfVars() const
{
   int Mi, i, k = 0;
   Constraint test;

   for(i = 0; i < numOfSets_; i++){
     test = constraints_[i];
     Mi   = test.getNumOfVars();
     k   += Mi; 
   }
   if( k!= 0 && k == Mi*numOfSets_ ) 
     return Mi;
   else
     return 0;
}

ColumnVector CompoundConstraint::getLower() const
{
   ColumnVector result(getNumOfCons());
   Constraint test;

   for(int i = 0; i < numOfSets_; i++){
     test = constraints_[i];
     ColumnVector temp =  test.getLower();
     if (i == 0)
        result = temp;
     else 
        result &= temp;
   }
   return result;
}

ColumnVector CompoundConstraint::getUpper() const
{
   ColumnVector result(getNumOfCons());
   Constraint test;

   for(int i = 0; i < numOfSets_; i++){
     test = constraints_[i];
     ColumnVector temp =  test.getUpper();
     if (i == 0)
        result = temp;
     else 
        result &= temp;
   }
   return result;

}

ColumnVector CompoundConstraint::getConstraintType() const
{
   ColumnVector result(getNumOfCons());
   Constraint test;

   for(int i = 0; i < numOfSets_; i++){
     test = constraints_[i];
     ColumnVector temp =  test.getConstraintType();
     if (i == 0)
        result = temp;
     else 
        result &= temp;
   }
   return result;
}

OptppArray<int> CompoundConstraint::getConstraintMappingIndices() const
{
   OptppArray<int> result;
   Constraint test;

   for(int i = 0; i < numOfSets_; i++){
     test = constraints_[i];
     OptppArray<int> temp =  test.getConstraintMappingIndices();
     for(int j = 0; j < temp.length(); j++)
        result.append(temp[j]);
   }
   return result;

}

ColumnVector CompoundConstraint::getConstraintValue() const 
{
   /* 
    * Returns the raw value of all the constraints. 
    */

   Constraint test;
   ColumnVector temp, value;

   for(int i = 0; i < numOfSets_; i++){
     test   = constraints_[i];
     temp   = test.getConstraintValue();

     if (i == 0)
        value = temp;
     else 
        value &= temp;
   }
   return value;
}

ColumnVector CompoundConstraint::getNLConstraintValue() const 
{
   /* 
    * Returns the raw value of the nonlinear equations and inequalities.  
    */

   int j = 0;
   Constraint test;
   ColumnVector temp, type, value, zero(1);

   // If there are no nonlinear constraints, initialize a return value of zero
   zero = 0;

   for(int i = 0; i < numOfSets_; i++){
     test   = constraints_[i];
     type   = test.getConstraintType();

     if(type(1) == NLeqn || type(1) == NLineq){
        temp   = test.getConstraintValue();
        if(j == 0)
          value = temp;
        else 
          value &= temp;
        j++;
     }
   }

   if( j == 0) value = zero;

   return value;
}

ColumnVector CompoundConstraint::getConstraintViolation() const 
{
   /* CPJW - Rethink! 
    * Returns the infeasibility with respect to the constraints.  
    */

   Constraint test;
   ColumnVector temp, value;

   for(int i = 0; i < numOfSets_; i++){
     test   = constraints_[i];
     temp   = test.getConstraintViolation();

     if(i == 0)
       value = temp;
     else 
       value &= temp;
   }
   return value;
}

void CompoundConstraint::computeFeasibleBounds(ColumnVector&xc, double epsilon) 
{
   /* 
    * Returns a point that is feasible with respect to the bound constraints 
    */

   int i, j, nvars;
   Constraint test;
   ColumnVector type, lower, upper;

   for(i = 0; i < numOfSets_; i++){
     test   = constraints_[i];
     type   = test.getConstraintType();

     if(type(1) == Bound){
       nvars = test.getNumOfVars();
       lower = test.getLower();
       upper = test.getUpper();

       for(j = 1; j < nvars; j++){
         if( xc(j) < lower(j) || xc(j) > upper(j)){
           if(lower(j) > -BIG_BND && upper(j) == MAX_BND)
              xc(j) = lower(j) + epsilon;
           else if(upper(j) < BIG_BND && lower(j) == MIN_BND)
              xc(j) = upper(j) + epsilon;
           else
              xc(j) = (lower(j) + upper(j))/2.0 + epsilon;
         }
       }
     }
   }
}

void CompoundConstraint::computeFeasibleInequalities(ColumnVector&xc, double ftol) 
{
   /* 
    * Returns a point that is feasible with respect to general inequalities 
    */

   int i, j, ncons;
   double alpha = 0.5;
   Constraint test;
   Matrix grad_c;
   ColumnVector g, gTg, type, v;

   for(i = 0; i < numOfSets_; i++){
     test   = constraints_[i];
     type   = test.getConstraintType();

     if(type(1) == Lineq || type(1) == NLineq){
       if(!test.amIFeasible(xc,ftol)){
         v      = test.getConstraintViolation();
         grad_c = test.evalGradient(xc);
         if(type(1) == Lineq || type(1) == NLineq){
           ncons = v.Nrows();
           gTg.ReSize(ncons);
           OptppArray<int> indices = test.getConstraintMappingIndices();
           for(j = 1; j <ncons; j++ ){
             if( abs(v(j)) > alpha ) {
               g = grad_c.Column(indices[j-1]);
               gTg(j) = Dot(g,g);
               xc += (-v(j)/gTg(j))*g;
             }
           }
         } 
       }
     }
   }
}


void CompoundConstraint::printConstraints() 
{
   /* 
    * Prints the raw value of the constraints.  
    */

   int i, j, index, ncons, nvars;
   Constraint test;
   ColumnVector lower, upper, type, value;
   OptppArray<int> mapping;
   char s[2];

   for(i = 0; i < numOfSets_; i++){
     test   = constraints_[i];
     type   = test.getConstraintType();

     value   = test.getConstraintValue();
     lower   = test.getLower();
     upper   = test.getUpper();

     if(type(1) == Bound){
          cout <<"\nBound Constraints: \n";
     }
     else if(type(1) == NLeqn || type(1) == NLineq)
          cout <<"\nNonlinear Constraints: \n";
     else if(type(1) == Leqn || type(1) == Lineq)
          cout <<"\nLinear Constraints: \n";

     if(type(1) != Bound){
          ncons   = test.getNumOfCons();
          mapping = test.getConstraintMappingIndices();
          cout << "Index  Type       Lower   \t Constraint \t Upper \n";
          for (j = 1; j<=ncons; j++) {
               index = mapping[j-1];
               if(type(index) == NLeqn ||  type(index) == Leqn)   strcpy(s,"E");
               if(type(index) == NLineq || type(index) == Lineq)  strcpy(s,"I");
               cout << d(index,5) << "\t"   << s << "\t" 
                    << e(lower(index),12,4) << "\t" 
                    << e(value(index),12,4) <<  "\t" 
                    << e(upper(index),12,4) << "\n";
          }
     }
     else{
          nvars = getNumOfVars();
          cout << "Index \t Lower \t\t\t X \t Upper \n";
          for (j = 1; j<=nvars; j++) {
               cout << d(j,5) << "\t"   <<  e(lower(j),12,4) << "\t" 
                    << e(value(j),12,4) <<  "\t" << e(upper(j),12,4)  
                    << "\n";
          }
     }
   }
}

void CompoundConstraint::evalCFGH(const ColumnVector & xc) const
{
   Constraint test;
   ColumnVector result(numOfSets_);

   for(int i = 0; i < numOfSets_; i++){
     test = constraints_[i];
     test.evalCFGH(xc);
   }
}

// Evaluation 
ColumnVector CompoundConstraint::evalResidual(const ColumnVector& xc ) const 
{
   Constraint test;
   ColumnVector result(numOfSets_);

   for(int i = 0; i < numOfSets_; i++){
     test = constraints_[i];
     ColumnVector temp =  test.evalResidual(xc);
     if (i == 0)
        result = temp;
     else 
        result &= temp;
   }
   return result;
}

Matrix CompoundConstraint::evalGradient(const ColumnVector& xc ) const 
{
   Matrix grad;
   Constraint test;

   for(int i = 0; i < numOfSets_; i++){
     test = constraints_[i];
     Matrix temp = test.evalGradient(xc);
       if(i == 0)
          grad = temp;
       else
   	  grad |= temp;
   }
   return grad;
}

SymmetricMatrix CompoundConstraint::evalHessian(ColumnVector& xc ) const 
{
  // Extremely adhoc.  Conceived on 12/07/2000.  Vertical Concatenation
   SymmetricMatrix hessian(xc.Nrows());
   hessian = 0;
   return hessian;
}

OptppArray<SymmetricMatrix> CompoundConstraint::evalHessian(ColumnVector& xc, int darg ) const 
{
  // Extremely adhoc.  Conceived on 12/07/2000.  Vertical Concatenation
   SymmetricMatrix hessianT(xc.Nrows());
   hessianT = 0;
   OptppArray<SymmetricMatrix> hessian(1);
   hessian[0] = hessianT;
   return hessian;
}

SymmetricMatrix CompoundConstraint::evalHessian(ColumnVector& xc,
                                                const ColumnVector& mult) const 
{
   int k, tncons;
   SymmetricMatrix hessian(xc.Nrows());
   OptppArray<SymmetricMatrix> temp;
   ColumnVector type;
   Constraint test;

   k       = 0;
   hessian = 0.0;

   for(int i = 0; i < numOfSets_; i++){
     test   = constraints_[i];
     type   = test.getConstraintType();
     tncons = test.getNumOfCons();
     k     += tncons;

     if(type(1) == NLeqn || type(1) == NLineq){
        temp   = test.evalHessian(xc, i);
	for(int j = 0; j < temp.length(); j++){
           hessian += mult(j+1)*temp[j];
        }
     }
   }
   return hessian;
}

bool CompoundConstraint::amIFeasible(const ColumnVector& xc, double epsilon ) const
{
   bool feasible = true;
   ColumnVector type;
   Constraint test;

   for(int i = 0; i < numOfSets_; i++){
     test = constraints_[i];
     type = test.getConstraintType();
     if(type(1) == Bound)
          feasible = test.amIFeasible(xc,epsilon);
     if(!feasible){
       //       OptppmathError("The current iterate is infeasible wrt to the bound and
       //                  linear constraints");
       break;
     }
   }
   return feasible;
}

void CompoundConstraint::insertSort(const OptppArray<Constraint>& constraints)
{
   Constraint ctemp;
   int dim = constraints.length();
   OptppArray<Constraint> sorted(dim);
   sorted = constraints;
   int i;

   if(dim == 1)
      constraints_ = sorted;
   else{
      for(int j = 1; j < dim; j++){
         ctemp = sorted[j];
	 i = j - 1;
	 while(i > -1 && compare(sorted[i],ctemp) > 0){
	    sorted[i+1] = sorted[i];
	    i--;
         }
	 sorted[i+1] = ctemp;
      }
      constraints_ = sorted;
  }
}

void CompoundConstraint::insertSort()
{
   Constraint ctemp;
   int dim = constraints_.length();
   int i;

   if(dim > 1){
      for(int j = 1; j < dim; j++){
         ctemp = constraints_[j];
	 i = j - 1;
	 while(i > -1 && compare(constraints_[i],ctemp) > 0){
	    constraints_[i+1] = constraints_[i];
	    i--;
         }
	 constraints_[i+1] = ctemp;
      }
   }
}

int CompoundConstraint::compare(const Constraint& c1, const Constraint& c2)
{
   ColumnVector ct1 = c1.getConstraintType();
   ColumnVector ct2 = c2.getConstraintType();
   
   if(ct1(1) < ct2(1))
      return -1;
   else if(ct1(1) > ct2(1))
     return 1;
   else
     return 0;
  
}

} // namespace OPTT
