//------------------------------------------------------------------------
// P.J. Williams
// Sandia National Laboratories
// pwillia@sandia.gov
// Last modified 03/20/2003
//------------------------------------------------------------------------

#include "Constraint.h"

using NEWMAT::ColumnVector;
using NEWMAT::Matrix;
using NEWMAT::SymmetricMatrix;

namespace OPTPP {

// Constructors 
Constraint::Constraint(): ptr_(0){;}
Constraint::Constraint(ConstraintBase* base): ptr_(base){;}

int Constraint::getNumOfCons() const
{
  int result = ptr_->getNumOfCons();
  return result;
}

int Constraint::getNumOfVars() const
{
  int result = ptr_->getNumOfVars();
  return result;
}

ColumnVector Constraint::getLower() const
{
  ColumnVector result = ptr_->getLower();
  return result;
}

ColumnVector Constraint::getUpper() const
{
  ColumnVector result = ptr_->getUpper();
  return result;
}

ColumnVector Constraint::getConstraintType() const
{
  ColumnVector result = ptr_->getConstraintType();
  return result;
}

ColumnVector Constraint::getConstraintValue() const
{
  ColumnVector result = ptr_->getConstraintValue();
  return result;
}

ColumnVector Constraint::getConstraintViolation() const
{
  ColumnVector result = ptr_->getConstraintViolation();
  return result;
}

OptppArray<int> Constraint::getConstraintMappingIndices() const
{
  OptppArray<int> result = ptr_->getConstraintMappingIndices();
  return result;
}

void Constraint::evalCFGH(const ColumnVector & xcurrent) const
{
   ptr_->evalCFGH(xcurrent);
}

ColumnVector Constraint::evalResidual(const ColumnVector& xcurrent) const 
{
   ColumnVector result;
   result = ptr_->evalResidual(xcurrent);
   return result;
}

Matrix Constraint::evalGradient(const ColumnVector& xcurrent) const 
{
   Matrix result;
   result = ptr_->evalGradient(xcurrent);
   return result;
}

SymmetricMatrix Constraint::evalHessian(ColumnVector& xcurrent) const 
{
   SymmetricMatrix result;
   result = ptr_->evalHessian(xcurrent);
   return result;
}

OptppArray<SymmetricMatrix> Constraint::evalHessian(ColumnVector& xcurrent, int darg) const 
{
   OptppArray<SymmetricMatrix> result;
   result = ptr_->evalHessian(xcurrent,darg);
   return result;
}

bool Constraint::amIFeasible(const ColumnVector& xcurrent,double epsilon) const 
{
   bool result;
   result =  ptr_->amIFeasible(xcurrent,epsilon);
   return result;
}

} // namespace OPTPP
