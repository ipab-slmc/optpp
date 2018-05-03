//------------------------------------------------------------------------
// P.J. Williams
// Sandia National Laboratories
// pwillia@sandia.gov
// Last modified 11/16/1999
//------------------------------------------------------------------------

#include "NLP.h"

using NEWMAT::ColumnVector;
using NEWMAT::Matrix;
using NEWMAT::SymmetricMatrix;

namespace OPTPP {

// Constructors 
NLP::NLP(): ptr_(0){;}
NLP::NLP(NLPBase* base): ptr_(base){;}

void NLP::setX(const int i, const real& x) 
{
   ptr_->setX(i,x);
}

void NLP::setX(const ColumnVector& x) 
{
   ptr_->setX(x);
}

void NLP::setF(const real& fx) 
{
   ptr_->setF(fx);
}

void NLP::setIsExpensive(const int e) 
{
   ptr_->setIsExpensive(e);
}

void NLP::setFcnAccrcy(const int i, const real& accrcy) 
{
   ptr_->setFcnAccrcy(i,accrcy);
}

void NLP::setFcnAccrcy(const ColumnVector& accrcy) 
{
   ptr_->setFcnAccrcy(accrcy);
}

int NLP::getDim() const
{
   int result = ptr_ -> getDim();
   return result;
}

int NLP::getFevals() const
{
   int result = ptr_ -> getFevals();
   return result;
}

int NLP::getIsExpensive() const
{
   int result = ptr_ -> getIsExpensive();
   return result;
}

real NLP::getF() const
{
   real result = ptr_ -> getF();
   return result;
}

ColumnVector NLP::getFcnAccrcy() const
{
   ColumnVector result = ptr_ -> getFcnAccrcy();
   return result;
}

ColumnVector NLP::getXc() const
{
   ColumnVector result = ptr_ -> getXc();
   return result;
}

real NLP::getFcnTime() const
{
   real result = ptr_ -> getFcnTime();
   return result;
}

int NLP::getNumOfCons() const
{
   int result = ptr_ -> getNumOfCons();
   return result;
}

int NLP::getNumOfNLCons() const
{
   int result = ptr_ -> getNumOfNLCons();
   return result;
}


bool NLP::hasConstraints() 
{
   bool result = ptr_ -> hasConstraints();
   return result;
}

void NLP::printConstraints() 
{
   ptr_ -> printConstraints();
}

void NLP::setDebug() 
{
   ptr_ -> setDebug();
}

bool NLP::getDebug() const 
{
   bool result = ptr_ -> getDebug();
   return result;
}

void NLP::reset()
{
   ptr_->reset();
}

void NLP::initFcn()
{
   ptr_->initFcn();
}

void NLP::eval()
{
   ptr_->eval();
}

real NLP::evalF()
{
   real result = ptr_->evalF();
   return result;
}

real NLP::evalF(const ColumnVector& x)
{
   real result = ptr_->evalF(x);
   return result;
}

ColumnVector NLP::evalG()
{
   ColumnVector result = ptr_->evalG();
   return result;
}

ColumnVector NLP::evalG(const ColumnVector& x)
{
   ColumnVector result = ptr_->evalG(x);
   return result;
}

SymmetricMatrix NLP::evalH()
{
   SymmetricMatrix result = ptr_->evalH();
   return result;
}

SymmetricMatrix NLP::evalH(ColumnVector& x)
{
   SymmetricMatrix result = ptr_->evalH(x);
   return result;
}

ColumnVector NLP::evalCF(const ColumnVector& x)
{
   ColumnVector result = ptr_->evalCF(x);
   return result;
}


Matrix NLP::evalCG(const ColumnVector& x)
{
   Matrix result = ptr_->evalCG(x);
   return result;
}

SymmetricMatrix NLP::evalCH(ColumnVector& x)
{
   SymmetricMatrix result = ptr_->evalCH(x);
   return result;
}

OptppArray<SymmetricMatrix> NLP::evalCH(ColumnVector& x, int darg)
{
   OptppArray<SymmetricMatrix> result = ptr_->evalCH(x,darg);
   return result;
}

void NLP::evalC(const ColumnVector& x)
{
  ptr_->evalC(x);
}

void NLP::printState(char* s)
{
   ptr_->printState(s);
}

void NLP::fPrintState(ostream *nlpout, char* s)
{
   ptr_->fPrintState(nlpout,s);
}

} // namespace OPTPP
