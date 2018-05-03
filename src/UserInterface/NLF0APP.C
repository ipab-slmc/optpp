//------------------------------------------------------------------------
// Copyright (C) 1993: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//------------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

using namespace std;

#include <iostream>
#ifdef HAVE_STD
#include <cfloat>
#include <cstring>
#else
#include <float.h>
#include <string.h>
#endif


#include "NLFAPP.h"
#include "NLP0.h"

using NEWMAT::ColumnVector;
using NEWMAT::Matrix;
using NEWMAT::SymmetricMatrix;

namespace OPTPP {

//-------------------------------------------------------------------------
// NLF0APP method routines.
//-------------------------------------------------------------------------
void NLF0APP::reset() // Reset parameter values
{
  init_flag = false;
  nfevals   = 0;
#ifdef WITH_MPI
  SpecFlag = Spec1; 
#else
  SpecFlag = NoSpec;
#endif
  application.reset();
}

void NLF0APP::initFcn() // Initialize Function
{
  if (init_flag == false)  {
    init_fcn(dim, mem_xc, launcher_);
    init_flag = true;
  }
  else  {
//    cout << "NLF0APP:initFcn: Warning - initialization called twice\n";
    init_fcn(dim, mem_xc, launcher_);
  }
}

double NLF0APP::evalF() // Evaluate Function
{
  int result = 0;
  double time0 = get_wall_clock_time();

  if (SpecFlag == NoSpec) {
    if (!application.getF(mem_xc,fvalue)) {
      fcn(dim, mem_xc, fvalue, result, launcher_);
      application.update(NLPFunction,dim,mem_xc,fvalue);
      nfevals++;
    }
  }
  else {
    SpecFlag = Spec1;
    (void) evalG();
    SpecFlag = Spec2;
  }

  function_time = get_wall_clock_time() - time0;
  return fvalue;
}

double NLF0APP::evalF(const ColumnVector& x) // Evaluate Function at x
{
  double fx;
  int result = 0;
  double time0 = get_wall_clock_time();

  if (SpecFlag == NoSpec) {
    if (!application.getF(x,fx)) {
      fcn(dim, x, fx, result, launcher_);
      application.update(NLPFunction,dim,x,fx);
      nfevals++;
    }
  }
  else {
    SpecFlag = Spec1;
    (void) evalG(x);
    fx = specF;
    SpecFlag = Spec2;
  }

  function_time = get_wall_clock_time() - time0;
  return fx;
}

ColumnVector NLF0APP::evalG() 
{
  ColumnVector grad(dim);
  ColumnVector sx(dim);
  sx = 1.0;

  // Since NLF0APP objects do not have analytic gradients supply
  // one by using finite differences

  grad = FDGrad(sx, mem_xc, fvalue, partial_grad);
  return grad;
}

ColumnVector NLF0APP::evalG(const ColumnVector& x) 
{
  ColumnVector gx(dim);
  ColumnVector sx(dim);
  sx = 1.0;

  // Since NLF0APP objects do not have analytic gradients supply
  // one by using finite differences

  if (SpecFlag == NoSpec) {
    int result = 0;
    if (!application.getF(x, specF)) {
      fcn(dim, x, specF, result, launcher_);
      nfevals++;
    }
  }

  gx = FDGrad(sx, x, specF, partial_grad);
  return gx;
}

SymmetricMatrix NLF0APP::evalH() 
{
// Since NLF0APP objects do not have analytic hessians supply
// one by using finite differences

  SymmetricMatrix hess(dim);
  hess = FD2Hessian(mem_xc);
  return hess;
}

SymmetricMatrix NLF0APP::evalH(ColumnVector& x) 
{
// Since NLF0APP objects do not have analytic hessians supply
// one by using finite differences

  SymmetricMatrix hess(dim);
  hess = FD2Hessian(x);
  return hess;
}

void NLF0APP::eval()
{
  (void) evalF();
}

real NLF0APP::evalLagrangian(const ColumnVector& xc , 
                          ColumnVector& multiplier,
                          const ColumnVector& type)
{
   real result = evalF(xc);
   if( hasConstraints()){
      ColumnVector resid = constraint_->evalResidual(xc);
      result  -=  Dot(resid, multiplier);
   }
   return result;
}

ColumnVector NLF0APP::evalLagrangianGradient(const ColumnVector& xc, 
                                          const ColumnVector& multiplier,
					  const ColumnVector& type) 
{
   ColumnVector grad  = evalG(xc);
   if(hasConstraints()){
      ColumnVector tmult = multiplier;
      for (int i = 1; i < getNumOfCons() + 1; i++){
         if(type(i) == NLineq || type(i) == Lineq)
            tmult(i)*= -1;
      }
      grad += constraint_->evalGradient(xc)*tmult;
   }
   return grad;
}


ColumnVector NLF0APP::evalCF(const ColumnVector& x) // Evaluate Nonlinear Constraint at x
{
  ColumnVector cfx(ncnln);
  int result = 0;

  double time0 = get_wall_clock_time();


  // *** CHANGE *** //
  if (!application.getCF(x,cfx)) {
    confcn(dim, ncnln, x, cfx, result, launcher_);
    application.constraint_update(NLPFunction,dim,ncnln,x,cfx);
   // nfevals++;
  }
  // *** CHANGE *** //
  function_time = get_wall_clock_time() - time0;
  
  setConstraintValue(cfx);
  return cfx;
}

Matrix NLF0APP::evalCG(const ColumnVector& x) // Evaluate Nonlinear Constraint Gradient at x
{
// Since NLF0APP objects do not have analytic gradients supply
// one by using finite differences

  Matrix grad(dim,ncnln);
  grad = CONFDGrad(x);
  return grad;
}

SymmetricMatrix NLF0APP::evalCH(ColumnVector& x) 
{
// CPJW - This is a placeholder routine.  NIPS is the only algorithm which
// supports nonlinear constraints and currently this routine is never accessed.
// The true evaluator will be implemented later.
// Since NLF0APP objects do not have analytic hessians supply
// one by using finite differences

  SymmetricMatrix hess(dim);
  hess = 0.0;
  return hess;
}

OptppArray<SymmetricMatrix> NLF0APP::evalCH(ColumnVector& x, int darg) 
{
// CPJW - This is a placeholder routine.  NIPS is the only algorithm which
// supports nonlinear constraints and currently this routine is never accessed.
// The true evaluator will be implemented later.
// Since NLF0APP objects do not have analytic hessians supply
// one by using finite differences

  OptppArray<SymmetricMatrix> hess(1);
  SymmetricMatrix hessT(dim);
  hessT = 0.0;
  hess[0] = hessT;
  return hess;
}

} // namespace OPTPP
