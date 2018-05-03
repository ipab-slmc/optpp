//------------------------------------------------------------------------
// Copyright (C) 1993: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//------------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#ifdef HAVE_STD
#include <cfloat>
#include <cstring>
#else
#include <float.h>
#include <string.h>
#endif

#include "NLF.h"

using namespace std;
using NEWMAT::ColumnVector;
using NEWMAT::Matrix;
using NEWMAT::SymmetricMatrix;


namespace OPTPP {

//-------------------------------------------------------------------------
// NLF0 method routines.
//-------------------------------------------------------------------------
void NLF0::reset() // Reset parameter values 
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

void NLF0::initFcn() // Initialize Function
{
  if (init_flag == false)  {
    init_fcn(dim, mem_xc);
    init_flag = true;
  }
  else  {
//    cerr << "NLF0:initFcn: Warning - initialization called twice\n";
    init_fcn(dim, mem_xc);
  }
}

double NLF0::evalF() // Evaluate Function
{
  int result = 0;
  double time0 = get_wall_clock_time();

  if (SpecFlag == NoSpec) {
    if (!application.getF(mem_xc,fvalue)) {
      fcn_v(dim, mem_xc, fvalue, result, vptr);
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

double NLF0::evalF(const ColumnVector& x) // Evaluate Function at x
{
  double fx;
  int result = 0;
  double time0 = get_wall_clock_time();

  if (SpecFlag == NoSpec) {
    if (!application.getF(x,fx)) {
      fcn_v(dim, x, fx, result, vptr);
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

ColumnVector NLF0::evalG() 
{
  ColumnVector grad(dim);
  ColumnVector sx(dim);
  sx = 1.0;

  // Since NLF0 objects do not have analytic gradients supply
  // one by using finite differences

  grad = FDGrad(sx, mem_xc, fvalue, partial_grad);
  return grad;
}

ColumnVector NLF0::evalG(const ColumnVector& x) 
{
  ColumnVector gx(dim);
  ColumnVector sx(dim);
  sx = 1.0;

  // Since NLF0 objects do not have analytic gradients supply
  // one by using finite differences

  if (SpecFlag == NoSpec) {
    int result = 0;
    if (!application.getF(x, specF)) {
      fcn_v(dim, x, specF, result, vptr);
      nfevals++;
    }
  }

  gx = FDGrad(sx, x, specF, partial_grad);
  return gx;
}

SymmetricMatrix NLF0::evalH() 
{
// Since NLF0 objects do not have analytic hessians supply
// one by using finite differences

  SymmetricMatrix hess(dim);

  hess = FD2Hessian(mem_xc);
  return hess;
}

SymmetricMatrix NLF0::evalH(ColumnVector& x) 
{
// Since NLF0 objects do not have analytic hessians supply
// one by using finite differences

  SymmetricMatrix hess(dim);

  hess = FD2Hessian(x);
  return hess;
}

void NLF0::eval()
{
  (void) evalF();
}

real NLF0::evalLagrangian(const ColumnVector& xc , 
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

ColumnVector NLF0::evalLagrangianGradient(const ColumnVector& xc, 
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


ColumnVector NLF0::evalCF(const ColumnVector& x) // Evaluate Nonlinear Constraint at x
{
  ColumnVector cfx(ncnln);
  int result = 0;

  double time0 = get_wall_clock_time();


  // *** CHANGE *** //
  if (!application.getCF(x,cfx)) {
    confcn(dim, x, cfx, result);
    application.constraint_update(NLPFunction,dim,ncnln,x,cfx);
   // nfevals++;
  }
  // *** CHANGE *** //
  function_time = get_wall_clock_time() - time0;
  
  setConstraintValue(cfx);
  return cfx;
}

Matrix NLF0::evalCG(const ColumnVector& x) // Evaluate Nonlinear Constraint Gradient at x
{
// Since NLF0 objects do not have analytic gradients supply
// one by using finite differences

  Matrix grad(dim,ncnln);
  grad = CONFDGrad(x);
  return grad;
}

SymmetricMatrix NLF0::evalCH(ColumnVector& x) 
{
// CPJW - This is a placeholder routine.  NIPS is the only algorithm which
// supports nonlinear constraints and currently this routine is never accessed.
// The true evaluator will be implemented later.
// Since NLF0 objects do not have analytic hessians supply
// one by using finite differences

  SymmetricMatrix hess(dim);
  hess = 0.0;
  return hess;
}

OptppArray<SymmetricMatrix> NLF0::evalCH(ColumnVector& x, int darg) 
{
// CPJW - This is a placeholder routine.  NIPS is the only algorithm which
// supports nonlinear constraints and currently this routine is never accessed.
// The true evaluator will be implemented later.
// Since NLF0 objects do not have analytic hessians supply
// one by using finite differences

  OptppArray<SymmetricMatrix> hess(1);
  SymmetricMatrix hessT(dim);
  hessT = 0.0;
  hess[0] = hessT;
  return hess;
}

void NLF0::evalC(const ColumnVector& x)
{
  (void) evalCF(x);
}

} // namespace OPTPP
