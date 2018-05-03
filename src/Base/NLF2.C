//------------------------------------------------------------------------
// Copyright (C) 1993: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//------------------------------------------------------------------------

#include "NLF.h"
#include "TOLS.h"
#include "cblas.h"

using namespace std;
using NEWMAT::ColumnVector;
using NEWMAT::Matrix;
using NEWMAT::SymmetricMatrix;

namespace OPTPP {

void NLF2::reset() // Reset parameter values  
{
  init_flag = false;    
  nfevals   = ngevals = nhevals = 0; 
#ifdef WITH_MPI
  SpecFlag = Spec1;
#else
  SpecFlag = NoSpec;
#endif
  application.reset();
}

void NLF2::initFcn() // Initialize Function
{
  if (init_flag == false) {
    init_fcn(dim, mem_xc);
    init_flag = true;
  }
  else {
//    cerr << "NLF2:initFcn: Warning - initialization called twice\n";
    init_fcn(dim, mem_xc);
  }
}

double NLF2::evalF() // Evaluate Function
{
  int  result = 0;
  ColumnVector gtmp(dim);
  SymmetricMatrix Htmp(dim);
  //cout << "NLF2:evalF \n";

  double time0 = get_wall_clock_time();
  // *** CHANGE *** //
  if (!application.getF(mem_xc,fvalue)) {
    fcn_v(NLPFunction, dim, mem_xc, fvalue, gtmp, Htmp,result,vptr);
    application.update(result,dim,mem_xc,fvalue,gtmp,Htmp);
    nfevals++;
  }
  // *** CHANGE *** //
  function_time = get_wall_clock_time() - time0;

  if (debug_)  cout << "NLF2::evalF()\n" 
                   << "nfevals       = " << nfevals   << "\n"
                   << "fvalue        = " << fvalue << "\n"
                   << "function time = " << function_time << "\n";
  return fvalue;
}

double NLF2::evalF(const ColumnVector& x) // Evaluate Function at x
{
  int    result = 0;
  double fx;
  ColumnVector gtmp(dim);
  SymmetricMatrix Htmp(dim);

  double time0 = get_wall_clock_time();
  // *** CHANGE *** //
  if (!application.getF(x,fx)) {
    fcn_v(NLPFunction, dim, x, fx, gtmp, Htmp,result,vptr);
    application.update(result,dim,x,fx,gtmp,Htmp);
    nfevals++;
  }
  // *** CHANGE *** //
  function_time = get_wall_clock_time() - time0;

  if (debug_) cout << "NLF2::evalF(x)\n" 
                  << "nfevals       = " << nfevals   << "\n"
                  << "fvalue        = " << fvalue << "\n"
                  << "function time = " << function_time << "\n";
  return fx;
}

ColumnVector NLF2::evalG() // Evaluate the gradient
{
  int    result = 0;
  double fx;
  SymmetricMatrix Htmp(dim);

  // *** CHANGE *** //
  if (!application.getGrad(mem_xc,mem_grad)) {
    fcn_v(NLPGradient, dim, mem_xc, fx, mem_grad, Htmp,result,vptr);
    application.update(result,dim,mem_xc,fx,mem_grad,Htmp);
    ngevals++;
  }
  // *** CHANGE *** //
  return mem_grad;
}

ColumnVector NLF2::evalG(const ColumnVector& x) // Evaluate the gradient at x
{
  int    result = 0;
  double fx;
  ColumnVector gx(dim);
  SymmetricMatrix Htmp(dim);

  // *** CHANGE *** //
  if (!application.getGrad(x,gx)) {
    fcn_v(NLPGradient, dim, x, fx, gx, Htmp, result,vptr);
    application.update(result,dim,x,fx,gx,Htmp);
    ngevals++;
  }
  // *** CHANGE *** //
  return gx;
}

SymmetricMatrix NLF2::evalH() // Evaluate the Hessian
{
  int    result = 0;
  double fx;
  ColumnVector gtmp(dim);

  // *** CHANGE *** //
  if (!application.getHess(mem_xc,Hessian)) {
    fcn_v(NLPHessian, dim, mem_xc, fx, gtmp, Hessian, result,vptr);
    application.update(result,dim,mem_xc,fx,gtmp,Hessian);
    nhevals++;
  }
  // *** CHANGE *** //
  return Hessian;
}

SymmetricMatrix NLF2::evalH(ColumnVector& x) // Evaluate the hessian at x
{
  int    result = 0;
  double fx;
  ColumnVector gx(dim);
  SymmetricMatrix Hx(dim);

  // *** CHANGE *** //
  if (!application.getHess(x,Hx)) {
    fcn_v(NLPHessian, dim, x, fx, gx, Hx, result,vptr);
    application.update(result,dim,x,fx,gx,Hx);
    nhevals++;
  }
  // *** CHANGE *** //
  return Hx;
}

void NLF2::eval() // Evaluate Function, Gradient, and Hessian
{
  int mode = NLPFunction | NLPGradient | NLPHessian, result = 0;

  double time0 = get_wall_clock_time();
  // *** CHANGE *** //
  if (!application.getF(mem_xc,fvalue) || !application.getGrad(mem_xc,mem_grad) ||
      !application.getHess(mem_xc,Hessian)) { 
    fcn_v(mode, dim, mem_xc, fvalue, mem_grad, Hessian,result,vptr);
    application.update(result,dim,mem_xc,fvalue,mem_grad,Hessian);
    nfevals++; ngevals++; nhevals++;
  }
  // *** CHANGE *** //
  function_time = get_wall_clock_time() - time0;

  if (debug_)  cout << "NLF2::eval()\n" 
                   << "mode          = " << mode   << "\n"
                   << "nfevals       = " << nfevals   << "\n"
	           << "fvalue        = " << fvalue << "\n"
	           << "function time = " << function_time << "\n";
  
}

real NLF2::evalLagrangian(const ColumnVector& xc , 
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

ColumnVector NLF2::evalLagrangianGradient(const ColumnVector& xc, 
                                          const ColumnVector& multiplier,
					  const ColumnVector& type) 
{
   ColumnVector grad  = evalG(xc);
   if(hasConstraints()){
      ColumnVector tmult = -multiplier;
      grad += constraint_->evalGradient(xc)*tmult;
   }
   return grad;
}

SymmetricMatrix NLF2::evalLagrangianHessian(ColumnVector& xc,
                                            const ColumnVector& multiplier,
                                            const ColumnVector& type) 
{
   SymmetricMatrix hessian = evalH(xc);
   if(hasConstraints()){
       SymmetricMatrix thessian = constraint_->evalHessian(xc);
       Print(thessian);
   }
   return hessian;
}

ColumnVector NLF2::evalCF(const ColumnVector& x) // Evaluate Function at x
{
  int    result = 0;
  ColumnVector cfx(ncnln);
  Matrix gtmp(dim,ncnln);
  OptppArray<SymmetricMatrix> Htmp(ncnln);

  double time0 = get_wall_clock_time();
  // *** CHANGE *** //
  if (!application.getCF(x,cfx)) {

    if(confcn1 != NULL){   
       confcn1(NLPFunction, dim, x, cfx, gtmp, result);
       application.constraint_update(result,dim,ncnln,x,cfx,gtmp);
    }
    else if(confcn2 != NULL){   
       confcn2(NLPFunction, dim, x, cfx, gtmp, Htmp,result);
       application.constraint_update(result,dim,ncnln,x,cfx,gtmp,Htmp);
    }
  }
  // *** CHANGE *** //
  function_time = get_wall_clock_time() - time0;

  if (debug_) cout << "NLF2::evalCF(x)\n" 
                  << "nfevals       = " << nfevals   << "\n"
                  << "fvalue(1)        = " << cfx(1) << "\n"
                  << "function time = " << function_time << "\n";
  return cfx;
}

Matrix NLF2::evalCG(const ColumnVector& x) // Evaluate the gradient at x
{
  int    result = 0;
  ColumnVector cfx(ncnln);
  Matrix cgx(dim,ncnln);
  OptppArray<SymmetricMatrix> Htmp(ncnln);

  // *** CHANGE *** //
  if (!application.getCGrad(x,cgx)) {
    if(confcn1 != NULL){
      confcn1(NLPGradient, dim, x, cfx, cgx, result);
      application.constraint_update(result,dim,ncnln,x,cfx,cgx);
    }
    if(confcn2 != NULL){
      confcn2(NLPGradient, dim, x, cfx, cgx, Htmp, result);
      application.constraint_update(result,dim,ncnln,x,cfx,cgx,Htmp);
    }
  }
  // *** CHANGE *** //
  return cgx;
}

SymmetricMatrix NLF2::evalCH(ColumnVector& x) // Evaluate the hessian at x
{
  ColumnVector cfx(ncnln);
  Matrix cgx(dim,ncnln);
  SymmetricMatrix cHx(dim);

  cHx = 0;
  return cHx;
}
OptppArray<SymmetricMatrix> NLF2::evalCH(ColumnVector& x, int darg) // Evaluate the hessian at x
{
  int    result = 0;
  ColumnVector cfx(ncnln);
  Matrix cgx(dim,ncnln);
  OptppArray<SymmetricMatrix> cHx(ncnln);

  // *** CHANGE *** //
  if (!application.getCHess(x,cHx)) {
    if(confcn2 != NULL){
       confcn2(NLPHessian, dim, x, cfx, cgx, cHx, result);
       application.constraint_update(result,dim,ncnln,x,cfx,cgx,cHx);
       nhevals++;
    }
  }
  // *** CHANGE *** //
  return cHx;
}

void NLF2::evalC(const ColumnVector& x)
{
  int mode1 = NLPFunction | NLPGradient;
  int mode2 = NLPFunction | NLPGradient | NLPHessian;
  int result = 0;
  ColumnVector cfx(ncnln);
  Matrix cgx(dim,ncnln);
  OptppArray<SymmetricMatrix> cHx(ncnln);

  double time0 = get_wall_clock_time();

  // *** CHANGE *** //
  if (!application.getCF(x, cfx) || !application.getCGrad(x, cgx) || !application.getCHess(x, cHx)) {
    if(confcn1 != NULL){
      confcn1(mode1, dim, x, cfx, cgx, result);
      application.constraint_update(result, dim, ncnln, x, cfx, cgx);
    }
    if(confcn2 != NULL){
       confcn2(mode2, dim, x, cfx, cgx, cHx, result);
       application.constraint_update(result, dim, ncnln, x, cfx, cgx, cHx);
       nhevals++;
    }
  }

  function_time = get_wall_clock_time() - time0;
}

} // namespace OPTPP
