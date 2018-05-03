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

void NLF1::reset() // Reset parameter values  
{
  init_flag = false;    
  nfevals   = ngevals = 0; 
#ifdef WITH_MPI
  SpecFlag = Spec1;
#else
  SpecFlag = NoSpec;
#endif
  application.reset();
}

void NLF1::initFcn() // Initialize Function
{
  if (init_flag == false)  {
      init_fcn(dim, mem_xc);
      init_flag = true;
  }
  else  {
//    cerr << "NLF1:initFcn: Warning - initialization called twice\n";
    init_fcn(dim, mem_xc);
  }
}

double NLF1::evalF() // Evaluate Function
{
  int result = 0;
  ColumnVector gtmp(dim);

  double time0 = get_wall_clock_time();
  // *** CHANGE *** //
  if (!application.getF(mem_xc,fvalue)) {
    fcn_v(NLPFunction, dim, mem_xc, fvalue, gtmp, result, vptr);
    application.update(result,dim,mem_xc,fvalue,gtmp);
    nfevals++;
  }
  // *** CHANGE *** //
  function_time = get_wall_clock_time() - time0;

  if (debug_)
  cout << "NLF1::evalF()\n" 
    << "nfevals       = " << nfevals << "\n"
    << "fvalue        = " << fvalue << "\n"
    << "function time = " << function_time << "\n";
  return fvalue;
}

double NLF1::evalF(const ColumnVector& x) // Evaluate Function at x
{
  int    result = 0;
  double fx;
  ColumnVector gtmp(dim);

  double time0 = get_wall_clock_time();
  // *** CHANGE *** //
  if (!application.getF(x,fx)) {
    fcn_v(NLPFunction, dim, x, fx, gtmp, result, vptr);
    application.update(result,dim,x,fx,gtmp);
    nfevals++;
  }
  // *** CHANGE *** //
  function_time = get_wall_clock_time() - time0;

  if (debug_)
  cout << "NLF1::evalF(x)\n" 
    << "nfevals       = " << nfevals << "\n"
    << "fvalue        = " << fx << "\n"
    << "function time = " << function_time << "\n";
  return fx;
}

ColumnVector NLF1::evalG() // Evaluate the gradient
{
  int    result = 0;
  double fx;

  // *** CHANGE *** //
  if (!application.getGrad(mem_xc,mem_grad)) {
    fcn_v(NLPGradient, dim, mem_xc, fx, mem_grad, result, vptr);
    application.update(result,dim,mem_xc,fx,mem_grad);
    ngevals++;
  }
  // *** CHANGE *** //
  return mem_grad;
}

ColumnVector NLF1::evalG(const ColumnVector& x) // Evaluate the gradient at x
{
  int    result = 0 ;
  double fx;
  ColumnVector gx(dim);

  // *** CHANGE *** //
  if (!application.getGrad(x,gx)) {
    fcn_v(NLPGradient, dim, x, fx, gx, result, vptr);
    application.update(result,dim,x,fx,gx);
    ngevals++;
  }
  // *** CHANGE *** //
  return gx;
}

SymmetricMatrix NLF1::evalH() // Evaluate the Hessian
{
  ColumnVector sx(dim);
  SymmetricMatrix Hessian(dim);

  sx = 1.0;
  Hessian = FDHessian(sx);
  return Hessian;
}

SymmetricMatrix NLF1::evalH(ColumnVector& x) // Evaluate the Hessian at x
{
  SymmetricMatrix Hessian(dim);

  Hessian = FDHessian(x);
  return Hessian;
}

void NLF1::eval() // Evaluate Function and Gradient
{
  int mode = NLPFunction | NLPGradient, result = 0;

  double time0 = get_wall_clock_time();
  // *** CHANGE *** //
  if (!application.getF(mem_xc,fvalue) || !application.getGrad(mem_xc,mem_grad)) { 
    fcn_v(mode, dim, mem_xc, fvalue, mem_grad, result, vptr);
    application.update(result,dim,mem_xc,fvalue,mem_grad);
    nfevals++; ngevals++;
  }
  // *** CHANGE *** //

  function_time = get_wall_clock_time() - time0;

  if (debug_)
  cout << "NLF1::eval()\n" 
    << "mode          = " << mode   << "\n"
    << "nfevals       = " << nfevals << "\n"
    << "fvalue        = " << fvalue << "\n"
    << "function time = " << function_time << "\n";
}

real NLF1::evalLagrangian(const ColumnVector& xc , 
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

ColumnVector NLF1::evalLagrangianGradient(const ColumnVector& xc, 
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


ColumnVector NLF1::evalCF(const ColumnVector& x) // Evaluate Function at x
{
  int    result = 0;
  ColumnVector cfx(ncnln);
  Matrix gtmp(dim,ncnln);

  double time0 = get_wall_clock_time();
  // *** CHANGE *** //
  if (!application.getCF(x,cfx)) {
    confcn(NLPFunction, dim, x, cfx, gtmp, result);
    application.constraint_update(result,dim,ncnln,x,cfx,gtmp);
  }
  // *** CHANGE *** //
  function_time = get_wall_clock_time() - time0;

  if (debug_)
    cout << "NLF1::evalCF(x)\n" 
         << "nfevals       = " << nfevals << "\n"
       //  << "fvalue        = " << cfx << "\n"
         << "function time = " << function_time << "\n";
  return cfx;
}

Matrix NLF1::evalCG(const ColumnVector& x) // Evaluate the gradient at x
{
  int    result = 0 ;
  ColumnVector cfx(ncnln);
  Matrix cgx(dim,ncnln);

  // *** CHANGE *** //
  if (!application.getCGrad(x,cgx)) {
    confcn(NLPGradient, dim, x, cfx, cgx, result);
    application.constraint_update(result,dim,ncnln,x,cfx,cgx);
  }
  // *** CHANGE *** //
  return cgx;
}

SymmetricMatrix NLF1::evalCH(ColumnVector& x) // Evaluate the Hessian at x
{
  SymmetricMatrix Hessian(dim);

  // PJW This is a dummy routine.  NIPS is the only algorithm which supports
  // nonlinear constraints and currently this routine is not being accessed.
  Hessian = 0.0;
  return Hessian;
}


OptppArray<SymmetricMatrix> NLF1::evalCH(ColumnVector& x, int darg) // Evaluate the Hessian at x
{
  OptppArray<SymmetricMatrix> Hessian(ncnln);
  Hessian = CONFDHessian(x);
  return Hessian;
}

void NLF1::evalC(const ColumnVector& x)
{
  int mode = NLPFunction | NLPGradient, result = 0;
  ColumnVector cfx(ncnln);
  Matrix cgx(dim, ncnln);

  double time0 = get_wall_clock_time();

  if (!application.getCF(x, cfx) || !application.getCGrad(x, cgx)) {
    confcn(mode, dim, x, cfx, cgx, result);
    application.constraint_update(result, dim, ncnln, x, cfx, cgx);
  }

  function_time = get_wall_clock_time() - time0;
}

} // namespace OPTPP
