//------------------------------------------------------------------------
// Copyright (C) 1993: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//------------------------------------------------------------------------

using namespace std;

#include "cblas.h"
#include "NLFAPP.h"
#include "NLP1.h"
#include "TOLS.h"
#include "ioformat.h"

using NEWMAT::ColumnVector;
using NEWMAT::Matrix;
using NEWMAT::SymmetricMatrix;

namespace OPTPP {

//-------------------------------------------------------------------------
// FDNLF1APP method routines.
//-------------------------------------------------------------------------

void FDNLF1APP::reset() // Reset parameter values
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

void FDNLF1APP::initFcn() // Initialize Function
{
  if (init_flag == false) {
      init_fcn(dim, mem_xc, launcher_);
      init_flag = true;
  }
  else  {
//    cout << "FDNLF1APP:initFcn: Warning - initialization called twice\n";
    init_fcn(dim, mem_xc, launcher_);
  }
}

void FDNLF1APP::eval() // Evaluate Function and Gradient
{
  (void) evalF();
  (void) evalG();
}

double FDNLF1APP::evalF() // Evaluate Function
{
  int result = 0;
  double time0 = get_wall_clock_time();

  if (SpecFlag == NoSpec) {
    if (!application.getF(mem_xc, fvalue)) {
      fcn(dim, mem_xc, fvalue, result, launcher_);
      function_time = get_wall_clock_time() - time0;
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

double FDNLF1APP::evalF(const ColumnVector& x) // Evaluate Function at x
{
  double fx;
  int result = 0;
  double time0 = get_wall_clock_time();

  if (SpecFlag == NoSpec) {
    if (!application.getF(x, fx)) {
      fcn(dim, x, fx, result, launcher_);
      function_time = get_wall_clock_time() - time0;
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

ColumnVector FDNLF1APP::evalG() // Evaluate the gradient
{ 
  ColumnVector sx(dim);
  sx = 1.0;
  ngevals++;

  if (finitediff == ForwardDiff)
    mem_grad =  FDGrad(sx, mem_xc, fvalue, partial_grad);
  else if (finitediff == BackwardDiff)
    mem_grad =  BDGrad(sx, mem_xc, fvalue, partial_grad);
  else if (finitediff == CentralDiff)
    mem_grad =  CDGrad(sx, mem_xc, fvalue, partial_grad);
  else {
//    cout << "FDNLF1APP::evalG: Unrecognized difference option\n";
//    cout << "FDNLF1APP::evalG: Using forward difference option\n";
    mem_grad =  FDGrad(sx, mem_xc, fvalue, partial_grad);
  }
 return mem_grad;
}

ColumnVector FDNLF1APP::evalG(const ColumnVector& x) // Evaluate the gradient at x
{
  ColumnVector gx(dim);
  ColumnVector sx(dim);
  sx = 1.0;
  ngevals++;

  if (SpecFlag == NoSpec) {
    int result = 0;
    if (!application.getF(x, specF)) {
      fcn(dim, x, specF, result, launcher_);
      nfevals++;
    }
  }

  if (finitediff == ForwardDiff) 
    gx = FDGrad(sx, x, specF, partial_grad);
  else if (finitediff == BackwardDiff)
    gx = BDGrad(sx, x, specF, partial_grad);
  else if (finitediff == CentralDiff)
    gx = CDGrad(sx, x, specF, partial_grad);
  else {
//    cout << "FDNLF1APP::evalG: Unrecognized difference option\n";
//    cout << "FDNLF1APP::evalG: Using forward difference option\n";
    mem_grad =  FDGrad(sx, x, specF, partial_grad);
  }
  
  return gx;
}

SymmetricMatrix FDNLF1APP::FDHessian(ColumnVector& sx) // Evaluate the Hessian
{
  SymmetricMatrix Hessian(dim);

  sx = 1.0;
  Hessian = FD2Hessian(sx);
  return Hessian;
}

SymmetricMatrix FDNLF1APP::evalH() // Evaluate the Hessian
{
  ColumnVector sx(dim);
  SymmetricMatrix Hessian(dim);

  sx = 1.0;
  Hessian = FD2Hessian(sx);
  return Hessian;
}

SymmetricMatrix FDNLF1APP::evalH(ColumnVector& x) // Evaluate the Hessian
{
  SymmetricMatrix Hessian(dim);

  Hessian = FD2Hessian(x);
  return Hessian;
}

void FDNLF1APP::printState(char * s) 
{ // Print out current state: x current, gradient and Function value
  cout <<"\n\n=========  " << s << "  ===========\n\n";
  cout <<"\n   i\t    xc \t\t grad \t\t fcn_accrcy \n";
  for (int i=1; i<=dim; i++) 
    cout << d(i,6) << e(mem_xc(i),12,4)<< "\t" << e(mem_grad(i),12,4) 
         << "\t"   << e(mem_fcn_accrcy(i),12,4) << "\n";
  cout <<"\nFunction Value     = " << e(fvalue,12,4) << "\n";
  double gnorm = Norm2(mem_grad);
  cout <<"Norm of gradient   = " << e(gnorm,12,4) << "\n";
  cout <<"Derivative Option  = " << finitediff << "\n\n";
}

void FDNLF1APP::fPrintState(ostream *nlpout, char * s) 
{ // Print out current state: x current, gradient and Function value
  (*nlpout) <<"\n\n=========  " << s << "  ===========\n\n";
  (*nlpout) <<"\n   i\t    xc \t\t grad \t\t fcn_accrcy \n";
  for (int i=1; i<=dim; i++) 
    (*nlpout) << d(i,6) << e(mem_xc(i),12,4)<< "\t" << e(mem_grad(i),12,4) 
         << "\t"   << e(mem_fcn_accrcy(i),12,4) << "\n";
  (*nlpout) <<"\nFunction Value     = " << e(fvalue,12,4) << "\n";
  double gnorm = Norm2(mem_grad);
  (*nlpout) <<"Norm of gradient   = " << e(gnorm,12,4) << "\n";
//  (*nlpout) <<"Function Accuracy  = " << e(mem_fcn_accrcy,12,4) << "\n";
  (*nlpout) <<"Derivative Option  = " << finitediff << "\n\n";
}


real FDNLF1APP::evalLagrangian(const ColumnVector& xc , 
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

ColumnVector FDNLF1APP::evalLagrangianGradient(const ColumnVector& xc, 
                                          const ColumnVector& multiplier,
                                          const ColumnVector& type) 
{
   mem_grad  = evalG(xc);
   ColumnVector grad  = mem_grad;
   if(hasConstraints())
      grad -= constraint_->evalGradient(xc)*multiplier;
   return grad;
}

ColumnVector FDNLF1APP::evalCF(const ColumnVector& x) // Evaluate Constraint Fcn at x
{
  int result = 0;
  ColumnVector cfx(ncnln);
  double time0 = get_wall_clock_time();
  confcn(dim, ncnln, x, cfx, result, launcher_);
  function_time = get_wall_clock_time() - time0;

  //nfevals++;
  return cfx;
}

Matrix FDNLF1APP::evalCG(const ColumnVector& x) // Evaluate the gradient at x
{
  ColumnVector sx(dim);
  sx = 1.0;
  ColumnVector xsave(dim);
  Matrix gx(dim, ncnln);
  xsave = getXc();
  setX(x);
  if (finitediff == ForwardDiff) 
    gx = CONFDGrad(sx);
  else if (finitediff == BackwardDiff)
    gx = CONBDGrad(sx);
  else if (finitediff == CentralDiff)
    gx = CONCDGrad(sx);
  else {
//    cout <<"FDNLF1APP::evalG: Unrecognized difference option\n";
  }
  
  setX(xsave);
  //ngevals++;
  return gx;
}

SymmetricMatrix FDNLF1APP::evalCH(ColumnVector& x) // Evaluate the Hessian
{
  SymmetricMatrix Hessian(dim);

  Hessian = FD2Hessian(x);
  return Hessian;
}

OptppArray<SymmetricMatrix> FDNLF1APP::evalCH( ColumnVector& x, int darg) // Evaluate the Hessian
{
  SymmetricMatrix Hessian(dim);

  Hessian = FD2Hessian(x);
  OptppArray<SymmetricMatrix> H(1);
  H[0] = Hessian;
  return H;
}

} // namespace OPTPP
