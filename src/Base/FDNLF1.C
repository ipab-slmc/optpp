//------------------------------------------------------------------------
// Copyright (C) 1993: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//------------------------------------------------------------------------

#include "NLF.h"
#include "TOLS.h"
#include "ioformat.h"
#include "cblas.h"

using namespace std;
using NEWMAT::ColumnVector;
using NEWMAT::Matrix;
using NEWMAT::SymmetricMatrix;


namespace OPTPP {

//-------------------------------------------------------------------------
// FDNLF1 method routines.
//-------------------------------------------------------------------------
void FDNLF1::reset() // Reset parameter values  
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

void FDNLF1::initFcn() // Initialize Function
{
  if (init_flag == false) {
      init_fcn(dim, mem_xc);
      init_flag = true;
  }
  else  {
    //cerr << "FDNLF1:initFcn: Warning - initialization called twice\n";
    init_fcn(dim, mem_xc);
  }
}

void FDNLF1::eval() // Evaluate Function and Gradient
{
  (void)evalF();
  (void)evalG();
}

double FDNLF1::evalF() // Evaluate Function
{
  int result = 0;
  double time0 = get_wall_clock_time();

  if (SpecFlag == NoSpec) {
    if (!application.getF(mem_xc, fvalue)) {
      fcn_v(dim, mem_xc, fvalue, result, vptr);
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

double FDNLF1::evalF(const ColumnVector& x) // Evaluate Function at x
{
  double fx;
  int result = 0;
  double time0 = get_wall_clock_time();

  if (SpecFlag == NoSpec) {
    if (!application.getF(x, fx)) {
      fcn_v(dim, x, fx, result, vptr);
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

ColumnVector FDNLF1::evalG() // Evaluate the gradient
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
//    cout << "FDNLF1::evalG: Unrecognized difference option\n";
//    cout << "FDNLF1::evalG: Using forward difference option\n";
    mem_grad =  FDGrad(sx, mem_xc, fvalue, partial_grad);
  }
 return mem_grad;
}

ColumnVector FDNLF1::evalG(const ColumnVector& x) // Evaluate the gradient at x
{
  ColumnVector gx(dim);
  ColumnVector sx(dim);
  sx = 1.0;
  ngevals++;

  if (SpecFlag == NoSpec) {
    int result = 0;
    if (!application.getF(x, specF)) {
      fcn_v(dim, x, specF, result, vptr);
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
//    cout << "FDNLF1::evalG: Unrecognized difference option\n";
//    cout << "FDNLF1::evalG: Using forward difference option\n";
    mem_grad =  FDGrad(sx, x, specF, partial_grad);
  }
  
  return gx;
}

SymmetricMatrix FDNLF1::FDHessian(ColumnVector& sx) // Evaluate the Hessian
{
  SymmetricMatrix Hessian(dim);

  sx = 1.0;
  Hessian = FD2Hessian(sx);
  return Hessian;
}

SymmetricMatrix FDNLF1::evalH() // Evaluate the Hessian
{
  ColumnVector sx(dim);
  SymmetricMatrix Hessian(dim);

  sx = 1.0;
  Hessian = FD2Hessian(sx);
  return Hessian;
}

SymmetricMatrix FDNLF1::evalH(ColumnVector& x) // Evaluate the Hessian
{
  SymmetricMatrix Hessian(dim);

  Hessian = FD2Hessian(x);
  return Hessian;
}

void FDNLF1::printState(char * s) 
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

void FDNLF1::fPrintState(ostream *nlpout, char * s) 
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


real FDNLF1::evalLagrangian(const ColumnVector& xc , 
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

ColumnVector FDNLF1::evalLagrangianGradient(const ColumnVector& xc, 
                                          const ColumnVector& multiplier,
                                          const ColumnVector& type) 
{
   mem_grad  = evalG(xc);
   ColumnVector grad  = mem_grad;
   if(hasConstraints())
      grad -= constraint_->evalGradient(xc)*multiplier;
   return grad;
}

ColumnVector FDNLF1::evalCF(const ColumnVector& x) // Evaluate Constraint Fcn at x
{
  int result = 0;
  ColumnVector cfx(ncnln);
  double time0 = get_wall_clock_time();
  confcn(dim, x, cfx,result);
  function_time = get_wall_clock_time() - time0;

  //nfevals++;
  return cfx;
}

Matrix FDNLF1::evalCG(const ColumnVector& x) // Evaluate the gradient at x
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
//    cout <<"FDNLF1::evalG: Unrecognized difference option\n";
  }
  
  setX(xsave);
  //ngevals++;
  return gx;
}

SymmetricMatrix FDNLF1::evalCH(ColumnVector& x) // Evaluate the Hessian
{
  SymmetricMatrix Hessian(dim);

  Hessian = FD2Hessian(x);
  return Hessian;
}

OptppArray<SymmetricMatrix> FDNLF1::evalCH( ColumnVector& x, int darg) // Evaluate the Hessian
{
  SymmetricMatrix Hessian(dim);

  Hessian = FD2Hessian(x);
  OptppArray<SymmetricMatrix> H(1);
  H[0] = Hessian;
  return H;
}

void FDNLF1::evalC(const ColumnVector& x) // Evaluate Function and Gradient
{
  (void) evalCF(x);
  (void) evalCG(x);
}

} // namespace OPTPP
