//------------------------------------------------------------------------
// Copyright (C) 1993: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//------------------------------------------------------------------------

#include "NLP1.h"
#include "TOLS.h"
#include "cblas.h"
#include "precisio.h"
#include "ioformat.h"

using namespace std;
using NEWMAT::ColumnVector;
using NEWMAT::Matrix;
using NEWMAT::SymmetricMatrix;
using NEWMAT::Real;
using NEWMAT::FloatingPointPrecision;

//------------------------------------------------------------------------
// external subroutines referenced by this module 
// Included to eliminate compilation errors when -ansi flag is used
//------------------------------------------------------------------------

#if !(defined(__GNUC__) && __GNUC__ >= 3)
extern "C" {
  double copysign(double, double);
}
#else
#ifdef CYGWIN
extern "C" {
  extern double copysign _PARAMS((double, double));
}
#endif
#endif


namespace OPTPP {

//----------------------------------------------------------------------------
// Evaluate the Hessian using finite differences
// Assume that analytical gradients are available 
//----------------------------------------------------------------------------

SymmetricMatrix NLP1::FDHessian(ColumnVector& sx) 
{
//  Tracer trace("NLP1::FDHessian");
  Real mcheps = FloatingPointPrecision::Epsilon();
  ColumnVector fcn_accrcy = getFcnAccrcy();

  int i;
  double hi, hieps;
  double xtmp;

  int nr = getDim();

  ColumnVector gx(nr), grad(nr), xc(nr);
  Matrix Htmp(nr,nr);
  SymmetricMatrix H(nr);
		     
  xc = getXc();
  gx = getGrad();


  for (i=1; i<=nr; i++) {

    hieps = sqrt(max(mcheps,fcn_accrcy(i) ));
    hi = hieps*max(fabs(xc(i)),sx(i));
    hi = copysign(hi,xc(i));
    xtmp = xc(i);
    xc(i) = xtmp + hi;
    grad = evalG(xc);
    Htmp.Column(i) << (grad - gx) / hi;
    xc(i) = xtmp;
 }

 H << (Htmp.t() + Htmp)/2.0;
 return H;
}

OptppArray<SymmetricMatrix> NLP1::CONFDHessian(ColumnVector& sx) 
{
//  Tracer trace("NLP1::FDHessian");
  Real mcheps = FloatingPointPrecision::Epsilon();
  ColumnVector fcn_accrcy = getFcnAccrcy();

  int i, counter;
  double hi, hieps;
  double xtmp;

  int nr = getDim();

  ColumnVector xc(nr);
  Matrix grad(nr, ncnln), gx(nr, ncnln), Htmp(nr,nr);
  SymmetricMatrix H(nr);

  OptppArray<SymmetricMatrix> Hessian(ncnln);
		     
  xc = getXc();
  gx = evalCG(xc);

  for (counter=0; counter< ncnln; counter++) {

    for (i=1; i<=nr; i++) {

      hieps = sqrt(max(mcheps,fcn_accrcy(i) ));
      hi = hieps*max(fabs(xc(i)),sx(i));
      hi = copysign(hi,xc(i));
      xtmp = xc(i);
      xc(i) = xtmp + hi;
      grad = evalCG(xc);
      Htmp.Column(i) << (grad - gx) / hi;
      xc(i) = xtmp;
   }

   H << (Htmp.t() + Htmp)/2.0;

   Hessian[counter] = H;

  }
 return Hessian;
}

//-------------------------------------------------------------------------
// Output Routines
//-------------------------------------------------------------------------

void NLP1::printState(char * s) 
{ // Print out current state: x current, gradient and Function value
  cout << "\n\n=========  " << s << "  ===========\n\n";
  cout << "\n    i\t    xc \t\t grad  \t\t fcn_accrcy \n";
  for (int i=1; i<=dim; i++) 
    cout << d(i,6) << e(mem_xc(i),12,4)<< "\t" << e(mem_grad(i),12,4) << "\t"
         << e(mem_fcn_accrcy(i),12,4) << "\n";
  cout <<"Function Value     = " << e(fvalue,12,4) << "\n";
  double gnorm = Norm2(mem_grad);
  cout <<"Norm of gradient   = " << e(gnorm,12,4) << "\n";
  cout <<"\n\n==============================================\n\n";
}

void NLP1::fPrintState(ostream *nlpout, char * s) 
{ // Print out current state: x current, gradient and Function value
  (*nlpout) << "\n\n=========  " << s << "  ===========\n\n";
  (*nlpout) << "\n    i\t    xc \t\t grad  \t\t fcn_accrcy \n";
  for (int i=1; i<=dim; i++) 
    (*nlpout) << d(i,6) << e(mem_xc(i),12,4)<< "\t" 
              << e(mem_grad(i),12,4) << "\t"
              << e(mem_fcn_accrcy(i),12,4) << "\n";
  (*nlpout) <<"Function Value     = " << e(fvalue,12,4) << "\n";
  double gnorm = Norm2(mem_grad);
  (*nlpout) <<"Norm of gradient   = " << e(gnorm,12,4) << "\n";
//  (*nlpout) <<"Function Accuracy  = " << e(mem_fcn_accrcy,12,4) << "\n";
  (*nlpout) <<"\n\n==============================================\n\n";
}

} // namespace OPTPP
