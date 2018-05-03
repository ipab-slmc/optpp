//------------------------------------------------------------------------
// P.J. Williams
// Sandia National Laboratories
// pwillia@sandia.gov
//------------------------------------------------------------------------

#include <iostream>
#include "NLP.h"

using NEWMAT::ColumnVector;
using NEWMAT::Matrix;
using NEWMAT::SymmetricMatrix;
using namespace OPTPP;

extern "C" {
void init_hs14(int ndim, ColumnVector& x)
{
  if (ndim != 2)
  {
    cerr << "Number of variables for Hock Problem 14 should be 2."
	 << "  The number of variables given is " << ndim << endl;
    exit (1);
  }

  x(1) = 2+(x(1) -1)* 1.1771243447;
  x(2) = 2+(x(2) -1)* 1.0885621722;
}
void hs14(int mode, int n, const ColumnVector& x, double& fx, ColumnVector& g, int& result)
{ // Hock and Schittkowski's Problem 14
  double f1, f2, x1, x2;

  x1 = x(1);
  x2 = x(2);
  f1 = x1 - 2.0;
  f2 = x2 - 1.0;

  if (mode & NLPFunction) {
    fx  = f1*f1 + f2*f2;
    result = NLPFunction;
  }
  if (mode & NLPGradient) {
    g(1) = 2*(x1 - 2.0);
    g(2) = 2*(x2 - 1.0);
    result = NLPGradient;
  }
}

void ineq_hs14(int mode, int n, const ColumnVector& x, ColumnVector& fx, Matrix& g, int& result)
{ // Hock and Schittkowski's Problem 14
  double f1, f2, x1, x2;

  x1 = x(1);
  x2 = x(2);
  f1 = x1*x1;
  f2 = x2*x2;

  if (mode & NLPFunction) {
    fx  = -.25*f1 - f2 + 1.0;
    result = NLPFunction;
  }
  if (mode & NLPGradient) {
    g(1,1) = -0.5*x1;
    g(2,1) = -2.0*x2;
    result = NLPGradient;
  }
}

void init_hs65(int ndim, ColumnVector& x)
{
  if (ndim != 3)
  {
    cerr << "Number of variables for Hock Problem 65 should be 3."
	 << "  The number of variables given is " << ndim << endl;
    exit (1);
  }

  x(1) = -5.0  - (x(1) - 1)*8.6505  ;
  x(2) =  5.0  + (x(2) - 1)*1.3495  ;
  x(3) =  0.0  - (x(3) - 1)*4.6204  ;
}
void hs65_2(int mode, int n, const ColumnVector& x, double& fx, ColumnVector& g, SymmetricMatrix& H, int& result)
{ // Hock and Schittkowski's Problem 65 (the objective fcn)
  double f1, f2, f3, x1, x2, x3;

  x1 = x(1);
  x2 = x(2);
  x3 = x(3);
  f1 = x1 - x2;
  f2 = x1 + x2 - 10.0;
  f3 = x3 - 5.0;

  if (mode & NLPFunction) {
    fx  = f1*f1+ (f2*f2)/9.0 +f3*f3;
    result = NLPFunction;
  }
  if (mode & NLPGradient) {
    g(1) =  2*f1 + (2.0/9.0)*f2;
    g(2) = -2*f1 + (2.0/9.0)*f2;
    g(3) =  2*f3;
    result = NLPGradient;
  }
  if (mode & NLPHessian) {
    H(1,1) =  2 + (2.0/9.0);

    H(2,1) = -2 + (2.0/9.0);
    H(2,2) =  2 + (2.0/9.0);

    H(3,1) = 0.0;
    H(3,2) = 0.0;
    H(3,3) = 2.0;
    result = NLPHessian;
  }
}

void ineq_hs65_2(int mode, int n, const ColumnVector& x, ColumnVector& fx, Matrix& g, OptppArray<SymmetricMatrix>& H, int& result)
{ // Hock and Schittkowski's Problem 65 
  double f1, f2, f3, x1, x2, x3;
  SymmetricMatrix Htmp(n);

  x1 = x(1);
  x2 = x(2);
  x3 = x(3);
  f1 = x1;
  f2 = x2;
  f3 = x3;

  if (mode & NLPFunction) {
    fx(1)  = 48 - f1*f1 - f2*f2 - f3*f3;
    result = NLPFunction;
  }
  if (mode & NLPGradient) {
    g(1,1) = -2*x1;
    g(2,1) = -2*x2;
    g(3,1) = -2*x3;
    result = NLPGradient;
  }
  if (mode & NLPHessian) {
    Htmp(1,1) = -2;
    Htmp(1,2) = 0.0;
    Htmp(1,3) = 0.0;
    Htmp(2,1) = 0.0;
    Htmp(2,2) = -2;
    Htmp(2,3) = 0.0;
    Htmp(3,1) = 0.0;
    Htmp(3,2) = 0.0;
    Htmp(3,3) = -2;

    H[0] = Htmp;
    result = NLPHessian;
  }
}
}
