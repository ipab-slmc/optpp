
#include <iostream>
#include "NLP.h"

using NEWMAT::ColumnVector;
using NEWMAT::SymmetricMatrix;

using namespace OPTPP;

/* Example file to demonstrate the calling sequence to a 
 * simple NLF1 function
 */

extern "C" {
void init_rosen (int ndim, ColumnVector& x)
{
  if (ndim != 2) {
    cerr << "Number of variables for Rosenbrock's function should be 2."
	 << "  The number of variables given is " << ndim << endl;
    exit (1);
  }
}
void rosen(int mode, int n, const ColumnVector& x, double& fx, ColumnVector& g, int& result)
{ // Rosenbrock's function
  double f1, f2, x1, x2;

  x1 = x(1);
  x2 = x(2);
  f1 = (x2 - x1 * x1);
  f2 = 1. - x1;
  
  if (mode & NLPFunction) {
    fx  = 100.* f1*f1 + f2*f2;
    result = NLPFunction;
  }
  if (mode & NLPGradient) {
    g(1) = -400.*f1*x1 - 2.*f2; 
    g(2) = 200.*f1;
    result = NLPGradient;
  }
}
void rosen0(int n, const ColumnVector& x, double& fx, int& result)
{ // Rosenbrock's function
  double f1, f2, x1, x2;

  x1 = x(1);
  x2 = x(2);
  f1 = (x2 - x1 * x1);
  f2 = 1. - x1;
  
  fx  = 100.* f1*f1 + f2*f2;
  result = NLPFunction;
}
void rosen2(int mode, int n, const ColumnVector& x, double& fx, 
	ColumnVector& g, SymmetricMatrix& H, int &result)
// Rosenbrock's function, n = 2 with first and second order derivatives

{ 
  double f1, f2, x1, x2;

  x1 = x(1);
  x2 = x(2);
  f1 = (x2 - x1 * x1);
  f2 = 1. - x1;
  
  if (mode & NLPFunction) {
    fx  = 100.* f1*f1 + f2*f2;
    result = NLPFunction;
  }
  if (mode & NLPGradient) {
    g(1) = -400.*f1*x1 - 2.*f2; 
    g(2) = 200.*f1;
    result = NLPGradient;
  }
  
  if (mode & NLPHessian) {
    f1 = (x2 - 3.0*x1*x1);
    H(1,1) = -400.0*f1 + 2.0;
    H(2,1) = -400.0*x1;
    H(2,2) = 200.0;
    result = NLPHessian;
  }
}
}
