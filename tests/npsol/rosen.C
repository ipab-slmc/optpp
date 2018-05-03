#include "NLP2.h"

using NEWMAT::ColumnVector;
using NEWMAT::SymmetricMatrix;
using std::cout;

using namespace OPTPP;

extern double alpha;
static int fcn_count = 0;
/* Example file to demonstrate the calling sequence to a 
 * simple NLF1 function
 */
void init_rosen (int ndim, ColumnVector& x)
{
  if (ndim != 2)
  {
    exit (1);
  }
  x(1) = -1.2;
  x(2) =  1.0;
}

void rosen(int mode, int n, const ColumnVector& x, double& fx, ColumnVector& g, int& result)
{ // Rosenbrock's function
  double f1, f2, x1, x2;
  
  if (n != 2) return;

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

  if (n != 2) return;

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
  fcn_count++;

  double f1, f2, x1, x2;

   cout << "\nrosen2: mode = " << mode
        << " count = " << fcn_count << "\n";
   for(int i=1; i<=n; i++) 
       cout << "x(" << i << ") = " << x(i) << "\n";

  if (n != 2) return;

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
