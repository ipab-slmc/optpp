/** \example sixtyfive.C
 * <em> minimize </em>
 * \f[ (x_1 - x_2)^2 + (1/9)(x_1 + x_2 - 10)^2 + (x_3 - 5)^2 \f]
 * <em> subject to </em> \f[ x_1^2 + x_2^2 + x_3^2 \le 48, \f]
 * <em> </em>  \f[-4.5 \le x_1 \le 4.5, \f]
 * <em> </em>  \f[-4.5 \le x_2 \le 4.5, \f]
 * <em> </em>  \f[ -5.0 \le x_3 \le 5.0 \f]
 */

//------------------------------------------------------------------------
// P.J. Williams
// Sandia National Laboratories
// pwillia@sandia.gov
//------------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#if defined(SGI) || defined(RS6K)
#define WANT_MATH
#else
#define WANT_STREAM
#define WANT_MATH
#endif

#ifdef HAVE_STD
#include <cstdio>
#include <cstdlib>
#include <cmath>
#else
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#endif

#include "newmatap.h"

#include "NLP.h"
#include "NLF.h"

using NEWMAT::ColumnVector;
using NEWMAT::Matrix;

using namespace OPTPP;

void init_hs65(int ndim, ColumnVector& x)
{
  if (ndim != 3)
  {
    exit (1);
  }
  double factor = 0.0;
  x(1) = -5.0 ;// - (factor - 1)*8.6505  ;
  x(2) =  5.0 ;// + (factor - 1)*1.3495  ;
  x(3) =  0.0 ;// - (factor - 1)*4.6204  ;
}

void hs65(int mode, int n, const ColumnVector& x, double& fx, ColumnVector& g, int& result)
{ // Hock and Schittkowski's Problem 65 (the objective fcn)
  double f1, f2, f3, x1, x2, x3;

  if (n != 3) return;

  x1 = x(1);
  x2 = x(2);
  x3 = x(3);
  f1 = x1 - x2;
  f2 = x1 + x2 - 10.0;
  f3 = x3 - 5.0;

  if (mode & NLPFunction) {
    fx  = f1*f1+ (f2*f2)*(1.0/9.0) +f3*f3;
    result = NLPFunction;
  }
  if (mode & NLPGradient) {
    g(1) =  2*f1 + (2.0/9.0)*f2;
    g(2) = -2*f1 + (2.0/9.0)*f2;
    g(3) =  2*f3;
    result = NLPGradient;
  }
}

void ineq_hs65(int mode, int n, const ColumnVector& x, ColumnVector& cfx, Matrix& g, int& result)
{ // Hock and Schittkowski's Problem 65 
      double f1, f2, f3, fx, x1, x2, x3;

  if (n != 3) return;

  x1 = x(1);
  x2 = x(2);
  x3 = x(3);
  f1 = x1;
  f2 = x2;
  f3 = x3;

  if (mode & NLPConstraint) {
    fx  = f1*f1 + f2*f2 + f3*f3;
    cfx = fx; 
    result = NLPConstraint;
  }

  if (mode & NLPCJacobian) {
    g(1,1) = 2*x1;
    g(1,2) = 2*x2;
    g(1,3) = 2*x3;
    result = NLPCJacobian;
  }
}
