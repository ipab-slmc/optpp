//------------------------------------------------------------------------
// Copyright (C) 1993, 1994: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//------------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
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

#if (defined(__sgi) || defined(__xlc__) || defined (__xlC__))
#define WANT_MATH
#else
#define WANT_STREAM
#define WANT_MATH
#endif

#include "newmatap.h"

#include "NLP.h"
#include "NLF.h"

using std::cerr;

using NEWMAT::ColumnVector;
using NEWMAT::SymmetricMatrix;

using namespace OPTPP;

static double square(double x) { return x*x; }

void init_rosen (int ndim, ColumnVector& x)
{
  if (ndim != 2)
  {
    exit (1);
  }
  x(1) = -1.2;
  x(2) =  1.0;
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

void rosen2(int mode, int n, const ColumnVector& x, double& fx, 
	ColumnVector& g, SymmetricMatrix& H, int &result)
// Rosenbrock's function, n = 2 with first and second order derivatives

{ 
  //fcn_count++;

  double f1, f2, x1, x2;


  // cout << "\nrosen2: mode = " << mode
  //      << " count = " << fcn_count << "\n";
  // for(int i=1; i<=n; i++) 
  //     cout << "x(" << i << ") = " << x(i) << "\n";

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

void erosen1(int mode, int ndim, const ColumnVector& x, double& fx, 
	ColumnVector& g, int &result)
// Extended Rosenbrock's function, with analytic  derivatives

{ 
  static int i;
  static double f1, f2, x1, x2;

  fx = 0.;

  for (i = 1; i <= ndim/2; ++i) {
    x1 = x(2*i-1);
    x2 = x(2*i);
    
    f1 = (x2 - x1 * x1);
    f2 = 1. - x1;
    
    fx += 100.*square(f1) + square(f2);

    if (mode & NLPGradient) {
      g(2*i-1) = -400.*f1*x1 - 2.*f2; 
      g(2*i) = 200.*f1;
    }
  }
  result = NLPFunction;
  if (mode & NLPGradient) result = NLPGradient;

}
void init_epowell (int ndim, ColumnVector& x)
{
  if (ndim % 4 != 0) 
  {
    printf ("init_epowell: ndim must be multiple of 4, ndim = %d\n", ndim);
    exit (1);
  }

  for (int i=1; i<=ndim/4; i++)
  {
    x(4*i - 3) =  3;
    x(4*i - 2) = -1;
    x(4*i - 1) =  0;
    x(4*i    ) =  1;
  }
}
void epowell(int ndim, const ColumnVector& x, double& fx, int& result)

{
  // Extended Powell Singular Function, n variable, multiple of 4
  // Standard initial guess  x0   = (psi1, psi2, psi3, psi4, ...)
  //                         psi1 =  3
  //                         psi2 = -1
  //                         psi3 =  0
  //                         psi4 =  1
  // Solution                fmin = 0 at origin

  double f1, f2, f3, f4;

  if (ndim % 4 != 0) {
    cerr << "epowell: ndim must be a multiple of 4, ndim = \n" << ndim;
    return;
  }

  fx = 0.0;
  for (int i=1; i<=ndim/4; i++) {
    f1 = x(4*i-3) + 10*x(4*i-2);
    f2 = sqrt(5.0) * ( x(4*i-1) - x(4*i));
    f3 = square(x(4*i-2) - 2*x(4*i-1));
    f4 = sqrt(10.0) * (x(4*i-3) - x(4*i));
    fx += square(f1) + square(f2) + square(f3) + square(f4);
  }  
  result = NLPFunction;
}

void init_trig (int ndim, ColumnVector& x)
{
  for (int i=1; i<=ndim; i++)
  {
    x(i) =  1/ double(ndim);
  }
}

void trig(int ndim, const ColumnVector& x, double& fx, int& result)

{ 
  // Trigonometric Function, ndim variable, m = ndim
  // Standard initial guess  x0   = (1/ndim, ..., 1/ndim)
  // Solution                fmin = 0

  double fi;

  fx = 0.0;

  for (int i=1; i<=ndim; ++i) {
    double sum = 0.0;
    for (int j=1; j<=ndim; ++j) {
      sum += cos(x(j)) + double(i)*(1.0-cos(x(i)) - sin(x(i)));
    }
    fi = ndim - sum;
    fx += fi*fi;
  }
  
  result = NLPFunction;
}
void init_erosen (int ndim, ColumnVector& x)
{

  if (ndim % 2 != 0) 
  {
    printf ("init_erosen: ndim must be even, ndim = %d\n", ndim);
    exit (1);
  }

  for (int i = 1; i <= ndim / 2; ++i)
  {
    x(2*i - 1) = -1.2;
    x(2*i)     =  1.0;
  }

}
void erosen(int ndim, const ColumnVector& x, double& fx, int& result)
{
//  erosen
//  Extended rosenbrock function for any even value of n. 
//                          2 
//     f    (x) = 10(x   - x    ) 
//      2i-1          2i    2i-1 
//
//     f  (x) = 1 - x 
//      2i           2i-1 
//                                           2         2 
//     f(x) = sum from i = 1 to n/2 (f    (x)  + f  (x) ). 
//                                    2i-1        2i 
//
//     x  = (epsilon ) where epsilon    = -1.2, epsilon  = 1 
//      0           i               2i-1               2i 
//
//     f(x ) = 0 at x  = (1,...,1) 
//        *          * 
//
//     parameters      
//
//     input 
//        n             dimension of the problem to be solved must be even
//        x             point at which the function is to be evaluated 
//     output 
//        f             function value at x 


  static int i;
  static double f1, f2, x1, x2;

  fx = 0.;

  for (i = 1; i <= ndim/2; ++i) {
    x1 = x(2*i-1);
    x2 = x(2*i);
    
    f1 = (x2 - x1 * x1);
    f2 = 1. - x1;
    
    fx += 100.*square(f1) + square(f2);
  }
  result = NLPFunction;
}

void init_penalty1 (int ndim, ColumnVector& x)
{

  for (int i = 1; i <= ndim; i++) x(i) = double(i);

}
void penalty1(int ndim, const ColumnVector& x, double& fx, int& result)
{
//  Penalty 1
//     parameters      
//
//     input 
//        ndim          dimension of the problem to be solved must be even
//        x             point at which the function is to be evaluated 
//     output 
//        f             function value at x 


  int i;
  double fi;
  double k1  = 1.e-5;
  double sum = 0.0;

  fx = 0.;

  for (i=1; i <= ndim; i++) sum += square(x(i));
  fx = (sum - .25)*(sum - .25);
  

  for (i = 1; i <= ndim; i++) {
    fi = (x(i) - 1.0);
    fx += k1*square(fi);
  }
  result = NLPFunction;
}
void init_penalty2 (int ndim, ColumnVector& x)
{

  for (int i = 1; i <= ndim; i++) x(i) = 0.5;

}
void penalty2(int ndim, const ColumnVector& x, double& fx, int& result)
{
//  Penalty 2
//     parameters      
//
//     input 
//        ndim          dimension of the problem to be solved must be even
//        x             point at which the function is to be evaluated 
//     output 
//        f             function value at x 


  int i;
  double f1, fi, yi;

  double k1  = 1.e-5;
  double sum = 0.0;

  fx = 0.0;

  // Take care of the two end cases first
  //

  f1  = x(1) - 0.2;
  fx  = square(f1);
  for (i=1; i <= ndim; i++) sum += double(ndim - i + 1)*square(x(i));
  sum = sum - 1.0;
  fx += square(sum);
  
  //
  // Now take care of all of the cases in between
  //
  for (i = 2; i <= ndim; i++) {
    yi = exp(double(i)/10.0) + exp(double(i-1)/10.0);
    fi = exp(x(i)/10.0) + exp(x(i-1)/10.0) - yi;
    fx += k1*square(fi);
  }
  for (i = ndim+1; i < 2*ndim; i++) {
    fi = exp(x(i-ndim+1)/10.0) - exp(double(-1)/10.0);
    fx += k1*square(fi);
  }
  result = NLPFunction;
}
void init_vardim (int ndim, ColumnVector& x)
{

  for (int i = 1; i <= ndim; i++) x(i) = 1.0 - (double(i)/double(ndim));

}
void vardim(int ndim, const ColumnVector& x, double& fx, int& result)
{
//  Variable dimensioned
//     parameters      
//
//     input 
//        ndim          dimension of the problem to be solved must be even
//        x             point at which the function is to be evaluated 
//     output 
//        f             function value at x 


  int i;
  double fn_plus_1, fn_plus_2, fi;

  double sum = 0.0;

  fx = 0.;

  // Take care of the two end cases first
  //
  for (i=1; i <= ndim; i++) sum += double(i)*(x(i) - 1.0);

  fn_plus_1  = sum;
  fn_plus_2  = square(sum);
  fx        += square(fn_plus_1) + square(fn_plus_2);
  
  //
  // Now take care of all of the cases in between
  //
  for (i = 1; i <= ndim; i++) {
    fi = x(i) - 1.0;
    fx += square(fi);
  }
  result = NLPFunction;
}
void init_broyden_tridiag (int ndim, ColumnVector& x)
{

  for (int i = 1; i <= ndim; i++) x(i)   = -1.0;

}
void broyden_tridiag(int ndim, const ColumnVector& x, double& fx, int& result)
{
//  Broyden Tridiagonal function
//     parameters      
//
//     input 
//        ndim          dimension of the problem to be solved must be even
//        x             point at which the function is to be evaluated 
//     output 
//        f             function value at x 


  int i;
  double f1, fi, fn;
  double x0 = 0.0, xn_plus_1 = 0.0;

  fx = 0.;

  f1 = (3.0 - 2.0*x(1))*x(1) - x0 - 2.0*x(2) + 1.0;
  fn = (3.0 - 2.0*x(ndim))*x(ndim) - x(ndim-1) - 2.0*xn_plus_1 + 1.0;

  fx += square(f1) + square(fn);
  //
  // Now take care of all of the cases in between
  //
  for (i = 2; i < ndim; i++) {
    fi  = (3.0 - 2.0*x(i))*x(i) - x(i-1) - 2.0*x(i+1) + 1.0;
    fx += square(fi);
  }
  result = NLPFunction;
}
//-----------------------------------------------------------------------------
// Problems from Conn, Gould, Toint
//
//-----------------------------------------------------------------------------
void init_gen_brown (int ndim, ColumnVector& x)
{

  for (int i = 1; i < ndim; i=i+2){
    x(i)   = -1.0;
    x(i+1) =  1.0;
  }

}
void gen_brown(int ndim, const ColumnVector& x, double& fx, int& result)
{
//  Generalization of a function due to A. Brown
//  Conn, Gould, Toint, Test function 21, pp. 427
//     parameters      
//
//     input 
//        ndim          dimension of the problem to be solved must be even
//        x             point at which the function is to be evaluated 
//     output 
//        f             function value at x 


  int i;
  double x0, x1, y1, y2;

  fx = 0.;

  //
  for (i = 1; i < ndim; i++) {
    x0 = x(i);
    x1 = x(i+1);
    y1 = square(x1) + 1;
    y2 = square(x0) + 1;
    fx += pow(square(x0),y1) + pow(square(x1),y2);
  }
  result = NLPFunction;
}
void init_chain_singular (int ndim, ColumnVector& x)
{

  if (ndim % 4 != 0) 
  {
    printf ("init_chain_singular: ndim must be multiple of 4, ndim = %d\n", ndim);
    exit (1);
  }

  for (int i=1; i<=ndim/4; i++)
  {
    x(4*i - 3) =  3;
    x(4*i - 2) = -1;
    x(4*i - 1) =  0;
    x(4*i    ) =  1;
  }

}
void chain_singular(int ndim, const ColumnVector& x, double& fx, int& result)
{
//  Chained Singular Function
//  Conn, Gould, Toint, Test function 21, pp. 422
//     parameters      
//
//     input 
//        ndim          dimension of the problem to be solved must be even
//        x             point at which the function is to be evaluated 
//     output 
//        f             function value at x 


  int i;
  double x0, x1, x2, x3;
  double t1, t2, t3, t4;

  fx = 0.;

  if (ndim % 4 != 0) 
  {
    printf ("chain_singular: ndim must be multiple of 4, ndim = %d\n", ndim);
    exit (1);
  }

  //
  for (i = 1; i <= ndim-3; i=i+2) {
    x0 = x(i);
    x1 = x(i+1);
    x2 = x(i+2);
    x3 = x(i+3);
    t1 = square(x0 + 10.0*x1);
    t2 = 5.0*square(x2 - x3);
    t3 = square(x1 - 2.0*x2)*square(x1 - 2.0*x2);
    t4 = 10.0*square(x0 - 10.0*x3)*square(x0 - 10.0*x3);
    fx += t1 + t2 + t3 + t4;
  }
  result = NLPFunction;
}
void init_gen_wood (int ndim, ColumnVector& x)
{

  if (ndim % 4 != 0) 
  {
    printf ("init_gen_wood: ndim must be multiple of 4, ndim = %d\n", ndim);
    exit (1);
  }

  x(1) =  -3;
  x(2) =  -1;
  x(3) =  -3;
  x(4) =  -1;

  for (int i=5; i<=ndim; i=i+2)
  {
    x(i)   = -2;
    x(i+1) =  0;
  }
  //  for (int i=1; i<=ndim; i++) x(i) = 1.1;
}
void gen_wood(int ndim, const ColumnVector& x, double& fx, int& result)
{
//  Chained Singular Function
//  Conn, Gould, Toint, Test function 21, pp. 422
//     parameters      
//
//     input 
//        ndim          dimension of the problem to be solved must be even
//        x             point at which the function is to be evaluated 
//     output 
//        f             function value at x 


  int i;
  double x0, x1, x2, x3;
  double t1, t2, t3, t4, t5, t6;

  fx = 0.0;
  if (ndim % 4 != 0) 
  {
    printf ("gen_wood: ndim must be multiple of 4, ndim = %d\n", ndim);
    exit (1);
  }

  //
  fx = 1.0;
  for (i = 1; i <= ndim-3; i=i+4) {
    x0 = x(i);
    x1 = x(i+1);
    x2 = x(i+2);
    x3 = x(i+3);

    t1 = 100.0*square(x1 - square(x0));
    t2 = square(1.0 - x0);
    t3 = 90.0*square(x3 - square(x2));
    t4 = square(1.0 - x2);
    t5 = 10.0*square(x1 + x3 - 2.0);
    t6 =  0.1*square(x1 - x3);
    fx += t1 + t2 + t3 + t4 + t5 + t6;
  }
  result = NLPFunction;
}
void init_chain_wood (int ndim, ColumnVector& x)
{

  if (ndim % 4 != 0) 
  {
    printf ("init_chain_wood: ndim must be multiple of 4, ndim = %d\n", ndim);
    exit (1);
  }

  x(1) =  -3;
  x(2) =  -1;
  x(3) =  -3;
  x(4) =  -1;

  for (int i=5; i<=ndim; i=i+2)
  {
    x(i)   = -2;
    x(i+1) =  0;
  }

}
void chain_wood(int ndim, const ColumnVector& x, double& fx, int& result)
{
//  Chained Singular Function
//  Conn, Gould, Toint, Test function 21, pp. 422
//     parameters      
//
//     input 
//        ndim          dimension of the problem to be solved must be even
//        x             point at which the function is to be evaluated 
//     output 
//        f             function value at x 


  int i;
  double x0, x1, x2, x3;
  double t1, t2, t3, t4, t5, t6;

  fx = 0.0;
  if (ndim % 4 != 0) 
  {
    printf ("chain_wood: ndim must be multiple of 4, ndim = %d\n", ndim);
    exit (1);
  }

  //
  fx = 1.0;
  for (i = 1; i <= ndim-3; i=i+2) {
    x0 = x(i);
    x1 = x(i+1);
    x2 = x(i+2);
    x3 = x(i+3);

    t1 = 100.0*square(x1 - square(x0));
    t2 = square(1.0 - x0);
    t3 = 90.0*square(x3 - square(x2));
    t4 = square(1.0 - x2);
    t5 = 10.0*square(x1 + x3 - 2.0);
    t6 =  0.1*square(x1 - x3);
    fx += t1 + t2 + t3 + t4 + t5 + t6;
  }
  result = NLPFunction;
}
void init_broyden (int ndim, ColumnVector& x)
{

  for (int i = 1; i <= ndim; i++) x(i)   = -1.0;

}
void broyden1a(int ndim, const ColumnVector& x, double& fx, int& result)
{
//  Broyden Tridiagonal function
//     parameters      
//
//     input 
//        ndim          dimension of the problem to be solved must be even
//        x             point at which the function is to be evaluated 
//     output 
//        f             function value at x 


  int i;
  double f1, fi, fn;
  double x0 = 0.0, xn_plus_1 = 0.0;
  double p  = double(7)/double(3);

  fx = 1.0;

  // First the end cases
  //
  f1 = fabs((3.0 - 2.0*x(1))*x(1) - x0 - 2.0*x(2) + 1.0);
  fn = fabs((3.0 - 2.0*x(ndim))*x(ndim) - x(ndim-1) - 2.0*xn_plus_1 + 1.0);

  fx += pow(f1,p) + pow(fn,p);
  //
  // Now take care of all of the cases in between
  //
  for (i = 2; i < ndim; i++) {
    fi  = fabs((3.0 - 2.0*x(i))*x(i) - x(i-1) - 2.0*x(i+1) + 1.0);
    fx += pow(fi,p);
  }
  result = NLPFunction;
}
void broyden1b(int ndim, const ColumnVector& x, double& fx, int& result)
{
//  Broyden Tridiagonal function
//     parameters      
//
//     input 
//        ndim          dimension of the problem to be solved must be even
//        x             point at which the function is to be evaluated 
//     output 
//        f             function value at x 


  int i;
  double f1, fi, fn;
  double x0  = 0.0, xn_plus_1 = 0.0;
  double p   = 2.0;

  fx = 1.0;

  // First the end cases
  //
  f1 = fabs((3.0 - 2.0*x(1))*x(1) - x0 - 2.0*x(2) + 1.0);
  fn = fabs((3.0 - 2.0*x(ndim))*x(ndim) - x(ndim-1) - 2.0*xn_plus_1 + 1.0);

  fx += pow(f1,p) + pow(fn,p);
  //
  // Now take care of all of the cases in between
  //
  for (i = 2; i < ndim; i++) {
    fi  = fabs((3.0 - 2.0*x(i))*x(i) - x(i-1) - 2.0*x(i+1) + 1.0);
    fx += pow(fi,p);
  }
  result = NLPFunction;
}
void broyden2a(int ndim, const ColumnVector& x, double& fx, int& result)
{
//  Generalization of Broyden Banded function, p = 7/3
//     parameters      
//
//     input 
//        ndim          dimension of the problem to be solved
//        x             point at which the function is to be evaluated 
//     output 
//        f             function value at x 


  int i, j;
  double fi;
  double sum;
  double p  = double(7)/double(3);
  int lower_bound, upper_bound;

  fx = 1.0;

  int ml = 5;
  int mu = 1;

  for (i = 1; i < ndim; i++) {

    lower_bound = max(1, i-ml);
    upper_bound = min(ndim, i+mu);

    sum = 0.0;
    for (j = 1; j < ndim; j++) {
      if ( (j >= lower_bound) & (j <= upper_bound))
	sum += x(j)*(1+x(j));
    }

    fi  = fabs((2.0 + 5.0*square(x(i)))*x(i) + 1.0 - sum);
    fx += pow(fi,p);
  }
  result = NLPFunction;
}
void broyden2b(int ndim, const ColumnVector& x, double& fx, int& result)
{
//  Generalization of Broyden Banded function, p = 2
//     parameters      
//
//     input 
//        ndim          dimension of the problem to be solved must be even
//        x             point at which the function is to be evaluated 
//     output 
//        f             function value at x 


  int i, j;
  double fi;
  double sum;
  double p  = 2.0;
  int lower_bound, upper_bound;

  fx = 1.0;

  int ml = 5;
  int mu = 1;

  for (i = 1; i < ndim; i++) {

    lower_bound = max(1, i-ml);
    upper_bound = min(ndim, i+mu);

    sum = 0.0;
    for (j = 1; j < ndim; j++) {
      if ( (j >= lower_bound) & (j <= upper_bound))
	sum += x(j)*(1+x(j));
    }

    fi  = fabs((2.0 + 5.0*square(x(i)))*x(i) + 1.0 - sum);
    fx += pow(fi,p);
  }
  result = NLPFunction;
}
void tointbroy(int ndim, const ColumnVector& x, double& fx, int& result)
{
//  Broyden Tridiagonal function
//     parameters      
//
//     input 
//        ndim          dimension of the problem to be solved must be even
//        x             point at which the function is to be evaluated 
//     output 
//        f             function value at x 


  int i;
  double f1, fi, fn;
  double x0 = 0.0, xn_plus_1 = 0.0;
  double p  = double(7)/double(3);

  if (ndim % 2 != 0) 
  {
    printf ("tointbroy: ndim must be multiple of 2, ndim = %d\n", ndim);
    exit (1);
  }

  fx = 1.0;

  // First the end cases
  //
  f1 = fabs((3.0 - 2.0*x(1))*x(1) - x0 - x(2) + 1.0);
  fn = fabs((3.0 - 2.0*x(ndim))*x(ndim) - x(ndim-1) - xn_plus_1 + 1.0);

  fx += pow(f1,p) + pow(fn,p);
  //
  // Now take care of all of the cases in between
  //
  for (i = 2; i < ndim; i++) {
    fi  = fabs((3.0 - 2.0*x(i))*x(i) - x(i-1) - x(i+1) + 1.0);
    fx += pow(fi,p);
  }
  //
  // Add extra term due to Toint
  //
  int mid = ndim/2;
  for (i = 1; i <= ndim/2; i++) {
    fi  = fabs(x(i) + x(i+mid));
    fx += pow(fi,p);
  }
  result = NLPFunction;
}
void init_gen_cragg_levy (int ndim, ColumnVector& x)
{

  if (ndim % 4 != 0) 
  {
    printf ("init_gen_wood: ndim must be multiple of 4, ndim = %d\n", ndim);
    exit (1);
  }

  x(1) =  1;
  for (int i=2; i<=ndim; i++) x(i)   = 2;
}
void gen_cragg_levy(int ndim, const ColumnVector& x, double& fx, int& result)
{
//  Generalized Cragg and Levy
//  Conn, Gould, Toint, Test function 17, pp. 422
//     parameters      
//
//     input 
//        ndim          dimension of the problem to be solved must be even
//        x             point at which the function is to be evaluated 
//     output 
//        f             function value at x 


  int i;
  double x0, x1, x2, x3;
  double t1, t2, t3, t4, t5;

  fx = 0.0;
  if (ndim % 4 != 0) 
  {
    printf ("gen_cragg_levy: ndim must be multiple of 4, ndim = %d\n", ndim);
    exit (1);
  }

  //
  for (i = 1; i <= ndim-3; i=i+4) {
    x0 = x(i);
    x1 = x(i+1);
    x2 = x(i+2);
    x3 = x(i+3);

    t1 = pow((exp(x0) - x1),4.0);
    t2 = 100.0*pow((x1 - x2),6.0);
    t3 = pow((tan(x2 - x3)),4.0);
    t4 = pow(x0,8.0);
    t5 = square(x3 - 1.0);
    fx += t1 + t2 + t3 + t4 + t5;
  }
  result = NLPFunction;
}
void init_toint_trig (int ndim, ColumnVector& x)
{

  for (int i = 1; i <= ndim; i++) x(i)   = 1.0;

}void toint_trig(int ndim, const ColumnVector& x, double& fx, int& result)
{
//  Toint Trig
//     parameters      
//
//     input 
//        ndim          dimension of the problem to be solved
//        x             point at which the function is to be evaluated 
//     output 
//        f             function value at x 


  int i, j;
  ColumnVector beta(ndim);
  double alpha_ij, c_ij, t;
  fx = 0.0;

  for (i = 1; i <= ndim; i++) {
    beta(i) = 1.0 + double(i)/10.0;
  }

  for (i = 1; i <= ndim; i++) {

    for (j = 1; j <= ndim; j++) {

      if ((abs(i-j))%4 == 0) {
	alpha_ij = 5.0*(1.0 + (i%5) + (j%5));
	c_ij     = double(i+j)/10.0;
	t        = beta(i)*x(i) + beta(j)*x(j) + c_ij;
	fx      += alpha_ij*sin(t);
      }
    }
  }
  result = NLPFunction;
}
void init_chebyquad (int ndim, ColumnVector& x)
{

  for (int i = 1; i <= ndim; i++) x(i)   = double(i)/double(ndim+1);

}void chebyquad(int ndim, const ColumnVector& x, double& fx, int& result)
{
//  Chebyquad
//     parameters      
//
//     input 
//        ndim          dimension of the problem to be solved
//        x             point at which the function is to be evaluated 
//     output 
//        f             function value at x 


  int i, i2, j;
  int m2   = ndim/2;
  int mdim = ndim;

  ColumnVector fi(ndim);
  double coeff = 1.0/double(ndim);
  double t, t0, t1, tj;

  fx = 0.0;
  fi = 0.0;

  for (i = 1; i <= m2; i++) {
    i2 = 2*i;
    fi(i2) = 1.0 / double((i2*i2) - 1.0);
  }

  for (j = 1; j <= ndim; j++) {
    t0 = 1.0;
    t1 = 2.0 * x(j) - 1.0;
    tj = t1;
    fi(1) = fi(1) + coeff*t1;

    for (i = 2; i <= mdim; i++) {
      t = 2.0 * tj * t1 - t0;
      fi(i) = fi(i) + coeff*t;
      t0 = t1;
      t1 = t;
    }
  }
  for (i = 1; i <= ndim; i++) {
    fx += fi(i)*fi(i);
  }

  result = NLPFunction;
}
void init_nelder (int ndim, ColumnVector& x)
{

  for (int i = 1; i <= ndim; i++) x(i)   = double(i)/double(ndim+1);

}
void nelder(int ndim, const ColumnVector& x, double& fx, int& result)
{
//  Nelder
//     parameters      
//
//     input 
//        ndim          dimension of the problem to be solved
//        x             point at which the function is to be evaluated 
//     output 
//        f             function value at x 

  double f1, f2, f3;

  if (ndim != 2)
  {
    exit (1);
  }

  f1 = x(4) - x(1)*x(2)*x(3);
  f2 = x(3) - x(1)*x(2);
  f3 = x(2) - x(1);

  fx = f1*f1 + f2*f2 + f3*f3 + x(1)*x(1);

  result = NLPFunction;
}
