//------------------------------------------------------------------------
// P.J. Williams
// Sandia National Laboratories
// pwillia@sandia.gov
//------------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#if (defined(__sgi) || defined(__xlC__) || defined (__xlc__))
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

#include "OptppArray.h"
#include "NLP.h"
#include "NLF.h"
#include "Constraint.h"
#include "CompoundConstraint.h"
#include "BoundConstraint.h"
#include "LinearEquation.h"
#include "LinearInequality.h"
#include "NonLinearEquation.h"
#include "NonLinearInequality.h"
#include "NonLinearConstraint.h"

using NEWMAT::ColumnVector;
using NEWMAT::Matrix;
using NEWMAT::SymmetricMatrix;

using namespace OPTPP;

void init_hs1(int ndim, ColumnVector& x)
{
  if (ndim != 2)
  {
    exit (1);
  }
  double factor = 0.0;
  x(1) = -2 - (factor - 1)*3.0;
  x(2) =  1 + (factor - 1)*0.0;
}

void hs1(int mode, int n, const ColumnVector& x, double& fx, ColumnVector& g, int& result)
{ // Hock and Schittkowski's Problem 1 (the objective fcn)
  double f1, f2, x1, x2;

  if (n != 2) return;

  x1 = x(1);
  x2 = x(2);
  f1 = x2 - x1*x1;
  f2 = 1  - x1 ;
  if (mode & NLPFunction) {
    fx     = 100*f1*f1+ f2*f2;
    result = NLPFunction;
  }

  if (mode & NLPGradient) {
    g(1)   = -400*x1*f1 - 2*f2;
    g(2)   = 200*f1;
    result = NLPGradient;
  }
}

CompoundConstraint* create_constraint_hs1(int n)
{ // Hock and Schittkowski's Problem 1 

  ColumnVector lower(n);
  lower << -1.0e+40 << -1.5;
  Constraint bc = new BoundConstraint(n,lower);
  CompoundConstraint* constraints  = new CompoundConstraint(bc);
  return constraints;
}

void init_hs2(int ndim, ColumnVector& x)
{
  if (ndim != 2)
  {
    exit (1);
  }

  double factor = 0.0;
  x(1) = -2.0000 - (factor - 1)*3.0;
  x(2) = 1.5000  + (factor - 1)*0.0;
}

void hs2(int mode, int n, const ColumnVector& x, double& fx, ColumnVector& g, int& result)
{ // Hock and Schittkowski's Problem 2 (the objective fcn)
  double f1, f2, x1, x2;

  if (n != 2) return;

  x1 = x(1);
  x2 = x(2);
  f1 = x2 - x1*x1;
  f2 = 1  - x1 ;

  if (mode & NLPFunction) {
    fx  = 100*f1*f1+ f2*f2;
    result = NLPFunction;
  }
  if (mode & NLPGradient) {
    g(1) = -400*x1*f1 - 2*f2;
    g(2) = 200*f1;
    result = NLPGradient;
  }
}

CompoundConstraint* create_constraint_hs2(int n)
{ // Hock and Schittkowski's Problem 2 

  ColumnVector lower(n);
  lower << -1.0e40 << 1.5;
  Constraint bc = new BoundConstraint(n,lower);
  CompoundConstraint* constraints = new CompoundConstraint( bc );
  return constraints;
}

void init_hs5(int ndim, ColumnVector& x)
{
  if (ndim != 2)
  {
    exit (1);
  }
  double factor = 0.0;
  x(1) = 0.0000 + (factor - 1) *(1.0472 - .5); //3.0
  x(2) = 0.0000 + (factor - 1) *(1.0472 + .5);  //2.0
}

void hs5(int mode, int n, const ColumnVector& x, double& fx, ColumnVector& g, int& result)
{ // Hock and Schittkowski's Problem 5 (the objective fcn)
  double f1, f2, x1, x2;

  if (n != 2) return;

  x1 = x(1);
  x2 = x(2);
  f1 = x1 + x2;
  f2 = x1 - x2;

  if (mode & NLPFunction) {
    fx  = sin(f1)+ f2*f2 -1.5*x1 + 2.5*x2 + 1;
    result = NLPFunction;
  }
  if (mode & NLPGradient) {
    g(1) =  cos(f1) + 2.0*f2 - 1.5;
    g(2) =  cos(f1) - 2.0*f2 + 2.5;
    result = NLPGradient;
  }
}

CompoundConstraint* create_constraint_hs5(int n)
{ // Hock and Schittkowski's Problem 5 

  ColumnVector lower(n);  
  lower << -1.5 << -3.0 ;
  ColumnVector upper(n);  
  upper << 4.0 << 3.0 ;
  CompoundConstraint* constraints = 
           new CompoundConstraint( new BoundConstraint(n, lower, upper) );
  return constraints;
}

void init_hs6(int ndim, ColumnVector& x)
{
  if (ndim != 2)
  {
    exit (1);
  }
  double factor = 0.0;
  x(1) = -1.2  - (factor - 1)*0.219999998951148e+01   ;
  x(2) = 1 + (factor - 1)*0.213770012802428e-07   ;
}

void hs6(int mode, int n, const ColumnVector& x, double& fx, ColumnVector& g, int& result)
{ // Hock and Schittkowski's Problem 6 (the objective fcn)
  double f1, x1;

  if (n != 2) return;

  x1 = x(1);
  f1 = 1.0 - x1;

  if (mode & NLPFunction) {
    fx  = f1*f1;
    result = NLPFunction;
  }
  if (mode & NLPGradient) {
    g(1) = -2*f1;
    g(2) = 0.0;
    result = NLPGradient;
  }
}

void eqn_hs6(int mode, int n, const ColumnVector& x, ColumnVector& fx, Matrix& g, int& result)
{ // Hock and Schittkowski's Problem 6 (the constraints - nonlinear equation)
  double f1, f2, x1, x2;

  if (n != 2) return;

  x1 = x(1);
  x2 = x(2);
  f1 = x2;
  f2 = x1*x1;

  if (mode & NLPFunction) {
    fx  = 10*(f1 - f2);
    result = NLPFunction;
  }
  if (mode & NLPGradient) {
    g(1,1) = -20*x1;
    g(2,1) =  10;
    result = NLPGradient;
  }
}

CompoundConstraint* create_constraint_hs6(int n)
{ // Hock and Schittkowski's Problem 10 
  // (the constraints - nonlinear inequality)

  // Construct the nonlinear equation
  NLP* chs6        = new NLP( new NLF1(n,1,eqn_hs6,init_hs6));
  Constraint nleqn = new NonLinearEquation(chs6);
  CompoundConstraint* constraints =  new CompoundConstraint( nleqn );
  return constraints;
}

void init_hs7(int ndim, ColumnVector& x)
{
  if (ndim != 2)
  {
    exit (1);
  }
  double factor = 0.0;
  x(1) = 2   + (factor - 1)*2;
  x(2) = 2   + (factor - 1)*(2 - sqrt(3.0));
}

void hs7(int mode, int n, const ColumnVector& x, double& fx, ColumnVector& g, int& result)
{ // Hock and Schittkowski's Problem 7 (the objective fcn)
  double f1, f2, x1, x2;

  if (n != 2) return;

  x1 = x(1);
  x2 = x(2);
  f1 = 1 + x1*x1;
  f2 = x2;

  if (mode & NLPFunction) {
    fx  = log(f1) - f2;
    result = NLPFunction;
  }
  if (mode & NLPGradient) {
    g(1) = (2*x1)/f1;
    g(2) = -1.0;
    result = NLPGradient;
  }
}

void eqn_hs7(int mode, int n, const ColumnVector& x, ColumnVector& fx, Matrix& g, int& result)
{ // Hock and Schittkowski's Problem 7 (the constraints - nonlinear equation)
  double f1, f2, x1, x2;

  if (n != 2) return;

  x1 = x(1);
  x2 = x(2);
  f1 = 1 + x1*x1;
  f2 = x2;

  if (mode & NLPFunction) {
    fx  = f1*f1 + f2*f2 - 4;
    result = NLPFunction;
  }
  if (mode & NLPGradient) {
    g(1,1) = 4*f1*x1;
    g(2,1) = 2*f2;
    result = NLPGradient;
  }
}

CompoundConstraint* create_constraint_hs7(int n)
{ // Hock and Schittkowski's Problem 70 

  // Construct the nonlinear equation
  NLP* chs7        = new NLP( new NLF1(n,1,eqn_hs7,init_hs7));
  Constraint nleqn = new NonLinearEquation(chs7);
  CompoundConstraint* constraints = new CompoundConstraint( nleqn );
  return constraints;
}

void init_hs10(int ndim, ColumnVector& x)
{
  if (ndim != 2)
  {
    exit (1);
  }
  double factor = 0.0;
  x(1) = -10.0 - (factor -1)*10.0;
  x(2) =  10.0 + (factor -1)*9.0;
}


void hs10(int mode, int n, const ColumnVector& x, double& fx, ColumnVector& g, int& result)
{ // Hock and Schittkowski's Problem 10
  double f1, f2, x1, x2;

  if (n != 2) return;

  x1 = x(1);
  x2 = x(2);
  f1 = x1;
  f2 = x2;

  if (mode & NLPFunction) {
    fx  = f1 - f2;
    result = NLPFunction;
  }
  if (mode & NLPGradient) {
    g(1) = 1;
    g(2) = -1;
    result = NLPGradient;
  }
}

void ineq_hs10(int mode, int n, const ColumnVector& x, ColumnVector& fx, Matrix& g, int& result)
{ // Hock and Schittkowski's Problem 10 (the constraints - nonlinear inequality)
  double f1, f2, x1, x2;

  if (n != 2) return;

  x1 = x(1);
  x2 = x(2);
  f1 = x1*x1;
  f2 = x2*x2;

  if (mode & NLPFunction) {
    fx  = -3.0*f1 + 2*x1*x2 - f2 + 1.0;
    result = NLPFunction;
  }
  if (mode & NLPGradient) {
    g(1,1) = -6.0*x1 + 2.0*x2;
    g(2,1) =  2.0*x1 - 2.0*x2;
    result = NLPGradient;
  }
}


CompoundConstraint* create_constraint_hs10(int n)
{ // Hock and Schittkowski's Problem 10 
  // (the constraints - nonlinear inequality)

  // Construct the nonlinear inequality
  NLP* chs10        = new NLP( new NLF1(n,1,ineq_hs10,init_hs10) );
  CompoundConstraint* constraints = 
         new CompoundConstraint( new NonLinearInequality(chs10) );
  return constraints;
}

void init_hs13(int ndim, ColumnVector& x)
{
  if (ndim != 2)
  {
    exit (1);
  }
  x(1) = -2.0; // -2.0
  x(2) = -2.0; // -2.0
}


void hs13(int mode, int n, const ColumnVector& x, double& fx, ColumnVector& g, int& result)
{ // Hock and Schittkowski's Problem 13
  double f1, f2, x1, x2;

  if (n != 2) return;

  x1 = x(1);
  x2 = x(2);
  f1 = x1 - 2.0;
  f2 = x2;

  if (mode & NLPFunction) {
    fx  = f1*f1 + f2*f2;
    result = NLPFunction;
  }
  if (mode & NLPGradient) {
    g(1) = 2*f1 ;
    g(2) = 2*f2;
    result = NLPGradient;
  }
}

void ineq_hs13(int mode, int n, const ColumnVector& x, ColumnVector& fx, Matrix& g, int& result)
{ // Hock and Schittkowski's Problem 13 (the constraints - nonlinear inequality)
  double f1, f2, x1, x2;

  if (n != 2) return;

  x1 = x(1);
  x2 = x(2);
  f1 = 1.0 - x1;
  f2 = x2;

  if (mode & NLPFunction) {
    fx  = f1*f1*f1 - f2;
    result = NLPFunction;
  }
  if (mode & NLPGradient) {
    g(1,1) = -3.0*f1*f1;
    g(2,1) = -1.0;
    result = NLPGradient;
  }
}

CompoundConstraint* create_constraint_hs13(int n)
{ // Hock and Schittkowski's Problem 13 

  // Construct the nonlinear inequality
  NLP* chs13        = new NLP( new NLF1(n,1,ineq_hs13,init_hs13) );
  Constraint nlineq = new NonLinearInequality(chs13);
  ColumnVector lower(2);
  lower = 0.0;
  Constraint bc = new BoundConstraint(n,lower);
  CompoundConstraint* constraints = new CompoundConstraint(nlineq,bc);
  return constraints;
}

void init_hs14(int ndim, ColumnVector& x)
{
  if (ndim != 2)
  {
    exit (1);
  }
  double factor = 0.0;
  x(1) = 2+(factor -1)* 1.1771243447;
  x(2) = 2+(factor -1)* 1.0885621722;
}
void hs14(int mode, int n, const ColumnVector& x, double& fx, ColumnVector& g, int& result)
{ // Hock and Schittkowski's Problem 14
  double f1, f2, x1, x2;

  if (n != 2) return;

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

  if (n != 2) return;

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

CompoundConstraint* create_constraint_hs14(int n)
{ // Hock and Schittkowski's Problem 14 
  // (the constraints - linear equation and nonlinear inequality)

  // Construct the linear equation
  Matrix A(1,n);
  A << 1.0 << -2.0 ;
  ColumnVector b(1);
  b << -1.0;
  Constraint leqn   = new LinearEquation(A,b);

  // Construct the nonlinear inequality
  NLP* chs14         = new NLP( new NLF1(n,1,ineq_hs14,init_hs14) );
  Constraint nlineq = new NonLinearInequality(chs14);
  CompoundConstraint* constraints = new CompoundConstraint(nlineq, leqn);
  return constraints;
}

void init_hsuncon(int ndim, ColumnVector& x)
{
  if (ndim != 2)
  {
    exit (1);
  }
  x(1) = 1.0;
  x(2) = 1.0;
}

void hsuncon(int mode, int n, const ColumnVector& x, double& fx, 
             ColumnVector& g, int& result)
{ // Hock and Schittkowski's Problem (Unconstrained)
  double f1, f2, x1, x2;

  if (n != 2) return;

  x1 = x(1);
  x2 = x(2);
  f1 = x1*x1;
  f2 = x2*x2;

  if (mode & NLPFunction) {
    fx  = f1 + f2 ;
    result = NLPFunction;
  }
  if (mode & NLPGradient) {
    g(1) =  2.0*x1;
    g(2) =  2.0*x2;
    result = NLPGradient;
  }
}

void hsuncon2(int mode, int n, const ColumnVector& x, double& fx, 
             ColumnVector& g, SymmetricMatrix& H, int& result)
{ // Hock and Schittkowski's Problem (Unconstrained)
  double f1, f2, x1, x2;

  if (n != 2) return;

  x1 = x(1);
  x2 = x(2);
  f1 = x1*x1;
  f2 = x2*x2;

  if (mode & NLPFunction) {
    fx  = f1 + f2 ;
    result = NLPFunction;
  }
  if (mode & NLPGradient) {
    g(1) =  2.0*x1;
    g(2) =  2.0*x2;
    result = NLPGradient;
  }
  if (mode & NLPHessian) {
    H(1,1) =  2.0;
    H(2,1) =  0.0;
    H(2,2) =  2.0;
    result = NLPHessian;
  }
}

void init_hs26(int ndim, ColumnVector& x)
{
  if (ndim != 3)
  {
    exit (1);
  }
  double factor = 0.0;
  x(1) = -2.6 - (factor - 1)*3.6; 
  x(2) =  2.0 + (factor - 1)*1.0;
  x(3) =  2.0 + (factor - 1)*.10;
}

void hs26(int mode, int n, const ColumnVector& x, double& fx, ColumnVector& g, int& result)
{ // Hock and Schittkowski's Problem 26 (the objective fcn)
  double f1, f2, x1, x2, x3;

  if (n != 3) return;

  x1 = x(1);
  x2 = x(2);
  x3 = x(3);
  f1 = x1 - x2;
  f2 = x2 - x3;

  if (mode & NLPFunction) {
    fx  =  f1*f1 + f2*f2*f2*f2;
    result = NLPFunction;
  }
  if (mode & NLPGradient) {
    g(1) =  2*f1;
    g(2) = -2*f1 +4*f2*f2*f2;
    g(3) = -4*f2*f2*f2;
    result = NLPGradient;
  }
}

void eqn_hs26(int mode, int n, const ColumnVector& x, ColumnVector& fx, Matrix& g, int& result)
{ // Hock and Schittkowski's Problem 26 (the constraints - nonlinear equation)
  double f1, f2, x1, x2, x3;

  if (n != 3) return;

  x1 = x(1);
  x2 = x(2);
  x3 = x(3);
  f1 = 1 + x2*x2;
  f2 = x3;

  if (mode & NLPFunction) {
    fx  = f1*x1 +f2*f2*f2*f2 - 3;
    result = NLPFunction;
  }
  if (mode & NLPGradient) {
    g(1,1) = f1;
    g(2,1) = 2*x1*x2;
    g(3,1) = 4*f2*f2*f2;
    result = NLPGradient;
  }
}

CompoundConstraint* create_constraint_hs26(int n)
{ // Hock and Schittkowski's Problem 26 
  // (the constraints - nonlinear inequality)

  // Construct the nonlinear equation
  NLP* chs26       = new NLP( new NLF1(n,1,eqn_hs26,init_hs26));
  Constraint nleqn = new NonLinearEquation(chs26);
  CompoundConstraint* constraints = new CompoundConstraint(nleqn);
  return constraints;
}

void init_hs28(int ndim, ColumnVector& x)
{
  if (ndim != 3)
  {
    exit (1);
  }
  x(1) = -4;
  x(2) = 1;
  x(3) = 1;
}

void hs28(int mode, int n, const ColumnVector& x, double& fx, ColumnVector& g, int& result)
{ // Hock and Schittkowski's Problem 28 (the objective fcn)
  double f1, f2, x1, x2, x3;

  if (n != 3) return;

  x1 = x(1);
  x2 = x(2);
  x3 = x(3);
  f1 = x1 + x2;
  f2 = x2 + x3 ;

  if (mode & NLPFunction) {
    fx  = f1*f1+ f2*f2;
    result = NLPFunction;
  }
  if (mode & NLPGradient) {
    g(1) = 2*f1;
    g(2) = 2*f1 + 2*f2;
    g(3) = 2*f2;
    result = NLPGradient;
  }
}

CompoundConstraint* create_constraint_hs28(int n)
{ // Hock and Schittkowski's Problem 28 

  Matrix A(1,n);
  ColumnVector b(A.Nrows());
  A <<  1.0 <<  2.0 <<  3.0;
  b <<  1.0; 
  Constraint c1    = new LinearEquation(A,b); 
  CompoundConstraint* constraints = new CompoundConstraint(c1);
  return constraints;
}

void init_hs35(int ndim, ColumnVector& x)
{
  if (ndim != 3)
  {
    exit (1);
  }
  double factor = 0.0;
  x(1) = 0.5 - (factor - 1)*.8333333333333;
  x(2) = 0.5 - (factor - 1)*.2777777777778;
  x(3) = 0.5 + (factor - 1)*.0555555555558;
  //  x(1) = 0.5 ;
  //  x(2) = 0.5 ;
  //  x(3) = 0.5 ;
}

void hs35(int mode, int n, const ColumnVector& x, double& fx, ColumnVector& g, int& result)
{ // Hock and Schittkowski's Problem 35 (the objective fcn)
  double f1, f2, x1, x2, x3;

  if (n != 3) return;

  x1 = x(1);
  x2 = x(2);
  x3 = x(3);
  f1 = x1*x2;
  f2 = x1*x3;

  if (mode & NLPFunction) {
    fx  =  9 - 8*x1 - 6*x2 - 4*x3 +2*x1*x1 + 2*x2*x2 + x3*x3 + 2*f1 + 2*f2;
    result = NLPFunction;
  }
  if (mode & NLPGradient) {
    g(1) =  -8 + 4*x1 + 2*x2 + 2*x3;
    g(2) =  -6 + 4*x2 + 2*x1;
    g(3) =  -4 + 2*x3 + 2*x1;
    result = NLPGradient;
  }
}

CompoundConstraint* create_constraint_hs35(int n)
{ // Hock and Schittkowski's Problem 35 

  // Construct the nonlinear equation
  Matrix A(1,n);
  ColumnVector b(A.Nrows()), lower(n);
  A <<  -1.0 << -1.0 << -2.0 ;
  b     = -3.0;
  lower =  0.0 ;
  Constraint c1    = new LinearInequality(A,b); 
  Constraint bc    = new BoundConstraint(n,lower); 
  CompoundConstraint* constraints = new CompoundConstraint(c1, bc);
  return constraints;
}

void init_hs65(int ndim, ColumnVector& x)
{
  if (ndim != 3)
  {
    exit (1);
  }
  double factor = 0.0;
  x(1) = -5.0  - (factor - 1)*8.6505  ;
  x(2) =  5.0  + (factor - 1)*1.3495  ;
  x(3) =  0.0  - (factor - 1)*4.6204  ;
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

void hs65_2(int mode, int n, const ColumnVector& x, double& fx, ColumnVector& g, SymmetricMatrix& H, int& result)
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

void ineq_hs65(int mode, int n, const ColumnVector& x, ColumnVector& fx, Matrix& g, int& result)
{ // Hock and Schittkowski's Problem 65 
  double f1, f2, f3, x1, x2, x3;

  if (n != 3) return;

  x1 = x(1);
  x2 = x(2);
  x3 = x(3);
  f1 = x1;
  f2 = x2;
  f3 = x3;

  if (mode & NLPFunction) {
    fx  = 48 - f1*f1 - f2*f2 - f3*f3;
    result = NLPFunction;
  }
  if (mode & NLPGradient) {
    g(1,1) = -2*x1;
    g(2,1) = -2*x2;
    g(3,1) = -2*x3;
    result = NLPGradient;
  }
}

void ineq_hs65_2(int mode, int n, const ColumnVector& x, ColumnVector& fx, Matrix& g, OptppArray<SymmetricMatrix>& H, int& result)
{ // Hock and Schittkowski's Problem 65 
  double f1, f2, f3, x1, x2, x3;
  SymmetricMatrix Htmp(n);

  if (n != 3) return;

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


CompoundConstraint* create_constraint_hs65(int n)
{ // Hock and Schittkowski's Problem 65 
  // (the constraints - nonlinear inequality)

  // Construct the nonlinear equation
  NLP* chs65        = new NLP( new NLF1(n,1,ineq_hs65,init_hs65));
  Constraint nleqn = new NonLinearInequality(chs65);
  ColumnVector lower(n); 
  lower << -4.5 << -4.5 << -5.0;
  ColumnVector upper(n); 
  upper <<  4.5 <<  4.5 <<  5.0 ;
  Constraint c1    = new BoundConstraint(n,lower,upper); 
  CompoundConstraint* constraints = new CompoundConstraint(nleqn,c1);
  return constraints;
}

CompoundConstraint* create_constraint_hs65_2(int n)
{ // Hock and Schittkowski's Problem 65 
  // (the constraints - nonlinear inequality)

  // Construct the nonlinear equation
  NLP* chs65        = new NLP( new NLF2(n,1,ineq_hs65_2,init_hs65));
  Constraint nleqn = new NonLinearInequality(chs65);
  ColumnVector lower(n); 
  lower << -4.5 << -4.5 << -5.0;
  ColumnVector upper(n); 
  upper <<  4.5 <<  4.5 <<  5.0 ;
  Constraint c1    = new BoundConstraint(n,lower,upper); 
  CompoundConstraint* constraints = new CompoundConstraint(nleqn,c1);
  return constraints;
}

void init_hs77(int ndim, ColumnVector& x)
{
  if (ndim != 5)
  {
    exit (1);
  }
  x(1) =  2.0; // +  (factor - 1)*0.895140986716219e+01;
  x(2) =  2.0; // +  (factor - 1)*0.803325789100716e+01;
  x(3) =  2.0; // +  (factor - 1)*0.464737743817106e+01;
  x(4) =  2.0;
  x(5) =  2.0;
}

void hs77(int mode, int n, const ColumnVector& x, double& fx, ColumnVector& g, int& result)
{ // Hock and Schittkowski's Problem 77 (the objective fcn)
  double f1, f2, f3, f4, f5, x1, x2, x3, x4, x5;

  if (n != 5) return;

  x1 = x(1);
  x2 = x(2);
  x3 = x(3);
  x4 = x(4);
  x5 = x(5);
  f1 = x1 - 1.0;
  f2 = x1 - x2;
  f3 = x3 - 1.0;
  f4 = x4 - 1.0;
  f5 = x5 - 1.0;

  if (mode & NLPFunction) {
    fx  = f1*f1+ f2*f2 +f3*f3 + pow(f4,4) + pow(f5,6);
    result = NLPFunction;
  }
  if (mode & NLPGradient) {
    g(1) =  2*f1 + 2.0*f2;
    g(2) = -2*f2;
    g(3) =  2*f3;
    g(4) =  4*f4*f4*f4;
    g(5) =  6*pow(f5,5);
    result = NLPGradient;
  }
}

void ineq_hs77(int mode, int n, const ColumnVector& x, ColumnVector& fx, Matrix& g, int& result)
{ // Hock and Schittkowski's Problem 77a 
  double f1, f2, x1, x2, x3, x4, x5;

  if (n != 5) return;

  x1 = x(1);
  x2 = x(2);
  x3 = x(3);
  x4 = x(4);
  x4 = x(4);
  x5 = x(5);
  f1 = x1*x1;
  f2 = x4 - x5;

  if (mode & NLPFunction) {
    fx(1)  = f1*x4 + sin(f2) - 2.0*sqrt(2.0);
    f1     = x2;
    f2     = pow(x3,4)*x4*x4;
    fx(2)  = f1 + f2 - 8 - sqrt(2.0);
    result = NLPFunction ;
  }
  if (mode & NLPGradient) {
    g(1,1) = 2*x1*x4;
    g(2,1) = 0.0;
    g(3,1) = 0.0;
    g(4,1) = f1 + cos(f2);
    g(5,1) = -cos(f2);
    g(1,2) = 0.0;
    g(2,2) = 1.0;
    g(3,2) = 4*pow(x3,3)*x4*x4;
    g(4,2) = 2*x4*pow(x3,4);
    g(5,2) = 0.0;
    result = NLPGradient;
  }
}


CompoundConstraint* create_constraint_hs77(int n)
{ // Hock and Schittkowski's Problem 77 
  // (the constraints - nonlinear inequality)

  // Construct the nonlinear equation
  NLP* chs77a        = new NLP( new NLF1(n,2,ineq_hs77,init_hs77));
  Constraint c1 = new NonLinearEquation(chs77a,2);
  CompoundConstraint* constraints = new CompoundConstraint(c1);
  return constraints;
}

void init_hs78(int ndim, ColumnVector& x)
{
  if (ndim != 5)
  {
    exit (1);
  }
  x(1) = -2.0;
  x(2) =  1.5;
  x(3) =  2.0;
  x(4) = -1.0;
  x(5) = -1.0;
}

void hs78(int mode, int n, const ColumnVector& x, double& fx, ColumnVector& g, int& result)
{ // Hock and Schittkowski's Problem 77 (the objective fcn)
  double x1, x2, x3, x4, x5;

  if (n != 5) return;

  x1 = x(1);
  x2 = x(2);
  x3 = x(3);
  x4 = x(4);
  x5 = x(5);

  if (mode & NLPFunction) {
    fx  = x1*x2*x3*x4*x5;
    result = NLPFunction;
  }
  if (mode & NLPGradient) {
    g(1) =  x2*x3*x4*x5;
    g(2) =  x1*x3*x4*x5;
    g(3) =  x2*x1*x4*x5;
    g(4) =  x2*x3*x1*x5;
    g(5) =  x2*x3*x4*x1;
    result = NLPGradient;
  }
}

void ineq_hs78(int mode, int n, const ColumnVector& x, ColumnVector& fx, Matrix& g, int& result)
{ // Hock and Schittkowski's Problem 78a 
  double x1, x2, x3, x4, x5;

  if (n != 5) return;

  x1 = x(1);
  x2 = x(2);
  x3 = x(3);
  x4 = x(4);
  x5 = x(5);

  if (mode & NLPFunction) {
    fx(1)  = x1*x1 + x2*x2 + x3*x3 + x4*x4 + x5*x5 -10;
    fx(2)  = x2*x3 - 5*x4*x5;
    fx(3)  = pow(x1,3) + pow(x2,3) + 1;
    result = NLPFunction;
  }
  if (mode & NLPGradient) {
    g(1,1) = 2*x1;
    g(2,1) = 2*x2;
    g(3,1) = 2*x3;
    g(4,1) = 2*x4;
    g(5,1) = 2*x5;
    g(1,2) = 0.0;
    g(2,2) = x3;
    g(3,2) = x2;
    g(4,2) = -5*x5;
    g(5,2) = -5*x4;
    g(1,3) = 3.0*x1*x1;
    g(2,3) = 3.0*x2*x2;
    g(3,3) = 0.0;
    g(4,3) = 0.0;
    g(5,3) = 0.0;
    result = NLPGradient;
  }
}


CompoundConstraint* create_constraint_hs78(int n)
{ // Hock and Schittkowski's Problem 78 
  // (the constraints - nonlinear inequality)

  // Construct the nonlinear equation
  NLP* chs78a   = new NLP( new NLF1(n,3,ineq_hs78,init_hs78));
  ColumnVector b(3);
  b = 0.0;
  Constraint c1 = new NonLinearEquation(chs78a, b, 3);
  CompoundConstraint* constraints = new CompoundConstraint(c1);
  return constraints;
}

void hs78_2(int mode, int n, const ColumnVector& x, double& fx, ColumnVector& g, SymmetricMatrix& H, int& result)
{ // Hock and Schittkowski's Problem 77 (the objective fcn)
  double x1, x2, x3, x4, x5;

  if (n != 5) return;

  x1 = x(1);
  x2 = x(2);
  x3 = x(3);
  x4 = x(4);
  x5 = x(5);

  if (mode & NLPFunction) {
    fx  = x1*x2*x3*x4*x5;
    result = NLPFunction;
  }
  if (mode & NLPGradient) {
    g(1) =  x2*x3*x4*x5;
    g(2) =  x1*x3*x4*x5;
    g(3) =  x2*x1*x4*x5;
    g(4) =  x2*x3*x1*x5;
    g(5) =  x2*x3*x4*x1;
    result = NLPGradient;
  }
  if(mode & NLPHessian){
    H(1,1) = 0.0;

    H(2,1) = x3*x4*x5;
    H(2,2) = 0.0;

    H(3,1) = x2*x4*x5;
    H(3,2) = x1*x4*x5;
    H(3,3) = 0.0;

    H(4,1) = x2*x3*x5;
    H(4,2) = x1*x3*x5;
    H(4,3) = x2*x1*x5;
    H(4,4) = 0.0;

    H(5,1) = x2*x3*x4;
    H(5,2) = x1*x3*x4;
    H(5,3) = x2*x1*x4;
    H(5,4) = x2*x3*x1;
    H(5,5) = 0.0;
    result = NLPHessian;
  }
}

void ineq_hs78_2(int mode, int n, const ColumnVector& x, ColumnVector& fx, Matrix& g, OptppArray<SymmetricMatrix>& H, int& result)
{ // Hock and Schittkowski's Problem 78a 
  double x1, x2, x3, x4, x5;
  SymmetricMatrix H1(n), H2(n), H3(n);

  if (n != 5) return;

  x1 = x(1);
  x2 = x(2);
  x3 = x(3);
  x4 = x(4);
  x5 = x(5);

  if (mode & NLPFunction) {
    fx(1)  = x1*x1 + x2*x2 + x3*x3 + x4*x4 + x5*x5 -10;
    fx(2)  = x2*x3 - 5*x4*x5;
    fx(3)  = pow(x1,3) + pow(x2,3) + 1;
    result = NLPFunction;
  }
  if (mode & NLPGradient) {
    g(1,1) = 2*x1;
    g(2,1) = 2*x2;
    g(3,1) = 2*x3;
    g(4,1) = 2*x4;
    g(5,1) = 2*x5;
    g(1,2) = 0.0;
    g(2,2) = x3;
    g(3,2) = x2;
    g(4,2) = -5*x5;
    g(5,2) = -5*x4;
    g(1,3) = 3.0*x1*x1;
    g(2,3) = 3.0*x2*x2;
    g(3,3) = 0.0;
    g(4,3) = 0.0;
    g(5,3) = 0.0;
    result = NLPGradient;
  }
  if (mode & NLPHessian) {
    H1(1,1) = 2.0; 

    H1(2,1) = 0.0; 
    H1(2,2) = 2.0; 

    H1(3,1) = 0.0; 
    H1(3,2) = 0.0; 
    H1(3,3) = 2.0; 

    H1(4,1) = 0.0; 
    H1(4,2) = 0.0; 
    H1(4,3) = 0.0; 
    H1(4,4) = 2.0; 

    H1(5,1) = 0.0; 
    H1(5,2) = 0.0; 
    H1(5,3) = 0.0; 
    H1(5,4) = 0.0; 
    H1(5,5) = 2.0; 

    H2(1,1) = 0.0; 

    H2(2,1) = 0.0; 
    H2(2,2) = 0.0; 

    H2(3,1) = 0.0; 
    H2(3,2) = 1.0; 
    H2(3,3) = 0.0; 

    H2(4,1) = 0.0; 
    H2(4,2) = 0.0; 
    H2(4,3) = 0.0; 
    H2(4,4) = 0.0; 

    H2(5,1) = 0.0; 
    H2(5,2) = 0.0; 
    H2(5,3) = 0.0; 
    H2(5,4) = -5.0; 
    H2(5,5) = 0.0; 

    H3(1,1) = 6.0*x1; 

    H3(2,1) = 0.0; 
    H3(2,2) = 6.0*x2; 

    H3(3,1) = 0.0; 
    H3(3,2) = 0.0; 
    H3(3,3) = 0.0; 

    H3(4,1) = 0.0; 
    H3(4,2) = 0.0; 
    H3(4,3) = 0.0; 
    H3(4,4) = 0.0; 

    H3(5,1) = 0.0; 
    H3(5,2) = 0.0; 
    H3(5,3) = 0.0; 
    H3(5,4) = 0.0; 
    H3(5,5) = 0.0; 

    H[0] = H1;
    H[1] = H2;
    H[2] = H3;
    result = NLPHessian;
 }
}


CompoundConstraint* create_constraint_hs78_2(int n)
{ // Hock and Schittkowski's Problem 78 
  // (the constraints - nonlinear inequality)

  // Construct the nonlinear equation
  NLP* chs78a   = new NLP( new NLF2(n,3,ineq_hs78_2,init_hs78));
  Constraint c1 = new NonLinearEquation(chs78a, 3);
  CompoundConstraint* constraints = new CompoundConstraint(c1);
  return constraints;
}
