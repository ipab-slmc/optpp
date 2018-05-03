//-----------------------------------------------------------------------
// Copyright (C) 1996:
// Opt++ group, Livermore
// Sandia National Laboratories
//------------------------------------------------------------------------

// 10/01/01 PJW On Solaris operating systems, we need to include
// iostream for compilation purposes.  Without this statement, 
// we get an error with regards to the limits includes.
 
#include<iostream>

#include "Appl_Data.h"

using NEWMAT::ColumnVector;
using NEWMAT::Matrix;
using NEWMAT::SymmetricMatrix;

namespace OPTPP {

//------------------------------------------------------------------------
// Constructor 
//------------------------------------------------------------------------
Appl_Data::Appl_Data()
{
  reset();
}

//------------------------------------------------------------------------
// destructor 
//------------------------------------------------------------------------
Appl_Data::~Appl_Data()
{
  if (xparm        != NULL) delete xparm;
  if (gradient     != NULL) delete gradient;
  if (Hessian      != NULL) delete Hessian;
  if (constraint_value    != NULL) delete constraint_value;
  if (constraint_gradient != NULL) delete constraint_gradient;
  if (constraint_Hessian  != NULL) delete constraint_Hessian;
  if (lsq_residuals       != NULL) delete lsq_residuals;
  if (lsq_jacobian        != NULL) delete lsq_jacobian;
}

void Appl_Data::reset()
{
  xparm = gradient = NULL; Hessian = NULL;
  constraint_value = NULL; constraint_gradient = NULL; 
  constraint_Hessian = NULL;
  lsq_residuals = NULL;    lsq_jacobian = NULL; 
  function_current = gradient_current = Hessian_current = false;
}

//------------------------------------------------------------------------
// compare two vectors 
//------------------------------------------------------------------------
bool Appl_Data::Compare(const ColumnVector &x)
{
  int    flag;
  double *x1, *x2;

  x1 = x.Store();
  x2 = (*xparm).Store();
  flag = memcmp(x1,x2,dimension*sizeof(double));
  if (flag == 0) return true; else return false;
} 

//------------------------------------------------------------------------
// get function value 
//------------------------------------------------------------------------
bool Appl_Data::getF(const ColumnVector &x, real &fvalue)
{
  if (function_current && Compare(x)) {
    fvalue = function_value; return true;
  } else return false;  
}

//------------------------------------------------------------------------
// get gradient value 
//------------------------------------------------------------------------
bool Appl_Data::getGrad(const ColumnVector &x,ColumnVector &g)
{
  if (gradient_current && Compare(x)) {
    g = (*gradient); return true;
  } else return false;  
}

//------------------------------------------------------------------------
// get Hessian 
//------------------------------------------------------------------------
bool Appl_Data::getHess(const ColumnVector &x, SymmetricMatrix &h)
{
  if (Hessian_current && Compare(x)) {
    h = (*Hessian); return true;
  } else return false;  
}

//------------------------------------------------------------------------
// get constraint function value 
//------------------------------------------------------------------------
bool Appl_Data::getCF(const ColumnVector &x, ColumnVector& cvalue)
{
  if (function_current && Compare(x)) {
    cvalue = (*constraint_value); return true;
  } else return false;  
}

//------------------------------------------------------------------------
// get constraint gradient value 
//------------------------------------------------------------------------
bool Appl_Data::getCGrad(const ColumnVector &x, Matrix &g)
{
  if (gradient_current && Compare(x)) {
    g = (*constraint_gradient); return true;
  } else return false;  
}

//------------------------------------------------------------------------
// get constraint Hessian 
//------------------------------------------------------------------------
bool Appl_Data::getCHess(const ColumnVector &x, OptppArray<SymmetricMatrix> &h)
{
  if (Hessian_current && Compare(x)) {
    h = (*constraint_Hessian); return true;
  } else return false;  
}

//------------------------------------------------------------------------
// get residuals of least squares objective function 
//------------------------------------------------------------------------
bool Appl_Data::getLSQF(const ColumnVector &x,ColumnVector &lsqf)
{
  if (function_current && Compare(x)) {
    lsqf = (*lsq_residuals); return true;
  } else return false;  
}

//------------------------------------------------------------------------
// get Jacobian of least squares objective function
//------------------------------------------------------------------------
bool Appl_Data::getLSQJac(const ColumnVector &x, Matrix &j)
{
  if (gradient_current && Compare(x)) {
    j = (*lsq_jacobian); return true;
  } else return false;  
}

//------------------------------------------------------------------------
// update the local data 
//------------------------------------------------------------------------
void Appl_Data::update(int mode,int dim, const ColumnVector &x, real fv)
{
  dimension = dim;
  if (xparm != NULL) delete xparm;
  xparm = new ColumnVector(dimension); (*xparm) = x;
  function_current = gradient_current = Hessian_current = false;
  if (mode & NLPFunction) {function_value = fv; function_current = true;}
}

//------------------------------------------------------------------------
// update the local data 
//------------------------------------------------------------------------
void Appl_Data::update(int mode,int dim, const ColumnVector &x,real fv,
                       ColumnVector &g)
{
  update(mode, dim, x, fv);
  if (mode & NLPGradient) {
    if (gradient != NULL) delete gradient;
    gradient = new ColumnVector(dimension); 
    (*gradient) = g; gradient_current = true;
  }
}

//------------------------------------------------------------------------
// update local data 
//------------------------------------------------------------------------
void Appl_Data::update (int mode, int dim, const ColumnVector & x, real fv,
                       ColumnVector &g, SymmetricMatrix &h)
{
  update(mode, dim, x, fv, g);
  if (mode & NLPHessian) {
    if (Hessian != NULL) delete Hessian;
    Hessian = new SymmetricMatrix(dimension);
    (*Hessian) = h; Hessian_current = true;
  }
}


//------------------------------------------------------------------------
// update the local constraint data 
//------------------------------------------------------------------------
void Appl_Data::constraint_update(int mode, int dim, int ncnln, 
                       const ColumnVector &x, ColumnVector& fv)
{
  dimension = dim;
  if (xparm != NULL) delete xparm;
  xparm = new ColumnVector(dimension); (*xparm) = x;
  function_current = gradient_current = Hessian_current = false;
  if (mode & NLPFunction) {
    if (constraint_value != NULL) delete constraint_value;
    constraint_value = new ColumnVector(ncnln); 
    (*constraint_value) = fv; function_current = true;
  }
}

//------------------------------------------------------------------------
// update the local constraint data 
//------------------------------------------------------------------------
void Appl_Data::constraint_update(int mode, int dim, int ncnln,
                       const ColumnVector &x, ColumnVector& fv, Matrix &g)
{
  constraint_update(mode, dim, ncnln, x, fv);
  if (mode & NLPGradient) {
    if (constraint_gradient != NULL) delete constraint_gradient;
    constraint_gradient = new Matrix(dimension,ncnln); 
    (*constraint_gradient) = g; gradient_current = true;
  }
}

//------------------------------------------------------------------------
// update local constraint data 
//------------------------------------------------------------------------
void Appl_Data::constraint_update (int mode, int dim, int ncnln,
                       const ColumnVector & x, 
                       ColumnVector& fv, Matrix &g, OptppArray<SymmetricMatrix> &h)
{
  constraint_update(mode, dim, ncnln, x, fv, g);
  if (mode & NLPHessian) {
    if (constraint_Hessian != NULL) delete constraint_Hessian;
    constraint_Hessian = new OptppArray<SymmetricMatrix>(ncnln);
    (*constraint_Hessian) = h; Hessian_current = true;
  }
}


//------------------------------------------------------------------------
// update the local least squares data 
//------------------------------------------------------------------------
void Appl_Data::lsq_update(int mode, int dim, int lsqterms, 
                       const ColumnVector &x, ColumnVector& lsqf)
{
  dimension = dim;
  if (xparm != NULL) delete xparm;
  xparm = new ColumnVector(dimension); (*xparm) = x;
  function_current = gradient_current = false;
  if (mode & NLPFunction) {
    if (lsq_residuals != NULL) delete lsq_residuals;
    lsq_residuals = new ColumnVector(lsqterms); 
    (*lsq_residuals) = lsqf; function_current = true;
  }
}

//------------------------------------------------------------------------
// update the local least squares data 
//------------------------------------------------------------------------
void Appl_Data::lsq_update(int mode, int dim, int lsqterms,
                       const ColumnVector &x, ColumnVector& lsqf, Matrix &lsqj)
{
  lsq_update(mode, dim, lsqterms, x, lsqf);
  if (mode & NLPGradient) {
    if (lsq_jacobian != NULL) delete lsq_jacobian;
    lsq_jacobian = new Matrix(lsqterms, dimension); 
    (*lsq_jacobian) = lsqj; gradient_current = true;
  }
}
} // namespace OPTPP
