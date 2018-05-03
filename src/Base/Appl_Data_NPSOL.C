//************************************************************************
// Copyright (C) 1996:
// Opt++ group, Livermore
// Sandia National Laboratories
//************************************************************************

#include <iostream>

using namespace std;

#include "Appl_Data_NPSOL.h"
#include "OptppFatalError.h"

using NEWMAT::ColumnVector;
using NEWMAT::Matrix;
using NEWMAT::SymmetricMatrix;

namespace OPTPP {

//************************************************************************
// Constructor 
//------------------------------------------------------------------------
Appl_Data_NPSOL::Appl_Data_NPSOL()
{
  buffer_len    = 1;

  x_buffer      = grad_buffer = constr_buffer = NULL;
  cjac_buffer   = NULL;
  fvalue_status = grad_status = false;
  constr_status = cjac_status = false;

  buffer_pointer = 0;
  dimension = ncnln = -1;
}

//************************************************************************
// Another Constructor 
//------------------------------------------------------------------------
Appl_Data_NPSOL::Appl_Data_NPSOL(int length)
{
  buffer_len    = length;
  if (buffer_len == 0) return;

  x_buffer      = grad_buffer = constr_buffer = NULL;
  cjac_buffer   = NULL;
  fvalue_status = grad_status = false;
  constr_status = cjac_status = false;

  buffer_pointer = 0;
  dimension = ncnln = -1;
}

//************************************************************************
// destructor 
//------------------------------------------------------------------------
Appl_Data_NPSOL::~Appl_Data_NPSOL()
{
  if (buffer_len == 0) return;
  if (x_buffer      != NULL) delete x_buffer;
  if (grad_buffer   != NULL) delete grad_buffer;
  if (constr_buffer != NULL) delete constr_buffer;
  if (cjac_buffer   != NULL) delete cjac_buffer;

}

//************************************************************************
// compare the incoming x with the local ones : return index if it is
// found, otherwise return -1
//------------------------------------------------------------------------
bool Appl_Data_NPSOL::Compare(ColumnVector &x)
{
  int    flag;
  double *x1, *x2;

  x1 = x.Store();
  x2 = (*x_buffer).Store();
  flag = memcmp(x1,x2,dimension*sizeof(double));
  if (flag == 0) return true; else return false;
} 

//************************************************************************
// get function value 
//------------------------------------------------------------------------
bool Appl_Data_NPSOL::getF(ColumnVector &x, real &fvalue)
{
  if (buffer_len == 0) return false;

  if (fvalue_status && Compare(x)){
    fvalue = fvalue_buffer; return true;
  } else return false;  
}

//************************************************************************
// get gradient value 
//------------------------------------------------------------------------
bool Appl_Data_NPSOL::getGrad(ColumnVector &x, ColumnVector &g)
{
  if (buffer_len == 0) return false;

  if (grad_status && Compare(x)){
     g = (*grad_buffer); return true;
  } else return false;  
}

//************************************************************************
// get constraint value 
//------------------------------------------------------------------------
bool Appl_Data_NPSOL::getConstraint(ColumnVector &x, ColumnVector &c)
{
  if (buffer_len == 0) return false;

  if (constr_status && Compare(x)){
     c = (*constr_buffer); return true;
  } else return false;  
}

//************************************************************************
// get constraint Jacobian 
//------------------------------------------------------------------------
bool Appl_Data_NPSOL::getCJacobian(ColumnVector &x, Matrix &cj)
{
  if (buffer_len == 0) return false;
  if (ncnln == 0)      return false;
    
  if (cjac_status && Compare(x)){
     cj = (*cjac_buffer); return true;
  } else return false;  
}

//************************************************************************
// update the local data 
//------------------------------------------------------------------------
void Appl_Data_NPSOL::update(int mode, int dim, ColumnVector &x, real fv)
{
  if (buffer_len == 0) return;

  if (dimension != -1 && dimension != dim) {
      OptppmathError("Dimensions are inconsistent."); 
  } else dimension = dim;

  if (x_buffer != NULL) delete x_buffer; 
  x_buffer = new ColumnVector(dimension); (*x_buffer) = x;
  grad_status = constr_status = cjac_status = false;
  if(mode & NLPFunction) { fvalue_buffer = fv; fvalue_status = true;}
}

//************************************************************************
// update the local data 
//------------------------------------------------------------------------
void Appl_Data_NPSOL::update(int mode, int dim,ColumnVector &x,ColumnVector &g)
{

  if (buffer_len == 0) return;

  if (dimension != -1 && dimension != dim) { 
      OptppmathError("Dimensions are inconsistent."); 
  } else dimension = dim;

  if (x_buffer != NULL) delete x_buffer; 
  x_buffer = new ColumnVector(dimension); (*x_buffer) = x;
  fvalue_status = constr_status = cjac_status = false;

  if(mode & NLPGradient){ 
     if (grad_buffer != NULL) delete grad_buffer; 
     grad_buffer = new ColumnVector(dimension); (*grad_buffer) = g;
     grad_status = true;
  }
}

//************************************************************************
// update the local data 
//------------------------------------------------------------------------
void Appl_Data_NPSOL::update(int dim,ColumnVector &x,int nc,ColumnVector &c)
{
  if (buffer_len == 0) return;

  if ((dimension != -1 && dimension != dim) || (ncnln != -1 && nc != ncnln)) { 
      OptppmathError("Dimensions are inconsistent."); 
  } else { dimension = dim; ncnln = nc; }

  if (x_buffer != NULL) delete x_buffer; 
  x_buffer = new ColumnVector(dimension); (*x_buffer) = x;
  fvalue_status = grad_status = cjac_status = false;

  if (constr_buffer != NULL) delete constr_buffer; 
  constr_buffer = new ColumnVector(ncnln); (*constr_buffer) = c;
  constr_status = true;
}

//************************************************************************
// update the local data 
//------------------------------------------------------------------------
void Appl_Data_NPSOL::update(int mode, int dim,ColumnVector &x,int nc,
		             ColumnVector &c, Matrix &cj)
{
  if (buffer_len == 0) return;

  if ((dimension != -1 && dimension != dim) || (ncnln != -1 && nc != ncnln)) { 
      OptppmathError("Dimensions are inconsistent."); 
  } else { dimension = dim; ncnln = nc; }

  update(dim, x, nc, c);

  if( mode & NLPCJacobian){
    if (cjac_buffer != NULL) delete cjac_buffer; 
    cjac_buffer = new Matrix(dimension,ncnln); (*cjac_buffer) = cj;
    cjac_status = true;
  }
}

//************************************************************************
// update the local data 
//------------------------------------------------------------------------
void Appl_Data_NPSOL::update(int mode, int dim,ColumnVector &x,double fx,int nc,
			     ColumnVector &c)
{
  if (buffer_len == 0) return;

  if ((dimension != -1 && dimension != dim) || (ncnln != -1 && nc != ncnln)) { 
      OptppmathError("Dimensions are inconsistent."); 
  } else { dimension = dim; ncnln = nc; }

  update(dim, x, nc, c);

  update(mode, dim, x, fx);
  if(mode & NLPFunction) { fvalue_buffer = fx; fvalue_status = true; }
}

//************************************************************************
// update the local data 
//------------------------------------------------------------------------
void Appl_Data_NPSOL::update(int mode, int dim,ColumnVector &x,double fx,int nc,
			     ColumnVector &c, Matrix &cj)
{
  if (buffer_len == 0) return;
  update(mode, dim, x, fx);
  update(mode, dim, x, nc, c, cj);
}

} // namespace OPTPP
