//------------------------------------------------------------------------
// Copyright (C) 1993: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//------------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#ifdef HAVE_STD
#include <cstring>
#include <cmath>
#else
#include <string.h>
#include <math.h>
#endif

#include "OptConstrQNewton.h"
#include "precisio.h"
#include "cblas.h"
#include "ioformat.h"

using NEWMAT::Real;
using NEWMAT::FloatingPointPrecision;
using NEWMAT::ColumnVector;
using NEWMAT::Matrix;
using NEWMAT::DiagonalMatrix;
using NEWMAT::SymmetricMatrix;


namespace OPTPP {

//------------------------------------------------------------------------
//
//   Quasi-Newton Method member functions
//   checkDeriv()
//   updateH()
//------------------------------------------------------------------------

// static char* class_name = "OptConstrQNewton";

// Check the analytic gradient with FD gradient
int OptConstrQNewton::checkDeriv()
{ return checkAnalyticFDGrad(); }


//---------------------------------------------------------------------------- 
//
// Update Hessian using a Quasi-Newton update
//
//---------------------------------------------------------------------------- 
SymmetricMatrix OptConstrQNewton::updateH(SymmetricMatrix& Hk, int k) 
{

  Real mcheps = FloatingPointPrecision::Epsilon();
  Real sqrteps = sqrt(mcheps);

  int i;

  NLP1* nlp = nlprob();
  int nr     = nlp->getDim();
  ColumnVector grad(nr), xc;
  xc     = nlp->getXc();
  grad   = nlp->getGrad();

  DiagonalMatrix D(nr);
// BFGS formula
  
  if (k == 0) { // Initial Hessian is set equal to the Identity Matrix
    Hessian = 0.0;
//    D = sx.AsDiagonal()*sx.AsDiagonal();
//    D = sfx.AsDiagonal();
    Real typx, xmax, gnorm;
//    Real gamma;
    gnorm = Norm2(grad);

    // Initialize xmax, typx and D to default values
    xmax   = -1.e30; typx   =  1.0; D      =  1.0;

    for (i=1; i <= nr; i++) xmax = max(xmax,fabs(xc(i)));
    if( xmax != 0.0) typx = xmax;
    if( gnorm!= 0.0) D    = gnorm/typx;
    if (debug_) {
      *optout << "UpdateH: gnorm0 = " << gnorm
	<< "typx = " << typx << "\n";
    }
    for (i=1; i <= nr; i++) Hessian(i,i) = D(i);
    return Hessian;
  }
  
  ColumnVector yk(nr), sk(nr), Bsk(nr);
  Matrix Htmp(nr,nr);
  
  yk = grad - gprev;
  sk = xc   - xprev;
  
  if (debug_) {
    Print(yk);
    Print(sk);
  }

  Real gts = Dot(gprev,sk);
  Real yts = Dot(yk,sk);
  
  Real snorm = Norm2(sk);
  Real ynorm = Norm2(yk);
  
  if (debug_) {
    *optout << "UpdateH: gts   = " << gts 
         << "  yts = " << yts << "\n";
    *optout << "UpdateH: snorm = " << snorm 
         << "  ynorm = " << ynorm << "\n";
  }

  if (yts <= sqrteps*snorm*ynorm) {
    if (debug_) {
      *optout << "UpdateH: <y,s> = " << e(yts,12,4) << " is too small\n";
      *optout << "UpdateH: The BFGS update is skipped\n";
    }
    Hessian = Hk; return Hk;
  }
  
  ColumnVector res(nr);
  res = yk - Hk*sk;
  if (res.NormInfinity() <= sqrteps) {
    if (debug_) {
      *optout << "UpdateH: <y,s> = " << e(yts,12,4) << " is too small\n";
      *optout << "UpdateH: The BFGS update is skipped\n";
    }
    Hessian = Hk; return Hk;
  }
  
  Bsk = Hk*sk;
  Real sBs = Dot(sk,Bsk);
  Real etol = 1.e-8;

  if (sBs <= etol*snorm*snorm) {
    if (debug_) {
      *optout << "UpdateH: <y,s> = " << e(yts,12,4) << " is too small\n";
      *optout << "UpdateH: The BFGS update is skipped\n";
    }
    D = sx.AsDiagonal()*sx.AsDiagonal();
    Hk = 0;
    for (i=1; i <= nr; i++) Hk(i,i) = D(i);
    Hessian = Hk; return Hk;
  }
  
// Otherwise update the Hessian approximation
  if (debug_) {
    //    *optout << "\nUpdateH: before update, k = " << k << "\n";
    //    FPrint(optout, Hk);
  }

  Htmp = - (Bsk * Bsk.t()) / sBs;
  Htmp = Htmp + (yk * yk.t()) / yts;
  Htmp = Hk + Htmp;
  Hk << Htmp;
  Htmp.Release(); 
  ColumnVector Bgk(nr);
  Bgk = Hk*grad;
  Real gBg = Dot(grad,Bgk);
  Real gg  = Dot(grad,grad);
  Real ckp1= gBg/gg;
  if (debug_) {
    //    *optout << "\nUpdateH: after update, k = " << k << "\n";
    //    FPrint(optout, Hk);
    *optout << "UpdateH: sBs  = " << sBs << "\n";
    *optout << "UpdateH: ckp1 = " << ckp1 << "\n";
  }
  Hessian = Hk;
  return Hk;
}

} // namespace OPTPP
