//------------------------------------------------------------------------
// P.J. Williams
// Sandia National Laboratories
// pwillia@sandia.gov
// Last Modified December 2000
//------------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#include <typeinfo>
#ifdef HAVE_STD
#include <cstring>
#include <ctime>
#else
#include <string.h>
#include <time.h>
#endif

#include "OptQNIPS.h"
#include "precisio.h"
#include "cblas.h"
#include "ioformat.h"

using NEWMAT::Real;
using NEWMAT::FloatingPointPrecision;
using NEWMAT::ColumnVector;
using NEWMAT::Matrix;
using NEWMAT::SymmetricMatrix;

namespace OPTPP {

int OptQNIPS::checkDeriv() // check the analytic gradient with FD gradient
{return GOOD;}

SymmetricMatrix OptQNIPS::updateH(SymmetricMatrix&Hk, int k)
{ 
  Real mcheps  = FloatingPointPrecision::Epsilon();
  Real sqrteps = sqrt(mcheps);

  NLP1* nlp1 = nlprob();
  int ndim  = nlp1->getDim();
  real yts, snrm, ynrm, sBs, maxres;
  ColumnVector xc, yk, sk, res, Bsk;
  ColumnVector gradl_curr, gradl_prev, yzmultiplier;
  Matrix Htmp(ndim,ndim);

  if(k == 0){
     initHessian();
     Hk = hessl;
     return Hk;
  }

  // Compute change in x and Lagrangian gradient 
  xc = nlp1->getXc();
  sk = xc - xprev; 
  yzmultiplier = y & z;
  gradl_curr   = getGradL();

  if( nlp->hasConstraints() )
      gradl_prev = gprev - constraintGradientPrev*yzmultiplier;
  else
      gradl_prev = gprev;

  yk = gradl_curr - gradl_prev; 

  // If yts too small, skip bfgs
  yts  = Dot(sk,yk);
  snrm = Norm2(sk);
  ynrm = Norm2(yk);

  if(yts <= sqrteps*snrm*ynrm){
    if (debug_) {
      *optout << "UpdateH: <y,s> = " << e(yts,12,4) << " is too small\n";
      *optout << "UpdateH: The BFGS update is skipped\n";
    }
    hessl = Hk;
    return Hk;
  }

  res = yk - Hk*sk;
  maxres = res.NormInfinity() ;
     
  if(maxres <= sqrteps){
    if (debug_) {
      *optout << "UpdateH: y-Hs = " << e(maxres,12,4) << " is too small\n";
      *optout << "UpdateH: The BFGS update is skipped\n";
    }
    hessl = Hk;
    return Hk;
  }

  // If s'Hs too small, skip bfgs
  Bsk = Hk*sk;
  sBs  = Dot(sk,Bsk);

  if(sBs <= sqrteps*snrm*snrm){
    if (debug_) {
      *optout << "UpdateH: The BFGS update is skipped\n";
    }
    hessl = Hk;
    return Hk;
  }

  // Perfom BFGS update 
  Htmp  = -(Bsk * Bsk.t()) / sBs;
  Htmp +=  (yk * yk.t()) / yts;
  Htmp  = Hk + Htmp;
  Hk << Htmp;
  Htmp.Release();
  if (debug_) {
    *optout << "\nUpdateH: after update, k = " << k << "\n";
    *optout << "UpdateH: sBs  = " << sBs << "\n";
  }
  hessl = Hk;
  return Hk;
}

} // namespace OPTPP
