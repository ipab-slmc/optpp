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

#include "OptFDNIPS.h"
#include "precisio.h"
#include "cblas.h"

using NEWMAT::Real;
using NEWMAT::FloatingPointPrecision;
using NEWMAT::ColumnVector;
using NEWMAT::Matrix;
using NEWMAT::SymmetricMatrix;

//------------------------------------------------------------------------
// external subroutines referenced by this module 
// Included to prevent compilation error when -ansi flag is used.
//------------------------------------------------------------------------

#if !(defined(__GNUC__) && __GNUC__ >= 3)
extern "C" {
  double copysign(double, double);
}
#else
#ifdef CYGWIN
extern "C" {
  extern double copysign _PARAMS((double, double));
}
#endif
#endif

namespace OPTPP {

int OptFDNIPS::checkDeriv() // check the analytic gradient with FD gradient

{return GOOD;}

SymmetricMatrix OptFDNIPS::updateH(SymmetricMatrix&Hk, int k)
{ 
  NLP1* nlp1 = nlprob();
  int ndim  = nlp1->getDim();
  ColumnVector xc, xtmp, gradtmp, yzmultiplier,grad;
  Matrix Htmp(ndim,ndim);
  double hi, hieps;

  Real mcheps = FloatingPointPrecision::Epsilon();
  ColumnVector fcn_accrcy = nlp1->getFcnAccrcy();

  Htmp   = 0;
  xc     = nlp1->getXc();
  yzmultiplier = y & z;

  // Get the gradient of the Lagrangian 
  grad  = getGradL();

  // Build approximation column by column 
  for(int j = 1; j <= ndim ; j++){

     hieps = sqrt(max(mcheps,fcn_accrcy(j)));
     hi     = hieps*max(fabs(xc(j)), sx(j));
     hi     = copysign(hi, xc(j));
     xtmp   = xc;
     xtmp(j)= xc(j) + hi;

     // Evaluate the gradient of the Lagrangian at the new point
     gradtmp  = nlp1->evalLagrangianGradient(xtmp,yzmultiplier,constrType); 

     /* Calculate jth column of Hessian using a first-order 
        finite-difference approximation  */
        Htmp.Column(j) << (gradtmp - grad)/hi;
  }

  // Symmetrize the Hessian Approximation
  Hk   << (Htmp + Htmp.t())/2.0;
  Htmp.Release();
  hessl = Hk;
  return Hk;
}

} // namespace OPTPP
