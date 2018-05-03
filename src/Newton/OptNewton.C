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
#else
#include <string.h>
#endif

#include "OptNewton.h"
#include "precisio.h"
#include "cblas.h"
#include "ioformat.h"
#include "globals.h"

using NEWMAT::Real;
using NEWMAT::FloatingPointPrecision;
using NEWMAT::ColumnVector;
using NEWMAT::SymmetricMatrix;

namespace OPTPP {

//------------------------------------------------------------------------
//
//   Newton Method member functions
//   checkDeriv()
//   initHessian()
//   printStatus()
//   stepTolNorm()
//   updateH()
//------------------------------------------------------------------------

int OptNewton::checkDeriv() // check the analytic gradient with FD gradient
{
  NLP2* nlp = nlprob2();
  int retcode = checkAnalyticFDGrad();

  Real mcheps = FloatingPointPrecision::Epsilon();
  Real third = 0.3333333;
  double gnorm = nlp->getGrad().NormInfinity();
  double eta   = pow(mcheps,third)*max(1.0,gnorm);
  *optout <<"\ncheck_Deriv: checking Hessian versus finite-differences\n";
  SymmetricMatrix Hess(dim), FDHess(dim), ErrH(dim);
  FDHess = nlp->FDHessian(sx); 
  Hess   = nlp->getHess();
  ErrH   = Hess - FDHess;
  Print(ErrH);
  Real maxerr = ErrH.NormInfinity();
  *optout << "maxerror = " << e(maxerr, 12,4) 
     << "tolerance =  " << e(eta, 12,4) << "\n";
  if (maxerr > eta) retcode = false;

  return retcode;
}

void OptNewton::initHessian()
{
  if (debug_) *optout << "OptNewton::initHessian: \n";
  NLP2* nlp = nlprob2();
  Hessian = nlp->getHess();
  return;
}

void OptNewton::printStatus(char *s) // set Message
{
  NLP1* nlp = nlprob();
  *optout << "\n\n=========  " << s << "  ===========\n\n";
  *optout << "Optimization method       = " << method << "\n";
  *optout << "Dimension of the problem  = " << nlp->getDim()   << "\n";
  *optout << "Return code               = " << ret_code << " ("
       << mesg << ")\n";
  *optout << "No. iterations taken      = " << iter_taken  << "\n";
  *optout << "No. function evaluations  = " << nlp->getFevals() << "\n";
  *optout << "No. gradient evaluations  = " << nlp->getGevals() << "\n";

  if (debug_) {
    *optout << "Hessian \n";
    Print(Hessian);
  }

  tol.printTol(optout);

  nlp->fPrintState(optout, s);
}


real OptNewton::stepTolNorm() const
{
  NLP1* nlp = nlprob();
  ColumnVector step(sx.AsDiagonal()*(nlp->getXc() - xprev));
  return Norm2(step);
}

SymmetricMatrix OptNewton::updateH(SymmetricMatrix&, int)
{
  return nlprob()->evalH();
}

} // namespace OPTPP
