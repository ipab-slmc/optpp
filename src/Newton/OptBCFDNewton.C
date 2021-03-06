//------------------------------------------------------------------------
// Copyright (C) 1993: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//------------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#if (defined(__sgi) || defined(__xlc__) || defined(__xlC__))
#define WANT_MATH
#else
#define WANT_STREAM
#define WANT_MATH
#endif

#ifdef HAVE_STD
#include <cstring>
#else
#include <string.h>
#endif

#include "OptBCFDNewton.h"
#include "precisio.h"
#include "cblas.h"

using NEWMAT::SymmetricMatrix;

namespace OPTPP {

//------------------------------------------------------------------------
//
//   Bound Constrained Finite-Difference Newton Method member functions
//   checkDeriv()
//   updateH()
//   reset()
//------------------------------------------------------------------------

static char* class_name = "OptBCFDNewton";

int OptBCFDNewton::checkDeriv() // check the analytic gradient with FD gradient
{return GOOD;}

//---------------------------------------------------------------------------- 
//
// Update Hessian using a Finite-Difference approximation
//
//---------------------------------------------------------------------------- 
SymmetricMatrix OptBCFDNewton::updateH(SymmetricMatrix&, int) 
{
  if (trace) *optout << class_name << ":UpdateH\n";
  return nlprob()->evalH();
}

void OptBCFDNewton::reset()
{
   NLP1* nlp = nlprob();
   int   n   = nlp->getDim();
   nlp->reset();
   OptimizeClass::defaultReset(n);
   nactive  = 0;
   work_set = false;
}

} // namespace OPTPP
