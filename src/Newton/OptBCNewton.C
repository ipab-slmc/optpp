//------------------------------------------------------------------------
// Copyright (C) 1996:
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

#include "OptBCNewton.h"
#include "precisio.h"
#include "cblas.h"
#include "ioformat.h"

using NEWMAT::Real;
using NEWMAT::FloatingPointPrecision;
using NEWMAT::ColumnVector;
using NEWMAT::SymmetricMatrix;
using NEWMAT::LowerTriangularMatrix;

namespace OPTPP {

static char* class_name = "OptBCNewton";

//------------------------------------------------------------------------
// BCNewton functions
// initHessian
// checkDeriv
// checkConvg
// initOpt
// printStatus
// stepTolNorm
// updateH
// computeSearch
// updateConstraints
// reset
//------------------------------------------------------------------------

void OptBCNewton::initHessian()
{ 
  if (debug_) *optout << class_name << "::initHessian: \n";
  NLP2* nlp = nlprob2();
  Hessian = nlp->getHess();
  return;
}

int OptBCNewton::checkConvg() // Check convergence
{
  NLP1* nlp = nlprob();
  ColumnVector xc(nlp->getXc());
  int   i, n = nlp->getDim();

// Test 1. step tolerance 

  double step_tol = tol.getStepTol();
  double snorm = stepTolNorm();
  double xnorm =  Norm2(xc);
  double stol  = step_tol*max(1.0,xnorm);
  if (snorm  <= stol) {
    strcpy(mesg,"Step tolerance test passed");
    *optout << "checkConvg: snorm = " << snorm 
      << "  stol = " << stol << "\n";
    return 1;
  }
  
// Test 2. function tolerance
  double ftol = tol.getFTol();
  double fvalue = nlp->getF();
  double rftol = ftol*max(1.0,fabs(fvalue));
  Real deltaf = fprev - fvalue;
  if (deltaf <= rftol) {
    strcpy(mesg,"Function tolerance test passed");
    *optout << "BCNewton::checkConvg: Delta f = " << e(deltaf,12,4) 
         << "  ftol = " << e(ftol,12,4) << "\n";
    return 2;
  }
  

// Test 3. gradient tolerance 

  ColumnVector grad(nlp->getGrad());
  double gtol = tol.getGTol();
  double rgtol = gtol*max(1.0,fabs(fvalue));
  for (i=1; i<=n; i++) if (work_set(i) == true) grad(i) = 0.0;
  double gnorm = Norm2(grad);
  if (gnorm <= rgtol) {
    strcpy(mesg,"Gradient tolerance test passed");
    *optout << "checkConvg: gnorm = " << e(gnorm,12,4) 
      << "  gtol = " << e(rgtol, 12,4) << "\n";
    return 3;
  }
  

// Test 4. absolute gradient tolerance 

  if (gnorm <= gtol) {
    strcpy(mesg,"Gradient tolerance test passed");
    *optout << "checkConvg: gnorm = " << e(gnorm,12,4) 
      << "  gtol = " << e(gtol, 12,4) << "\n";
    return 4;
  }
  
  // Nothing to report 

  return 0;

}

int OptBCNewton::checkDeriv() // Check the analytic gradient with FD gradient
{
  NLP2* nlp = nlprob2();
  SymmetricMatrix Hess(dim), FDHess(dim), ErrH(dim);

  int retcode   = checkAnalyticFDGrad();
  Real mcheps   = FloatingPointPrecision::Epsilon();
  Real third   = 0.3333333;
  double gnorm = nlp->getGrad().NormInfinity();
  double eta   = pow(mcheps,third)*max(1.0,gnorm);
  FDHess       = nlp->FDHessian(sx); 
  Hess         = nlp->getHess();
  ErrH         = Hess - FDHess;
  Real maxerr = ErrH.NormInfinity();
  if(debug_){
     *optout <<"\nCheck_Deriv: Checking Hessian versus finite-differences\n";
     *optout << "maxerror = " << e(maxerr, 12,4) 
            << "tolerance =  " << e(eta, 12,4) << "\n";
  }
  if (maxerr > eta) retcode = BAD;
  return retcode;
}


void OptBCNewton::initOpt()
{ 
  OptBCNewtonLike::initOpt();
  updateConstraints(0);
  return;
}

void OptBCNewton::printStatus(char *s) // set Message
{
  NLP2* nlp = nlprob2();
  *optout << "\n\n=========  " << s << "  ===========\n\n";
  *optout << "Optimization method       = " << method << "\n";
  *optout << "Dimension of the problem  = " << nlp->getDim()   << "\n";
  *optout << "No. of bound constraints  = " << nlp->getDim()   << "\n";
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

real OptBCNewton::stepTolNorm() const
{
  NLP2* nlp = nlprob2();
  ColumnVector step(sx.AsDiagonal()*(nlp->getXc() - xprev));
  return Norm2(step);
}

SymmetricMatrix OptBCNewton::updateH(SymmetricMatrix&, int)
{
  return nlprob()->evalH();
}

//---------------------------------------------------------------------------- 
// 
// Compute the maximum step allowed along the search direction sk
// before we hit a constraint
//
//--------------------------------------------------------------------------- 
double OptBCNewton::computeMaxStep(ColumnVector &sk)
{
  NLP2* nlp = nlprob2();
  int i, n = nlp->getDim();
  double gamma=FLT_MAX, delta;
  ColumnVector lower = nlp->getConstraints()->getLower();
  ColumnVector upper = nlp->getConstraints()->getUpper();
  ColumnVector xc    = nlp->getXc();

  double snorm = Norm2(sk);
  double feas_tol = 1.e-3;

  for (i=1; i<=n; i++) {
    if (work_set(i) == false) {
      if      (sk(i) > 0.0e0) {
        delta = (upper(i)-xc(i)) / sk(i);
        if (delta <= feas_tol) {
	  if (debug_)
	    *optout << "Hit an upper constraint for variable " << i << "\n";
	}
      }
      else if (sk(i) < 0.0e0) {
        delta = (lower(i)-xc(i)) / sk(i);
        if (delta <= feas_tol) {
	  if (debug_)
	    *optout << "Hit a  lower constraint for variable " << i << "\n";
	}
      }
      gamma = min(gamma,delta);
    }
  }
  if (debug_)
    *optout << "computeMaxStep: maximum step allowed = " << gamma*snorm << "\n";
  return gamma*snorm;
}

ColumnVector OptBCNewton::computeSearch(SymmetricMatrix &H)
{
  NLP2*	nlp = nlprob2();
  int   i, j, ncnt=0, *index_array, n = nlp->getDim();

  ColumnVector          gg(n), sk2(n), sk(n);
  SymmetricMatrix       H1;
  LowerTriangularMatrix L;

  // set up index_array to count the number of free variables

  index_array = new int[n+1];
  for (i=1; i<=n; i++) index_array[i] = 0;
  for (i=1; i<=n; i++) 
    if (work_set(i) == false) index_array[i] = ++ncnt;
  if (ncnt != (n-nactive)) {
    *optout << "Number of fixed and free variables do not correspond. \n";
    exit(-1);
  }

  // Form the projected Hessian

  H1.ReSize(ncnt);
  for (i=1; i<=n; i++) {
     for (j=1; j<=n; j++) 
       if (index_array[i] != 0 && index_array[j] != 0)  
	  H1(index_array[i],index_array[j]) = H(i,j);
  }

  // Form the projected gradient

  gg.ReSize(ncnt,1);
  for (i=1; i<=n; i++) 
    if (index_array[i] != 0) gg(index_array[i]) = gprev(i); 

  // Solve (H1 * sk2 = - gg) for projected search direction sk2 

  L.ReSize(ncnt);
  sk2.ReSize(ncnt,1);
  if (ncnt == 1) sk2(1) = - gg(1) / H1(1,1);
  else {
    L   = MCholesky(H1);
    sk2 = -(L.t().i()*(L.i()*gg));
  }

  // Form search direction sk from from projected search direction sk2 

  for (i=1; i<=n; i++) sk(i) = 0.0;
  for (i=1; i<=n; i++) 
    if (index_array[i] != 0) sk(i) = sk2(index_array[i]);

  // Sanitation and return

  delete [] index_array; 
  return sk;
}

int OptBCNewton::updateConstraints(int step_type)
{
  NLP2*       	nlp = nlprob2();
  int          	n = nlp->getDim(), ret_flag=0;
  int          	i, j, j2, k, *new_active, actcnt=0, notnew;
  double       	reduced_grad_norm, feas_tol=1.0e-12;
  ColumnVector 	lower(n), upper(n), xc(n), gg(n);

  // initialization

  lower      = nlp->getConstraints()->getLower();
  upper      = nlp->getConstraints()->getUpper();
  xc         = nlp->getXc();
  new_active = new int[n];

// cpjw - Use current gradient info or rather gradient at xc
  gg = nlp->evalG(xc);

  // Add variables to the working set

  for (i=1; i<=n; i++) {
    if ((fabs(upper(i)-xc(i))<feas_tol) || (fabs(lower(i)-xc(i))<feas_tol)) { 
      if (work_set(i) == false) {
        new_active[actcnt++] = i; work_set(i) = true; nactive++;
        *optout << "OptBCNewton : variable added to working set : " << i << "\n";
      }
    } 
  }

  // Delete variables from the active set 
  // First compute the norm of the reduced gradient

  int    jdel=0;
  double ggdel=0;

//gg = nlp->getGrad();
  for (i=1; i<=n; i++) if(work_set(i) == true) gg(i) = 0.0;
  reduced_grad_norm = Norm2(gg);
  if (m_nconvgd > 0 || step_type < 0) {
    gg = nlp->getGrad();
    ret_flag = -1;
    *optout << "OptBCNewton : reduced_grad_norm = " << reduced_grad_norm << "\n";
    for (i=1; i<=n; i++) {
      notnew = true;
      for (j=0; j<actcnt; j++) if (new_active[j] == i) notnew = false;
      if (work_set(i) == true && notnew) 
        if (((fabs(upper(i)-xc(i))<feas_tol) && gg(i)>feas_tol) ||
            ((fabs(lower(i)-xc(i))<feas_tol) && gg(i)<(-feas_tol))) {
	 if (fabs(gg(i)) > ggdel) {jdel = i; ggdel = fabs(gg(i)); }
	}
    }
    if (jdel != 0) {
      work_set(jdel) = false; nactive--;
      *optout << "OptBCNewton : variable deleted from working set : " << jdel << "\n";
      ret_flag = 1;
    }
  }
  if (nactive > 0) *optout << "OptBCNewton: Current working set  \n";
  k = 1;
  for (i=1; i<=nactive; i+=10) {
    *optout << " ----- variable index: ";
    j2 = min(i*10,nactive);
    for (j=(i-1)*10+1; j<=j2; j++) {
      while (work_set(k) == false) k++;
      *optout << d(k,6)  << "\t" << xc(k); k++;
    }
    *optout << "\n ";
  }
  return ret_flag;
}

void OptBCNewton::reset()
{
   NLP1* nlp = nlprob();
   int   n   = nlp->getDim();
   nlp->reset();
   OptimizeClass::defaultReset(n);
   nactive  = 0;
   work_set = false;
}

} // namespace OPTPP
