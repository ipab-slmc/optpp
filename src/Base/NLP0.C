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
#include <cfloat>
#include <cstring>
#else
#include <float.h>
#include <string.h>
#endif

#ifdef WITH_MPI
#include "mpi.h"
#endif

#include "NLP0.h"
#include "TOLS.h"
#include "cblas.h"
#include "ioformat.h"

using namespace std;

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

//----------------------------------------------------------------------------
// Evaluate the Hessian using finite differences
// No analytical gradients available so use function values
//----------------------------------------------------------------------------

SymmetricMatrix NLP0::FD2Hessian(ColumnVector & sx) 
{
  Real mcheps = FloatingPointPrecision::Epsilon();
  ColumnVector fcn_accrcy = getFcnAccrcy();
  double hieps, eta;
  int i;
  int nr = getDim();

  double xtmpi, xtmpj;
  double fii, fij, fx;
  
  ColumnVector fhi(nr), step(nr);
  SymmetricMatrix H(nr);

  // do we need this??? Dougm xc = getXc();
  fx = getF();
  

  for (i=1; i<=nr; i++) {
    hieps = max(mcheps,fcn_accrcy(i));
    eta   = pow(hieps,0.333333);
    step(i) = eta*max(fabs(mem_xc(i)),sx(i));
    step(i) = copysign(step(i),mem_xc(i));
    xtmpi = mem_xc(i);
    mem_xc(i) = xtmpi + step(i);
    fhi(i) = evalF(mem_xc);
    mem_xc(i) = xtmpi;
  }
  
  for (i=1; i<=nr; i++) {
    xtmpi = mem_xc(i);
    mem_xc(i) = mem_xc(i) + step(i)*2.0;
    fii = evalF(mem_xc); 
    H(i,i) = ((fx - fhi(i)) + (fii - fhi(i))) / (step(i)*step(i));
    mem_xc(i) = xtmpi + step(i);
    for (int j=i+1; j<=nr; ++j) {
      xtmpj = mem_xc(j);
      mem_xc(j) = mem_xc(j) + step(j);
      fij = evalF(mem_xc);
      H(i,j) = ((fx - fhi(i)) + (fij - fhi(j))) / (step(i)*step(j));
      mem_xc(j) = xtmpj;
    }
    mem_xc(i) = xtmpi;
  } 
  return H;
}

// Compute gradient using backward finite differences

ColumnVector NLP0::BDGrad(const ColumnVector& sx, const ColumnVector& x,
			  double& fx, ColumnVector& grad)
{
  int i, j, gradStart, gradEnd, nBcasts;
  double xtmp, fminus, hi, hieps;

  int me = 0;
  int nprocs = 1;
  int ndim = getDim();
  const int tmpSize = (int) ceil((double) ndim/nprocs);
  double *tmpGradMinus = new double[tmpSize];
  ColumnVector xcurrent = x;
  ColumnVector fcn_accrcy = getFcnAccrcy();
  Real mcheps = FloatingPointPrecision::Epsilon();
  SpecOption SpecPass = getSpecOption();

#ifdef WITH_MPI

  int error, resultlen, flag;
  char buffer[MPI_MAX_ERROR_STRING];

  // Check to see if MPI has been initialized.

  error = MPI_Initialized(&flag);
  if (error != MPI_SUCCESS) {
    MPI_Error_string(error, buffer, &resultlen);
//    cerr << "NLP0::BDGrad: MPI Error - " << buffer << endl;
  }

  // If it has, reset me and nprocs accordingly.

  if (flag == 1) {
    error = MPI_Comm_rank(MPI_COMM_WORLD, &me);
    if (error != MPI_SUCCESS) {
      MPI_Error_string(error, buffer, &resultlen);
//      cerr << "NLP0::BDGrad: MPI Error - " << buffer << endl;
    }
    error = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if (error != MPI_SUCCESS) {
      MPI_Error_string(error, buffer, &resultlen);
//      cerr << "NLP0::BDGrad: MPI Error - " << buffer << endl;
    }
  }

#endif

  // Set loop endpoints, f, and x according to which pass of
  // speculative gradient evaluation this is.

  if (SpecPass == Spec1) {
    if (me == nprocs-1) {
      setSpecOption(NoSpec);
      fx = evalF(xcurrent);
      setSpecOption(SpecPass);
#ifdef WITH_MPI
      if (nprocs > 1)
	MPI_Bcast(&fx, 1, MPI_DOUBLE, me, MPI_COMM_WORLD);
#endif
    }
    gradStart = 1;
    gradEnd = min(ndim, nprocs-1);
    nBcasts = min(ndim, nprocs-1);
  }
  else if (SpecPass == Spec2) {
    gradStart = nprocs;
    gradEnd = ndim;
    nBcasts = min(gradEnd-gradStart+1, nprocs);
  }
  else {
    gradStart = 1;
    gradEnd = ndim;
    nBcasts = min(ndim, nprocs);
    if (SpecPass != NoSpec) {
//    cerr << "NLP0::BDGrad: Invalid speculative gradient option - "
//	 << "SpecFlag = " << SpecPass << "\n"
//	 << "Assuming NoSpec..." << endl;
    }
  }

  // Compute my piece of the gradient.

  for (i=me+gradStart; i<=gradEnd; i+=nprocs) {

    hieps = sqrt(max(mcheps, fcn_accrcy(i)));
    hi = hieps * max(fabs(xcurrent(i)), sx(i));
    hi = copysign(hi, xcurrent(i));
    xtmp = xcurrent(i);
    xcurrent(i) = xtmp - hi;

    setSpecOption(NoSpec);
    fminus = evalF(xcurrent);
    setSpecOption(SpecPass);
#ifdef WITH_MPI
    if (SpecPass == Spec1)
      MPI_Bcast(&fx, 1, MPI_DOUBLE, nprocs-1, MPI_COMM_WORLD);
#endif
    grad(i) = (fx - fminus) / hi;
    xcurrent(i) = xtmp;
  }

  // Share my piece of the gradient with everyone else, and
  // incorporate their pieces.

  if (nprocs > 1) {

    for (i=0; i<nBcasts; i++) {

      for (j=me+gradStart; j<=gradEnd; j+=nprocs)
	tmpGradMinus[(j-me-gradStart)/nprocs] = grad(j);

#ifdef WITH_MPI
      MPI_Bcast(tmpGradMinus, tmpSize, MPI_DOUBLE, i, MPI_COMM_WORLD);
#endif

      for (j=i+gradStart; j<=gradEnd; j+=nprocs)
	grad(j) = tmpGradMinus[(j-i-gradStart)/nprocs];
    }
  }

  if (tmpGradMinus != NULL)
    delete[] tmpGradMinus;

  return grad;
}

// Compute gradient using forward finite differences

ColumnVector NLP0::FDGrad(const ColumnVector& sx, const ColumnVector& x,
			  double& fx, ColumnVector& grad) 
{
  int i, j, gradStart, gradEnd, nBcasts;
  double xtmp, fplus, hi, hieps;

  int me = 0;
  int nprocs = 1;
  int ndim = getDim();
  const int tmpSize = (int) ceil((double) ndim/nprocs);
  double *tmpGradPlus = new double[tmpSize];
  ColumnVector xcurrent = x;
  ColumnVector fcn_accrcy = getFcnAccrcy();
  Real mcheps = FloatingPointPrecision::Epsilon();
  SpecOption SpecPass = getSpecOption();

#ifdef WITH_MPI

  int error, resultlen, flag;
  char buffer[MPI_MAX_ERROR_STRING];

  // Check to see if MPI has been initialized.

  error = MPI_Initialized(&flag);
  if (error != MPI_SUCCESS) {
    MPI_Error_string(error, buffer, &resultlen);
//    cerr << "NLP0::FDGrad: MPI Error - " << buffer << endl;
  }

  // If it has, reset me and nprocs accordingly.

  if (flag == 1) {
    error = MPI_Comm_rank(MPI_COMM_WORLD, &me);
    if (error != MPI_SUCCESS) {
      MPI_Error_string(error, buffer, &resultlen);
//      cerr << "NLP0::FDGrad: MPI Error - " << buffer << endl;
    }
    error = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if (error != MPI_SUCCESS) {
      MPI_Error_string(error, buffer, &resultlen);
//      cerr << "NLP0::FDGrad: MPI Error - " << buffer << endl;
    }
  }

#endif
  
  // Set loop endpoints, f, and x according to which pass of
  // speculative gradient evaluation this is.

  if (SpecPass == Spec1) {
    if (me == nprocs-1) {
      setSpecOption(NoSpec);
      fx = evalF(xcurrent);
      setSpecOption(SpecPass);
#ifdef WITH_MPI
      if (nprocs > 1)
	MPI_Bcast(&fx, 1, MPI_DOUBLE, me, MPI_COMM_WORLD);
#endif
    }
    gradStart = 1;
    gradEnd = min(ndim, nprocs-1);
    nBcasts = min(ndim, nprocs-1);
  }
  else if (SpecPass == Spec2) {
    gradStart = nprocs;
    gradEnd = ndim;
    nBcasts = min(gradEnd-gradStart+1, nprocs);
  }
  else {
    gradStart = 1;
    gradEnd = ndim;
    nBcasts = min(ndim, nprocs);
    if (SpecPass != NoSpec) {
//    cerr << "NLP0::FDGrad: Invalid speculative gradient option - "
//	 << "SpecFlag = " << SpecPass << "\n"
//	 << "Assuming NoSpec..." << endl;
    }
  }

  // Compute only my piece of the gradient.

  for (i=me+gradStart; i<=gradEnd; i+=nprocs) {

    hieps = sqrt(max(mcheps, fcn_accrcy(i)));
    hi = hieps * max(fabs(xcurrent(i)), sx(i));
    hi = copysign(hi, xcurrent(i));
    xtmp = xcurrent(i);
    xcurrent(i) = xtmp + hi;
    setSpecOption(NoSpec);
    fplus = evalF(xcurrent);
    setSpecOption(SpecPass);
#ifdef WITH_MPI
    if (SpecPass == Spec1)
      MPI_Bcast(&fx, 1, MPI_DOUBLE, nprocs-1, MPI_COMM_WORLD);
#endif
    grad(i) = (fplus - fx) / hi;
    xcurrent(i) = xtmp;
  }

  // Share my piece of the gradient with everyone else, and
  // incorporate their pieces.

  if (nprocs > 1) {

    for (i=0; i<nBcasts; i++) {

      for (j=me+gradStart; j<=gradEnd; j+=nprocs)
	tmpGradPlus[(j-me-gradStart)/nprocs] = grad(j);

#ifdef WITH_MPI
      MPI_Bcast(tmpGradPlus, tmpSize, MPI_DOUBLE, i, MPI_COMM_WORLD);
#endif

      for (j=i+gradStart; j<=gradEnd; j+=nprocs)
	grad(j) = tmpGradPlus[(j-i-gradStart)/nprocs];

    }
  }

  if (tmpGradPlus != NULL)
    delete[] tmpGradPlus;

  return grad;
}

// Compute gradient using central differences

ColumnVector NLP0::CDGrad(const ColumnVector& sx, const ColumnVector& x,
			  double& fx, ColumnVector& grad) 
{
  int i, gradStart, gradEnd, myStart, inc, nBcasts;
  double xtmp, fplus, fminus, hi, hieps;
  int j, tmpSize;

  int me = 0;
  int nprocs = 1;
  int ndim = getDim();
  ColumnVector xcurrent = x;
  ColumnVector fcn_accrcy = getFcnAccrcy();
  Real mcheps = FloatingPointPrecision::Epsilon();
  SpecOption SpecPass = getSpecOption();

#ifdef WITH_MPI

  char buffer[MPI_MAX_ERROR_STRING];
  int error, resultlen, flag;

  // Check to see if MPI has been initialized.

  error = MPI_Initialized(&flag);
  if (error != MPI_SUCCESS) {
    MPI_Error_string(error, buffer, &resultlen);
//    cerr << "NLP0::CDGrad: MPI Error - " << buffer << endl;
  }

  // If it has, reset me and nprocs accordingly.

  if (flag == 1) {
    error = MPI_Comm_rank(MPI_COMM_WORLD, &me);
    if (error != MPI_SUCCESS) {
      MPI_Error_string(error, buffer, &resultlen);
//      cerr << "NLP0::CDGrad: MPI Error - " << buffer << endl;
    }
    error = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if (error != MPI_SUCCESS) {
      MPI_Error_string(error, buffer, &resultlen);
//      cerr << "NLP0::CDGrad: MPI Error - " << buffer << endl;
    }
  }

#endif
  
  // Set loop endpoints, f, and x according to which pass of
  // speculative gradient evaluation this is.

  if (SpecPass == Spec1) {
    if (me == nprocs-1) {
      setSpecOption(NoSpec);
      fx = evalF(xcurrent);
      setSpecOption(SpecPass);
    }
    gradStart = 1;
    gradEnd = min(ndim, (int) floor((double) (nprocs-1)/2));
    if (nprocs > 1)
      inc = (int) floor((double) (nprocs-1)/2);
    else
      inc = 1;
    nBcasts = min(ndim, (int) floor((double) (nprocs-1)/2));
  }
  else if (SpecPass == Spec2) {
    gradStart = (int) ceil((double) nprocs/2);
    gradEnd = ndim;
    if (nprocs > 1)
      inc = (int) floor((double) nprocs/2);
    else
      inc = 1;
    nBcasts = min(gradEnd-gradStart+1, (int) floor((double) nprocs/2));
  }
  else {
    gradStart = 1;
    gradEnd = ndim;
    if (nprocs > 1)
      inc = (int) floor((double) nprocs/2);
    else
      inc = 1;
    nBcasts = min(ndim, (int) floor((double) nprocs/2));
    if (SpecPass != NoSpec) {
//    cerr << "NLP0::FDGrad: Invalid speculative gradient option - "
//	 << "SpecFlag = " << SpecPass << "\n"
//	 << "Assuming NoSpec..." << endl;
    }
  }

  // Compute only my piece of the gradient.

  myStart = (int) floor((double) me/2) + gradStart;

  for (i=myStart; i<=gradEnd; i+=inc) {

    hieps = max(mcheps, fcn_accrcy(i));
    hieps = pow(hieps, 0.333333);
    hi = hieps*max(fabs(xcurrent(i)), sx(i));
    hi = copysign(hi, xcurrent(i));
    xtmp = xcurrent(i);

#ifdef WITH_MPI
    if (nprocs > 1) {

      // For multiple processors, even processors look forward, and
      // odd look backward.

      if (me%2 == 0)
	xcurrent(i) = xtmp + hi;
      else
	xcurrent(i) = xtmp - hi;

      setSpecOption(NoSpec);
      grad(i) = evalF(xcurrent)/(2*hi);
      setSpecOption(SpecPass);
    }
    else {
      // Otherwise, do the same as in the serial case.
#endif
      xcurrent(i) = xtmp + hi;
      setSpecOption(NoSpec);
      fplus = evalF(xcurrent);
      setSpecOption(SpecPass);

      xcurrent(i) = xtmp - hi;
      setSpecOption(NoSpec);
      fminus = evalF(xcurrent);
      setSpecOption(SpecPass);

      grad(i)= (fplus - fminus) / (2*hi);
#ifdef WITH_MPI
    }
#endif
    xcurrent(i) = xtmp;
  }

  if (nprocs > 1) {

    if (nprocs%2 == 0)
      tmpSize = (int) ceil((double) (2*ndim)/nprocs);
    else
      // If there are an odd number of processors, the last one doesn't
      // count.
      tmpSize = (int) ceil((double) (2*ndim)/(nprocs-1));

    double *tmpGradPlus = new double[tmpSize];
    double *tmpGradMinus = new double[tmpSize];

    for (i=0; i<nBcasts; i++) {

      for (j=myStart; j<=gradEnd; j+=inc) {
	if (me%2 == 0)
	  tmpGradPlus[(j-myStart)/inc] = grad(j);
	else
	  tmpGradMinus[(j-myStart)/inc] = grad(j);
      }

#ifdef WITH_MPI
      MPI_Bcast(tmpGradPlus, tmpSize, MPI_DOUBLE, 2*i, MPI_COMM_WORLD);
      MPI_Bcast(tmpGradMinus, tmpSize, MPI_DOUBLE, (2*i)+1, MPI_COMM_WORLD);
#endif

      for (j=i+gradStart; j<=gradEnd; j+=inc)
	grad(j) = tmpGradPlus[(j-i-gradStart)/inc] - 
	                     tmpGradMinus[(j-i-gradStart)/inc];
    }
    if (tmpGradPlus != NULL)
      delete[] tmpGradPlus;
    if (tmpGradMinus != NULL)
      delete[] tmpGradMinus;
  }
  return grad;
}

// Compute gradient of nonlinear constraints using backward finite differences
Matrix NLP0::CONBDGrad(const ColumnVector& sx) 
{
  Real mcheps = FloatingPointPrecision::Epsilon();
  ColumnVector fcn_accrcy = getFcnAccrcy();
  int i, n;
  double xtmp, hi, hieps;
  ColumnVector fminus, fx;
  
  n = dim;
  Matrix grad(n,ncnln), gtmp(ncnln,n);
  fx = evalCF(mem_xc);
  //fx = getConstraintValue();

  for (i=1; i<=n; i++) {
    hieps = sqrt(max(mcheps,fcn_accrcy(i) ));
    hi = hieps*max(fabs(mem_xc(i)),sx(i));
    hi = copysign(hi,mem_xc(i));
    xtmp = mem_xc(i);
    mem_xc(i) = xtmp - hi;
    fminus = evalCF(mem_xc);
    gtmp.Column(i) = (fx - fminus) / hi;
    mem_xc(i) = xtmp;
  }
  grad = gtmp.t();
  return grad;
}

// Compute gradient of nonlinear constraints using forward finite differences
Matrix NLP0::CONFDGrad(const ColumnVector& sx) 
{
  Real mcheps = FloatingPointPrecision::Epsilon();
  ColumnVector  fcn_accrcy = getFcnAccrcy();
  int i, n;
  double xtmp, hi, hieps;
  ColumnVector fx, fplus;
  
  n = dim;
  ColumnVector xcurrent(n);
  Matrix grad(n,ncnln), gtmp(ncnln,n);
  xcurrent = getXc();
  fx = evalCF(xcurrent);
  //fx = getConstraintValue();

  for (i=1; i<=n; i++) {
    hieps = sqrt(max(mcheps,fcn_accrcy(i) ));
    hi = hieps*max(fabs(xcurrent(i)),sx(i));
    hi = copysign(hi,xcurrent(i));
    xtmp = xcurrent(i);
    xcurrent(i) = xtmp + hi;
    fplus = evalCF(xcurrent);
    gtmp.Column(i) = (fplus - fx) / hi;
    xcurrent(i) = xtmp;
  }
  grad = gtmp.t();
  return grad;
}

// Compute gradient of nonlinear constraints using central differences
Matrix NLP0::CONCDGrad(const ColumnVector& sx) 
{
  Real mcheps = FloatingPointPrecision::Epsilon();
  ColumnVector fcn_accrcy = getFcnAccrcy();
  int i, n;
  double xtmp, hi, hieps; 
  ColumnVector fplus, fminus;
  
  n = dim;
  Matrix grad(n, ncnln), gtmp(ncnln,n);

  for (i=1; i<=n; i++) {

    hieps = max(mcheps,fcn_accrcy(i) );
    hieps = pow(hieps,0.333333);

    hi = hieps*max(fabs(mem_xc(i)),sx(i));
    hi = copysign(hi,mem_xc(i));

    xtmp   = mem_xc(i);
    mem_xc(i)  = xtmp + hi;
    fplus  = evalCF(mem_xc);

    mem_xc(i)  = xtmp - hi;
    fminus = evalCF(mem_xc);

    gtmp.Column(i)= (fplus - fminus) / (2*hi);
    mem_xc(i) = xtmp;
  }
  grad = gtmp.t();
  return grad;
}

//-------------------------------------------------------------------------
// Output Routines
//-------------------------------------------------------------------------

void NLP0::printState(char * s) 
{ // Print out current state: x current, gradient and Function value
  cout << "\n\n=========  " << s << "  ===========\n\n";
  cout << "\n    i\t   x  \t      grad   \t\t fcn_accrcy \n\n";
  for (int i=1; i<=dim; i++) 
    cout << d(i,5) << "\t" << e(mem_xc(i),12,4)<< "\t\t"
         << e(mem_fcn_accrcy(i),12,4) << "\n";
  cout <<"Function Value     = " << e(fvalue,12,4) << "\n";
  //cout <<"Function Accuracy  = " << e(mem_fcn_accrcy,12,4) << "\n";
  cout <<"\n\n===================================================\n\n";
}

void NLP0::fPrintState(ostream *nlpout, char * s) 
{ // Print out current state: x current, gradient and Function value
  (*nlpout) << "\n\n=========  " << s << "  ===========\n\n";
  (*nlpout) << "\n    i\t   x  \t      grad   \t\t fcn_accrcy \n\n";
  for (int i=1; i<=dim; i++) 
    (*nlpout) << d(i,5) << "\t" << e(mem_xc(i),12,4) << "\t\t"
              << e(mem_fcn_accrcy(i),12,4) << "\n";
  (*nlpout) <<"Function Value     = " << e(fvalue,12,4) << "\n";
 // (*nlpout) <<"Function Accuracy  = " << e(mem_fcn_accrcy,12,4) << "\n";
  (*nlpout) <<"\n\n===================================================\n\n";
}

void NLP0::saveState() 
{ // Save current state: x current, gradient and Function value
  cout << dim << "\n";
  for (int i=1; i<=dim; i++) 
     cout << e(mem_xc(i),24,16) << "\t" << e(mem_fcn_accrcy(i),24,16) << "\n";
  cout << e(fvalue,24,16) << "\n" 
	<< nlp_name << "\n"
	<< nfevals << "\n"
	<< is_expensive << "\n"
	<< debug_ << "\n"
	<< e(function_time,24,16) << "\n";
}

bool NLP0::hasConstraints()
{
  bool nonempty = false;
  if( constraint_ )
     nonempty = true;
  return nonempty;
}

void NLP0::printConstraints()
{
   
  constraint_->printConstraints();
  cout <<"\n\n===================================================\n\n";

}

} // namespace OPTPP
