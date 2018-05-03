
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

#include "LSQNLF.h"

using namespace std;
using NEWMAT::ColumnVector;
using NEWMAT::Matrix;
using NEWMAT::SymmetricMatrix;
using NEWMAT::Real;
using NEWMAT::FloatingPointPrecision;

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

//-------------------------------------------------------------------------
// LSQNLF method routines.
// Derived from the NLP2 class.
//-------------------------------------------------------------------------
void LSQNLF::reset() // Reset parameter values  
{
  init_flag = false;    
  nfevals   = ngevals = 0; 
#ifdef WITH_MPI
  SpecFlag  = Spec1;
#else
  SpecFlag  = NoSpec;
#endif
  application.reset();
}

void LSQNLF::initFcn() // Initialize Function
{
  if (init_flag == false)  {
    init_fcn(dim, mem_xc);
    init_flag = true;
  }
  else  {
//    cout << "LSQNLF:initFcn: Warning - initialization called twice\n";
    init_fcn(dim, mem_xc);
  }
}

double LSQNLF::evalF() // Evaluate Function
{
  int result   = 0;
  double time0 = get_wall_clock_time();

  if(fcn0_v != NULL){
     if(SpecFlag == NoSpec) {
        if (!application.getLSQF(mem_xc,fvector)) {
           fcn0_v(dim, mem_xc, fvector, result, vptr);
	   application.lsq_update(NLPFunction,dim,lsqterms_,mem_xc,fvector);
           nfevals++;
           Jacobian_current = false;
        }
     } 
     else {
       SpecFlag = Spec1;
       (void) evalG();
       SpecFlag = Spec2;
     }
  }
  else if(fcn1_v != NULL){
     Matrix jac(lsqterms_,dim); 
     if (!application.getLSQF(mem_xc,fvector)) {
        fcn1_v(NLPFunction, dim, mem_xc, fvector, jac, result, vptr);
        application.lsq_update(result,dim,lsqterms_,mem_xc,fvector,jac);
        nfevals++;
        Jacobian_current = false;
     }
  } 
  else{
    //cerr << "Error: A function has not been declared. \n";
    exit(1);
  }

  fvalue = Dot(fvector,fvector);  

  setFcnResidual(fvector);

  function_time = get_wall_clock_time() - time0;
  if(debug_)
//  cout << "LSQNLF::evalF( )\n"
//       << "nfevals       = " << nfevals << "\n"
//       << "fvalue        = " << fvalue << "\n"
//       << "function time = " << function_time << "\n";
  return fvalue;

}

double LSQNLF::evalF(const ColumnVector& x) // Evaluate Function at x
{
  int result = 0;
  ColumnVector fx(lsqterms_);
  double ftmp, time0 = get_wall_clock_time();

  if (fcn0_v != NULL){
     if (SpecFlag == NoSpec) {
         if (!application.getLSQF(x,fx)) {
            fcn0_v(dim, x, fx, result, vptr);
	    application.lsq_update(NLPFunction,dim,lsqterms_,x,fx);
            nfevals++;
            Jacobian_current = false;
         }
     } 
     else {
       SpecFlag = Spec1;
       (void) evalG(x);
       fx = specLSQF;
       SpecFlag = Spec2;
     }
  } 
  else if(fcn1_v != NULL){
     Matrix gx(lsqterms_,dim);
     if (!application.getLSQF(x,fx)) {
        fcn1_v(NLPFunction, dim, x, fx, gx, result, vptr);
        application.lsq_update(result,dim,lsqterms_,x,fx,gx);
        nfevals++;
        Jacobian_current = false;
     }
  } 
  else{
    //cerr << "Error: A function has not been declared. \n";
    exit(1);
  }

  ftmp = Dot(fx,fx);  
  /*
   * In the search strategy routines, an evalF(xplus) is
   * followed with a setF(xplus) and evalG().  The setF
   * call updates the function value.  A similiar routine
   * is needed for the function residuals otherwise we would have
   * to recompute the function residuals.  My kluge is
   * to update setFcnResidual each time the function is evaluated.
   */
  setFcnResidual(fx);

  function_time = get_wall_clock_time() - time0;

  if(debug_)
//  cout << "LSQNLF::evalF(x)\n"
//       << "nfevals       = " << nfevals << "\n"
//       << "fvalue        = " << ftmp << "\n"
//       << "function time = " << function_time << "\n";
  return ftmp;

}

ColumnVector LSQNLF::evalG() 
{

  int result;

  if(fcn0_v != NULL){
    ColumnVector sx(dim);
    sx = 1.0;

    if (!application.getLSQF(mem_xc,fvector)) {
        fcn0_v(dim, mem_xc, fvector, result, vptr);
	application.lsq_update(NLPFunction,dim,lsqterms_,mem_xc,fvector);
	nfevals++;
    }
    else
	fvector = getFcnResidual();

    if(finitediff == ForwardDiff)
       Jacobian_ = LSQFDJac(sx, mem_xc, fvector, partial_jac);
    else if(finitediff == BackwardDiff)
       Jacobian_ = LSQBDJac(sx, mem_xc, fvector, partial_jac);
    else if(finitediff == CentralDiff)
       Jacobian_ = LSQCDJac(sx, mem_xc, fvector, partial_jac);
    else{
//       cout << "LSQNLF::evalG: Unrecognized difference option\n";
//       cout << "LSQNLF::evalG: Using forward difference option\n";
       Jacobian_ = LSQFDJac(sx, mem_xc, fvector, partial_jac);
    }
    mem_grad = 2*Jacobian_.t()*fvector;

  }
  else if(fcn1_v != NULL){
    if (!application.getLSQF(mem_xc,fvector) || !application.getLSQJac(mem_xc,Jacobian_)){
       int mode = NLPGradient;
       if (!application.getLSQF(mem_xc,fvector)) {
	 mode = NLPFunction | NLPGradient;
	 nfevals++;
       }
       fcn1_v(mode, dim, mem_xc, fvector, Jacobian_, result, vptr);
       application.lsq_update(result,dim,lsqterms_,mem_xc,fvector,Jacobian_);
       mem_grad = 2*Jacobian_.t()*fvector;  
       ngevals++;
    }
    else {
       /*
        * If an external package overrode the mode
        * and computed the function and gradient simultaneously
        * the if block wouldn't be activated and hence the
        * gradient would never change.  We grab the stored
        * residual values.
        */ 
      mem_grad = 2*Jacobian_.t()*getFcnResidual();
    }
  }

  Jacobian_current = true;
  return mem_grad;

}

ColumnVector LSQNLF::evalG(const ColumnVector& x) 
{
  int result = 0;
  ColumnVector fx(lsqterms_), gtmp(dim);
  Matrix gx(lsqterms_,dim);

  if(fcn0_v != NULL){
    ColumnVector sx(dim);
    sx = 1.0;

    if(SpecFlag == NoSpec){
       if(!application.getLSQF(x,specLSQF)){
         fcn0_v(dim, x, specLSQF, result, vptr);
         nfevals++;
       }
    }

    if(finitediff == ForwardDiff)
       gx = LSQFDJac(sx, x, specLSQF, partial_jac);
    else if(finitediff == BackwardDiff)
       gx = LSQBDJac(sx, x, specLSQF, partial_jac);
    else if(finitediff == CentralDiff)
       gx = LSQCDJac(sx, x, specLSQF, partial_jac);
    else{
//       cout << "LSQNLF::evalG: Unrecognized difference option\n";
//       cout << "LSQNLF::evalG: Using forward difference option\n";
       gx =  LSQFDJac(sx, x, specLSQF, partial_jac);
    }
    gtmp    = 2*gx.t()*specLSQF;
    Hessian <<  ( gx.t()*gx ) * 2.0;
  }
  else if(fcn1_v != NULL){
    if (!application.getLSQF(x,specLSQF) || !application.getLSQJac(x,gx)){
       int mode = NLPGradient;
       if (!application.getLSQF(x,specLSQF)) {
	 mode = NLPFunction | NLPGradient;
	 nfevals++;
       }
       fcn1_v(mode, dim, x, fx, gx, result, vptr);
       application.lsq_update(result,dim,lsqterms_,x,fx,gx);
       gtmp        = 2*gx.t()*fx;
       ngevals++;
    }
    else {
       /*
        * If an external package overrode the mode
        * and computed the function and gradient simultaneously
        * the if block wouldn't be activated and hence the
        * gradient would never change.  We grab the stored
        * residual values.
        */
      gtmp = 2*gx.t()*getFcnResidual();
    }
    Hessian <<  ( gx.t()*gx ) * 2.0;
  }

  Jacobian_current = true;
  return gtmp;
}

SymmetricMatrix LSQNLF::evalH() 
{
  if (!application.getLSQJac(mem_xc,Jacobian_))
    (void) evalG();
  Hessian << (Jacobian_.t()*Jacobian_)*2.0;
  return Hessian;
}

SymmetricMatrix LSQNLF::evalH(ColumnVector& x) 
{
  Matrix gx(lsqterms_,dim);

  if (!application.getLSQJac(x,gx))
    (void) evalG(x);
  return Hessian;
}

void LSQNLF::eval()
{
  (void) evalG();

  fvalue = Dot(fvector,fvector);  
  setFcnResidual(fvector);

  Hessian << (Jacobian_.t()*Jacobian_)*2.0;
}

real LSQNLF::evalLagrangian(const ColumnVector& xc , 
                          ColumnVector& multiplier,
                          const ColumnVector& type)
{
   real result = evalF(xc);
   if( hasConstraints()){
      ColumnVector resid = constraint_->evalResidual(xc);
      result  -=  Dot(resid, multiplier);
   }
   return result;
}

ColumnVector LSQNLF::evalLagrangianGradient(const ColumnVector& xc, 
                                          const ColumnVector& multiplier,
					  const ColumnVector& type) 
{
   ColumnVector grad = evalG(xc);
   if(hasConstraints()){
      ColumnVector tmult = multiplier;
      for (int i = 1; i < getNumOfCons() + 1; i++){
         if(type(i) == NLineq || type(i) == Lineq)
            tmult(i)*= -1;
      }
      grad += constraint_->evalGradient(xc)*tmult;
   }
   return grad;
}

ColumnVector LSQNLF::evalCF(const ColumnVector& x) // Evaluate Function at x
{
 
//  cerr << "Error: OPT++ does not support the Gauss Newton operator \n"
//       << "for nonlinear constraints.  Please select a different   \n"
//       << "NLF object, say an FDNLF.  "
//       << endl;
  exit(1);
  return(x);

}

Matrix LSQNLF::evalCG(const ColumnVector& x) // Evaluate the gradient at x
{

//  cerr << "Error: OPT++ does not support the Gauss Newton operator \n"
//       << "for nonlinear constraints.  Please select a different   \n"
//       << "NLF object, say an FDNLF.  "
//       << endl;
  exit(1);
  return(x);

}

SymmetricMatrix LSQNLF::evalCH(ColumnVector& x) // Evaluate the Hessian at x
{
  
//    cerr << "Error: OPT++ does not support the Gauss Newton operator \n"
//         << "for nonlinear constraints.  Please select a different   \n"
//         << "NLF object, say an FDNLF.  "
//         << endl;
  
    exit(1);
    SymmetricMatrix H(dim);
    H = 0;
    return(H);

}

OptppArray<SymmetricMatrix> LSQNLF::evalCH(ColumnVector& x, int darg) // Evaluate the Hessian at x
{
//    cerr << "Error: OPT++ does not support the Gauss Newton operator \n"
//       << "for nonlinear constraints.  Please select a different   \n"
//       << "NLF object, say an FDNLF.  "
//       << endl;
    exit(1);
    SymmetricMatrix H(dim);
    H = 0;
    OptppArray<SymmetricMatrix> HH;
    HH.append(H);
    return(HH);
  
}

void LSQNLF::evalC(const ColumnVector& x)
{
//  cerr << "Error: OPT++ does not support the Gauss Newton operator \n"
//       << "for nonlinear constraints.  Please select a different   \n"
//       << "NLF object, say an FDNLF.  "
//       << endl;
  exit(1);
}

// Compute Jacobian of function vector using backward finite differences
Matrix LSQNLF::LSQBDJac(const ColumnVector& sx, const ColumnVector& xc,
	                ColumnVector& fx, Matrix& jac)
{
  int i, j, k, jacStart, jacEnd, nBcasts;
  double xtmp, hi, hieps;
  ColumnVector fminus(lsqterms_); 

  int me = 0;
  int nprocs = 1;
  int n = getDim(), result = 0;
  const int tmpSize = (int) ceil((double) n/nprocs);
  double *tmpJacMinus = new double[tmpSize*lsqterms_];
  double *tmpF = new double[lsqterms_];
  ColumnVector fcn_accrcy = getFcnAccrcy();
  ColumnVector xcurrent   = xc;
  Real mcheps = FloatingPointPrecision::Epsilon();
  SpecOption SpecPass = getSpecOption();

#ifdef WITH_MPI

  int error, resultlen, flag;
  char buffer[MPI_MAX_ERROR_STRING];

  // Check to see if MPI has been initialized.

  error = MPI_Initialized(&flag);
  if (error != MPI_SUCCESS) {
    MPI_Error_string(error, buffer, &resultlen);
//    cerr << "LSQNLF::LSQBDJac: MPI Error - " << buffer << endl;
  }

  // If it has, reset me and nprocs accordingly.

  if (flag == 1) {
    error = MPI_Comm_rank(MPI_COMM_WORLD, &me);
    if (error != MPI_SUCCESS) {
      MPI_Error_string(error, buffer, &resultlen);
//      cerr << "LSQNLF::LSQBDJac: MPI Error - " << buffer << endl;
    }
    error = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if (error != MPI_SUCCESS) {
      MPI_Error_string(error, buffer, &resultlen);
//      cerr << "LSQNLF::LSQBDJac: MPI Error - " << buffer << endl;
    }
  }

#endif

  // Set loop endpoints, f, and x according to which pass of
  // speculative Jacobian evaluation this is.

  if (SpecPass == Spec1) {
    if (me == nprocs-1) {
      fcn0_v(n, xcurrent, fx, result, vptr);
#ifdef WITH_MPI
      if (nprocs > 1) {
	for (i=0; i<lsqterms_; i++)
	  tmpF[i] = fx(i+1);
	MPI_Bcast(tmpF, lsqterms_, MPI_DOUBLE, me, MPI_COMM_WORLD);
      }
#endif
    }
    jacStart = 1;
    jacEnd = min(n, nprocs-1);
    nBcasts = min(n, nprocs-1);
  }
  else if (SpecPass == Spec2) {
    jacStart = nprocs;
    jacEnd = n;
    nBcasts = min(jacEnd-jacStart+1, nprocs);
  }
  else {
    jacStart = 1;
    jacEnd = n;
    nBcasts = min(n, nprocs);
    if (SpecPass != NoSpec) {
//    cerr << "LSQNLF::LSQBDJac: Invalid speculative Jacobian option - "
//	 << "SpecFlag = " << SpecPass << "\n"
//	 << "Assuming NoSpec..." << endl;
    }
  }

  // Compute only my piece of the Jacobian.

  for (i=me+jacStart; i<=jacEnd; i+=nprocs) {

    hieps = sqrt(max(mcheps,fcn_accrcy(i)));
    hi    = hieps*max(fabs(xcurrent(i)),sx(i));
    hi    = copysign(hi,xcurrent(i));

    xtmp          = xcurrent(i);
    xcurrent(i)   = xtmp - hi;

    fcn0_v(n, xcurrent, fminus, result, vptr);
#ifdef WITH_MPI
    if (SpecPass == Spec1) {
      MPI_Bcast(tmpF, lsqterms_, MPI_DOUBLE, nprocs-1, MPI_COMM_WORLD);
      for (j=0; j<lsqterms_; j++)
	fx(j+1) = tmpF[j];
    }
#endif
    jac.Column(i) = (fx - fminus) / hi;
    xcurrent(i)   = xtmp;
  }

  // Share my piece of the Jacobian with everyone else, and
  // incorporate their pieces.

  if (nprocs > 1) {

    for (i=0; i<nBcasts; i++) {

      for (j=me+jacStart; j<=jacEnd; j+=nprocs) {
	for (k=0; k<lsqterms_; k++)
	  tmpJacMinus[k+lsqterms_*((j-me-jacStart)/nprocs)] = jac(k+1,j);
      }

#ifdef WITH_MPI
      MPI_Bcast(tmpJacMinus, tmpSize*lsqterms_, MPI_DOUBLE, i, MPI_COMM_WORLD);
#endif

      for (j=i+jacStart; j<=jacEnd; j+=nprocs) {
	for (k=0; k<lsqterms_; k++)
	  jac(k+1,j) = tmpJacMinus[k+lsqterms_*((j-i-jacStart)/nprocs)];
      }
    }
  }

  if (tmpJacMinus != NULL)
    delete[] tmpJacMinus;
  if (tmpF != NULL)
    delete[] tmpF;

  return jac;
}

// Compute Jacobian of function vector using forward finite differences
Matrix LSQNLF::LSQFDJac(const ColumnVector& sx, const ColumnVector& xc,
	                ColumnVector& fx, Matrix& jac)
{
  int i, j, k, jacStart, jacEnd, nBcasts;
  double xtmp, hi, hieps;
  ColumnVector fplus(lsqterms_); 

  int me = 0;
  int nprocs = 1;
  int n = getDim(), result = 0;
  const int tmpSize = (int) ceil((double) n/nprocs);
  double *tmpJacPlus = new double[tmpSize*lsqterms_];
  double *tmpF = new double[lsqterms_];
  ColumnVector fcn_accrcy = getFcnAccrcy();
  ColumnVector xcurrent   = xc;
  Real mcheps = FloatingPointPrecision::Epsilon();
  SpecOption SpecPass = getSpecOption();

#ifdef WITH_MPI

  int error, resultlen, flag;
  char buffer[MPI_MAX_ERROR_STRING];

  // Check to see if MPI has been initialized.

  error = MPI_Initialized(&flag);
  if (error != MPI_SUCCESS) {
    MPI_Error_string(error, buffer, &resultlen);
//    cerr << "LSQNLF::LSQFDJac: MPI Error - " << buffer << endl;
  }

  // If it has, reset me and nprocs accordingly.

  if (flag == 1) {
    error = MPI_Comm_rank(MPI_COMM_WORLD, &me);
    if (error != MPI_SUCCESS) {
      MPI_Error_string(error, buffer, &resultlen);
//      cerr << "LSQNLF::LSQFDJac: MPI Error - " << buffer << endl;
    }
    error = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if (error != MPI_SUCCESS) {
      MPI_Error_string(error, buffer, &resultlen);
//      cerr << "LSQNLF::LSQFDJac: MPI Error - " << buffer << endl;
    }
  }

#endif

  // Set loop endpoints, f, and x according to which pass of
  // speculative Jacobian evaluation this is.

  if (SpecPass == Spec1) {
    if (me == nprocs-1) {
      fcn0_v(n, xcurrent, fx, result, vptr);
#ifdef WITH_MPI
      if (nprocs > 1) {
	for (i=0; i<lsqterms_; i++)
	  tmpF[i] = fx(i+1);
	MPI_Bcast(tmpF, lsqterms_, MPI_DOUBLE, me, MPI_COMM_WORLD);
      }
#endif
    }
    jacStart = 1;
    jacEnd = min(n, nprocs-1);
    nBcasts = min(n, nprocs-1);
  }
  else if (SpecPass == Spec2) {
    jacStart = nprocs;
    jacEnd = n;
    nBcasts = min(jacEnd-jacStart+1, nprocs);
  }
  else {
    jacStart = 1;
    jacEnd = n;
    nBcasts = min(n, nprocs);
    if (SpecPass != NoSpec) {
//    cerr << "LSQNLF::LSQFDJac: Invalid speculative Jacobian option - "
//	 << "SpecFlag = " << SpecPass << "\n"
//	 << "Assuming NoSpec..." << endl;
    }
  }

  // Compute only my piece of the Jacobian.

  for (i=me+jacStart; i<=jacEnd; i+=nprocs) {

    hieps = sqrt(max(mcheps,fcn_accrcy(i)));
    hi    = hieps*max(fabs(xcurrent(i)),sx(i));
    hi    = copysign(hi,xcurrent(i));

    xtmp          = xcurrent(i);
    xcurrent(i)   = xtmp + hi;

    fcn0_v(n, xcurrent, fplus, result, vptr);
#ifdef WITH_MPI
    if (SpecPass == Spec1) {
      MPI_Bcast(tmpF, lsqterms_, MPI_DOUBLE, nprocs-1, MPI_COMM_WORLD);
      for (j=0; j<lsqterms_; j++)
	fx(j+1) = tmpF[j];
    }
#endif
    jac.Column(i) = (fplus - fx) / hi;
    xcurrent(i)   = xtmp;
  }

  // Share my piece of the Jacobian with everyone else, and
  // incorporate their pieces.

  if (nprocs > 1) {

    for (i=0; i<nBcasts; i++) {

      for (j=me+jacStart; j<=jacEnd; j+=nprocs) {
	for (k=0; k<lsqterms_; k++)
	  tmpJacPlus[k+lsqterms_*((j-me-jacStart)/nprocs)] = jac(k+1,j);
      }

#ifdef WITH_MPI
      MPI_Bcast(tmpJacPlus, tmpSize*lsqterms_, MPI_DOUBLE, i, MPI_COMM_WORLD);
#endif

      for (j=i+jacStart; j<=jacEnd; j+=nprocs) {
	for (k=0; k<lsqterms_; k++)
	  jac(k+1,j) = tmpJacPlus[k+lsqterms_*((j-i-jacStart)/nprocs)];
      }
    }
  }

  if (tmpJacPlus != NULL)
    delete[] tmpJacPlus;
  if (tmpF != NULL)
    delete[] tmpF;

  return jac;
}

// Compute Jacobian of function vector using central differences
Matrix LSQNLF::LSQCDJac(const ColumnVector& sx, const ColumnVector& xc,
	                ColumnVector& fx, Matrix& jac)
{
  int i, j, k, tmpSize, jacStart, jacEnd, myStart, inc, nBcasts;
  double xtmp, hi, hieps;
  ColumnVector fplus(lsqterms_), fminus(lsqterms_); 

  int me = 0;
  int nprocs = 1;
  int n = getDim(), result = 0;
  ColumnVector fcn_accrcy = getFcnAccrcy();
  ColumnVector xcurrent   = xc;
  Real mcheps = FloatingPointPrecision::Epsilon();
  SpecOption SpecPass = getSpecOption();
#ifdef WITH_MPI

  int error, resultlen, flag;
  char buffer[MPI_MAX_ERROR_STRING];

  // Check to see if MPI has been initialized.

  error = MPI_Initialized(&flag);
  if (error != MPI_SUCCESS) {
    MPI_Error_string(error, buffer, &resultlen);
//    cerr << "LSQNLF::LSQCDJac: MPI Error - " << buffer << endl;
  }

  // If it has, reset me and nprocs accordingly.

  if (flag == 1) {
    error = MPI_Comm_rank(MPI_COMM_WORLD, &me);
    if (error != MPI_SUCCESS) {
      MPI_Error_string(error, buffer, &resultlen);
//      cerr << "LSQNLF::LSQCDJac: MPI Error - " << buffer << endl;
    }
    error = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if (error != MPI_SUCCESS) {
      MPI_Error_string(error, buffer, &resultlen);
//      cerr << "LSQNLF::LSQCDJac: MPI Error - " << buffer << endl;
    }
  }

#endif

  // Set loop endpoints, f, and x according to which pass of
  // speculative Jacobian evaluation this is.

  if (SpecPass == Spec1) {
    if (me == nprocs-1) {
      fcn0_v(n, xcurrent, fx, result, vptr);
    }
    jacStart = 1;
    jacEnd = min(n, (int) floor((double) (nprocs-1)/2));
    if (nprocs > 1)
      inc = (int) floor((double) (nprocs-1)/2);
    else
      inc = 1;
    nBcasts = min(n, (int) floor((double) (nprocs-1)/2));
  }
  else if (SpecPass == Spec2) {
    jacStart = (int) ceil((double) nprocs/2);
    jacEnd = n;
    if (nprocs > 1)
      inc = (int) floor((double) nprocs/2);
    else
      inc = 1;
    nBcasts = min(jacEnd-jacStart+1, (int) floor((double) nprocs/2));
  }
  else {
    jacStart = 1;
    jacEnd = n;
    if (nprocs > 1)
      inc = (int) floor((double) nprocs/2);
    else
      inc = 1;
    nBcasts = min(n, (int) floor((double) nprocs/2));
    if (SpecPass != NoSpec) {
//    cerr << "LSQNLF::LSQCDJac: Invalid speculative Jacobian option - "
//	 << "SpecFlag = " << SpecPass << "\n"
//	 << "Assuming NoSpec..." << endl;
    }
  }

  // Compute only my piece of the Jacobian.

  myStart = (int) floor((double) me/2) + jacStart;

  for (i=myStart; i<=jacEnd; i+=inc) {

    hieps = max(mcheps,fcn_accrcy(i) );
    hieps = pow(hieps,0.333333);
    hi    = hieps*max(fabs(xcurrent(i)),sx(i));
    hi    = copysign(hi,xcurrent(i));

    xtmp  = xcurrent(i);

#ifdef WITH_MPI
    if (nprocs > 1) {

      // For multiple processors, even processors look forward, and
      // odd look backward.

      if (me%2 == 0)
	xcurrent(i) = xtmp + hi;
      else
	xcurrent(i) = xtmp - hi;

      fcn0_v(n, xcurrent, fplus, result, vptr);
      jac.Column(i) = fplus/(2*hi);
    }
    else {
      // Otherwise, do the same as in the serial case.
#endif

      xcurrent(i) = xtmp + hi;
      fcn0_v(n, xcurrent, fplus, result, vptr);

      xcurrent(i) = xtmp - hi;
      fcn0_v(n, xcurrent, fminus, result, vptr);

      jac.Column(i) = (fplus - fminus) / (2*hi);
#ifdef WITH_MPI
    }
#endif
    xcurrent(i) = xtmp;
  }

  if (nprocs > 1) {

    if (nprocs%2 == 0)
      tmpSize = (int) ceil((double) (2*n)/nprocs);
    else
      // If there are an odd number of processors, the last one doesn't
      // count.
      tmpSize = (int) ceil((double) (2*n)/(nprocs-1));

    double *tmpJacPlus = new double[tmpSize*lsqterms_];
    double *tmpJacMinus = new double[tmpSize*lsqterms_];

    for (i=0; i<nBcasts; i++) {

      for (j=myStart; j<=jacEnd; j+=inc) {
	for (k=0; k<lsqterms_; k++) {
	  if (me%2 == 0)
	    tmpJacPlus[k+lsqterms_*((j-myStart)/inc)] = jac(k+1,j);
	  else
	    tmpJacMinus[k+lsqterms_*((j-myStart)/inc)] = jac(k+1,j);
	}
      }

#ifdef WITH_MPI
      MPI_Bcast(tmpJacPlus, tmpSize*lsqterms_, MPI_DOUBLE, 2*i, MPI_COMM_WORLD);
      MPI_Bcast(tmpJacMinus, tmpSize*lsqterms_, MPI_DOUBLE, (2*i)+1, MPI_COMM_WORLD);
#endif

      for (j=i+jacStart; j<=jacEnd; j+=inc) {
	for (k=0; k<lsqterms_; k++)
	  jac(k+1,j) = tmpJacPlus[k+lsqterms_*((j-i-jacStart)/inc)] - 
	               tmpJacMinus[k+lsqterms_*((j-i-jacStart)/inc)];
      }
    }

    if (tmpJacPlus != NULL)
      delete[] tmpJacPlus;
    if (tmpJacMinus != NULL)
      delete[] tmpJacMinus;
  }

  return jac;
}

} // namespace OPTPP
