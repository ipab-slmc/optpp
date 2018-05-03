//************************************************************************
//   Object Interface to NPSOL software package 
//************************************************************************
//------------------------------------------------------------------------
// system libraries 
//------------------------------------------------------------------------

#include "OptNPSOL.h"

#ifdef HAVE_STD
#include <cstdio>
#include <cstring>
#include <ctime>
#include <sstream>
#else
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <strstream>
#endif

#include <string>
#include <iomanip>

using namespace std;

#include "OptppFatalError.h"

using NEWMAT::ColumnVector;
using NEWMAT::Matrix;

#ifndef __xlC__
#define npsol  npsol_
#endif

//------------------------------------------------------------------------
// external subroutines referenced by this module 
//------------------------------------------------------------------------

extern "C" {
 void npsol_setup(char* string);
 void npsol(int& n, int& nclin, int& ncnln, int& nrowa, int& nrowj, int& nrowr,
            double* a, double* bl, double* bu,
            void (*funcon)(int& npsol_mode, int& ncnln, int& n, int& nrowj,
                           int* needc, double* x, double* c,
                           double* cjac, int& nstate),
            void (*funobj)(int& npsol_mode, int& n, double* x, double& objf,
                        double* objgrd, int& nstate),
            int& inform, int& iter, int* istate,
            double* c, double* cjac, double* clambda, double& objf,
            double* grad, double* r, double* x, int* iw, int& leniw,
            double* w, int& lenw );
}

namespace OPTPP {

//------------------------------------------------------------------------
// global links 
//------------------------------------------------------------------------

static char* class_name = "OptNPSOL";
static USERFCN0        fcn0;
static USERFCN1        fcn1;
static USERNLNCON0     cfcn0;
static USERNLNCON1     cfcn1;
static int             dlevel;
static int             feval_cnt=0;
static int             geval_cnt=0;
static Appl_Data_NPSOL *app;

//------------------------------------------------------------------------
// local subroutines defined at the end of this module 
//------------------------------------------------------------------------

void get_F(int&,int&,double*,double&,double*,int&);
void get_C(int&,int&,int&,int&,int*,double*,double*,double*,int&); 

//************************************************************************
// Local functions and subroutines :
//    initOpt
//    optimize
//    printStatus
//    readOptInput
//    setUpper
//    setLower
//    setX
//    get_F
//    get_C
//    setDifferenceInterval
//************************************************************************

//************************************************************************
// Destructor
//------------------------------------------------------------------------
OptNPSOL:: ~OptNPSOL()
{
    delete [] xvalue;
    delete [] grad;
    delete [] A;
    delete [] cfcn;
    delete [] cjac;
    delete [] clambda;
    delete [] hessian;
    delete [] lowerbounds;
    delete [] upperbounds;
    delete [] istate;
    delete [] iwork;
    delete [] work;
}

//************************************************************************
// Allocate workspace to various arrays
//------------------------------------------------------------------------
void OptNPSOL::allocate(int& n, int& nclin, int& ncnln,
	                int& nrowA, int& nrowJac)
{
  int ntot = n + nclin + ncnln;

  // initialize length of work arrays
  liwork = 3*n + nclin + 2*ncnln;
  lwork  = 20*n + 11*nclin + 21*ncnln +2*n*n + n*nclin + 2*n*ncnln;

  // allocate space for gradient information
  if(grad == NULL) grad = new double[n];

  // allocate space for linear constraint coefficients
  if(A == NULL) A = new double[nrowA*n];

  // allocate space for nonlinear constraint value
  if(cfcn == NULL) cfcn = new double[n];

  // allocate space for nonlinear constraint Jacobian
  if(cjac == NULL) cjac = new double[nrowJac*n];

  // allocate space for Lagrange multiplier information
  if(clambda == NULL) clambda = new double[ntot];

  // allocate space for upper triangular Cholesky factor of Hessian 
  if(hessian == NULL) hessian = new double[n*n];

  // allocate istate space
  if(istate == NULL) istate = new int[ntot];

  // allocate integer workspace
  if(iwork == NULL) iwork = new int[liwork];

  // allocate real workspace
  if(work == NULL) work = new double[lwork];

}

//************************************************************************
// Initialization before performing the actual optimization 
//------------------------------------------------------------------------
void OptNPSOL::initOpt()
{
  time_t       t;
  int          i;
  char         *c;
  ColumnVector xcol(npsol_n);

  // get date and print out header

  t = time(NULL);
  c = asctime(localtime(&t));
  *optout << "**********************************************************\n";
  *optout << "OPT++ version " << OPT_GLOBALS::OPT_VERSION << "\n";
  *optout << "Job run at " << c << "\n";
  copyright();
  *optout << "**********************************************************\n";

  // Read in OPT++ input file if it exists

  readOptInput();

  ret_code = 0;

  // set up the link so that get_F and get_C can access them

  fcn0        = usrfcn0; 
  fcn1        = usrfcn1; 
  cfcn0       = confcn0; 
  cfcn1       = confcn1; 
  dlevel      = deriv_level;

  // instantiate a storage object to improve efficiency

  if (buf_len >= 0) {
    if (buf_len == 1 && deriv_level == 0) buf_len = 6;
    application = new Appl_Data_NPSOL(buf_len);
    app = application;
  } else {
    printf("OptNPSOL : buffer length not valid (setBufferLength) \n");
    exit(EXIT_FAILURE);
  }

  // call user-provided initialization subroutine to get initial x
  if(!setXFlag){
    initfcn(npsol_n,xcol);
    for (i=0; i<npsol_n; i++) xvalue[i] = xcol(i+1);
  }

  iter_taken  = ret_code = 0;

  // set the NPSOL verification of gradient parameter 
  int v_level = -1;
  char verify_string[72];
  #ifdef HAVE_STD
    ostringstream verify_stream(verify_string);
  #else
    ostrstream verify_stream(verify_string, 72);
  #endif

  verify_stream <<"Verify Level                = " << setiosflags(ios::left)
		<< setw(41) << v_level << ends;

  #ifdef HAVE_STD 
    for (int is=0;is<verify_stream.str().size();is++) 
      verify_string[is] = verify_stream.str().c_str()[is]; 
  #endif 

  npsol_setup(verify_string);
}

//************************************************************************
// call NPSOL to optimize
//------------------------------------------------------------------------ 
void OptNPSOL::optimize()
{
  int i, inform = 0, ntot = npsol_n + npsol_nclin + npsol_ncnln;

  if (trace) *optout << class_name << ": Optimize\n";

  // if x has not been set, set them to 0.0
  if (xvalue == NULL) {
    xvalue = new double[npsol_n];
    for (i=0; i<npsol_n; i++) xvalue[i] = 0.0;
  }

  // if lowerbounds has not been set, set them to -1.0E12
  if (lowerbounds == NULL) {
    lowerbounds = new double[ntot];
    for (i=0; i<ntot; i++) lowerbounds[i] = -1.0E12;
  }

  // if upperbounds has not been set, set them to 1.0E12
  if (upperbounds == NULL) {
    upperbounds = new double[ntot];
    for (i=0; i<ntot; i++) upperbounds[i] = 1.0E12;
  }

  if (npsol_nclin == 0) lda    = 1;
  if (npsol_ncnln == 0) ldcjac = 1;

  // allocate workspace 
  allocate(npsol_n,npsol_nclin,npsol_ncnln,lda,ldcjac);

  // setup optimization parameters 
  initOpt();

  if (ret_code == 0) {
    // call npsol solution method
    npsol(npsol_n, npsol_nclin, npsol_ncnln, lda,
	  ldcjac, npsol_n, A, lowerbounds, upperbounds,
	  get_C, get_F, inform, iter_taken, istate, cfcn, cjac, clambda,
	  fvalue, grad, hessian, xvalue, iwork, liwork, work, lwork);

    setReturnCode(inform);

    // debug information
    *optout << "\n\n========= Output Summary ===========\n\n";
    *optout << "Function value    = " << fvalue << "\n";
    for (i=0; i<npsol_n; i++) 
      *optout << "  gradient " << i << " = " << grad[i] << "\n";
    for (i=0; i<npsol_n+npsol_nclin+npsol_ncnln; i++) 
      *optout << "  lower and upper bounds = " << lowerbounds[i] << " " 
	      << upperbounds[i] << "\n";
    *optout << "\n\n=========== End Summary =============\n\n";
  }
}

//************************************************************************
// A VERY simple routine for reading the optimization parameters
// We should really make this more general, but as a first pass this
// will have to do.
// 
// The input file should be of the form keyword = value
// where keyword is one of the following
// 
// max_iter      = 100
// max_step      = 100.0
// fcn_accrcy    = 1.e-9
//------------------------------------------------------------------------ 
void OptNPSOL::readOptInput() 
{
  int      p_level, v_level;
  char     token[80], ignore[80], equals[1];
  string keyword;
  string cderiv_level("deriv_level");
  string cdiff_interval("diff_interval");
  string cfcn_accrcy("fcn_accrcy");
  string cmaxiter("maxiter");
  string cbuffer_length("buf_len");
  string cprint_level("print_level");
  string cverify_level("verify_level");

  int keyword_count = 0;

  char *opt_input  = {"opt.input"};

  // Open opt.input file and check to see if we succeeded

  ifstream optin(opt_input);
  if (!optin.rdbuf()->is_open()) {
    *optout << "ReadOptInput: No opt.input file found\n";
    *optout << "ReadOptInput: Default values will be used\n";
    return;
  }

  *optout << "ReadOptInput: Reading opt.input file\n";

  *optout << "\n\n======  Summary of input file  ======\n\n";

  optin >> token;

  while (!optin.eof()) {

    keyword = token;
    keyword_count++;

    // check for difference interval 

    if (keyword == cdiff_interval) {
      optin >> equals >> diff_interval;
      char diff_string[72];
      #ifdef HAVE_STD
        ostringstream diff_stream(diff_string);
      #else
        ostrstream diff_stream(diff_string, 72);
      #endif
      diff_stream << "Difference Interval         = " << setiosflags(ios::left)
                  << setw(41) << diff_interval << ends;
      #ifdef HAVE_STD 
        for (int is=0;is<diff_stream.str().size();is++) 
          diff_string[is] = diff_stream.str().c_str()[is]; 
      #endif 
     
      npsol_setup(diff_string);

      diff_stream << "Central Difference Interval = " << setiosflags(ios::left)
              << setw(41) << diff_interval << ends;
      #ifdef HAVE_STD 
        for (int is=0;is<diff_stream.str().size();is++) 
          diff_string[is] = diff_stream.str().c_str()[is]; 
      #endif 
      npsol_setup(diff_string);
      *optout << cdiff_interval << " = " << diff_interval << "\n";
    }    

    // check for function precision  

    else if (keyword == cfcn_accrcy) {
      optin >> equals >> fcn_accrcy;
      char ftol_string[72];
      #ifdef HAVE_STD 
        ostringstream ftol_stream(ftol_string);
      #else
        ostrstream ftol_stream(ftol_string, 72);
      #endif
      ftol_stream <<   "Function Precision          = " << setiosflags(ios::left)
                  << setw(26) << fcn_accrcy << "               " << ends;
      #ifdef HAVE_STD 
        for (int is=0;is<ftol_stream.str().size();is++) 
          ftol_string[is] = ftol_stream.str().c_str()[is]; 
      #endif 
      npsol_setup(ftol_string);
      *optout << cfcn_accrcy  << " = " << fcn_accrcy << "\n";
    }    

    // check for iteration limit 

    else if (keyword == cmaxiter) {
      optin >> equals >> maxiter;
      char iter_string[72];
      #ifdef HAVE_STD 
        ostringstream iter_stream(iter_string);
      #else
        ostrstream iter_stream(iter_string, 72);
      #endif
      iter_stream << "Major Iteration Limit       = " << setiosflags(ios::left)
                  << setw(41) << maxiter << ends;
      npsol_setup(iter_string);
      *optout << cmaxiter << " = " << maxiter << "\n";
    }    

    // set print level 

    else if (keyword == cprint_level) {
      optin >> equals >> p_level;
      char plev_string[72];
      #ifdef HAVE_STD 
        ostringstream plev_stream(plev_string);
      #else
        ostrstream plev_stream(plev_string, 72);
      #endif
      plev_stream << "Print Level                 = " << setiosflags(ios::left)
                  << setw(41) << p_level << ends;
      npsol_setup(plev_string);
      *optout << cprint_level  << " = " << p_level << "\n";
    }    

    // set verify level 

    else if (keyword == cverify_level) {
      optin >> equals >> v_level;
      char verify_string[72];
      #ifdef HAVE_STD 
        ostringstream verify_stream(verify_string);
      #else
        ostrstream verify_stream(verify_string, 72);
      #endif
      verify_stream <<"Verify Level                = " << setiosflags(ios::left)
                    << setw(41) << v_level << ends;
      npsol_setup(verify_string);
      *optout << cverify_level  << " = " << v_level << "\n";
    }    

    // set derivative level 

    else if (keyword == cderiv_level) {
      optin >> equals >> dlevel;
      deriv_level = dlevel;
      char dlevel_string[72];
      #ifdef HAVE_STD 
        ostringstream dlevel_stream(dlevel_string);
      #else
        ostrstream dlevel_stream(dlevel_string, 72);
      #endif
      dlevel_stream <<"Derivative Level            = " << setiosflags(ios::left)
                    << setw(41) << dlevel << ends;
      npsol_setup(dlevel_string);
      *optout << cderiv_level  << " = " << dlevel << "\n";
    }    

    // input buffer length in order to improve efficiency

    else if (keyword == cbuffer_length) {
      optin >> equals >> buf_len;
    }    

    else {
      *optout << "Unrecognized keyword '" << keyword << "'. "
	<< "Skipping the rest of this line\n";
      optin.getline(ignore, sizeof(ignore));
    }

    optin >> token;
  }
  *optout << "ReadOptInput: opt.input file read.\n";
}

//************************************************************************
// output module information
//------------------------------------------------------------------------ 
void OptNPSOL::printStatus(char *s) 
{
  int i;

  *optout << "\n\n*************  " << s << "  *************\n\n";
  *optout << "Optimization method          = NPSOL methods" << "\n";
  *optout << "Dimension of the problem     = " << npsol_n << "\n";
  *optout << "No. of linear constraints    = " << npsol_nclin << "\n";
  *optout << "No. of nonlinear constraints = " << npsol_ncnln << "\n";
  *optout << "Return code                  = " << ret_code << "\n";
  *optout << "No. iterations taken         = " << iter_taken  << "\n";
  *optout << "Final solutions : \n";
  for (i=0; i<npsol_n; i++) 
    *optout << "  x[" << i << "] = " << xvalue[i] << "\n";
  *optout << "      Function value    = " << fvalue << "\n";
  for (i=0; i<npsol_n; i++) 
    *optout << "      grad[" << i << "] = " << grad[i] << "\n";
  for (i=0; i<npsol_n; i++) 
    *optout << "      lambda[" << i << "] = " << clambda[i] << "\n";
  *optout << "\n\n*******************************************\n\n";

  printf("\n\n************* %s *************\n\n",s);
  printf("Optimization method          = NPSOL methods\n");
  printf("Dimension of the problem     = %d \n", npsol_n);
  printf("No. of linear constraints    = %d \n", npsol_nclin);
  printf("No. of nonlinear constraints = %d \n", npsol_ncnln);
  printf("Return code                  = %d \n", ret_code);
  printf("No. of iterations taken      = %d \n", iter_taken);
  for (i=0; i<npsol_n; i++) 
    printf("  x[%2d] = %16.8f \n", i, xvalue[i]);
  printf("      Function value = %16.8f \n", fvalue);
  for (i=0; i<npsol_n; i++) 
    printf("      grad[%2d]     = %16.8f \n", i, grad[i]);
  for (i=0; i<npsol_n; i++) 
    printf("      lambda[%2d]   = %16.8f \n", i, clambda[i]);
  printf("\n\n*******************************************\n\n");
}

//************************************************************************
// set initial guess
//------------------------------------------------------------------------ 
void OptNPSOL::setX(const ColumnVector& x) 
{
  if (xvalue != NULL) delete [] xvalue;
  xvalue = new double[npsol_n];
  for (int i=0; i<npsol_n; i++) xvalue[i] = x(i+1);
  setXFlag = true;
}

//************************************************************************
// set lower bound for the control and constraints
//----------------------------------------------------------------------- 
void OptNPSOL::setLower(const ColumnVector& bounds) 
{
  int i, ntot = npsol_n + npsol_nclin + npsol_ncnln;

  if (lowerbounds != NULL) delete [] lowerbounds;
  lowerbounds = new double[ntot];
  for (i=0; i<ntot; i++) lowerbounds[i] = bounds(i+1);
}

//************************************************************************
// set upper bound for the control and constraints
//------------------------------------------------------------------------ 
void OptNPSOL::setUpper(const ColumnVector& bounds) 
{
  int i, ntot = npsol_n + npsol_nclin + npsol_ncnln;

  if (upperbounds != NULL) delete [] upperbounds;
  upperbounds = new double[ntot];
  for (i=0; i<ntot; i++) upperbounds[i] = bounds(i+1);
}

//************************************************************************
// set linear constraint coefficient matrix
//------------------------------------------------------------------------ 
void OptNPSOL::setMatrix(const Matrix& CoefficientMatrix) 
{
  int i, j, index = 0;

  if (A != NULL) delete [] A;
  A = new double[lda*npsol_n];
  for (i=0; i<lda; i++) 
     for (j=0; j<npsol_n; j++) 
         A[index++] = CoefficientMatrix(i+1,j+1);
}

//************************************************************************
// set derivative level 
//------------------------------------------------------------------------ 
void OptNPSOL::setDerivativeLevel(int& dlevel) 
{
  deriv_level = dlevel;
  char dlevel_string[72];
  #ifdef HAVE_STD
    ostringstream dlevel_stream(dlevel_string);
  #else
    ostrstream dlevel_stream(dlevel_string, 72);
  #endif
  dlevel_stream <<"Derivative Level            = " << setiosflags(ios::left)
                << setw(41) << dlevel << ends;
  #ifdef HAVE_STD 
    for (int is=0;is<dlevel_stream.str().size();is++) 
      dlevel_string[is] = dlevel_stream.str().c_str()[is]; 
  #endif 
  npsol_setup(dlevel_string);
  *optout << "Derivative Level = " << dlevel << "\n";
}


//************************************************************************
// set a fixed difference interval
//------------------------------------------------------------------------ 
void OptNPSOL::setDifferenceInterval(double &di) 
{
  diff_interval = di;
  char diff_string[72];
  #ifdef HAVE_STD
    ostringstream diff_stream(diff_string);
  #else
    ostrstream diff_stream(diff_string, 72);
  #endif
  diff_stream << "Difference Interval         = " << setiosflags(ios::left)
              << setw(41) << diff_interval << ends;
  #ifdef HAVE_STD 
    for (int is=0;is<diff_stream.str().size();is++) 
      diff_string[is] = diff_stream.str().c_str()[is]; 
  #endif 
  npsol_setup(diff_string);
  diff_stream << "Central Difference Interval = " << setiosflags(ios::left)
              << setw(41) << diff_interval << ends;
  npsol_setup(diff_string);
  *optout << "Difference Interval = " << diff_interval << "\n";
}    

//************************************************************************
// set function accuracy 
//------------------------------------------------------------------------ 
void OptNPSOL::setFcnAccrcy(double &fa) 
{
  fcn_accrcy = fa;
  char ftol_string[72];
  #ifdef HAVE_STD
    ostringstream ftol_stream(ftol_string);
  #else
    ostrstream ftol_stream(ftol_string, 72);
  #endif
      ftol_stream <<   "Function Precision          = " << setiosflags(ios::left)
                  << setw(26) << fa << "               " << ends;
  #ifdef HAVE_STD 
    for (int is=0;is<ftol_stream.str().size();is++) 
      ftol_string[is] = ftol_stream.str().c_str()[is]; 
  #endif 

  npsol_setup(ftol_string);
  *optout << "Function accuracy = " << fcn_accrcy << "\n";
}

//************************************************************************
// external function - compute function value
//------------------------------------------------------------------------ 
void get_F(int& npsol_mode,int& n,double *x,double& objf,
	double *grad,int& nstate) 
{
  int          i, optpp_mode, nn=n, result;
  ColumnVector xx(n), ggrad(n);

  // first convert control variable into ColumnVector format
  for (i=0; i<nn; i++) xx(i+1) = x[i];

  // for the case when the user function does not provide gradients
  if (fcn0 != NULL) {

    // if the function value is not found, do evaluation
    if (app != NULL){ 
      if (!app->getF(xx,objf)) { 
        fcn0(nn, xx, objf, result);
        app->update(result,nn,xx,objf);
        feval_cnt++;
      }
    } 

  // for the case when some gradient information is provided
  } else if (fcn1 != NULL) {
    optpp_mode = NLPFunction;
    ggrad = -1.0E14;

    // compute function value 
    if (npsol_mode == 0 || npsol_mode == 2) {
      if (app != NULL){ 
         if (!app->getF(xx,objf) ) { 
            fcn1(NLPFunction, nn, xx, objf, ggrad, result);
            app->update(result,nn,xx,objf);
            feval_cnt++;
	 }
      } 
    } 
    // compute both function value and gradient
    if (npsol_mode == 1 || npsol_mode == 2) {
      optpp_mode = optpp_mode | NLPGradient;
      if (app != NULL){ 
        if ( !app->getGrad(xx,ggrad) ) { 
          fcn1(NLPGradient, nn, xx, objf, ggrad, result);
          app->update(result,nn,xx,ggrad);
          geval_cnt++;
        }
      } 
    }

    // if certain gradients are not specified by users, do not modify
    // the array sent back to NPSOL
    if (optpp_mode & NLPGradient) {
      for (i=0; i<nn; i++) {
        if (ggrad(i+1) != -1.0E14) grad[i] = ggrad(i+1);
        else printf("NOTE : grad[%2d] not initialized.\n",i+1);
      }
      //for (i=0; i<nn; i++) printf("** grad[%2d] = %e \n",i+1,grad[i]);
    }
  } else {
    printf("get_F : ERROR - no user function available.\n");
    exit(EXIT_FAILURE);
  }
}

//************************************************************************
// external functions - compute constraint values 
//------------------------------------------------------------------------ 
void get_C(int& npsol_mode,int& ncnln,int& n,int& nrowj,int *needc,double *x,
	   double *c,double *cjac,int& nstate) 
{
  int          i, j, icnt, result=0, optpp_mode;
  ColumnVector xx(n), cc(ncnln), xx2(n);
  Matrix       ccjac(nrowj,n);

  // re-format the control variables before calling user function
  for (i=0; i<n; i++) xx(i+1) = x[i];

  // if only the constraint functions are given, call it and update
  if (cfcn0 != NULL) {

    optpp_mode = NLPConstraint;
    if (app != NULL) {
      // if the constraints are not found, do evaluation
      if ( !app->getConstraint(xx,cc) ){
        cfcn0(n, xx, cc, result);

        // if no constraint values are returned, flag error
        if ((result & NLPConstraint) == 0) {
          printf("get_C : ERROR - data requested but not given.\n");
          exit(EXIT_FAILURE);
        } 
        // if the appl_data structure is installed, save data
        app->update(n,xx,ncnln,cc);
       }
    }

    // re-format the constraint data to be returned to NPSOL 
    for (i=0; i<ncnln; i++) c[i] = cc(i+1);

  } else if (cfcn1 != NULL) {

  // if both constraints and Jacobian may be given, call it and update
    optpp_mode = NLPConstraint;

    // error checking (safeguard)
    if (nrowj != ncnln) OptppmathError("Error in get_C : nrowj != ncnln");

    ccjac = -1.0E14;
    xx2 = xx;

    // If only the constraint values are needed 
    if ( npsol_mode == 0 || npsol_mode == 2) {

      if (app != NULL){ 
        if (!app->getConstraint(xx,cc)){ 
	   cfcn1(NLPConstraint, n, xx2, cc, ccjac, result);
           app->update(n,xx,ncnln,cc);
	}
      }
      // re-format the constraint data 
      for (i=0; i<ncnln; i++)  c[i] = cc(i+1); 
    } 
    if ( npsol_mode == 1 || npsol_mode == 2) {

      optpp_mode = optpp_mode | NLPCJacobian;

      if (app != NULL){ 
        if (!app->getCJacobian(xx,ccjac)){ 
	   cfcn1(NLPCJacobian, n, xx2, cc, ccjac, result);
           app->update(result,n,xx,ncnln,cc,ccjac);
	}
      }
    }

    // re-format the Jacobian for NPSOL 
    if (optpp_mode & NLPCJacobian ) {
      icnt = 0;
      for (i=0; i<n; i++) 
        for (j=0; j<ncnln; j++) {
          if (ccjac(j+1,i+1) != -1.0E14) cjac[icnt] = ccjac(j+1,i+1);
            icnt++;
          }
      }
   }
}
} // namespace OPTPP
