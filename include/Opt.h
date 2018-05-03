#ifndef Opt_h
#define Opt_h

/*---------------------------------------------------------------------
 Copyright (c) 2001, Sandia Corporation.  Under the terms of Contract
 DE-AC04-94AL85000, there is a non-exclusive license for use of this
 work by or on behalf of the U.S. Government.

 J.C. Meza, Sandia National Laboratories, meza@ca.sandia.gov
 ----------------------------------------------------------------------*/

#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#include <iostream>
#include <fstream>
#ifdef HAVE_STD
#include <cstring>
#else
#include <string.h>
#endif

#include "globals.h"

#include "include.h"
#include "newmatap.h"

#include "NLP.h"
#include "NLF.h"
#include "TOLS.h"

using std::cerr;
using std::cout;
using std::ifstream;
using std::ostream;
using std::endl;
using std::filebuf;
using std::flush;

namespace OPTPP {

//-------------------------------------------------------------------------
// Various Optimization methods and related support routines
//-------------------------------------------------------------------------

inline void abort_handler(int code)
{
//  if (code > 1) // code = 2 (Cntl-C signal), 0 (normal), & -1/1 (abnormal)
//    cout << "Signal Caught!" << endl;
 
  // Clean up
  //cout << flush; // flush cout or ofstream redirection
  //cerr << flush; // flush cerr or ofstream redirection
  exit(code);
}

inline void opt_default_update_model(int k, int dim, NEWMAT::ColumnVector x) {
  /*  clog << "opt_default_update_model: " 
       << "Iter =    "   << k 
       << ", dim  =    " << dim
       << ", x(1) =    " << x(1) << "\n"; */
}

int trustregion(NLP1*, ostream*, NEWMAT::SymmetricMatrix&, NEWMAT::ColumnVector&, 
		NEWMAT::ColumnVector&, real&, real&, real stpmax = 1.e3,
		real stpmin = 1.e-9);

int trustpds(NLP1*, ostream*, NEWMAT::SymmetricMatrix&, NEWMAT::ColumnVector&, 
		NEWMAT::ColumnVector&, real&, real&, real stpmax = 1.e3,
		real stpmin = 1.e-9, int searchSize = 64);

int dogleg(NLP1*, ostream*, NEWMAT::SymmetricMatrix&, NEWMAT::ColumnVector&, NEWMAT::ColumnVector&,
           NEWMAT::ColumnVector&, real&, real&, real);

int pdsstep(NLP1*, ostream*, NEWMAT::SymmetricMatrix&, NEWMAT::ColumnVector&, NEWMAT::ColumnVector&, 
	   NEWMAT::ColumnVector&, real&, real&, real, double&, bool, int);

int linesearch(NLP1*, ostream*, NEWMAT::ColumnVector&, NEWMAT::ColumnVector&, real *,
	       real stpmax = 1.e3, real stpmin = 1.e-9,
 	       int itnmax = 5, real ftol = 1.e-4, real xtol = 2.2e-16, 
	       real gtol = 0.9);

int backtrack(NLP1*, ostream*, NEWMAT::ColumnVector&, NEWMAT::ColumnVector&, real *,
	      int itnmax = 5, real ftol = 1.e-4, 
              real stpmax = 1.e3, real stpmin = 1.e-9);

int mcsrch(NLP1*, NEWMAT::ColumnVector&, ostream*, real *,
	   int itnmax = 5, real ftol = 1.e-4, real xtol = 2.2e-16, 
	   real gtol = 0.9, real stpmax = 1.e3, 
	   real stpmin = 1.e-9);

int mcstep(real *, real *, real *, real *, real *, 
	   real *, real *, real  , real  , bool *, 
	   real, real, int *);

NEWMAT::ReturnMatrix PertChol(NEWMAT::SymmetricMatrix&, NEWMAT::Real, 
                              NEWMAT::Real&);
NEWMAT::ReturnMatrix MCholesky(NEWMAT::SymmetricMatrix&);

/**
 *
 * Opt is the Base Optimization Class
 * All other Optimization classes are derived from this one
 *
 */

class OptimizeClass {

private:
  int x_optout_fd;

protected:
  /// Dimension of the problem
  int  dim;			
  /// Various tolerances assoc. with the problem
  TOLS tol;			
  /// Diagonal Scaling Matrix for x
  NEWMAT::ColumnVector sx;       	
  /// Diagonal Scaling Matrix for f
  NEWMAT::ColumnVector sfx;      	
  /// Previous iterate
  NEWMAT::ColumnVector xprev;	 	
  /// Objective function value at xprev
  real         fprev;	 	
  /// Current step direction
  NEWMAT::ColumnVector mem_step;	
  /// Length of step direction
  real         step_length;	
  virtual real stepTolNorm() const { return step_length; }

  /// What method is being used
  char method[80];  
  /// Optional message
  char mesg[80];    
  /// Return code from Optimization class
  int  ret_code;    
  /// Number of iterations taken 
  int  iter_taken;  
  /// Number of function evaluations taken
  int  fcn_evals;   
  /// Number of bactracks in linesearch  
  int  backtracks;  
  /// Print debug statements
  bool  debug_;	    
  int  trace;
  /// Compute time per iteration
  double iter_time;
  /// Total compute time 
  double total_time;

  /// User defined function to call after each nonlinear iteration
  UPDATEFCN  update_fcn;  

  filebuf file_buffer;
  /// Output file 
  ostream *optout;
  /// Output file success
  int     optout_fd;


/**
 * Provide default implementation of reset 
 */
  virtual void defaultReset(int n) 
  { 
     sfx.ReSize(n);
     sx.ReSize(n);
     xprev.ReSize(n);
     sx    = 1.0; 
     sfx   = 1.0; 
     xprev = 0.0;
     fcn_evals = backtracks = 0; 
  }

/**
 * Provide default implementation of AcceptStep 
 */
  virtual void defaultAcceptStep(int, int) 
    { if (debug_) *optout << "Optimize: AcceptStep not implemented yet.\n";};

/**
 * Provide default implementation of UpdateModel
 */
  virtual void defaultUpdateModel(int k, int ndim, NEWMAT::ColumnVector x) 
    {update_fcn(k, ndim, x);};

/**
 * Write copyright information to the screen
 */
  void copyright()
  {
    //pdh  did a quick fix here...path should not be relative
    //     for now, took out cerr statement and put else block
    //     around the rest.  
     char str[255];
     ifstream in("../../include/abbrev_copyright.h");
     if (!in){
       //        cerr << "Cannot open input file.\n";
     }
     else {
       while (in) {
	 in.getline(str,255);
	 if(in) *optout << str<< endl;
       }
       in.close();
     }
  }

public:
/**
 * Default Constructor
 * @see OptimizeClass(int n)
 * @see OptimizeClass(TOLS t)
 * @see OptimizeClass(int n, TOLS t)
 */
  OptimizeClass(): x_optout_fd(-1), dim(0), debug_(0), trace(0) {
    optout = new ostream(&file_buffer);
    file_buffer.open("OPT_DEFAULT.out", std::ios::out);
    if (!file_buffer.is_open() || !optout->good()) {
//      cout << "OptimizeClass:: Can't open default output file\n";
      optout_fd = 0;
    }
    update_fcn = &opt_default_update_model;
    tol.setDefaultTol();
  }

/**
 * @param n an integer argument
 */
  OptimizeClass(int n): x_optout_fd(-1), dim(n), sx(n), sfx(n), xprev(n),
    fcn_evals(0), backtracks(0), debug_(0), trace(0)      {
    optout = new ostream(&file_buffer);
    file_buffer.open("OPT_DEFAULT.out", std::ios::out);
    if (!file_buffer.is_open() || !optout->good()) {
//      cout << "OptimizeClass:: Can't open default output file\n";
      optout_fd = 0;
    }
    update_fcn = &opt_default_update_model;
    sx  = 1.0; sfx = 1.0; xprev = 0.0; 
    tol.setDefaultTol(); 
  }
  
/**
 * @param t a TOLS object
 */
  OptimizeClass(TOLS t): x_optout_fd(-1), dim(0), tol(t), debug_(0), trace(0){
    optout = new ostream(&file_buffer);
    file_buffer.open("OPT_DEFAULT.out", std::ios::out);
    if (!file_buffer.is_open() || !optout->good()) {
//      cout << "OptimizeClass:: Can't open default output file\n";
      optout_fd = 0;
    }
    update_fcn = &opt_default_update_model;
    sx  = 1.0; sfx = 1.0; xprev = 0.0; 
  }
  
/**
 * @param n an integer argument
 * @param t a TOLS object
 */
  OptimizeClass(int n, TOLS t): x_optout_fd(-1), dim(n), tol(t), sx(n),sfx(n),
      xprev(n), fcn_evals(0), backtracks(0), debug_(0), trace(0){
    optout = new ostream(&file_buffer);
    file_buffer.open("OPT_DEFAULT.out", std::ios::out);
    if (!file_buffer.is_open() || !optout->good()) {
//      cout << "OptimizeClass:: Can't open default output file\n";
      optout_fd = 0;
    }
      update_fcn = &opt_default_update_model;
      sx  = 1.0; sfx = 1.0; xprev = 0.0;
    }

  virtual ~OptimizeClass() { cleanup(); if (optout != NULL) delete optout;}
  void  cleanup() {optout->flush();};

// set various properties

  /// Set message 
  void setMesg(const char *s)   {strcpy(mesg,s);}     
  /// Set method of choice
  void setMethod(const char *s) {strcpy(method,s);}   
  /// Set maximum steplength
  void setMaxStep(real x) {tol.setMaxStep(x);}  
  /// Set minimum steplength
  void setMinStep(real x) {tol.setMinStep(x);}  
  /// Set step tolerance
  void setStepTol(real x) {tol.setStepTol(x);}	
  /// Set function tolerance
  void setFcnTol(real x)  {tol.setFTol(x);}	
  /// Set constraint tolerance
  void setConTol(real x)  {tol.setCTol(x);}	
  /// Set gradient tolerance
  void setGradTol(real x) {tol.setGTol(x);}	
  /// Set linesearch tolerance
  void setLineSearchTol(real x) {tol.setLSTol(x);} 
  /// Set maximum outer iterations
  void setMaxIter(int k)  {tol.setMaxIter(k);}	
  /// Set maximum backtrack iterations
  void setMaxBacktrackIter(int k)  {tol.setMaxBacktrackIter(k);} 
  /// Set maximum fcn evaluations
  void setMaxFeval(int k) {tol.setMaxFeval(k);} 
  /// Set update model 
  void setUpdateModel(UPDATEFCN u) {update_fcn = u;}

  /// Set step scale
  void setXScale(NEWMAT::ColumnVector x)  {sx = x;}  	
  /// Set function scale
  void setFcnScale(NEWMAT::ColumnVector x) {sfx = x;}  	
  /// Set function scale
  void setFcnScale(NEWMAT::Real x) {sfx = x;}		

  /// Set return code
  void setReturnCode(int val) {ret_code = val;}		

  /**
   * @return Problem dimension
   */
  int          getDim()      const {return dim;}	
  /**
   * @return Iterations 
   */
  int          getIter()     const {return iter_taken;}	
  /**
   * @return Return Code
   */
  int          getReturnCode()  const {return ret_code;}	
  /**
   * @return Previous iterate
   */
  NEWMAT::ColumnVector getXPrev()    const {return xprev;}
  /**
   * @return Scaling vector used for x
   */
  NEWMAT::ColumnVector getXScale()   const {return sx;}
  /**
   * @return Scaling vector used for f(x)
   */
  NEWMAT::ColumnVector getFcnScale() const {return sfx;}
  /**
   * @return Output Filename 
   */
  ostream*    getOutputFile() {return optout;};

  int setOutputFile(const char *filename, int append) { 

    if (x_optout_fd == -1) {  // Change the default output file
      file_buffer.close();
      if (append)
         file_buffer.open(filename, std::ios::out|std::ios::app);
      else
         file_buffer.open(filename, std::ios::out);
      if (!file_buffer.is_open() || !optout->good()) {
//	cout << "OptimizeClass::setOutputFile: Can't open " << filename << endl;
	optout_fd = 0;
      }
      else
	optout_fd = 1;
    }
    else {
//      cout << "OptimizeClass::setOutputFile: File already attached\n";
      optout_fd = 1;
    }
    return optout_fd;
  }

  int setOutputFile(int FileDescriptor) { 

    optout_fd = FileDescriptor;
//    cerr << "setOutputFile(int FileDescriptor) no longer supported.\n"
//	 << "Please use setOutputFile(const char *filename, int append)"
//	 << "or setOutputFile(ostream& fout)."
//	 << endl;
    optout_fd = 0;

    return optout_fd;
  }

  int setOutputFile(ostream& fout) { 

    optout->rdbuf(fout.rdbuf());
    if (!optout->good()) {
//      cout << "OptimizeClass::setOutputFile: Can't open file." << endl;
      optout_fd = 0;
    }
    else
      optout_fd = 1;

    return optout_fd;
  }

 /// Set debug flag to true
  void setDebug()         {debug_ = true;}   
/**
 * @return Debug parameter
 */
  bool Debug()            {return debug_;}

/**
 * @note Pure virtual functions 
 * @note Each derived class must define these functions for themselves
 */
  
  virtual void  acceptStep(int, int )  = 0;
  virtual int   checkConvg()           = 0;
  virtual NEWMAT::ColumnVector computeSearch(NEWMAT::SymmetricMatrix& ) = 0;
  virtual void  optimize()             = 0;
  virtual void  reset()                = 0;
  virtual void  readOptInput()         = 0;
  virtual void  printStatus(char *)    = 0;
  virtual void  updateModel(int, int, NEWMAT::ColumnVector)       = 0;
};

} // namespace OPTPP

#endif
