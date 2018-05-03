#ifndef OptNPSOL_h
#define OptNPSOL_h

/*----------------------------------------------------------------------
  Copyright (c) 2001, Sandia Corporation.   Under the terms of Contract 
  DE-AC04-94AL85000, there is a non-exclusive license for use of this 
  work by or on behalf of the U.S. Government.
 ----------------------------------------------------------------------*/

#include "Opt.h"
#include "NLF.h"
#include "Appl_Data_NPSOL.h"

namespace OPTPP {

/**
 * The interface to NPSOL software package for nonlinear programming. 
 * For more information, see P. Gill, W. Murray, M. Saunders, and M. Wright,
 * "User's Guide for NPSOL(Version 4.0): A Fortran Package
 * for Nonlinear Programming", TR SOL 86-2, Department of Operations
 * Research, Stanford, CA. 
 *
 * @note Modified by P.J. Williams 
 */

class OptNPSOL: public OptimizeClass {

  int             npsol_n;        ///< dimension of control vector
  int             npsol_nclin;    ///< number of linear constraints
  int             npsol_ncnln;    ///< number of nonlinear constraints
  int             buf_len;        ///< buffer length for Appl_Data_NPSOL
  int             deriv_level;    ///< NPSOL input variable
  int		  lda; 	  	  ///< NPSOL input variable
  int		  ldcjac; 	  ///< NPSOL input variable
  int             liwork;         ///< NPSOL input variable
  int             lwork;          ///< NPSOL input variable
  int             maxiter;        ///< NPSOL input variable
  int             iter_taken;     ///< NPSOL output variable
  int             ret_code;       ///< NPSOL output variable
  int		  *istate;	  ///< NPSOL input array
  int		  *iwork;	  ///< NPSOL input array
  double          diff_interval;  ///< a user-specified difference interval
  double          fcn_accrcy;     ///< a user-specified function accuracy
  double          *xvalue;        ///< control variables
  double          fvalue;         ///< function value
  double          *grad;          ///< gradient of the objective function
  double	  *A;		  ///< linear constraint coefficient matrix
  double	  *cfcn;          ///< nonlinear constraint vector matrix
  double          *clambda;       ///< Lagrange multiplier
  double          *cjac;          ///< Jacobian of constraints 
  double	  *hessian;	  ///< upper triangular Chol factor of Hessian 
  double          *lowerbounds;   ///< lower bound for the constraints
  double          *upperbounds;   ///< upper bound for the constraints
  double	  *work;	  ///< NPSOL input array 
  bool		  setXFlag;	  ///< Flag which allows user to bypass initfcn
  INITFCN         initfcn;        ///< pointer to initialization function
  USERFCN0        usrfcn0;        ///< function evaluator
  USERFCN1        usrfcn1;        ///< function evaluator (with gradient)
  USERNLNCON0     confcn0;        ///< constraint evaluator
  USERNLNCON1     confcn1;        ///< constraint evaluator (with Jacobian)
  Appl_Data_NPSOL *application;   ///< data store for efficiency

public:

  //----------------------------------------------------------------------
  // Constructors & Destructors
  //----------------------------------------------------------------------

  OptNPSOL(int n, int nclin, int ncnln, USERFCN0 f,INITFCN i) : npsol_n(n), 
           npsol_nclin(nclin), npsol_ncnln(ncnln), buf_len(1), 
	   deriv_level(0), lda(nclin), ldcjac(ncnln), fvalue(1.0e30),
           setXFlag(false)  
           {initfcn = i; usrfcn0 = f; usrfcn1 = NULL; 
	    confcn0 = NULL; confcn1 = NULL; application = NULL;
	    xvalue=grad=A=cfcn=clambda=cjac=lowerbounds=upperbounds=NULL;
	    hessian=work=NULL; istate=iwork=NULL;
           }

  OptNPSOL(int n, int nclin, int ncnln, USERFCN1 f,INITFCN i) : npsol_n(n), 
           npsol_nclin(nclin), npsol_ncnln(ncnln), buf_len(1), 
	   deriv_level(0), lda(nclin), ldcjac(ncnln), fvalue(1.0e30),
           setXFlag(false) 
           {initfcn = i; usrfcn0 = NULL; usrfcn1 = f; 
	    confcn0 = NULL; confcn1 = NULL; application = NULL;
	    xvalue=grad=A=cfcn=clambda=cjac=lowerbounds=upperbounds=NULL;
	    hessian=work=NULL; istate=iwork=NULL;
	   }

  OptNPSOL(int n,int nclin,int ncnln,USERFCN0 f,INITFCN i,USERNLNCON0 c) : 
	   npsol_n(n), npsol_nclin(nclin), npsol_ncnln(ncnln), buf_len(1),
	   deriv_level(0), lda(nclin), ldcjac(ncnln), fvalue(1.0e30),
           setXFlag(false) 
	   {initfcn = i; usrfcn0 = f; usrfcn1 = NULL; 
	    confcn0 = c; confcn1 = NULL; application = NULL;
	    xvalue=grad=A=cfcn=clambda=cjac=lowerbounds=upperbounds=NULL;
	    hessian=work=NULL; istate=iwork=NULL;
	   }

  OptNPSOL(int n,int nclin,int ncnln,USERFCN1 f,INITFCN i,USERNLNCON0 c) : 
	   npsol_n(n), npsol_nclin(nclin), npsol_ncnln(ncnln), buf_len(1),
	   deriv_level(0), lda(nclin), ldcjac(ncnln), fvalue(1.0e30),
           setXFlag(false) 
	   {initfcn = i; usrfcn0 = NULL; usrfcn1 = f; 
            confcn0 = c; confcn1 = NULL; application = NULL;
	    xvalue=grad=A=cfcn=clambda=cjac=lowerbounds=upperbounds=NULL;
	    hessian=work=NULL; istate=iwork=NULL;
	   }

  OptNPSOL(int n,int nclin,int ncnln,USERFCN1 f,INITFCN i,USERNLNCON1 c) : 
	   npsol_n(n), npsol_nclin(nclin), npsol_ncnln(ncnln), buf_len(1),
	   deriv_level(0), lda(nclin), ldcjac(ncnln), fvalue(1.0e30), 
           setXFlag(false) 
	   { initfcn = i; usrfcn0 = NULL; usrfcn1 = f; 
            confcn0 = NULL; confcn1 = c; application = NULL;
	    xvalue=grad=A=cfcn=clambda=cjac=lowerbounds=upperbounds=NULL;
	    hessian=work=NULL; istate=iwork=NULL;
	   }

  ~OptNPSOL(); ///< Destructor

  //----------------------------------------------------------------------
  // internal functions
  //----------------------------------------------------------------------

  void initOpt();			///< Initialize algorithmic parameters 
  void optimize();			///< Call the optimization method 
  void allocate(int&, int&, int&,	///< Allocate storage space for arrays 
	    int&, int&);
  void readOptInput();			///< Read user-specifed input options
  void printStatus(char *);		///< Print status of opt. method
  void setX(const NEWMAT::ColumnVector&);	///< Set the control vector
  void setLower(const NEWMAT::ColumnVector&);	///< Set lower bounds on the constraints
  void setUpper(const NEWMAT::ColumnVector&);	///< Set upper bounds on the constraints
  void setMatrix(const NEWMAT::Matrix&);	///< Set linear constraint matrix
  void setDerivativeLevel(int &); 	///< Set user-specifed der interval
  void setDifferenceInterval(double &); ///< Set user-specifed diff interval
  void setFcnAccrcy(double &); 		///< Set user-specifed function accuracy
  

  //----------------------------------------------------------------------
  // other functions that need to be defined to get around the pure
  // virtual function complaints from the compiler
  //----------------------------------------------------------------------

  virtual void  acceptStep(int, int ) { }
  virtual int   checkConvg()          {return 0;}
  virtual NEWMAT::ColumnVector computeSearch(NEWMAT::SymmetricMatrix& ) 
				      {return (NEWMAT::ColumnVector) 0;}
  virtual void  updateModel(int k, int ndim, NEWMAT::ColumnVector x)
    {OptimizeClass::defaultUpdateModel(k, ndim, x);}
  virtual void reset(){;}

};

} // namespace OPTPP

#endif

