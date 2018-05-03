#ifndef OptBCNewtonLike_h
#define OptBCNewtonLike_h

/*----------------------------------------------------------------------
 Copyright (c) 2001, Sandia Corporation.
 J.C. Meza, Sandia National Laboratories, meza@ca.sandia.gov
 ----------------------------------------------------------------------*/

#include "Opt.h"
#include "NLF.h"

using std::ostream;

namespace OPTPP {

/**
 * Bound Constrained Newton abstract data classes
 * OptBCNewtonLike
 * OptBCNewton1Deriv
 * OptBCNewton2Deriv
 * OptBCNewtonLike provides common data and functionality for
 * the OptBCQNewton, OptBCFDNewton, and OptBCNewton methods
 *
 * OptBCNewtonlike implements a active set algorithm for bound
 * constrained optimization.
 *
 * @author J.C. Meza, Sandia National Laboratories, meza@ca.sandia.gov
 * @note Modified by P.J. Williams, pwillia@sandia.gov
 * @date Last modified 02/2004
 */

class OptBCNewtonLike: public OptimizeClass {
protected:
  virtual NLP1   *nlprob()   const = 0; ///< returns an NLP1 pointer

  NEWMAT::ColumnVector gprev;		///< Gradient at previous iteration
  NEWMAT::SymmetricMatrix Hessian;	///< Current Hessian 
  int grad_evals;			///< Number of gradient evaluations
  SearchStrategy strategy;		///< User-specified globalization strategy
  DerivOption finitediff;		///< User-specified derivative option  
  real TR_size;				///< Trust region radius
  real gradMult;			///< Gradient multiplier to compute TR_size
  int searchSize;               	///< Search pattern size for TRPDS
  int m_nconvgd;			///< Syncs fcn & constraint convergence
  bool WarmStart;

  void defaultAcceptStep(int, int);
  NEWMAT::ColumnVector defaultComputeSearch(NEWMAT::SymmetricMatrix& );

public:

/**
 * Default Constructor
 * @see OptBCNewtonLike(int n)
 * @see OptBCNewtonLike(int n, UPDATEFCN u)
 * @see OptBCNewtonLike(int n, TOLS t)
 */

  OptBCNewtonLike() {;}
/**
 * @param n an integer argument
 * @see OptBCNewtonLike(int n, UPDATEFCN u)
 * @see OptBCNewtonLike(int n, TOLS t)
 */
  OptBCNewtonLike(int n): 
    OptimizeClass(n), gprev(n), Hessian(n), grad_evals(0),
    strategy(LineSearch), finitediff(ForwardDiff), TR_size(0.0),
    gradMult(0.1), searchSize(64), WarmStart(false){;}
/**
 * @param n an integer argument
 * @param u a function pointer.
 * @see OptBCNewtonLike(int n)
 * @see OptBCNewtonLike(int n, UPDATEFCN u)
 */
  OptBCNewtonLike(int n, UPDATEFCN u): 
    OptimizeClass(n), gprev(n), Hessian(n),grad_evals(0),
    strategy(LineSearch), finitediff(ForwardDiff),TR_size(0.0),
    gradMult(0.1), searchSize(64), WarmStart(false)
      {update_fcn = u;}
/**
 * @param n an integer argument
 * @param t tolerance class reference.
 * @see OptBCNewtonLike(int n)
 * @see OptBCNewtonLike(int n, UPDATEFCN u)
 */
  OptBCNewtonLike(int n, TOLS t): 
    OptimizeClass(n,t), gprev(n), Hessian(n), grad_evals(0),
    strategy(LineSearch), finitediff(ForwardDiff),TR_size(0.0),
    gradMult(0.1), searchSize(64), WarmStart(false){;}
  
/**
 * Destructor 
 */
  virtual ~OptBCNewtonLike(){;}

//-------------------------------------------
// Virtual functions 
//-------------------------------------------

//-----------------------------------------
// These have default values
//-----------------------------------------
  virtual void acceptStep(int k, int step_type)
    {defaultAcceptStep(k, step_type);}

  /// Compute search direction 
  virtual NEWMAT::ColumnVector computeSearch(NEWMAT::SymmetricMatrix& H )
    {return defaultComputeSearch(H);}

  virtual void updateModel(int k, int ndim, NEWMAT::ColumnVector x)
    {OptimizeClass::defaultUpdateModel(k, ndim, x);}

  virtual void reset()
    {int ndim = nlprob()->getDim(); OptimizeClass::defaultReset(ndim);}

  /// Add and remove variables from the working set
  virtual int  updateConstraints(int ) {return 0;}

//-----------------------------------------
// These have to be defined by other classes
//-----------------------------------------
  /// Check to see if algorithm satisfies the convergence criterion 
  virtual int    checkConvg();

  /// Compare the analytic gradient with the finite difference gradient
  virtual int    checkDeriv();

  /// Compute the maximum step allowed along the search direction 
  /// before we hit a constraint
  virtual double computeMaxStep(NEWMAT::ColumnVector&) {return FLT_MAX;}

  /// Compute the step length along the search direction
  virtual int    computeStep(NEWMAT::ColumnVector sk);

  /// Compute the Hessian or its approximation at the initial point
  virtual void   initHessian();

  /// Initialize algorithmic parameters 
  virtual void   initOpt();

  /// Initialize the size of the trust-region 
  virtual double initTrustRegionSize() const;

  /// Invoke a bound constrained Newton's method 
  virtual void   optimize();

  /// Read user-specified input options from a file  
  virtual void   readOptInput();

  /// Compute the Hessian or its approximation at the current point
  virtual NEWMAT::SymmetricMatrix updateH(NEWMAT::SymmetricMatrix& H, int k) = 0;

//-----------------------------------------
// Non-virtual functions
//-----------------------------------------
  /// Check analytic gradients against finite-difference gradients
  int checkAnalyticFDGrad() ;

  /// Print status of the bound constrained Newton's method
  void 	 printStatus(char *);

  /**
   * @return The number of function evaluations
   */
  int getFevals() const {return fcn_evals;}

  /**
   * @return The number of gradient evaluations
   */
  int getGevals() const {return grad_evals;}

  /**
   * @return Type of finite difference approximation.
   */
  DerivOption getDerivOption() const {return finitediff;}
  /// Set the type of finite difference routine
  void setDerivOption(DerivOption d) {finitediff = d;}

  /**
   * @return Radius of trust-region 
   */
  real getTRSize() const {return TR_size;}
  /// Set radius of trust-region 
  void setTRSize(real delta) {TR_size = delta;}

  /**
   * @return Gradient multiplier to compute TR_size 
   */
  real getGradMult() const {return gradMult;}
  /// Set gradient multiplier which is used to compute trust-region radius 
  void setGradMult(real tau) {gradMult = tau;}

  /**
   * @return Number of points in search scheme 
   *  which is used in trustpds search strategy
   */
  int getSearchSize() const {return searchSize;}
  /// Set number of points in search scheme for trust-pds search strategy 
  void setSearchSize(int sss) {searchSize = sss;}

  bool getWarmStart() const {return WarmStart;}
  void useWarmStart(NEWMAT::SymmetricMatrix& H) {Hessian = H; WarmStart = true;}

  /**
   * @return Globalization strategy for optimization algorithm 
   */
  SearchStrategy getSearchStrategy() const {return strategy;}
  /// Set globalization strategy for optimization algorithms
  void setSearchStrategy(SearchStrategy s) {strategy = s;}

  /**
   * @return Hessian matrix 
   */
  NEWMAT::SymmetricMatrix getHessian() const {return Hessian;}
  /// Store the current Hessian matrix 
  void setHessian(NEWMAT::SymmetricMatrix& H) {Hessian = H;}

  friend int trustregion(NLP1*, ostream*, NEWMAT::SymmetricMatrix&,
			 NEWMAT::ColumnVector&, NEWMAT::ColumnVector&,
			 real&, real&, real stpmax, real stpmin);

  friend int trustpds(NLP1*, ostream*, NEWMAT::SymmetricMatrix&,
		      NEWMAT::ColumnVector&, NEWMAT::ColumnVector&,
		      real&, real&, real stpmax, real stpmin, int);
};

/**
 * Bound constrained Newton class that will take either an NLP1 or NLP2.
 */

class OptBCNewton1Deriv: public OptBCNewtonLike {
public:
  OptBCNewton1Deriv(): 
    OptBCNewtonLike(), mem_nlp(0){;}

  OptBCNewton1Deriv(NLP1* p): 
    OptBCNewtonLike(p->getDim()), mem_nlp(p){;}

  OptBCNewton1Deriv(NLP1* p, UPDATEFCN u): 
    OptBCNewtonLike(p->getDim(),u), mem_nlp(p)
    {update_fcn = u;}

  OptBCNewton1Deriv(NLP1* p, TOLS t): 
    OptBCNewtonLike(p->getDim(),t),  mem_nlp(p) {;}
  
  virtual ~OptBCNewton1Deriv(){;}


protected:
 /// returns an NLP1 pointer
  NLP1*   nlprob() const { return mem_nlp; }

private:
  NLP1*   mem_nlp;  
};

/**
 * Bound constrained Newton class that requires an NLP2.
 */

class OptBCNewton2Deriv: public OptBCNewtonLike {
public:
  OptBCNewton2Deriv(): 
    OptBCNewtonLike(), mem_nlp(0) {;}

  OptBCNewton2Deriv(NLP2* p): 
    OptBCNewtonLike(p->getDim()), mem_nlp(p){;}

  OptBCNewton2Deriv(NLP2* p, UPDATEFCN u): 
    OptBCNewtonLike(p->getDim(),u), mem_nlp(p){ update_fcn = u; }

  OptBCNewton2Deriv(NLP2* p, TOLS t): 
    OptBCNewtonLike(p->getDim(),t), mem_nlp(p) {}

protected:
  /// returns an NLP1 pointer
  NLP1*   nlprob() const {return mem_nlp;}   
  /// returns an NLP2 pointer
  NLP2*   nlprob2() const {return mem_nlp; }
  
private:
  NLP2*   mem_nlp;  
  
};

} // namespace OPTPP

#endif
