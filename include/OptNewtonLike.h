#ifndef OptNewtonLike_h
#define OptNewtonLike_h

/*--------------------------------------------------------------------
 Copyright (c) 2001, Sandia Corporation.
 J.C. Meza, Sandia National Laboratories, meza@ca.sandia.gov
 ---------------------------------------------------------------------*/


#include "Opt.h"

using std::ostream;

namespace OPTPP {

/**
 * OptNewtonLike is the base class for Newton Methods. 
 * OptNewtonLike provides common data and functionality for
 * OptFDNewton, OptQNewton, and OptNewton methods. 
 *
 * @note Modified by P.J. Williams, pwillia@sandia.gov
 * @date Last modified 03/2007
 */

class OptNewtonLike: public OptimizeClass {
protected:
  virtual NLP1 *nlprob() const = 0;  ///< Returns pointer to an NLP1 object
  NEWMAT::ColumnVector gprev;		///< Gradient at the prev. iteration
  NEWMAT::SymmetricMatrix Hessian;	///< Current Hessian
  int grad_evals;		///< Number of gradient evaluations
  SearchStrategy strategy;	///< User-specified globalization Strategy 
  DerivOption finitediff;	///< User-specified derivative option 
  real TR_size;			///< Trust region radius
  real gradMult;		///< Gradient multiplier to compute TR_size
  int searchSize;               ///< Search pattern size for TRPDS
  void defaultAcceptStep(int, int);
  NEWMAT::ColumnVector defaultComputeSearch(NEWMAT::SymmetricMatrix& );
  bool WarmStart;

public:

 /**
  * Default Constructor
  * @see OptNewtonLike(int n)
  * @see OptNewtonLike(int n, UPDATEFCN u)
  * @see OptNewtonLike(int n, TOLS t)
  */
  OptNewtonLike() {}

 /**
  * @param n an integer argument.
  */
  OptNewtonLike(int n): 
    OptimizeClass(n), gprev(n), Hessian(n), grad_evals(0),
    strategy(TrustRegion), finitediff(ForwardDiff), TR_size(0.0),
    gradMult(0.1), searchSize(64), WarmStart(false){}

 /**
  * @param n an integer argument.
  * @param u a function pointer.
  */
  OptNewtonLike(int n, UPDATEFCN u): 
    OptimizeClass(n), gprev(n), Hessian(n),grad_evals(0),
    strategy(TrustRegion), finitediff(ForwardDiff),TR_size(0.0),
    gradMult(0.1), searchSize(64), WarmStart(false){update_fcn = u;}
 /**
  * @param n an integer argument.
  * @param t tolerance class reference.
  */
  OptNewtonLike(int n, TOLS t): 
    OptimizeClass(n,t), gprev(n), Hessian(n), grad_evals(0),
    strategy(TrustRegion), finitediff(ForwardDiff),TR_size(0.0),
    gradMult(0.1), searchSize(64), WarmStart(false){}
  
 /**
  * Destructor
  */
  virtual ~OptNewtonLike(){}

//-----------------------------------------------------------------
// Virtual functions 
//-----------------------------------------------------------------

//-----------------------------------------------------------------
// These have default values
//-----------------------------------------------------------------
  virtual void acceptStep(int k, int step_type)
    {defaultAcceptStep(k, step_type);}

  /// Compute Newton direction
  virtual NEWMAT::ColumnVector computeSearch(NEWMAT::SymmetricMatrix& H )
    {return defaultComputeSearch(H);}

  virtual void updateModel(int k, int ndim, NEWMAT::ColumnVector x)
    {OptimizeClass::defaultUpdateModel(k, ndim, x);}

//-----------------------------------------------------------------
// These have to be defined by other classes
//-----------------------------------------------------------------

  /// Check to see if algorithm satisfies the convergence criterion 
  virtual int  checkConvg();

  /// Compare the analytic gradient with the finite difference gradient
  virtual int  checkDeriv();


  /// Compute the step length along the Newton direction
  virtual int  computeStep(NEWMAT::ColumnVector sk);

  /// Initialize algorithmic parameters
  virtual void initOpt();

  /// Compute the Hessian or its approximation at the initial point
  virtual void initHessian();

  /// Initialize the size of the trust-region.  Only relevant when either 
  /// the trustregion or trustpds globalization strategies are selected 
  virtual double initTrustRegionSize() const;
 
  /// Invoke Newton's method on an unconstrained problem
  virtual void optimize();

  /// Read user-specified input options from a file
  virtual void readOptInput();

  /// Reset parameters
  virtual void reset();

//-----------------------------------------------------------------
// These are used by all derived classes
//-----------------------------------------------------------------

  /// Compare the analytic gradient with the finite difference gradient
  int checkAnalyticFDGrad();

  /**
   * @return The number of function evaluations
   */
  int getFevals() const {return fcn_evals;}

  /**
   * @return The number of gradient evaluations
   */
  int getGevals() const {return grad_evals;}

  /**
   * @return Radius of trust-region 
   */
  real getTRSize() const {return TR_size;}
  /// Set trust-region radius 
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
  void UseWarmStart(NEWMAT::SymmetricMatrix& H) {Hessian = H; WarmStart = true;}

  /**
   * @return Globalization strategy for optimization algorithms 
   */
  void setSearchStrategy(SearchStrategy s) {strategy = s;}
  /// Set globalization strategy for optimization algorithms
  SearchStrategy getSearchStrategy() const {return strategy;}

  
 /**
   * @return Type of finite difference approximation.
   */
  void setDerivOption(DerivOption d) {finitediff = d;}
  /// Set the type of finite difference routine
  DerivOption getDerivOption() const {return finitediff;}

  /**
   * @return Hessian matrix 
   */
  NEWMAT::SymmetricMatrix getHessian() const {return Hessian;}
  /// Store the current Hessian matrix
  void setHessian(NEWMAT::SymmetricMatrix& H) {Hessian = H;}
  
  /// Compute the Hessian of the objective function or its approximation at the current point
  virtual NEWMAT::SymmetricMatrix updateH(NEWMAT::SymmetricMatrix& H, int k) = 0;

  /// Print status of unconstrained Newton's method
  void printStatus(char *);

  friend int trustregion(NLP1*, ostream*, NEWMAT::SymmetricMatrix&,
			 NEWMAT::ColumnVector&, NEWMAT::ColumnVector&,
			 real&, real&, real stpmax, real stpmin);

  friend int trustpds(NLP1*, ostream*, NEWMAT::SymmetricMatrix&,
		      NEWMAT::ColumnVector&, NEWMAT::ColumnVector&,
		      real&, real&, real stpmax, real stpmin, int);

};

/**
 *  Unconstrained Newton class that accepts either an NLP or NLP2.
 */

class OptNewton1Deriv: public OptNewtonLike {
public:
  OptNewton1Deriv() {}

  OptNewton1Deriv(NLP1* p): 
    OptNewtonLike(p->getDim()), mem_nlp(p){}

  OptNewton1Deriv(NLP1* p, UPDATEFCN u): 
    OptNewtonLike(p->getDim(),u), mem_nlp(p){}

  OptNewton1Deriv(NLP1* p, TOLS t): 
    OptNewtonLike(p->getDim(),t), mem_nlp(p){}
  
  virtual ~OptNewton1Deriv(){}

private:
  NLP1* mem_nlp;
  
protected:
  /// returns an NLP1 pointer
  NLP1* nlprob() const { return mem_nlp; }
};

/**
 *  Unconstrained Newton class that requires an NLP2.
 */

class OptNewton2Deriv: public OptNewtonLike {
public:
  OptNewton2Deriv() {}

  OptNewton2Deriv(NLP2* p): 
    OptNewtonLike(p->getDim()), mem_nlp(p){}

  OptNewton2Deriv(NLP2* p, UPDATEFCN u): 
    OptNewtonLike(p->getDim(),u), mem_nlp(p){}

  OptNewton2Deriv(NLP2* p, TOLS t): 
    OptNewtonLike(p->getDim(),t), mem_nlp(p){}
  
  virtual ~OptNewton2Deriv(){}

private:
  NLP2* mem_nlp;
  
protected:
  /// returns an NLP1 pointer
  NLP1* nlprob() const { return mem_nlp; }
  /// returns an NLP2 pointer
  NLP2* nlprob2() const { return mem_nlp; }
};

} // namespace OPTPP

#endif
