#ifndef OptConstrNewtonLike_h
#define OptConstrNewtonLike_h

/*--------------------------------------------------------------------
  Copyright (c) 2001, Sandia Corporation.
  J.C. Meza, Sandia National Laboratories, meza@ca.sandia.gov
 -----------------------------------------------------------------------*/

#include "Opt.h"

using std::ostream;

namespace OPTPP {

/**
 * Constrained Newton abstract data classes
 * OptConstrNewtonLike
 * OptConstrNewton1Deriv
 * OptConstrNewton2Deriv
 *
 * OptConstrNewtonLike provides common data and functionality for
 * OptConstrFDNewton, OptConstrQNewton, OptConstrNewton, OptFDNIPS,
 * OptQNIPS, and OptNIPS methods.
 *
 * @author J.C. Meza, Lawrence Berkeley National Laboratory
 * @note Modified by P.J. Williams
 * @date 03/27/2007 
 */


class OptConstrNewtonLike: public OptimizeClass {
protected:
  virtual NLP1* nlprob() const = 0; ///< returns an NLP1 pointer
  int me;			///< Number of equality constraints
  int mi;			///< Number of inequality constraints
  int grad_evals;		///< Number of gradient evaluations 
  NEWMAT::ColumnVector gprev;		///< Gradient at previous iteration
  NEWMAT::ColumnVector z;		///< Lagrange mult. corr. to ineq constraints 
  NEWMAT::ColumnVector y;		///< Lagrange mult. corr. to eq constraints
  NEWMAT::ColumnVector s;		///< Slack variables corr. to ineq. constraints 
  NEWMAT::ColumnVector constrType;	///< Vector which contains the type of constraint
  NEWMAT::ColumnVector constraintResidual; ///< Constraint residual at xc
  NEWMAT::ColumnVector gradl;		///< Current gradient of the lagrangian
  NEWMAT::ColumnVector gradlprev;	///< Previous gradient of the lagrangian
  NEWMAT::Matrix constraintGradient;	///< Current constraint gradient 
  NEWMAT::Matrix constraintGradientPrev; ///< Previous constraint gradient
  NEWMAT::SymmetricMatrix Hessian;	///< Current Hessian
  NEWMAT::SymmetricMatrix hessl;	///< Current Hessian of the lagrangian
  SearchStrategy strategy;	///< User-specified globalization strategy
  DerivOption finitediff;	///< User-specified derivative option
  MeritFcn    mfcn;		///< Merit function option
  real TR_size;			///< Size of the trust region radius 
  real gradMult;		///< Gradient multiplier to compute TR_size
  int searchSize;               ///< Search pattern size for TRPDS
  real cost;			///< Value of the merit function
  void defaultAcceptStep(int, int);
  NEWMAT::ColumnVector defaultComputeSearch(NEWMAT::SymmetricMatrix& );
  bool WarmStart;
  bool feas_flag;		///< Switch to turn on the feasibility recovery method
  int max_feas_iter;            ///< Maximize number of iterations to achieve feasibility 

public:
 /**
  * Default Constructor
  * @see OptConstrNewtonLike(int n)
  * @see OptConstrNewtonLike(int n, UPDATEFCN u)
  * @see OptConstrNewtonLike(int n, TOLS t)
  */

  OptConstrNewtonLike() {}
 /**
  * @param n an integer argument
  * @see OptConstrNewtonLike(int n, UPDATEFCN u)
  * @see OptConstrNewtonLike(int n, TOLS t)
  */
  OptConstrNewtonLike(int n): 
    OptimizeClass(n), me(0), mi(0), grad_evals(0),   gprev(n), 
    z(n), y(n), s(n), constrType(n), constraintResidual(n), 
    gradl(n), gradlprev(n),constraintGradient(n,n), constraintGradientPrev(n,n),
    Hessian(n), hessl(n), strategy(TrustRegion), finitediff(ForwardDiff), 
    mfcn(ArgaezTapia), TR_size(0.0), 
    gradMult(0.1), searchSize(64), cost(0.0), WarmStart(false),
    feas_flag(false), max_feas_iter(3)
    {z = 0; y = 0; s = 0;}
 /**
  * @param n an integer argument
  * @param u a function pointer.
  * @see OptConstrNewtonLike(int n)
  * @see OptConstrNewtonLike(int n, TOLS t)
  */
  OptConstrNewtonLike(int n, UPDATEFCN u): 
    OptimizeClass(n), me(0), mi(0), grad_evals(0),   gprev(n), 
    z(n), y(n), s(n), constrType(n), constraintResidual(n), 
    gradl(n), gradlprev(n),constraintGradient(n,n), constraintGradientPrev(n,n),
    Hessian(n), hessl(n), strategy(TrustRegion), finitediff(ForwardDiff), 
    mfcn(ArgaezTapia), TR_size(0.0), 
    gradMult(0.1), searchSize(64), cost(0.0), WarmStart(false),
    feas_flag(false), max_feas_iter(3)
    {update_fcn = u; z = 0; y = 0; s = 0;}
 /**
  * @param n an integer argument
  * @param t tolerance class reference.
  * @see OptConstrNewtonLike(int n)
  * @see OptConstrNewtonLike(int n, UPDATEFCN u)
  */
  OptConstrNewtonLike(int n, TOLS t): 
    OptimizeClass(n,t), me(0), mi(0), grad_evals(0),   gprev(n), 
    z(n), y(n), s(n), constrType(n), constraintResidual(n), 
    gradl(n), gradlprev(n),constraintGradient(n,n), constraintGradientPrev(n,n),
    Hessian(n), hessl(n), strategy(TrustRegion), finitediff(ForwardDiff), 
    mfcn(ArgaezTapia), TR_size(0.0), 
    gradMult(0.1), searchSize(64), cost(0.0), WarmStart(false),
    feas_flag(false), max_feas_iter(3)
    {z = 0; y = 0; s = 0;}
  
 /**
  * Destructor
  */
  virtual ~OptConstrNewtonLike(){}

 //---------------------------------------
 // Virtual functions 
 //---------------------------------------

// These have default values
  /// Accept this step and update the nonlinear model 
  virtual void acceptStep(int k, int step_type)
    {defaultAcceptStep(k, step_type);}

  /// Solve for the Newton direction
  virtual NEWMAT::ColumnVector computeSearch(NEWMAT::SymmetricMatrix& H )
    {return defaultComputeSearch(H);}

  virtual void updateModel(int k, int ndim, NEWMAT::ColumnVector x)
    {OptimizeClass::defaultUpdateModel(k, ndim, x);}

  virtual void reset();

 // These have to be defined by other classes

  /// Check to see if algorithm satisfies the convergence criterion
  virtual int  checkConvg();

  /// Compare the analytic gradient with the finite-difference approximation 
  virtual int  checkDeriv();

  /// Compute the steplength along the search direction.  
  /// If an acceptable step not found, returns an error code = -1. 
  virtual int  computeStep(NEWMAT::ColumnVector sk);

  virtual real computeMaxStep(NEWMAT::ColumnVector sk);

  /// Initialize algorithmic parameters 
  virtual void initOpt();

  /// Compute the Hessian or its approximation at the initial point
  virtual void initHessian();

  /// Initialize the radius of the trust-region 
  virtual double initTrustRegionSize() const;

  /// Invoke a constrained Newton's method 
  virtual void optimize();

  /// Read user-specified input options from a file
  virtual void readOptInput();

// These are used by all derived classes

  /// Check finite difference approximation with analytic derivatives
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
   * @return Trust region radius 
   */
  real getTRSize() const {return TR_size;}
  /// Set radius of trust-region
  void setTRSize(real delta) {TR_size = delta;}

  /**
   * @return Gradient multiplier to compute radius of trust-region 
   */
  real getGradMult() const {return gradMult;}
  /// Set gradient multiplier which is used to compute radius of trust-region 
  void setGradMult(real tau) {gradMult = tau;}

  /**
   * @return Number of points in search scheme 
   * (only relevant when trustpds search strategy is selected)
   */
  int getSearchSize() const {return searchSize;}
  /// Set number of points in search scheme for trustpds strategy
  void setSearchSize(int sss) {searchSize = sss;}

  bool getWarmStart() const {return WarmStart;}
  void UseWarmStart(NEWMAT::SymmetricMatrix& H) {Hessian = H; WarmStart = true;}

  /**
   * @return Globalization strategy for optimization algorithm 
   */
  void setSearchStrategy(SearchStrategy search) {strategy = search;}
  /// Set globalization strategy for optimization algorithm
  SearchStrategy getSearchStrategy() const {return strategy;}

  /**
   * @return Type of finite difference approximation
   */
  void setDerivOption(DerivOption d) {finitediff = d;}
  /// Set type of finite difference approximation (forward, backward, or central)
  DerivOption getDerivOption() const {return finitediff;}

  /**
   * @return Lagrangian multipliers corresponding to equality constraints 
   */
  NEWMAT::ColumnVector getEqualityMultiplier() const { return y;}
  /// Store the current Lagrangian multipliers corresponding to equality constraints
  void setEqualityMultiplier(const NEWMAT::ColumnVector& ymult) { y  = ymult;}

  /**
   * @return Lagrangian multipliers corresponding to inequality constraints 
   */
  NEWMAT::ColumnVector getInequalityMultiplier() const { return z;}
  /// Store the current Lagrangian multipliers corresponding to inequality constraints
  void setInequalityMultiplier(const NEWMAT::ColumnVector& zmult){ z = zmult;}

  /**
   * @return Slack variables associated with inequality constraints 
   */
  NEWMAT::ColumnVector getSlacks() const { return s;}
  /// Store the current slack vector 
  void setSlacks(const NEWMAT::ColumnVector& slackVar) { s = slackVar;}

  /**
   * @return Merit function 
   */
  MeritFcn getMeritFcn() const  { return mfcn;}
  /// Specify the merit function to used in step acceptance test
  virtual void setMeritFcn(MeritFcn option) { mfcn = option;}

  /**
   * @return Gradient of the Lagrangian at the current iteration 
   */
  NEWMAT::ColumnVector getGradL() const  { return gradl;}
  /// Store the gradient of the Lagrangian at the current iteration
  virtual void setGradL(NEWMAT::ColumnVector gradl_value) { gradl = gradl_value;}

  /**
   * @return Gradient of the Lagrangian at the previous iteration 
   */
  NEWMAT::ColumnVector getGradLPrev() const  { return gradlprev;}
  /// Store the gradient of the Lagrangian at the previous iteration
  virtual void setGradLPrev(NEWMAT::ColumnVector gradl_value) { gradlprev = gradl_value;}

  /**
   * @return Residuals of the constraints at the current iteration 
   */
  NEWMAT::ColumnVector getConstraintResidual() const  { return constraintResidual;}
  /// Store the residuals of the constraints at the current iteration
  virtual void setConstraintResidual(const NEWMAT::ColumnVector& constraint_value) 
                       { constraintResidual = constraint_value;}

  /**
   * @return Gradient of the constraints at the current iteration 
   */
  NEWMAT::Matrix getConstraintGradient() const  { return constraintGradient;}
  /// Store the current gradients of the constraints 
  virtual void setConstraintGradient(const NEWMAT::Matrix& constraint_grad) 
                       { constraintGradient = constraint_grad;}

  /**
   * @return Value of merit function 
   */
  real getCost() const  { return cost;}
  /// Store current value of merit function
  void setCost(real value) { cost = value;}

  /**
   * @return Switch to turn on feasibility recovery method 
   */
  bool getFeasibilityRecoveryFlag() const {return feas_flag;}
  /// Set switch to turn on feasibility recovery method 
  void setFeasibilityRecoveryFlag(bool flag) {feas_flag = flag;}

  /**
   * @return Maximum number of iterations for feasibility recovery method 
   */
  int getMaxFeasIter() const {return max_feas_iter;}
  /// Set maximum number of iterations for feasibility recovery method of trust-region
  void setMaxFeasIter(int k) {max_feas_iter = k;}

  /**
   * @return Hessian matrix 
   */
  NEWMAT::SymmetricMatrix getHessian() const {return Hessian;}
  /// Store current Hessian matrix 
  void setHessian(NEWMAT::SymmetricMatrix& H) {Hessian = H;}
  
  /// Compute the Hessian of the Lagrangrian or its approximation at iteration k
  virtual NEWMAT::SymmetricMatrix updateH(NEWMAT::SymmetricMatrix& H, int k) = 0;

  /// Compute the active set according to Facchinei, Fischer, and Kanzow indicator 
  NEWMAT::ColumnVector computeFFK1Ind(const NEWMAT::ColumnVector& xc);
  /// Compute the active set according to Facchinei, Fischer, and Kanzow indicator 
  NEWMAT::ColumnVector computeFFK2Ind(const NEWMAT::ColumnVector& xc);
  /// Compute the Tapia indicators for the constrained problem
  NEWMAT::ColumnVector computeTapiaInd(const NEWMAT::ColumnVector& step);

  /// Print status of the constrained Newton's method
  void printStatus(char *);
  /// Output the Lagrangian multipliers to the screen 
  void printMultipliers(char *);
  /// Print the Lagrangian multipliers to a file 
  void fPrintMultipliers(ostream *nlpout, char *);
  /// Print second order sufficiency information to a file 
  void fPrintSecSuff(ostream *nlpout, NEWMAT::ColumnVector& info);

  friend int trustregion(NLP1*, ostream*, NEWMAT::SymmetricMatrix&,
			 NEWMAT::ColumnVector&, NEWMAT::ColumnVector&,
			 real&, real&, real stpmax, real stpmin);

  friend int trustpds(NLP1*, ostream*, NEWMAT::SymmetricMatrix&,
		      NEWMAT::ColumnVector&, NEWMAT::ColumnVector&,
		      real&, real&, real stpmax, real stpmin, int);
};

/**
 * Constrained Newton classes that will accept either an NLP1 or NLP2
 */

class OptConstrNewton1Deriv: public OptConstrNewtonLike {
public:

  OptConstrNewton1Deriv() {}

  OptConstrNewton1Deriv(NLP1* p): 
    OptConstrNewtonLike(p->getDim()), mem_nlp(p) {;}
   
  OptConstrNewton1Deriv(NLP1* p, UPDATEFCN u): 
    OptConstrNewtonLike(p->getDim(),u), mem_nlp(p) {;}

  OptConstrNewton1Deriv(NLP1* p, TOLS t): 
    OptConstrNewtonLike(p->getDim(),t), mem_nlp(p) {;}
  
  virtual ~OptConstrNewton1Deriv(){}

private:
  NLP1* mem_nlp;
  
protected:
  /**
   * @ returns a pointer to an NLP1
   */
  NLP1* nlprob() const { return mem_nlp; }
};

/**
 * Constrained Newton classes that require an NLP2
 */

class OptConstrNewton2Deriv: public OptConstrNewtonLike {
public:

  OptConstrNewton2Deriv() {}

  OptConstrNewton2Deriv(NLP2* p): 
    OptConstrNewtonLike(p->getDim()), mem_nlp(p){;}

  OptConstrNewton2Deriv(NLP2* p, UPDATEFCN u): 
    OptConstrNewtonLike(p->getDim(),u), mem_nlp(p){;}

  OptConstrNewton2Deriv(NLP2* p, TOLS t): 
    OptConstrNewtonLike(p->getDim(),t), mem_nlp(p){;}
  
  virtual ~OptConstrNewton2Deriv(){}

private:
  NLP2* mem_nlp;
  
protected:
  /**
   * @ returns a pointer to an NLP1
   */
  NLP1* nlprob() const { return mem_nlp;}
  /**
   * @ returns a pointer to an NLP2
   */
  NLP2* nlprob2() const { return mem_nlp;}
};

} // namespace OPTPP

#endif
