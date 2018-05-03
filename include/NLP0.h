#ifndef NLP0_h
#define NLP0_h

#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#include "NLPBase.h"
#include "Appl_Data.h"
#include "CompoundConstraint.h"

using std::ostream;

namespace OPTPP {

/**
 *
 * Base Class for NonLinear Programming Problem
 * For NLP0 the only assumption on the objective function is
 * that it be continuous. No derivative information is available.
 * Note that NLP0-2 are abstract data types
 *
 * @author J.C. Meza, Sandia National Laboratories, meza@ca.sandia.gov
 *
 */


class NLP0: public NLPBase{ 

protected:
  int          dim;          		///< Dimension of the problem
  NEWMAT::ColumnVector mem_xc;  	///< Current point
  real         fvalue;  		///< Objective function value at mem_xc
  char         nlp_name[80]; 		///< Optional identifying name
  NEWMAT::ColumnVector mem_fcn_accrcy;	///< Accuracy available from function 
  int          nfevals;      		///< Number of function evaluations
  int          is_expensive; 		///< Is this an expensive function?
  bool         debug_;			///< Print debug statements
  bool         modeOverride;	        
  double       function_time;		///< Function compute time
  CompoundConstraint* constraint_;  	///< Pointer to constraints
  NEWMAT::ColumnVector  constraint_value;///< Constraint residual 
  int          ncnln;      		///< Number of nonlinear constraints   
  Appl_Data    application;
  DerivOption finitediff;		///< User-specified derivative option
  SpecOption SpecFlag;          	///< Speculative gradient information
  NEWMAT::ColumnVector partial_grad;
  double specF;

public:
#ifdef WITH_MPI

// Constructors
 /**
  * Default Constructor
  * @see NLP0(int dim)
  * @see NLP0(int dim, int nlncons)
  * @see NLP0(int dim, CompoundConstraint* constraint)
  */
  NLP0():
    dim(0),mem_xc(0),fvalue(1.0e30), mem_fcn_accrcy(0),
    nfevals(0),is_expensive(0),debug_(0), modeOverride(0), function_time(0.0), 
    constraint_(0), constraint_value(0), ncnln(0), partial_grad(0)
    {mem_xc = 0; mem_fcn_accrcy = DBL_EPSILON; finitediff = ForwardDiff;
    SpecFlag = Spec1;}
 /**
  * @param ndim an int
  * @see NLP0()
  * @see NLP0(int dim, int nlncons)
  * @see NLP0(int dim, CompoundConstraint* constraint)
  */
  NLP0(int ndim):
    dim(ndim),mem_xc(ndim),fvalue(1.0e30), mem_fcn_accrcy(ndim),
    nfevals(0),is_expensive(0),debug_(0), modeOverride(0), function_time(0.0), 
    constraint_(0), constraint_value(0), ncnln(0), partial_grad(ndim)
    {mem_xc = 0; mem_fcn_accrcy = DBL_EPSILON; finitediff = ForwardDiff;
    SpecFlag = Spec1;}
 /**
  * @param ndim an int
  * @param nlncons an int
  * @see NLP0()
  * @see NLP0(int dim)
  * @see NLP0(int dim, CompoundConstraint* constraint)
  */
  NLP0(int ndim, int nlncons):
    dim(ndim),mem_xc(ndim),fvalue(1.0e30), mem_fcn_accrcy(ndim),
    nfevals(0),is_expensive(0),debug_(0), modeOverride(0), function_time(0.0), 
    constraint_(0), constraint_value(nlncons), ncnln(nlncons),
    partial_grad(ndim)
    {mem_xc = 0; mem_fcn_accrcy = DBL_EPSILON; finitediff = ForwardDiff;
    SpecFlag = Spec1; constraint_value = 0;}
 /**
  * @param ndim an int
  * @param constraint pointer to a CompoundConstraint
  * @see NLP0()
  * @see NLP0(int dim)
  * @see NLP0(int dim, int nlncons)
  */
  NLP0(int ndim, CompoundConstraint* constraint):
    dim(ndim),mem_xc(ndim),fvalue(1.0e30), mem_fcn_accrcy(ndim),
    nfevals(0),is_expensive(0),debug_(0), modeOverride(0), function_time(0.0), 
    constraint_(constraint), constraint_value(0), ncnln(0), partial_grad(ndim) 
    {mem_xc = 0; mem_fcn_accrcy = DBL_EPSILON; finitediff = ForwardDiff;
    SpecFlag = Spec1;}

#else
      
// Constructors
 /**
  * Default Constructor
  * @see NLP0(int dim)
  * @see NLP0(int dim, int nlncons)
  * @see NLP0(int dim, CompoundConstraint* constraint)
  */
  NLP0():
    dim(0),mem_xc(0),fvalue(1.0e30), mem_fcn_accrcy(0),
    nfevals(0),is_expensive(0),debug_(0), modeOverride(0), function_time(0.0), 
    constraint_(0), constraint_value(0), ncnln(0), partial_grad(0) 
    {mem_xc = 0; mem_fcn_accrcy = DBL_EPSILON; finitediff = ForwardDiff;
    SpecFlag = NoSpec;}
 /**
  * @param ndim an int
  * @see NLP0()
  * @see NLP0(int dim, int nlncons)
  * @see NLP0(int dim, CompoundConstraint* constraint)
  */
  NLP0(int ndim):
    dim(ndim),mem_xc(ndim),fvalue(1.0e30), mem_fcn_accrcy(ndim),
    nfevals(0),is_expensive(0),debug_(0), modeOverride(0), function_time(0.0), 
    constraint_(0), constraint_value(0), ncnln(0), partial_grad(ndim)
    {mem_xc = 0; mem_fcn_accrcy = DBL_EPSILON; finitediff = ForwardDiff;
    SpecFlag = NoSpec;}
 /**
  * @param ndim an int
  * @param nlncons an int
  * @see NLP0()
  * @see NLP0(int dim)
  * @see NLP0(int dim, CompoundConstraint* constraint)
  */
  NLP0(int ndim, int nlncons):
    dim(ndim),mem_xc(ndim),fvalue(1.0e30), mem_fcn_accrcy(ndim),
    nfevals(0),is_expensive(0),debug_(0), modeOverride(0), function_time(0.0), 
    constraint_(0), constraint_value(nlncons), ncnln(nlncons),
    partial_grad(ndim)
    {mem_xc = 0; mem_fcn_accrcy = DBL_EPSILON; finitediff = ForwardDiff;
    SpecFlag = NoSpec; constraint_value = 0;}
 /**
  * @param ndim an int
  * @param constraint pointer to a CompoundConstraint
  * @see NLP0()
  * @see NLP0(int dim)
  * @see NLP0(int dim, int nlncons)
  */
  NLP0(int ndim, CompoundConstraint* constraint):
    dim(ndim),mem_xc(ndim),fvalue(1.0e30), mem_fcn_accrcy(ndim),
    nfevals(0),is_expensive(0),debug_(0), modeOverride(0), function_time(0.0), 
    constraint_(constraint), constraint_value(0), ncnln(0),
    partial_grad(ndim) 
    {mem_xc = 0; mem_fcn_accrcy = DBL_EPSILON; finitediff = ForwardDiff;
    SpecFlag = NoSpec;}

#endif

/// Functions for setting various properties of this NLP problem
  virtual void setX(const int i, const real& x)  {mem_xc(i) = x;}
  virtual void setX(const NEWMAT::ColumnVector& x)       {mem_xc = x;}

  virtual void setF(const real& fx)              {fvalue = fx;}

  virtual void setIsExpensive(const int e) 	 {is_expensive = e;}
  virtual void setFcnAccrcy(const int i, const real& accrcy)  
                                                   {mem_fcn_accrcy(i)=accrcy;}
  virtual void setFcnAccrcy(const NEWMAT::ColumnVector& accrcy) {mem_fcn_accrcy=accrcy;}

  /**
   * @return Problem Dimension
   */
  virtual int  getDim()         	const {return dim;}
  /**
   * @return Number of Function Evaluations 
   */
  virtual int  getFevals()      	const {return nfevals;}
  /**
   * @return Is the function expensive? 
   */
  virtual int  getIsExpensive() 	const {return is_expensive;}

  virtual void setModeOverride(bool override_mode) {modeOverride = override_mode;}
   virtual bool getModeOverride() const {return modeOverride;}

  /**
   * @return The function value 
   */
  virtual real getF()           	const {return fvalue;}
  /**
   * @return NEWMAT::ColumnVector of the function accuracy with respect to each variable
   */
  virtual NEWMAT::ColumnVector getFcnAccrcy()   const {return mem_fcn_accrcy;}
  /**
   * @return Current point
   */
  virtual NEWMAT::ColumnVector getXc()  	const {return mem_xc;}
  /**
   * @return Function compute time 
   */
  virtual real getFcnTime()     	const {return function_time;}

// Debugging tools
  virtual void setDebug()   	{debug_ = true;}
  virtual bool getDebug()       const {return debug_;}

  // Set and get the finite differencing options.

  void setDerivOption(DerivOption d) {finitediff = d;}
  DerivOption getDerivOption() const {return finitediff;}

  // Set and get the speculative gradient options.

  void setSpecOption(SpecOption SpecEval) {SpecFlag = SpecEval;}
  SpecOption getSpecOption() const {return SpecFlag;}

// Function to reset parameter values 
  virtual void reset() = 0;   

// Functions for evaluating the objective fcn 
  virtual void initFcn() = 0;   
  virtual void eval()    = 0;
  virtual real evalF()   = 0;
  virtual real evalF(const NEWMAT::ColumnVector& x) = 0;

  // Constraint helper functions
  /**
   * @return Total number of constraints 
   */
  int  getNumOfCons()  const { return constraint_->getNumOfCons();}
  /**
   * @return Number of nonlinear constraints 
   */
  int  getNumOfNLCons()  const { return constraint_->getNumOfNLCons();}
  /**
   * @return Constraint value 
   */
  NEWMAT::ColumnVector  getConstraintValue()  const { return constraint_value;}
  /// Set the constraint value
  void  setConstraintValue(const NEWMAT::ColumnVector& cfx)  { constraint_value = cfx;}
  /**
   * @return 1 - if problem has constraints, 0 - otherwise
   */
  bool hasConstraints();
  /**
   * Print constraints to standard output
   */
  void printConstraints();
  /**
   * @return Pointer to a CompoundConstraint object 
   */
  CompoundConstraint* getConstraints() { return constraint_;}
  /// Set a pointer to a CompoundConstraint object 
  void setConstraints(CompoundConstraint* constraintSet) 
                { constraint_ = constraintSet;}

/// Destructor
  virtual ~NLP0() {;}        

//------------------------------------------------------------------------
// These are defined elsewhere
//------------------------------------------------------------------------

/// Finite-difference gradient and Hessian
  NEWMAT::ColumnVector FDGrad(const NEWMAT::ColumnVector &, 
    const NEWMAT::ColumnVector &, double &, NEWMAT::ColumnVector &);        
  NEWMAT::ColumnVector BDGrad(const NEWMAT::ColumnVector &, 
    const NEWMAT::ColumnVector &, double &, NEWMAT::ColumnVector &);
  NEWMAT::ColumnVector CDGrad(const NEWMAT::ColumnVector &,
    const NEWMAT::ColumnVector &, double &, NEWMAT::ColumnVector &);
  NEWMAT::SymmetricMatrix FD2Hessian(NEWMAT::ColumnVector &);

/// Finite-difference gradient and Hessian of nonlinear constraints
  NEWMAT::Matrix CONFDGrad(const NEWMAT::ColumnVector &);        
  NEWMAT::Matrix CONBDGrad(const NEWMAT::ColumnVector &);        
  NEWMAT::Matrix CONCDGrad(const NEWMAT::ColumnVector &);        

/// Evaluate a finite-difference gradient and Hessian 
  virtual NEWMAT::ColumnVector evalG() = 0;
  virtual NEWMAT::ColumnVector evalG(const NEWMAT::ColumnVector& x) = 0;
  virtual NEWMAT::SymmetricMatrix evalH() = 0;
  virtual NEWMAT::SymmetricMatrix evalH(NEWMAT::ColumnVector& x) = 0;

/// Evaluate the Lagrangian
  virtual real evalLagrangian(const NEWMAT::ColumnVector& x, NEWMAT::ColumnVector& mult,
    const NEWMAT::ColumnVector& type) = 0;
/// Evaluate the Lagrangian gradient 
  virtual NEWMAT::ColumnVector evalLagrangianGradient(const NEWMAT::ColumnVector& x,
    const NEWMAT::ColumnVector& mult, const NEWMAT::ColumnVector& type) = 0;

/// Evaluate the constraint at x
  virtual NEWMAT::ColumnVector evalCF(const NEWMAT::ColumnVector& x) = 0;
/// Evaluate the constraint gradient at x
  virtual NEWMAT::Matrix  evalCG(const NEWMAT::ColumnVector& x) = 0;
  virtual NEWMAT::SymmetricMatrix evalCH(NEWMAT::ColumnVector& x) = 0;
/// Evaluate the constraint Hessian at x
  virtual OptppArray<NEWMAT::SymmetricMatrix> evalCH(NEWMAT::ColumnVector& x, int darg) = 0;
  virtual void evalC(const NEWMAT::ColumnVector& x) = 0;

  /// Print the function
  virtual void printState(char *); 
  /// Print state of the function to a file
  virtual void fPrintState(ostream *, char *); 
  /// Save current state of the function
  void saveState(); 
};

} // namespace OPTPP
#endif
