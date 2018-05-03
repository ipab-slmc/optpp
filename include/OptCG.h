#ifndef OptCG_h
#define OptCG_h

/*----------------------------------------------------------------------
  Copyright (c) 2001, Sandia Corporation.
  J.C. Meza, Sandia National Laboratories, meza@ca.sandia.gov
 ----------------------------------------------------------------------*/

#ifndef Opt_h
#include "Opt.h"
#endif

namespace OPTPP {

 /**
  * CG-Like Methods
  * OptCG is a derived class from OptCGLike, which implements a
  * nonlinear conjugate gradient method. This version uses the
  * Polak-Ribiere formula and a linesearch routine due to More
  * and Thuente as implemented in the routines mcsrch and mcstep.
  *
  * @author J.C. Meza, Lawrence Berkeley National Laboratory
  */

class OptCGLike: public OptimizeClass {
 protected:
  virtual NLP1 *nlprob() const = 0;
  /// Gradient computed at previous iteration
  NEWMAT::ColumnVector gprev;		
  /// Number of gradient evaluations 
  int grad_evals;		
  /// User-specified globalization strategy
  SearchStrategy strategy;	

 public:
  /**
   * Default Constructor
   * @see OptCGLike(int n)
   * @see OptCGLike(int n, TOLS t)
   */
  OptCGLike(){}
  /**
   * @param n an integer argument
   * @see OptCGLike(int n, TOLS t)
   */
  OptCGLike(int n): OptimizeClass(n), 
                    gprev(n), grad_evals(0), strategy(LineSearch){}
  /**
   * @param n an integer argument
   * @param t a TOLS object
   * @see OptCGLike(int n)
   */
  OptCGLike(int n, TOLS t): OptimizeClass(n,t),
                            gprev(n),grad_evals(0), strategy(LineSearch){}

  /**
   * Destructor
   */
  virtual ~OptCGLike(){}

  /// Set the user-specified globalization strategy
  void setSearchStrategy(SearchStrategy s) {strategy = s;}

   /**
    * @return User-specified globalization strategy
    */
  SearchStrategy getSearchStrategy() const {return strategy;}

  virtual void acceptStep(int, int) = 0;
  virtual int checkConvg();
  virtual int checkDeriv();
  virtual void optimize() {} 
  virtual void readOptInput() {}
  virtual void updateModel(int, int, NEWMAT::ColumnVector) = 0;

};

 /**
  * CG-Like Methods
  * OptCG is a derived class from OptCGLike, which implements a
  * nonlinear conjugate gradient method. This version uses the
  * Polak-Ribiere formula and a linesearch routine due to More
  * and Thuente as implemented in the routines mcsrch and mcstep.
  *
  * @author J.C. Meza, Sandia National Laboratories, meza@ca.sandia.gov
  */

class OptCG: public OptCGLike {
private:
  NLP1* nlp;	/// Pointer to an NLP1 object
  
protected:
  /**
   * @return Pointer to an NLP1 object
   */
  NLP1* nlprob() const { return nlp; }

public:
 /**
  * Default Constructor
  * @see OptCG(NLP1* p)
  * @see OptCG(NLP1* p, TOLS t)
  */
  OptCG(){strcpy(method,"Nonlinear CG");}
 /**
  * @param p a pointer to an NLP1 object
  * @see OptCG(NLP1* p, TOLS t)
  */
  OptCG(NLP1* p): OptCGLike(p->getDim()),nlp(p)
    {strcpy(method,"Nonlinear CG");}
 /**
  * @param p a pointer to an NLP1 object
  * @param t a TOLS object
  * @see OptCG(NLP1* p)
  */
  OptCG(NLP1* p, TOLS t): OptCGLike(p->getDim(), t),nlp(p)
    {strcpy(method,"Nonlinear CG");}

 /**
  * Destructor
  */
  virtual ~OptCG(){}

  virtual NEWMAT::ColumnVector computeSearch(NEWMAT::SymmetricMatrix& ) 
    {return NEWMAT::ColumnVector();}

  virtual void acceptStep(int k, int step_type)
    {OptimizeClass::defaultAcceptStep(k, step_type);}

  virtual void updateModel(int k, int ndim, NEWMAT::ColumnVector x)
    {OptimizeClass::defaultUpdateModel(k, ndim, x);}

// These are defined elsewhere
  /// Compare the analytic gradient with the finite-difference gradient
  virtual int checkDeriv();
  /// Compute the step direction based upon specified globalization strategy 
  virtual int computeStep(NEWMAT::ColumnVector sk);
  /// Reset the parameters 
  virtual void reset();
  /// Initialize the optimization method
  virtual void initOpt();
  /// Run the optimization method
  virtual void optimize();
  /// Compute steplength
  virtual real stepTolNorm() const;
  /// Print the status to the optimization method at the current iteration
  virtual void printStatus(char *);
};

} // namespace OPTPP
#endif
