#ifndef OptLBFGS_h
#define OptLBFGS_h

/*----------------------------------------------------------------------
  Copyright (c) 2003
 ----------------------------------------------------------------------*/

#ifndef Opt_h
#include "Opt.h"
#endif

namespace OPTPP {

 /**
  * LBFGS-Like Methods
  * OptLBFGS is a derived class of OptLBFGSLike that 
  * implements the LBFGS method of J. Nocedal.
  *
  * @author R.A.Oliva, Lawrence Berkely National Laboratories, raoliva@lbl.gov
  */

class OptLBFGSLike: public OptimizeClass {
 protected:
  virtual NLP1 *nlprob() const = 0;

  NEWMAT::ColumnVector gprev;		///< Gradient computed at previous iteration

  int grad_evals;		///< Number of gradient evaluations 

  SearchStrategy strategy;	///< User-specified globalization strategy

 public:
  /**
   * Default Constructor
   * @see OptLBFGSLike(int n)
   * @see OptLBFGSLike(int n, TOLS t)
   */
  OptLBFGSLike(){}
  /**
   * @param n an integer argument
   * @see OptLBFGSLike(int n, TOLS t)
   */
  OptLBFGSLike(int n): OptimizeClass(n), 
                    gprev(n), grad_evals(0), strategy(LineSearch){}
  /**
   * @param n an integer argument
   * @param t a TOLS object
   * @see OptLBFGSLike(int n)
   */
  OptLBFGSLike(int n, TOLS t): OptimizeClass(n,t),
                            gprev(n),grad_evals(0), strategy(LineSearch){}

  /**
   * Destructor
   */
  virtual ~OptLBFGSLike(){}

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
   * The Limited Memory BFGS Method for Large Scale Optimization
   *
   * Solves the unconstrained minimization problem
   * 
   *      min F(x)    x = (x_1, x_2, ..., x_N),
   *
   * using the limited memory BFGS method, where N can be large.
   *
   * The inverse Hessian approximation Hk is computed via BFGS
   * updates to a diagonal matrix H0.  The number of updates
   * depend on the previous m steps kept in memory, as set by the
   * user. H0 can be any symmetric positive definite matrix 
   * specified by the user (else a default one is constructed).
   *
   * References:
   *
   * D. Liu and J. Nocedal, 
   * "On the limited memory BFGS method for large scale optimization"
   * Mathematical Programming B 45 (1989), 503-528
   *
   * @author R.A.Oliva, Lawrence Berkely National Laboratories, raoliva@lbl.gov
   */

class OptLBFGS: public OptLBFGSLike {
private:
  NLP1* nlp;	///< Pointer to an NLP1 object  

  int memM;     ///< number of memory vectors kept during optimization iteration

  bool printXs; ///< controls if final point is printed by printStatus()

protected:
  /**
   * @return Pointer to an NLP1 object
   */
  NLP1* nlprob() const { return nlp; }

  void initMem(int n) {

    if (n>=30) {memM = 15; return;}
    if (n>= 2) {memM =  2; return;}
    if (n>= 5) {memM =  3; return;}
    if (n>=10) {memM =  5; return;}
    if (n>=20) {memM = 10; return;}
    memM = 1;
  }

public:

 /**
  * Default Constructor
  * @see OptLBFGS(NLP1* p)
  * @see OptLBFGS(NLP1* p, TOLS t)
  */
  OptLBFGS(): memM(15), printXs(false) {strcpy(method,"Limited Memory BFGS method");}

 /**
  * @param p a pointer to an NLP1 object
  * @see OptLBFGS(NLP1* p, TOLS t)
  */
  OptLBFGS(NLP1* p): OptLBFGSLike(p->getDim()), nlp(p), printXs(false)
    {strcpy(method,"Limited Memory BFGS method"); initMem(p->getDim());}

 /**
  * @param p a pointer to an NLP1 object
  * @param m integer specifying number of memory vectors
  * @see OptLBFGS(NLP1* p)
  * @see OptLBFGS(NLP1* p, TOLS t)
  * @see OptLBFGS(NLP1* p, TOLS t, int m)
  */
  OptLBFGS(NLP1* p, int m): OptLBFGSLike(p->getDim()), nlp(p), printXs(false) {
      strcpy(method,"Limited Memory BFGS method");
      memM = (m <= p->getDim())? m : p->getDim();
  }

 /**
  * @param p a pointer to an NLP1 object
  * @param t a TOLS object
  * @see OptLBFGS(NLP1* p, TOLS t, int m)
  */
  OptLBFGS(NLP1* p, TOLS t): OptLBFGSLike(p->getDim(),t), nlp(p), printXs(false)
    {strcpy(method,"Limited Memory BFGS method"); initMem(p->getDim());}

 /**
  * @param p a pointer to an NLP1 object
  * @param t a TOLS object
  * @param m integer specifying number of memory vectors
  * @see OptLBFGS(NLP1* p, TOLS t)
  */
  OptLBFGS(NLP1* p, TOLS t, int m): OptLBFGSLike(p->getDim(),t), nlp(p), printXs(false) {
      strcpy(method,"Limited Memory BFGS method");
      memM = (m <= p->getDim())? m : p->getDim();
  }

 /**
  * Destructor
  */
  virtual ~OptLBFGS(){}

  virtual NEWMAT::ColumnVector 
    computeSearch(NEWMAT::SymmetricMatrix& ) {return NEWMAT::ColumnVector();}

  virtual void acceptStep(int k, int step_type)
    {OptimizeClass::defaultAcceptStep(k, step_type);}

  virtual void updateModel(int k, int ndim, NEWMAT::ColumnVector x)
    {OptimizeClass::defaultUpdateModel(k, ndim, x);}

  void setPrintFinalX(bool b) { printXs = b;}

  //--
  // Defined in OptLBFGS.C:
  //--
  virtual int checkDeriv();                    /// Compare analytic vs finite-difference gradient
  virtual int computeStep(NEWMAT::ColumnVector& sk, double stp=1.0);   /// Compute the step direction
  virtual void reset();                       /// Reset the parameters 
  virtual void initOpt();                     /// Initialize the optimization 
  virtual void optimize();                    /// Run the optimization 
  virtual real stepTolNorm() const;  
  virtual void printStatus(char *c);  
  virtual void printIter(int, double, double, double, double, int);
};

} // namespace OPTPP
#endif
