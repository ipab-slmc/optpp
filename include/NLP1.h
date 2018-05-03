#ifndef NLP1_h
#define NLP1_h

#include "NLP0.h"

using std::ostream;

/**
 * NLP1: NLP0 + First derivatives
 * NLP1 is derived from the base class NLP0
 * by adding the necessary information to compute
 * and store the gradient.
 *
 * @author J.C. Meza, Sandia National Laboratories, meza@ca.sandia.gov
 *
 * @note Modified by P. J. Williams to incorporate namespaces
 * @note Sandia National Laboratories, pwillia@sandia.gov
 */

namespace OPTPP {

class NLP1: public NLP0 {
protected:
  /// Gradient of objective function at mem_xc
  NEWMAT::ColumnVector mem_grad;     	
  /// Number of gradient evaluations
  int          ngevals;		
  /// Is an analytic gradient available?
  int          analytic_grad;

public:
// Constructors
/**
 * Default Constructor
 * @see NLP1(int dim)
 * @see NLP1(int dim, int nlncons)
 * @see NLP1(int dim, CompoundConstraint* constraint)
 */
  NLP1(): 
    NLP0(), mem_grad(0), ngevals(0) {;}
/**
 * @param ndim an integer argument
 * @see NLP1(int dim, int nlncons)
 * @see NLP1(int dim, CompoundConstraint* constraint)
 */
  NLP1(int ndim): 
    NLP0(ndim), mem_grad(ndim), ngevals(0) {;}
/**
 * @param ndim an integer argument
 * @param nlncons an integer argument
 * @see NLP1(int dim)
 * @see NLP1(int dim, CompoundConstraint* constraint)
 */
  NLP1(int ndim,int nlncons): 
    NLP0(ndim, nlncons), mem_grad(ndim), ngevals(0) {;}
/**
 * @param ndim an integer argument
 * @param constraint pointer to a CompoundConstraint object 
 * @see NLP1(int dim)
 * @see NLP1(int dim, int nlncons)
 */
  NLP1(int ndim, CompoundConstraint* constraint): 
    NLP0(ndim, constraint), mem_grad(ndim), ngevals(0) {;}

/**
 * Destructor
 */
  virtual ~NLP1() {;}                     

  /**
   * @return Gradient of objective function.
   */
  void setGrad(const NEWMAT::ColumnVector& set_grad) {mem_grad = set_grad;}

 /**
  * @return Gradient of objective function.
  */
  NEWMAT::ColumnVector getGrad() const {return mem_grad;}

 /**
  * @return Number of gradient evaluations.
  */
  int getGevals()        const {return ngevals;}

 /**
  * @return Is there analytic gradient information available?
  */
  int AnalyticGrad()     const {return analytic_grad;}

/// Reset the parameter values 
  virtual void reset() = 0;   

/// Evaluate the function
  virtual void initFcn() = 0;   
  virtual real evalF() = 0;
  virtual real evalF(const NEWMAT::ColumnVector& x) = 0;
  virtual void eval()  = 0;

/// Evaluate the gradient
  virtual NEWMAT::ColumnVector evalG() = 0;
  virtual NEWMAT::ColumnVector evalG(const NEWMAT::ColumnVector& x) = 0;

/// Evaluate a Finite-difference Hessian
  virtual NEWMAT::SymmetricMatrix evalH() = 0;
  virtual NEWMAT::SymmetricMatrix evalH(NEWMAT::ColumnVector& x) = 0;
  virtual NEWMAT::SymmetricMatrix FDHessian(NEWMAT::ColumnVector& x);


/// Evaluate the Lagrangian, its gradient and Hessian
  virtual real evalLagrangian(const NEWMAT::ColumnVector& x, NEWMAT::ColumnVector& mult,
    const NEWMAT::ColumnVector& type) = 0;
  virtual NEWMAT::ColumnVector evalLagrangianGradient(const NEWMAT::ColumnVector& x,
    const NEWMAT::ColumnVector& mult, const NEWMAT::ColumnVector& type) = 0;

/// Evaluate the constraints
  virtual NEWMAT::ColumnVector evalCF(const NEWMAT::ColumnVector& x) = 0;
  virtual NEWMAT::Matrix evalCG(const NEWMAT::ColumnVector& x) = 0;
  virtual NEWMAT::SymmetricMatrix evalCH(NEWMAT::ColumnVector& x) = 0;
  virtual OptppArray<NEWMAT::SymmetricMatrix> evalCH(NEWMAT::ColumnVector& x, int darg) = 0;
  virtual void evalC(const NEWMAT::ColumnVector& x) = 0;

/// Evaluate a Finite-difference Hessian for the nonlinear constraints
  virtual OptppArray<NEWMAT::SymmetricMatrix> CONFDHessian(NEWMAT::ColumnVector& x);

/// Print the function
  virtual void printState(char* s);   
  virtual void fPrintState(ostream *nlpout, char* s);   

};

} // namespace OPTPP
#endif
