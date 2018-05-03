#ifndef NLP2_h
#define NLP2_h

#include "NLP1.h"

using std::ostream;

/**
 * NLP2: NLP1 + Second derivatives
 * NLP2 is derived from NLP1 by adding the necessary
 * information to compute and store the Hessian.
 * @author J.C. Meza, Lawrence Berkeley National Laboratory, 
 *
 * @note Modified by P. J. Williams to incorporate namespaces
 * @date 02/2006 
 */

namespace OPTPP {

class NLP2: public NLP1 {
protected:
  /// Hessian of objective fcn (or approx.)
  NEWMAT::SymmetricMatrix Hessian; 	
  /// Number of Hessian evaluations
  int     nhevals;		

public:
 /**
  * Default Constructor
  * @see NLP2(int dim)
  * @see NLP2(int dim, int nlncons)
  * @see NLP2(int dim, CompoundConstraint* constraint)
  */
  NLP2(): 
     NLP1(), Hessian(0), nhevals(0) {;}
 /**
  * @param ndim an int
  * @see NLP2()
  * @see NLP2(int dim, int nlncons)
  * @see NLP2(int dim, CompoundConstraint* constraint)
  */
  NLP2(int ndim): 
     NLP1(ndim), Hessian(ndim), nhevals(0) {;}
 /**
  * @param ndim  problem dimension
  * @param nlncons number of nonlinear constraints
  * @see NLP2()
  * @see NLP2(int dim)
  * @see NLP2(int dim, CompoundConstraint* constraint)
  */
  NLP2(int ndim, int nlncons): 
     NLP1(ndim, nlncons), Hessian(ndim), nhevals(0) {;}
 /**
  * @param ndim  problem dimension
  * @param constraint pointer to a CompoundConstraint
  * @see NLP2()
  * @see NLP2(int dim)
  * @see NLP2(int dim, int nlncons)
  */
  NLP2(int ndim, CompoundConstraint* constraint): 
     NLP1(ndim,constraint), Hessian(ndim), nhevals(0) {;}

 /**
  * Destructor
  */
  virtual ~NLP2() {}                     

 /**
  * @return Hessian of the objective function
  */
  NEWMAT::SymmetricMatrix getHess() const {return Hessian;}

 /**
  * @return Number of Hessian evaluations
  */
  int getHevals()           const {return nhevals;}

/// Reset the parameter values 
  virtual void reset() = 0;   
/// Initialize the function
  virtual void initFcn() = 0;   
/// Evaluate the function at the current point
  virtual real evalF() = 0;
/// Evaluate the function at x
  virtual real evalF(const NEWMAT::ColumnVector& x) = 0;
/// Evaluate the function, gradient, and Hessian at the current point
  virtual void eval() = 0;   

/// Evaluate the gradient at the current point
  virtual NEWMAT::ColumnVector evalG() = 0;
/// Evaluate the gradient at x
  virtual NEWMAT::ColumnVector evalG(const NEWMAT::ColumnVector& x) = 0;

/// Evaluate the Hessian at the current point
  virtual NEWMAT::SymmetricMatrix evalH() = 0;
/// Evaluate the Hessian at x
  virtual NEWMAT::SymmetricMatrix evalH(NEWMAT::ColumnVector& x) = 0;


/// Evaluate the Lagrangian at the x 
  virtual real evalLagrangian(const NEWMAT::ColumnVector& x, NEWMAT::ColumnVector& mult,
                              const NEWMAT::ColumnVector& type) = 0;
/// Evaluate the gradient of the Lagrangian at x
  virtual NEWMAT::ColumnVector evalLagrangianGradient(const NEWMAT::ColumnVector& x,
                                              const NEWMAT::ColumnVector& mult,
					      const NEWMAT::ColumnVector& type) = 0;

/// Evaluate the constraint at x
  virtual NEWMAT::ColumnVector evalCF(const NEWMAT::ColumnVector& x) = 0;
/// Evaluate the constraint gradient at x
  virtual NEWMAT::Matrix evalCG(const NEWMAT::ColumnVector& x) = 0;
  virtual NEWMAT::SymmetricMatrix evalCH(NEWMAT::ColumnVector& x) = 0;
/// Evaluate the constraint Hessian at x
  virtual OptppArray<NEWMAT::SymmetricMatrix> evalCH(NEWMAT::ColumnVector& x, int darg) = 0;
  virtual void evalC(const NEWMAT::ColumnVector& x) = 0;


/// Print state of NLP2 object to screen 
  virtual void printState(char* s) ;
/// Print state of NLP2 object to file
  virtual void fPrintState(ostream *nlpout, char* s ) ;

};

/**
 * Print various quantities in newmat classes
 */

void Print(const NEWMAT::Matrix& X);
void Print(const NEWMAT::UpperTriangularMatrix& X);
void Print(const NEWMAT::DiagonalMatrix& X);
void Print(const NEWMAT::SymmetricMatrix& X);
void Print(const NEWMAT::LowerTriangularMatrix& X);

void FPrint(ostream *fout, const NEWMAT::Matrix& X);
void FPrint(ostream *fout, const NEWMAT::SymmetricMatrix& X);
void FPrint(ostream *fout, const NEWMAT::DiagonalMatrix& X);

} // namespace OPTPP

#endif
