#ifndef OptQNIPS_h
#define OptQNIPS_h


#ifndef OptNIPSLike_h
#include "OptNIPSLike.h"
#endif

namespace OPTPP {

/**
 *  OptQNIPS is a derived class of OptNIPSLike.
 *  This class implements a quasi-Newton nonlinear interior-point method.
 *  The Hessian of the Lagrangrian is approximated by an BFGS update.
 *  @author P.J. Williams
 */

class OptQNIPS: public OptNIPSLike {
 private:
  NLP1* nlp;  /// pointer to an NLP1

 protected:
  NLP1* nlprob() const {return nlp;}

 public:
  /**
   * Default Constructor
   * @ see OptQNIPS(NLP1* p)
   * @ see OptQNIPS(NLP1* p, UPDATEFCN u)
   * @ see OptQNIPS(NLP1* p, TOLS t)
   */
  OptQNIPS(): OptNIPSLike(),  nlp(0)
    {strcpy(method,"Nonlinear Interior-Point Method");}

  /**
   * @ param p a pointer to an NLP1
   */
  OptQNIPS(NLP1* p): OptNIPSLike(p->getDim()),  nlp(p)
    {strcpy(method,"Nonlinear Interior-Point Method");}

  /**
   * @ param p a pointer to an NLP1.
   * @ param u a function pointer.
   */
  OptQNIPS(NLP1* p, UPDATEFCN u): OptNIPSLike(p->getDim(),u), nlp(p) 
    {strcpy(method,"Nonlinear Interior-Point Method"); }

  /**
   * @ param p a pointer to an NLP1.
   * @ param t tolerance class reference.
   */
  OptQNIPS(NLP1* p, TOLS t): OptNIPSLike(p->getDim(),t),  nlp(p) 
    {strcpy(method,"Nonlinear Interior-Point Method"); }

  /**
   * Destructors
   */
  virtual ~OptQNIPS() {;}

  /// Compare the analytic gradient with the finite-difference approximation
  virtual int checkDeriv();
  /// Compute BFGS approximation to the Hessian of the Lagrangian 
  virtual NEWMAT::SymmetricMatrix updateH(NEWMAT::SymmetricMatrix& H, int k);

};
} // namespace OPTPP
#endif
