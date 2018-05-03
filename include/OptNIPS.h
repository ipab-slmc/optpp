#ifndef OptNIPS_h
#define OptNIPS_h

#ifndef OptNIPSLike_h
#include "OptNIPSLike.h"
#endif

namespace OPTPP {

/**
 *  OptNIPS is a derived class of OptNIPSLike.
 *  This class implements a Newton nonlinear interior-point method
 *  with analytic Hessian information.
 *  @author P.J. Williams
 */

class OptNIPS: public OptNIPSLike {
 private:
  NLP2*		nlp; ///< a pointer to an NLP2 object

 protected:
  /**
   * @return Pointer to an NLP2 object
   */
  NLP2* nlprob2() const {return nlp;}
  /**
   * @return Pointer to an NLP1 object
   */
  NLP1* nlprob()  const {return nlp;}

 public:

 /**
  * Default Constructor
  * @see OptNIPS(NLP2* p)
  * @see OptNIPS(NLP2* p, UPDATEFCN u)
  * @see OptNIPS(NLP2* p, TOLS t)
  */
  OptNIPS(): OptNIPSLike(), nlp(0)
    {strcpy(method,"Nonlinear Interior-Point Method");}
 /**
  * @param p a pointer to an NLP2.
  */
  OptNIPS(NLP2* p): OptNIPSLike(p->getDim()), nlp(p) 
    {strcpy(method,"Nonlinear Interior-Point Method");}
 /**
  * @param p a pointer to an NLP2.
  * @param u a function pointer.
  */
  OptNIPS(NLP2* p, UPDATEFCN u): OptNIPSLike(p->getDim(),u), nlp(p) 
    {strcpy(method,"Nonlinear Interior-Point Method"); }
 /**
  * @param p a pointer to an NLP2.
  * @param t tolerance class reference.
  */
  OptNIPS(NLP2* p, TOLS t): OptNIPSLike(p->getDim(),t), nlp(p)
    {strcpy(method,"Nonlinear Interior-Point Method"); }

 /**
  * Destructor
  */
  virtual ~OptNIPS(){}

  /// Initialize Hessian of the Lagrangian
  virtual void initHessian();
  /// Compute analytic Hessian of the Lagrangian
  virtual NEWMAT::SymmetricMatrix updateH(NEWMAT::SymmetricMatrix& H, int k);
  /// Print status of the NIPS method
  virtual void printStatus(char *s);
};

} // namespace OPTPP
#endif
