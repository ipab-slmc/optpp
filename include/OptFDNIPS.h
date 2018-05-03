#ifndef OptFDNIPS_h
#define OptFDNIPS_h


#ifndef OptNIPSLike_h
#include "OptNIPSLike.h"
#endif

namespace OPTPP {

/**
 *  OptFDNIPS is a derived class of OptNIPSLike.
 *  This class implements a Newton nonlinear interior-point method
 *  with a finite difference approximation to the Hessian. 
 *  @author P.J. Williams, Sandia National Laboratories, pwillia@sandia.gov
 */

class OptFDNIPS: public OptNIPSLike {
 private:
  NLP1* nlp;  ///< a pointer to an NLP1 object.

 protected:
  NLP1* nlprob() const {return nlp;}

 public:
 /**
  * Default Constructor
  * @see OptFDNIPS(NLP1* p)
  * @see OptFDNIPS(NLP1* p, UPDATEFCN u)
  * @see OptFDNIPS(NLP1* p, TOLS t)
  */
  OptFDNIPS(): OptNIPSLike(),  nlp(0)
    {strcpy(method,"Nonlinear Interior-Point Method");}
 /**
  * @param p a pointer to an NLP1.
  */
  OptFDNIPS(NLP1* p): OptNIPSLike(p->getDim()),  nlp(p)
    {strcpy(method,"Nonlinear Interior-Point Method");}
 /**
  * @param p a pointer to an NLP1.
  * @param u a function pointer.
  */
  OptFDNIPS(NLP1* p, UPDATEFCN u): OptNIPSLike(p->getDim(),u), nlp(p) 
    {strcpy(method,"Nonlinear Interior-Point Method"); }
 /**
  * @param p a pointer to an NLP1.
  * @param t tolerance class reference.
  */
  OptFDNIPS(NLP1* p, TOLS t): OptNIPSLike(p->getDim(),t),  nlp(p) 
    {strcpy(method,"Nonlinear Interior-Point Method"); }

  /// Destructor
  virtual ~OptFDNIPS() {;}

  /// Compare the analytic gradient to the finite-difference gradient
  virtual int checkDeriv();
  /// Computes finite-difference Hessian of the Lagrangian
  virtual NEWMAT::SymmetricMatrix updateH(NEWMAT::SymmetricMatrix& H, int k);

};

} // namespace OPTPP
#endif
