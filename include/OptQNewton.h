#ifndef OptQNewton_h
#define OptQNewton_h

#ifndef OptNewtonLike_h
#include "OptNewtonLike.h"
#endif


namespace OPTPP {

/**
 * OptQNewton is a derived class of OptNewtonLike.
 * This class implements an unconstrained Quasi-Newton Method
 * with BFGS approximation to the Hessian.  The user can select
 * from the following globalization strategies: linesearch, trust-region,
 * and trustpds.
 *
 * @author J.C. Meza, Sandia National Laboratories,meza@ca.sandia.gov
 * @note Modified by P.J. Williams, pwillia@sandia.gov 
 */

class OptQNewton: public OptNewton1Deriv {
 public:
  /**
   * Default Constructor
   * @see OptQNewton(NLP1* p)
   * @see OptQNewton(NLP1* p, UPDATEFCN u)
   * @see OptQNewton(NLP1* p, TOLS t)
   */
  OptQNewton(){strcpy(method,"Quasi-Newton");}
  /**
   * @param p a pointer to an NLP1.
   */
  OptQNewton(NLP1* p): OptNewton1Deriv(p)
    {strcpy(method,"Quasi-Newton");}
  /**
   * @param p a pointer to an NLP1.
   * @param u a function pointer.
   */
  OptQNewton(NLP1* p, UPDATEFCN u): OptNewton1Deriv(p, u)
    {strcpy(method,"Quasi-Newton"); }
  /**
   * @param p a pointer to an NLP1.
   * @param t tolerance class reference.
   */
  OptQNewton(NLP1* p, TOLS t): OptNewton1Deriv(p, t)
    {strcpy(method,"Quasi-Newton"); }

  /**
   * Destructors
   */
  virtual ~OptQNewton(){}

//------------------------------------------------
// These are defined elsewhere
//------------------------------------------------

 /// Compute BFGS approximation to the Hessian of the objective function
  NEWMAT::SymmetricMatrix updateH(NEWMAT::SymmetricMatrix& H, int k);

  /// Compare the analytic gradient with the finite difference gradient
  int checkDeriv();

 // virtual double initTrustRegionSize() const;
};

} // namespace OPTPP
#endif
