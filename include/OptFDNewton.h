#ifndef OptFDNewton_h
#define OptFDNewton_h

/*---------------------------------------------------------------------
 Copyright (c) 2001, Sandia Corporation
 J.C. Meza, Sandia National Laboratories, meza@ca.sandia.gov
 ---------------------------------------------------------------------- */

#ifndef OptNewtonLike_h
#include "OptNewtonLike.h"
#endif

namespace OPTPP {

/**
 * OptFDNewton is a derived class of OptNewtonLike.
 * This class implements an unconstrained Newton's Method
 * with a finite-difference approximation to the Hessian.  
 * The user can select from the following globalization strategies: 
 * linesearch, trust-region, and trustpds.
 *
 * @author J.C. Meza, Sandia National Laboratories, meza@ca.sandia.gov
 * @note Modified by P.J. Williams 
 */

class OptFDNewton: public OptNewton1Deriv {
 public:

 /**
  * Default Constructor
  * @see OptFDNewton(NLP1* p)
  * @see OptFDNewton(NLP1* p, UPDATEFCN u)
  * @see OptFDNewton(NLP1* p, TOLS t)
  */
  OptFDNewton(){strcpy(method,"Finite-Difference Newton");}
 /**
  * @param p a pointer to an NLP1.
  */
  OptFDNewton(NLP1* p): OptNewton1Deriv(p)
    {strcpy(method,"Finite-Difference Newton");}
 /**
  * @param p a pointer to an NLP1.
  * @param u a function pointer.
  */
  OptFDNewton(NLP1* p, UPDATEFCN u): OptNewton1Deriv(p, u)
    {strcpy(method,"Finite-Difference Newton");}
 /**
  * @param p a pointer to an NLP1.
  * @param t tolerance class reference.
  */
  OptFDNewton(NLP1* p, TOLS t): OptNewton1Deriv(p,t)
    {strcpy(method,"Finite-Difference Newton");}

 /**
  * Destructors
  */ 
  virtual ~OptFDNewton(){}

  /// Compute finite-difference approximation to Hessian of the objective function 
  NEWMAT::SymmetricMatrix updateH(NEWMAT::SymmetricMatrix& H, int k);

  /// Compare the analytic gradient with the finite difference gradient
  int checkDeriv();
};

} // namespace OPTPP

#endif
