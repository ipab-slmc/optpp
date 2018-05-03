#ifndef OptNewton_h
#define OptNewton_h

#ifndef OptNewtonLike_h
#include "OptNewtonLike.h"
#endif

namespace OPTPP {

/**
 * OptNewton is a derived class of OptNewtonLike. 
 * This class implements an unconstrained Newton Method
 * with analytic Hessian information.  The user can select
 * from the following globalization strategies: Linesearch, 
 * trust-region, and trustpds.
 *
 * Copyright (c) 2001, Sandia Corporation.
 * @author J.C. Meza, Sandia National Laboratories,meza@ca.sandia.gov
 * @note Modified by P.J. Williams 
 */

class OptNewton: public OptNewton2Deriv {
public:
/**
 * Default Constructor
 * @see OptNewton(NLP2* p)
 * @see OptNewton(NLP2* p, UPDATEFCN u)
 * @see OptNewton(NLP2* p, TOLS t)
 */

  OptNewton(): OptNewton2Deriv() 
    {strcpy(method,"Newton");}
/**
 * @param p a pointer to an NLP2.
 */
  OptNewton(NLP2* p): OptNewton2Deriv(p)
    {strcpy(method,"Newton");}
/**
 * @param p a pointer to an NLP2.
 * @param u a function pointer.
 */
  OptNewton(NLP2* p, UPDATEFCN u): OptNewton2Deriv(p,u)
    {strcpy(method,"Newton");}
/**
 * @param p a pointer to an NLP2.
 * @param t tolerance class reference.
 */
  OptNewton(NLP2* p, TOLS t): OptNewton2Deriv(p,t)
    {strcpy(method,"Newton");}
  
/**
 * Destructor
 */
  virtual ~OptNewton(){}

//-----------------------------------------------------------------
// These are defined elsewhere
//-----------------------------------------------------------------

  /// Compute the analytic Hessian at the initial point
  void initHessian();

  /// Returns the analytic Hessian
  NEWMAT::SymmetricMatrix updateH(NEWMAT::SymmetricMatrix& H, int k);

  /// Compare the analytic gradient with the finite difference gradient
  int checkDeriv();

  /// Print status of Newton's method
  void printStatus(char *);

  /// Compute length of the step direction 
  real stepTolNorm() const;
};

} // namespace OPTPP

#endif
