#ifndef OptConstrNewton_h
#define OptConstrNewton_h

/*----------------------------------------------------------------------
 Copyright (c) 2001, Sandia Corporation.
 J.C. Meza, Sandia National Laboratories, meza@ca.sandia.gov
 ----------------------------------------------------------------------*/

#ifndef OptConstrNewtonLike_h
#include "OptConstrNewtonLike.h"
#endif

namespace OPTPP {

/**
 * OptConstrNewton is a derived class of OptConstrNewtonLike.
 * This class implements an constrained Newton Method
 * with analytic Hessian information. 
 * The user can select from the following globalization strategies: 
 * Linesearch, trust-region, and trustpds.
 *
 * @author J.C. Meza, Sandia National Laboratories, meza@ca.sandia.gov
 * @note Modified by P.J. Williams, pwillia@sandia.gov
 */

class OptConstrNewton: public OptConstrNewton2Deriv {
public:
 /**
  * Default Constructor
  * @ see OptConstrNewton(NLP2* p)
  * @ see OptConstrNewton(NLP2* p, UPDATEFCN u)
  * @ see OptConstrNewton(NLP2* p, TOLS t)
  */
  OptConstrNewton() {strcpy(method,"Constrained Newton");}

 /**
  * @ param p a pointer to an NLP2.
  * @ see OptConstrNewton(NLP2* p, UPDATEFCN u)
  * @ see OptConstrNewton(NLP2* p, TOLS t)
  */
  OptConstrNewton(NLP2* p): OptConstrNewton2Deriv(p)
  {strcpy(method,"Constrained Newton");}

 /**
  * @ param p a pointer to an NLP2.
  * @ param u a function pointer.
  * @ see OptConstrNewton(NLP2* p)
  * @ see OptConstrNewton(NLP2* p, TOLS t)
  */
  OptConstrNewton(NLP2* p, UPDATEFCN u): OptConstrNewton2Deriv(p,u)
      {strcpy(method,"Constrained Newton");}
 /**
  * @ param p a pointer to an NLP2.
  * @ param t tolerance class reference.
  * @ see OptConstrNewton(NLP2* p)
  * @ see OptConstrNewton(NLP2* p, UPDATEFCN u)
  */
  OptConstrNewton(NLP2* p, TOLS t): OptConstrNewton2Deriv(p,t)
    {strcpy(method,"Constrained Newton");}
  
 /**
  * Destructor
  */
  virtual ~OptConstrNewton(){}

// These are defined elsewhere

  /// Initialie the Hessian 
  void initHessian();

  /// Compute the analytic Hessian 
  NEWMAT::SymmetricMatrix updateH(NEWMAT::SymmetricMatrix& H, int k);

  /// Compare the analytic gradient with finite-difference gradient
  int checkDeriv();

  /// Print status of constrained Newton's method
  void printStatus(char *);

  /// Compute length of step direction
  real stepTolNorm() const;
};

} // namespace OPTPP
#endif
