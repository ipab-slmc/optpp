#ifndef OptConstrQNewton_h
#define OptConstrQNewton_h

/*----------------------------------------------------------------------
 Copyright (c) 2001, Sandia Corporation.
 J.C. Meza, Sandia National Laboratories, meza@ca.sandia.gov
 ----------------------------------------------------------------------*/

#ifndef OptConstrNewtonLike_h
#include "OptConstrNewtonLike.h"
#endif

namespace OPTPP {

/**
 * OptConstrQNewton is a derived class of OptConstrNewtonLike.
 * This class implements a Constrained Quasi-Newton Method
 * with BFGS approximation to the Hessian.
 *
 * @author J.C. Meza, Lawrence Berkeley National Laboratory
 * @note Modified by P.J. Williams, pwillia@sandia.gov
 * @date 11/2005 
 */

class OptConstrQNewton: public OptConstrNewton1Deriv {
 public:
 /**
  * Default Constructor
  * @see OptConstrQNewton(NLP1* p)
  * @see OptConstrQNewton(NLP1* p, UPDATEFCN u)
  * @see OptConstrQNewton(NLP1* p, TOLS t)
  */

  OptConstrQNewton(){strcpy(method,"Constrained Quasi-Newton");}

 /**
  * @param p a pointer to an NLP1.
  */
  OptConstrQNewton(NLP1* p): OptConstrNewton1Deriv(p)
    {strcpy(method,"Constrained Quasi-Newton");}
 /**
  * @param p a pointer to an NLP1.
  * @param u a function pointer.
  */
  OptConstrQNewton(NLP1* p, UPDATEFCN u): OptConstrNewton1Deriv(p, u)
    {strcpy(method,"Constrained Quasi-Newton"); }
 /**
  * @param p a pointer to an NLP1.
  * @param t tolerance class reference.
  */
  OptConstrQNewton(NLP1* p, TOLS t): OptConstrNewton1Deriv(p, t)
    {strcpy(method,"Constrained Quasi-Newton"); }

 /**
  * Destructor
  */
  virtual ~OptConstrQNewton(){}

//----------------------------------
// These are defined elsewhere
//----------------------------------

  /// Compute BFGS appoximation to the Hessian 
  NEWMAT::SymmetricMatrix updateH(NEWMAT::SymmetricMatrix& H, int k);
  /// Compare the analytic gradient with the finite difference gradient
  int checkDeriv();

//  virtual double initTrustRegionSize() const;
};

} // namespace OPTPP
#endif
