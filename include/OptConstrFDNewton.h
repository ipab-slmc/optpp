#ifndef OptConstrFDNewton_h
#define OptConstrFDNewton_h

/*----------------------------------------------------------------------
 Copyright (c) 2001, Sandia Corporation.
 J.C. Meza, Sandia National Laboratories, meza@ca.sandia.gov
 ----------------------------------------------------------------------*/

#ifndef OptConstrNewtonLike_h
#include "OptConstrNewtonLike.h"
#endif

namespace OPTPP {

/**
 * OptConstrFDNewton is a derived class of OptConstrNewtonLike.
 * This class implements a Newton Method
 * with finite-difference approximation to the Hessian.
 *
 * @author J.C. Meza, Sandia National Laboratories, meza@ca.sandia.gov
 * @note Modified by P.J. Williams, pwillia@sandia.gov
 * @date Last modified 9/2001
 */

class OptConstrFDNewton: public OptConstrNewton1Deriv {
 public:
 /**
  * Default Constructor
  * @see OptConstrFDNewton(NLP1* p)
  * @see OptConstrFDNewton(NLP1* p, UPDATEFCN u)
  * @see OptConstrFDNewton(NLP1* p, TOLS t)
  */
  OptConstrFDNewton(){strcpy(method,"Constrained Finite-Difference Newton");}
 /**
  * @param p a pointer to an NLP1.
  * @see OptConstrFDNewton(NLP1* p, UPDATEFCN u)
  * @see OptConstrFDNewton(NLP1* p, TOLS t)
  */
  OptConstrFDNewton(NLP1* p): OptConstrNewton1Deriv(p)
    {strcpy(method,"Constrained Finite-Difference Newton");}
 /**
  * @param p a pointer to an NLP1.
  * @param u a function pointer.
  * @see OptConstrFDNewton(NLP1* p)
  * @see OptConstrFDNewton(NLP1* p, TOLS t)
  */
  OptConstrFDNewton(NLP1* p, UPDATEFCN u): OptConstrNewton1Deriv(p, u)
    {strcpy(method,"Constrained Finite-Difference Newton");}
 /**
  * @param p a pointer to an NLP1.
  * @param t tolerance class reference.
  * @see OptConstrFDNewton(NLP1* p)
  * @see OptConstrFDNewton(NLP1* p, UPDATEFCN u)
  */
  OptConstrFDNewton(NLP1* p, TOLS t): OptConstrNewton1Deriv(p,t)
    {strcpy(method,"Constrained Finite-Difference Newton");}

 /**
  * Destructor
  */
  virtual ~OptConstrFDNewton(){}

// These are defined elsewhere

/// Compute finite-difference Hessian 
  NEWMAT::SymmetricMatrix updateH(NEWMAT::SymmetricMatrix& H, int k);

/// Compare the analytic gradient with the finite-difference gradient
  int checkDeriv();
};

} // namespace OPTPP

#endif
