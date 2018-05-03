#ifndef optbcnewton1_h
#define optbcnewton1_h

/*----------------------------------------------------------------------
 Copyright (c) 2001, Sandia Corporation.
 J.C. Meza, Sandia National Laboratories, meza@ca.sandia.gov
 ----------------------------------------------------------------------*/


#include "OptBCNewtonLike.h"

namespace OPTPP {

/**
 * OptBaQNewton implements a barrier Quasi-Newton method.  The only allowable
 * constraints are simple bounds.  Barrier methods transform a constrained 
 * into an unconstrained problem via a barrier term.  In this implementation, 
 * the unconstrained subproblem's objective function consists of the original 
 * objective function and a logarithmic barrier term that prevents iterates 
 * from violating the bound constraints.
 *
 * @author J.C. Meza, Lawrence Berkeley National Laboratory
 * @note   Modified by P.J. Williams, pwillia@sandia.gov
 * @date Last modified 11/2005
 */

class OptBaQNewton: public OptBCNewton1Deriv {

  /// Barrier parameter 
  double          mu;
  /// Value of barrier function
  double          fvalue_barrier, fprev_barrier, fprev_outer;
  /// Gradient of barrier function
  NEWMAT::ColumnVector    grad_barrier, gprev_barrier; 

 public:
  /**
   * Default Constructor
   * @see OptBaQNewton(NLP1* p)
   * @see OptBaQNewton(NLP1* p, UPDATEFCN u)
   * @see OptBaQNewton(NLP1* p, TOLS t)
   */
  OptBaQNewton(){strcpy(method,"Bound constrained QNewton with barrier");}
  /**
   * @param p a pointer to an NLP1.
   * @see OptBaQNewton(NLP1* p, UPDATEFCN u)
   * @see OptBaQNewton(NLP1* p, TOLS t)
   */
  OptBaQNewton(NLP1* p): OptBCNewton1Deriv(p)
    {strcpy(method,"Bound constrained QNewton with barrier");}
  /**
   * @param p a pointer to an NLP1.
   * @param u a function pointer.
   * @see OptBaQNewton(NLP1* p)
   * @see OptBaQNewton(NLP1* p, TOLS t)
   */
  OptBaQNewton(NLP1* p, UPDATEFCN u): OptBCNewton1Deriv(p, u)
    {strcpy(method,"Bound constrained QNewton with barrier"); }
  /**
   * @param p a pointer to an NLP1.
   * @param t tolerance class reference.
   * @see OptBaQNewton(NLP1* p)
   * @see OptBaQNewton(NLP1* p, UPDATEFCN u)
   */
  OptBaQNewton(NLP1* p, TOLS t): OptBCNewton1Deriv(p, t)
    {strcpy(method,"Bound constrained QNewton with barrier"); }

  /**
   * Destructor
   */
  virtual ~OptBaQNewton(){}

//-----------------------------------------
// These are defined elsewhere
//-----------------------------------------
  void            initOpt();
  int             checkInnerConvg(int);
  int             checkConvg();
  void            acceptStep(int,int);
  void            optimize();
  void            updateBarrierMultiplier();
  NEWMAT::SymmetricMatrix updateH(NEWMAT::SymmetricMatrix& H, int k);
  double          compute_Barrier_Fvalue(double,NEWMAT::ColumnVector&);
  NEWMAT::ColumnVector    compute_Barrier_Gradient(NEWMAT::ColumnVector&,
    NEWMAT::ColumnVector&);
  NEWMAT::SymmetricMatrix compute_Barrier_Hessian(NEWMAT::SymmetricMatrix&,
    NEWMAT::ColumnVector&);
  NEWMAT::ColumnVector    computeSearch2(NEWMAT::SymmetricMatrix&, 
    NEWMAT::ColumnVector&);
  int             computeStep(NEWMAT::ColumnVector );
  double          computeMaxStep(NEWMAT::ColumnVector &);
  double          scalarNewton(double,double,double,double,double);
  void            setAsideCurrentVariables();
};

} // namespace OPTPP

#endif
