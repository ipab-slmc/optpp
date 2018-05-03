#ifndef OptBCQNewton_h
#define OptBCQNewton_h

/*----------------------------------------------------------------------
 Copyright (c) 2001, Sandia Corporation.
 J.C. Meza, Sandia National Laboratories, meza@ca.sandia.gov
 ----------------------------------------------------------------------*/

#include "OptBCNewtonLike.h"

using std::cerr;

namespace OPTPP {

/**
 * OptBCQNewton is a derived class of OptBCNewtonLike.
 * OptBCQNewton implements a bound constrained Quasi-Newton method 
 * These methods will use the active set method.
 *
 * @author J.C. Meza, Lawrence Berkeley National Laboratory
 * @note Modified by P.J. Williams, pwillia@sandia.gov
 * @date Last modified 03/2007
 */

class OptBCQNewton: public OptBCNewton1Deriv {
 protected:
  int 			nactive; 	///< Number of variables in the active set
  NEWMAT::ColumnVector	work_set;	///< Working set 

 public:
 /**
  * Default Constructor
  * @see OptBCQNewton(NLP1* p)
  * @see OptBCQNewton(NLP1* p, UPDATEFCN u)
  * @see OptBCQNewton(NLP1* p, TOLS t)
  */
  OptBCQNewton(): 
    OptBCNewton1Deriv(), nactive(0), work_set(0) 
    { //cerr << "OptBCQNewton :: instantiation \n";
      strcpy(method,"Bound constrained Quasi-Newton"); work_set = false; }
 /**
  * @param p a pointer to an NLP1.
  * @see OptBCQNewton(NLP1* p, UPDATEFCN u)
  * @see OptBCQNewton(NLP1* p, TOLS t)
  */
  OptBCQNewton(NLP1* p): 
    OptBCNewton1Deriv(p), nactive(0), work_set(p->getDim())
    { strcpy(method,"Bound constrained Quasi-Newton"); work_set = false; }
 /**
  * @param p a pointer to an NLP1.
  * @param u a function pointer.
  * @see OptBCQNewton(NLP1* p)
  * @see OptBCQNewton(NLP1* p, TOLS t)
  */
  OptBCQNewton(NLP1* p, UPDATEFCN u): 
    OptBCNewton1Deriv(p, u), nactive(0), work_set(p->getDim())
    { strcpy(method,"Bound constrained Quasi-Newton"); work_set = false; }
 /**
  * @param p a pointer to an NLP1.
  * @param t tolerance class reference.
  * @see OptBCQNewton(NLP1* p)
  * @see OptBCQNewton(NLP1* p, TOLS t)
  */
  OptBCQNewton(NLP1* p, TOLS t): 
    OptBCNewton1Deriv(p, t), nactive(0), work_set(p->getDim())
    { strcpy(method,"Bound constrained Quasi-Newton"); work_set = false; }

 /**
  * Destructor
  */
  virtual ~OptBCQNewton(){;}

  //-------------------------------------------
  // These are defined elsewhere
  //-------------------------------------------

  virtual int             checkConvg();
  virtual int             checkDeriv();
  virtual void            initHessian();
  virtual void            initOpt();
  virtual void            printStatus(char *);
  virtual real            stepTolNorm() const;
  NEWMAT::SymmetricMatrix updateH(NEWMAT::SymmetricMatrix& H, int k);
  double          computeMaxStep(NEWMAT::ColumnVector&);
  NEWMAT::ColumnVector    computeSearch(NEWMAT::SymmetricMatrix &);
  int             updateConstraints(int);
  virtual void            reset();

};

} // namespace OPTPP
#endif
