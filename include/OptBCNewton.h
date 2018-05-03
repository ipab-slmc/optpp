#ifndef OptBCNewton_h
#define OptBCNewton_h

/*----------------------------------------------------------------------
 Copyright (c) 2001, Sandia Corporation.
 J.C. Meza, Sandia National Laboratories, meza@ca.sandia.gov
 ----------------------------------------------------------------------*/

#include "OptBCNewtonLike.h"

using std::cerr;

namespace OPTPP {

/**
 * OptBCNewton is a derived class of OptBCNewtonLike.
 * OptBCNewton implements a bound constrained Newton method. 
 * These methods will use the active set method.
 *
 * @author J.C. Meza, Lawrence Berkeley National Laboratory
 * @note Modified by P.J. Williams, pwillia@sandia.gov
 * @date Last modified 03/2007
 */

class OptBCNewton: public OptBCNewton2Deriv {
protected:
  /// Number of variables in active set
  int 		nactive;  	
  /// Working set of variables
  NEWMAT::ColumnVector  work_set;  	

 public:
 /**
  * Default Constructor
  * @see OptBCNewton(NLP2* p)
  * @see OptBCNewton(NLP2* p, UPDATEFCN u)
  * @see OptBCNewton(NLP2* p, TOLS t)
  */
  OptBCNewton(): 
    OptBCNewton2Deriv(), nactive(0), work_set(0)
    { //cerr << "OptBCNewton :: instantiation \n";
      strcpy(method,"Bound constrained Newton");
    }

 /**
  * @param p a pointer to an NLP1.
  * @see OptBCNewton(NLP2* p, UPDATEFCN u)
  * @see OptBCNewton(NLP2* p, TOLS t)
  */
  OptBCNewton(NLP2* p): 
    OptBCNewton2Deriv(p), nactive(0), work_set(p->getDim())
    { strcpy(method,"Bound constrained Newton"); work_set = false; }

 /**
  * @param p a pointer to an NLP1.
  * @param u a function pointer.
  * @see OptBCNewton(NLP2* p)
  * @see OptBCNewton(NLP2* p, TOLS t)
  */
  OptBCNewton(NLP2* p, UPDATEFCN u): 
    OptBCNewton2Deriv(p, u), nactive(0), work_set(p->getDim())
    { strcpy(method,"Bound constrained Newton"); work_set = false; }

 /**
  * @param p a pointer to an NLP1.
  * @param t tolerance class reference.
  * @see OptBCNewton(NLP2* p)
  * @see OptBCNewton(NLP2* p, UPDATEFCN u)
  */
  OptBCNewton(NLP2* p, TOLS t): 
    OptBCNewton2Deriv(p, t), nactive(0), work_set(p->getDim())
    { strcpy(method,"Bound constrained Newton"); work_set = false; }

 /*
  * Destructor
  */
  virtual ~OptBCNewton(){;}

//----------------------------------
// These are defined elsewhere
//----------------------------------
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
  virtual void	          reset();

protected:
//  NLP2*   nlprob2() const {return nlp; }
};

} // namespace OPTPP

#endif
