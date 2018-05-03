#ifndef OptDHNIPS_h
#define OptDHNIPS_h

#ifndef OptNIPSLike_h
#include "OptNIPSLike.h"
#endif

#include "OptppArray.h"

namespace OPTPP {

/**
 *  OptDHNIPS is a derived class of OptNIPSLike.
 *  This class implements a disaggregated Hessian approximation 
 *  nonlinear interior-point method with either an NLF2 or least 
 *  squares function operator and Quasi-Newton approximations to 
 *  constraint Hessians.
 *  @author P.J. Williams
 *  @date Last modified 03/2007
 */

class OptDHNIPS: public OptNIPSLike {
 private:
  NLP2*		nlp; 			///< a pointer to an NLP2 object

 protected:
  OptppArray<NEWMAT::SymmetricMatrix> HCk_;  	///< Array of constraint Hessians
  OptppArray<int> indices;  	        ///< Indices of nonlinear constraints 

  /**
   * @return Pointer to an NLP2 object
   */
  NLP2* nlprob2() const {return nlp;}
  /**
   * @return Pointer to an NLP1 object
   */
  NLP1* nlprob()  const {return nlp;}


 public:

 /**
  * Default Constructor
  * @see OptDHNIPS(NLP2* p)
  * @see OptDHNIPS(NLP2* p, UPDATEFCN u)
  * @see OptDHNIPS(NLP2* p, TOLS t)
  */
  OptDHNIPS(): OptNIPSLike(), nlp(0), HCk_(0), indices(0)
    {strcpy(method,"Nonlinear Interior-Point Method w/ Disaggregated Hessian");}
 /**
  * @param p a pointer to an NLP2.
  */
  OptDHNIPS(NLP2* p): OptNIPSLike(p->getDim()), nlp(p), 
    HCk_(0), indices(0)
    {strcpy(method,"Nonlinear Interior-Point Method w/ Disaggregated Hessian");}
 /**
  * @param p a pointer to an NLP2.
  * @param u a function pointer.
  */
  OptDHNIPS(NLP2* p, UPDATEFCN u): OptNIPSLike(p->getDim(),u), nlp(p), 
    HCk_(0), indices(0)
    {strcpy(method,"Nonlinear Interior-Point Method w/ Disaggregated Hessian");}
 /**
  * @param p a pointer to an NLP2.
  * @param t tolerance class reference.
  */
  OptDHNIPS(NLP2* p, TOLS t): OptNIPSLike(p->getDim(),t), nlp(p), 
    HCk_(0), indices(0)
    {strcpy(method,"Nonlinear Interior-Point Method w/ Disaggregated Hessian");}

 /**
  * Destructor
  */
  virtual ~OptDHNIPS(){}

  /**
   *    @return Array which contains quasi-Newton Hessians of the constraints
   */ 
  OptppArray<NEWMAT::SymmetricMatrix> getConstraintHessian() const { return HCk_;}

  /**
   *    @return Array which contains indices of nonlinear constraints
   */ 

  void nonLinearConstraintIndices(const NEWMAT::ColumnVector& types);

  /// Reset parameters 
  virtual void reset();
  /// Initialize Hessian of the Lagrangian
  virtual void initHessian();
  /// Compute disaggregated approximation to the Hessian of the Lagrangian 
  virtual NEWMAT::SymmetricMatrix updateH(NEWMAT::SymmetricMatrix& H, int k);
  /// Print status of the disaggregated Hessian nonlinear interior-point method
  virtual void printStatus(char *s);
};

} // namespace OPTPP
#endif
