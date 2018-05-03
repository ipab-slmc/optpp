#ifndef OptNIPSLike_h
#define OptNIPSLike_h

/*----------------------------------------------------------------------
 Copyright (c) 2001, Sandia Corporation.   Under the terms of Contract 
 DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 work by or on behalf of the U.S. Government.

 P.J. Williams, Sandia National Laboratories, pwillia@ca.sandia.gov
 ----------------------------------------------------------------------*/

#ifndef OptConstrNewtonLike_h
#include "OptConstrNewtonLike.h"
#endif

namespace OPTPP {

/**
 * @class OptNIPSLike 
 * OptNIPSLike is a derived class of OptConstrNewtonLike.
 * OptNIPSLike provides common data and functionality for OptFDNIPS,
 * OptQNIPS, and OptNIPS.
 *
 * The OptNIPS algorithm is a C++ implementation of NIPSM, a nonlinear
 * interior-point code developed under MATLAB by Amr El-Bakry at Rice
 * University and NIPSF, a Fortran implementation of the same code written
 * by Frederik Saaf.  Additional features include the merit functions 
 * proposed by Miguel Argaez and Richard Tapia in "Global Convergence of a
 * Primal-Dual Newton Interior-Point Method for Nonlinear Programming
 * Using a Modified Augmented Lagrange Function" as well as 
 * Robert Vanderbei and David Shanno in "An Interior-Point Algorithm For
 * Nonconvex Nonlinear Programming".
 *
 * @author P.J. Williams, Sandia National Laboratories, pwillia@sandia.gov
 */

class OptNIPSLike: public OptConstrNewtonLike {
 protected:
  virtual NLP1* nlprob() const= 0; ///< pointer to NLP1
  real          beta_;  ///< penalty parameter for merit function 3
  real          dirder_; ///< directional derivative of a merit function
  real          mu_; 	///< pertubation parameter 
  real          penalty_; ///< penalty parameter for merit function 2
  real		sigmin_; ///< centering parameter 
  real		taumin_; ///< percentage of steplength to boundary
  const real	rho_;    ///< constant set to .5 
  const real	sw_;	///<  constant

 public:
 /**
  * Default Constructor
  * @see OptNIPSLike(int n)
  * @see OptNIPSLike(int n, UPDATEFCN u)
  * @see OptNIPSLike(int n, TOLS t)
  */
  OptNIPSLike(): OptConstrNewtonLike(), beta_(0.0e0), 
    dirder_(0.0e0), mu_(0.0e0), penalty_(0.0e0),
    sigmin_(0.0e0), taumin_(0.0e0), rho_(0.0e0), sw_(0.0e0)
    {strcpy(method,"Nonlinear Interior-Point Method");}
 /**
  * @param n an integer argument.
  */
  OptNIPSLike(int n): OptConstrNewtonLike(n), beta_(0.0e0), 
    dirder_(0.0e0), mu_(0.0e0), penalty_(1.0e2),
    sigmin_(0.1e0), taumin_(0.95e0), rho_(5.0e-1), sw_(1.0e2)
    { strcpy(method,"Nonlinear Interior-Point Method"); }
 /**
  * @param n an integer argument.
  * @param u a function pointer.
  */
  OptNIPSLike(int n, UPDATEFCN u): OptConstrNewtonLike(n,u), beta_(0.0e0), 
    dirder_(0.0e0), mu_(0.0e0), penalty_(1.0e2),
    sigmin_(0.1e0), taumin_(0.95e0), rho_(5.0e-1), sw_(1.0e2)
    { strcpy(method,"Nonlinear Interior-Point Method"); }
 /**
  * @param n an integer argument.
  * @param t tolerance class reference.
  */
  OptNIPSLike(int n, TOLS t): OptConstrNewtonLike(n,t), beta_(0.0e0), 
    dirder_(0.0e0), mu_(0.0e0), penalty_(1.0e2),
    sigmin_(0.1e0), taumin_(0.95e0), rho_(5.0e-1), sw_(1.0e2)
    { strcpy(method,"Nonlinear Interior-Point Method"); }

 /**
  * Destructor
  */
  virtual ~OptNIPSLike(){}

//-------------------------------------------------------------------
// Virtual Functions 
//-------------------------------------------------------------------

  virtual void setMeritFcn(MeritFcn option);
  virtual NEWMAT::Matrix setupMatrix(const NEWMAT::ColumnVector& xc);
  virtual real merit(int flag, const NEWMAT::ColumnVector& xc, const NEWMAT::ColumnVector& yc,
                     NEWMAT::ColumnVector& zc, NEWMAT::ColumnVector& sc);
  virtual NEWMAT::ColumnVector setupRHS(const NEWMAT::ColumnVector& xc, real mu);
  virtual NEWMAT::ColumnVector setupRHS(const NEWMAT::ColumnVector& xplus,
                                const NEWMAT::ColumnVector& yplus,
				const NEWMAT::ColumnVector& zplus,
				const NEWMAT::ColumnVector& splus, real mu);

  virtual NEWMAT::SymmetricMatrix updateH(NEWMAT::SymmetricMatrix& H, int k) = 0;

  virtual int checkConvg();
  virtual int checkDeriv();
  virtual int computeStep(NEWMAT::ColumnVector step);

  virtual void initOpt();		///< Initialize algorithmic parameters
  virtual void initHessian();		///< Initialize Hessian of Lagrangian 
  virtual void optimize();		///< Call the interior-point method
  virtual void printStatus(char *s);	///< Print status of opt. method
  virtual void readOptInput();		///< Read user-specified input options

//-------------------------------------------------------------------
// Accessor Methods
//-------------------------------------------------------------------

/**
 * @return The value of mu_, the pertubation parameter
 */
  real getMu() const { return mu_;}
/**
 * Set the value of the perturbation parameter
 */
  void setMu(real newMu) { mu_ = newMu;}

/**
 * Sets the value of the centering parameter.
 */
  void setCenteringParameter(real newSigma) { sigmin_ = newSigma;}

/**
 * Sets the percentage of step taken towards the boundary 
 */
  void setStepLengthToBdry(real newTau) { taumin_ = newTau;}

//-------------------------------------------------------------------
// These are used by the derived classes 
//-------------------------------------------------------------------

 /**
  * Takes zero arguments with void return.
  * Resets parameter values.
  */
  virtual void reset();

  void recoverFeasibility(NEWMAT::ColumnVector xinit, CompoundConstraint* constraints, 
                          double ftol);
  NEWMAT::ColumnVector computeSearch2(NEWMAT::Matrix& Jacobian, const NEWMAT::ColumnVector& rhs);
  /**
   * Takes two arguments and returns a NEWMAT::ColumnVector.
   * @param df a NEWMAT::ColumnVector - gradient of obj. function
   * @param dcon a Matrix  - gradient of constraints
   * @return The initial value of the Lagrange multiplier z 
   */
  NEWMAT::ColumnVector initMultipliers(const NEWMAT::ColumnVector& df, NEWMAT::Matrix& dcon);

  /**
   * Takes one arguments and updates the perturbation parameter
   * @param k an integer - iteration counter 
   */
  void updateMu(int k);

  /**
   * Takes five arguments and returns a real value.
   * @param flag an integer argument
   * @param xc a NEWMAT::ColumnVector 
   * @param yc a NEWMAT::ColumnVector of Lagrange multipliers 
   * @param zc a NEWMAT::ColumnVector of Lagrange multipliers 
   * @param sc a NEWMAT::ColumnVector of slack variables
   * @see merit(flag,xc,yc,zc,sc)
   * @return The value of the Argaez-Tapia merit function.
   */
  real merit2(int flag, const NEWMAT::ColumnVector& xc, const NEWMAT::ColumnVector& yc,
              NEWMAT::ColumnVector& zc, NEWMAT::ColumnVector& sc);
  /**
   * Takes four arguments and returns a real value.
   * @param flag an integer argument
   * @param xc a NEWMAT::ColumnVector 
   * @param zc a NEWMAT::ColumnVector of Lagrange multipliers 
   * @param sc a NEWMAT::ColumnVector of slack variables
   * @see merit(flag,xc,yc,zc,sc)
   * @return The value of the Vanderbei et al merit function.
   */
  real merit3(int flag, const NEWMAT::ColumnVector& xc, NEWMAT::ColumnVector& zc,
              NEWMAT::ColumnVector& sc);

  /**
   * Takes three arguments and void return.
   * @param sk a NEWMAT::ColumnVector with contains the search direction 
   * @param xc a NEWMAT::ColumnVector of current point 
   * @param derivative a NEWMAT::ColumnVector of derivative of cost function 
   */
  void computeDirDeriv(NEWMAT::ColumnVector& sk, const NEWMAT::ColumnVector& xc,
                       NEWMAT::ColumnVector& derivative);
  /**
   * Takes one arguments and returns a real value.
   * @param step a NEWMAT::ColumnVector with contains the search direction 
   */
  double dampenStep(NEWMAT::ColumnVector& step);

};

} // namespace OPTPP

#endif
