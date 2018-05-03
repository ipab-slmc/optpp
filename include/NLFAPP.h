#ifndef NLF0APP_h
#define NLF0APP_h

#include "NLF.h"
#include "AppLauncher.h"

using std::ostream;

namespace OPTPP {

/** 
 * These classes and typedefs are used for Application Launching where 
 * an AppLauncher Object is also required so that the launcher specific 
 * data can be used for the function evaluation.
 */


typedef void (*USERFCN0APP)(int, const NEWMAT::ColumnVector&, real&, int&, 
  AppLauncher * launcher);
typedef void (*USERNLNCON0APP)(int, int, const NEWMAT::ColumnVector&, 
  NEWMAT::ColumnVector&, int&, AppLauncher * launcher);
typedef void (*INITFCNAPP)(int, NEWMAT::ColumnVector&, AppLauncher * launcher);

/** 
 * These classes and typedefs are used for Application Launching where 
 * an AppLauncher Object is also required so that the launcher specific 
 * data can be used for the function evaluation.
 */

class NLF0APP: public NLF0 {
protected:
USERFCN0APP fcn;		///< User-defined objective function
USERNLNCON0APP confcn;		///< User-defined nonlinear constraints 
INITFCNAPP init_fcn;		///< Initializes the objective function
INITCONFCN init_confcn;		///< Initializes the constraints 
AppLauncher * launcher_;  	///< holds info for Launching a black box application
bool init_flag;			///< Has the function been initialized?

public:
// Constructors
NLF0APP(): 
  NLF0(){;}
NLF0APP(int ndim): 
  NLF0(ndim){;}
NLF0APP(int ndim, USERFCN0APP f): 
  NLF0(ndim) {fcn = f;}
NLF0APP(int ndim, USERFCN0APP f, INITFCNAPP i, AppLauncher * launcher, CompoundConstraint* constraint = 0):
  NLF0(ndim, NULL, NULL, constraint) {fcn = f; init_fcn = i; init_flag = false; launcher_ = launcher;}
NLF0APP(int ndim, USERFCN0APP f, INITFCNAPP i, INITCONFCN c):
  NLF0(ndim)
  {fcn = f; init_fcn = i; init_confcn = c; init_flag = false; 
   constraint_ = init_confcn(ndim);}				
NLF0APP(int ndim,int nlncons, USERNLNCON0APP f, INITFCNAPP i, AppLauncher *launcher):
  NLF0(ndim, nlncons, NULL, NULL) {confcn = f; init_fcn = i; init_flag = false; launcher_ = launcher;}

// Destructor
virtual ~NLF0APP() {;}               

/// Reset the parameter values 
virtual void reset(); 			

/// Initialize selected function
virtual void initFcn(); 			

/// Evaluate objective function, gradient, and Hessian 
virtual void eval(); 				

/// Evaluate objective function
virtual real evalF();               		

/// Evaluate objective function at x
virtual real evalF(const NEWMAT::ColumnVector& x); 	

/// Evaluate the nonlinear constraints at x
virtual NEWMAT::ColumnVector evalCF(const NEWMAT::ColumnVector& x); 	

// Default Gradient function for NLF0 is to do a finite-difference
/// Evaluate the finite-difference gradient
virtual NEWMAT::ColumnVector evalG();

/// Evaluate the Lagrangian at x
virtual real evalLagrangian(const NEWMAT::ColumnVector& x, NEWMAT::ColumnVector& mult, 
  const NEWMAT::ColumnVector& type) ;

/// Evaluate the gradient of the Lagrangian at x
virtual NEWMAT::ColumnVector evalLagrangianGradient(const NEWMAT::ColumnVector& x,
						const NEWMAT::ColumnVector& mult,
						const NEWMAT::ColumnVector& type) ;

private:
/// Evaluate gradient at x
virtual NEWMAT::ColumnVector evalG(const NEWMAT::ColumnVector& x);  
/// Evaluate the Hessian 
virtual NEWMAT::SymmetricMatrix evalH();		
/// Evaluate the Hessian at x
virtual NEWMAT::SymmetricMatrix evalH(NEWMAT::ColumnVector &x);
/// Evaluate the constraint gradient at x
virtual NEWMAT::Matrix evalCG(const NEWMAT::ColumnVector& x);  	
/// Evaluate constraint hessian at x
virtual NEWMAT::SymmetricMatrix evalCH(NEWMAT::ColumnVector &x);	
/// Evaluate constraint hessian at x
virtual OptppArray<NEWMAT::SymmetricMatrix> evalCH(NEWMAT::ColumnVector &x, int darg);
virtual void evalC(const NEWMAT::ColumnVector& x); 	
};

/** 
 * These classes and typedefs are used for Application Launching where 
 * an AppLauncher Object is also required so that the launcher specific 
 * data can be used for the function evaluation.
 */

class FDNLF1APP: public FDNLF1 {
protected:
USERFCN0APP fcn;		///< User-defined objective function
USERNLNCON0APP confcn;		///< User-defined nonlinear constraints
INITFCNAPP init_fcn;		///< Initializes the objective function
INITCONFCN init_confcn;		///< Initializes the constraints
AppLauncher * launcher_;
bool init_flag;			///< Has the function been initialized?

public:
// Constructor
FDNLF1APP(): 
  FDNLF1(){;}
FDNLF1APP(int ndim): 
  FDNLF1(ndim){;}
FDNLF1APP(int ndim, USERFCN0APP f, INITFCNAPP i, AppLauncher * launcher, CompoundConstraint* constraint = 0): 
  FDNLF1(ndim, NULL, NULL, constraint)
  { fcn = f; init_fcn = i; init_flag = false; analytic_grad = 0; launcher_ = launcher;}
FDNLF1APP(int ndim, int nlncons, USERNLNCON0APP f, INITFCNAPP i): 
  FDNLF1(ndim, nlncons, NULL, NULL)
  { confcn = f; init_fcn = i; init_flag = false; analytic_grad = 0;}
FDNLF1APP(int ndim, USERFCN0APP f, INITFCNAPP i, INITCONFCN c): 
  FDNLF1(ndim){ fcn = f; init_fcn = i; init_confcn = c; 
  init_flag = false; analytic_grad = 0; constraint_ = init_confcn(ndim);}
FDNLF1APP(int ndim, int nlncons, USERNLNCON0APP f, INITFCNAPP i, AppLauncher *launcher):
  FDNLF1(ndim, nlncons, NULL, NULL) {confcn = f; init_fcn = i; init_flag = false; 
  analytic_grad = 0; launcher_ = launcher;}
			
// Destructor
virtual ~FDNLF1APP() {;}                     

virtual void reset(); 				///< Reset the parameter values 
virtual void initFcn();                       	///< Initialize function
virtual void eval();                       	///< Evaluate everything 
virtual real evalF();                         	///< Evaluate f 
virtual real evalF(const NEWMAT::ColumnVector& x);    	///< Evaluate f at x 
virtual NEWMAT::ColumnVector evalG();              	///< Evaluate gradient 
/// Evaluate gradient at x 
virtual NEWMAT::ColumnVector evalG(const NEWMAT::ColumnVector& x);  	
virtual NEWMAT::SymmetricMatrix evalH();              		///< Evaluate hessian 
NEWMAT::SymmetricMatrix FDHessian(NEWMAT::ColumnVector& x);     ///< Evaluate Hessian

// Print state
virtual void printState(char *);    
virtual void fPrintState(ostream *, char *);    

/// Evaluate the Lagrangian, its gradient and Hessian
virtual real evalLagrangian(const NEWMAT::ColumnVector& x, NEWMAT::ColumnVector& mult,
  const NEWMAT::ColumnVector& type) ;
virtual NEWMAT::ColumnVector evalLagrangianGradient(const NEWMAT::ColumnVector& x,
						const NEWMAT::ColumnVector& mult,
						const NEWMAT::ColumnVector& type) ;
/// Evaluate constraint at x
virtual NEWMAT::ColumnVector evalCF(const NEWMAT::ColumnVector& x);  	
/// Evaluate constraint gradient at x
virtual NEWMAT::Matrix evalCG(const NEWMAT::ColumnVector& x);  	

private:
/// Evaluate hessian at x
virtual NEWMAT::SymmetricMatrix evalH(NEWMAT::ColumnVector& x);
/// Evaluate constraint hessian at x
virtual NEWMAT::SymmetricMatrix evalCH(NEWMAT::ColumnVector &x);
/// Evaluate constraint hessian at x
virtual OptppArray<NEWMAT::SymmetricMatrix> evalCH(NEWMAT::ColumnVector &x, int darg);
virtual void evalC(const NEWMAT::ColumnVector& x); 	
};

} // namespace OPTPP

#endif
