#ifndef Appl_Data_h
#define Appl_Data_h

/*----------------------------------------------------------------------
  Copyright (c) 2001, Sandia Corporation.   Under the terms of Contract 
  DE-AC04-94AL85000, there is a non-exclusive license for use of this 
  work by or on behalf of the U.S. Government.
 ----------------------------------------------------------------------*/

#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#ifdef HAVE_STD
#include <cstring>
#else
#include <string.h>
#endif

#include "globals.h"
#include "OptppArray.h"

/**
 * @author J. C. Meza, Sandia National Laboratories, meza@ca.sandia.gov
 * @note Modified by P.J. Williams 02/2006
 */

namespace OPTPP {

class Appl_Data {
private :
  /// Dimension of the problem
  int             	dimension;		
  /// Current point
  NEWMAT::ColumnVector    *xparm;		
  /// Objective function value 
  double          function_value;		
  /// Gradient of the objective function 
  NEWMAT::ColumnVector    *gradient;		
  /// Hessian of the objective function 
  NEWMAT::SymmetricMatrix *Hessian;		
  /// Constraint value 
  NEWMAT::ColumnVector    *constraint_value;	
  /// Gradient of the constraints 
  NEWMAT::Matrix          *constraint_gradient;	
  /// Hessian of the constraints 
  OptppArray<NEWMAT::SymmetricMatrix> *constraint_Hessian; 
  /// Residuals of the least square objective function 
  NEWMAT::ColumnVector    *lsq_residuals;	
  /// Jacobian of the least square objective function 
  NEWMAT::Matrix          *lsq_jacobian;        
  /// Is the function value current? 
  bool            function_current;		
  /// Is the gradient current? 
  bool            gradient_current;		
  /// Is the Hessian current? 
  bool            Hessian_current;		

public:
  /**
   * Default Constructor
   */
  Appl_Data();

  /**
   * Destructor
   */
  ~Appl_Data();

  void reset();

  bool Compare(const NEWMAT::ColumnVector&);
  bool getF(const NEWMAT::ColumnVector&, real&);
  bool getGrad(const NEWMAT::ColumnVector&, NEWMAT::ColumnVector&);
  bool getHess(const NEWMAT::ColumnVector&, NEWMAT::SymmetricMatrix&);

  bool getCF(const NEWMAT::ColumnVector&, NEWMAT::ColumnVector&);
  bool getCGrad(const NEWMAT::ColumnVector&, NEWMAT::Matrix&);
  bool getCHess(const NEWMAT::ColumnVector&, OptppArray<NEWMAT::SymmetricMatrix>&);

  bool getLSQF(const NEWMAT::ColumnVector&, NEWMAT::ColumnVector&);
  bool getLSQJac(const NEWMAT::ColumnVector&, NEWMAT::Matrix&);
  

  /// Update the objective function value
  void update(int,int,const NEWMAT::ColumnVector&,real);
  /// Update the objective function and gradient
  void update(int,int,const NEWMAT::ColumnVector&,real,NEWMAT::ColumnVector&);
  /// Update the objective function, gradient, and Hessian
  void update(int,int,const NEWMAT::ColumnVector&,real,
    NEWMAT::ColumnVector&,NEWMAT::SymmetricMatrix&);

  /// Update the nonlinear constraint functions
  void constraint_update(int,int,int,const NEWMAT::ColumnVector&,
    NEWMAT::ColumnVector&);
  /// Update the nonlinear constraint functions and Jacobian
  void constraint_update(int,int,int,const NEWMAT::ColumnVector&,
    NEWMAT::ColumnVector&,NEWMAT::Matrix&);
  /// Update the nonlinear constraint functions, Jacobian, and Hessians
  void constraint_update(int,int,int,const NEWMAT::ColumnVector&,
    NEWMAT::ColumnVector&,NEWMAT::Matrix&,OptppArray<NEWMAT::SymmetricMatrix>&);

  /// Update the least square residuals 
  void lsq_update(int,int,int,const NEWMAT::ColumnVector&,NEWMAT::ColumnVector&);

  /// Update the least square residuals and Jacobian
  void lsq_update(int,int,int,const NEWMAT::ColumnVector&,
    NEWMAT::ColumnVector&,NEWMAT::Matrix&);
};

} // namespace OPTPP

#endif
