#ifndef Appl_Data_npsol_h
#define Appl_Data_npsol_h

/*----------------------------------------------------------------------
  Copyright (c) 2001, Sandia Corporation.   Under the terms of Contract 
  DE-AC04-94AL85000, there is a non-exclusive license for use of this 
  work by or on behalf of the U.S. Government.
 --------------------------------------------------------------------- */

#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#ifdef HAVE_STD
#include <cstring>
#else
#include <string.h>
#endif

#include "globals.h"

#include "newmatap.h"

namespace OPTPP {

/** 
 * @note Modified by P.J. Williams on 02/2006 
 */

class Appl_Data_NPSOL {
private :
  int     	buffer_len;	///< buffer length of data 
  int     	buffer_pointer; ///< pointer to the buffer
  int     	dimension;	///< dimension of control vector 
  int     	ncnln;		///< number of nonlinear constraints
  double  	fvalue_buffer;	///< buffer for the function value
  NEWMAT::ColumnVector  *x_buffer;      ///< buffer for the control vector
  NEWMAT::ColumnVector  *grad_buffer;	///< buffer for the gradient
  NEWMAT::ColumnVector  *constr_buffer;	///< buffer for the constraints
  NEWMAT::Matrix  	*cjac_buffer;	///< buffer for the constraint jacobian
  bool    	fvalue_status;	///< status of the function
  bool    	grad_status;	///< status of the gradient of the obj. function
  bool    	constr_status;	///< status of the constraints
  bool    	cjac_status;	///< status of the constraint jacobian

public:
  Appl_Data_NPSOL();
  Appl_Data_NPSOL(int length);
  ~Appl_Data_NPSOL();
  bool Compare(NEWMAT::ColumnVector&);
  bool getF(NEWMAT::ColumnVector&, real&);
  bool getGrad(NEWMAT::ColumnVector&, NEWMAT::ColumnVector&);
  bool getConstraint(NEWMAT::ColumnVector&, NEWMAT::ColumnVector&);
  bool getCJacobian(NEWMAT::ColumnVector&, NEWMAT::Matrix&);

  /// Update the objective function value
  void update(int, int, NEWMAT::ColumnVector&, real); 
  /// Update the objective function gradient
  void update(int, int, NEWMAT::ColumnVector&, NEWMAT::ColumnVector&); 
  /// Update the nonlinear constraint functions
  void update(int, NEWMAT::ColumnVector&, int, NEWMAT::ColumnVector&); 
  /// Update the nonlinear constraint functions and constraint Jacobian
  void update(int, int, NEWMAT::ColumnVector&, int, NEWMAT::ColumnVector&, 
    NEWMAT::Matrix&); 
  /// Update the objective function value and nonlinear constraint functions
  void update(int, int, NEWMAT::ColumnVector&, double, int, 
    NEWMAT::ColumnVector&); 
  /// Update the objective function value, constraint functions, and Jacobian
  void update(int, int, NEWMAT::ColumnVector&, double, int, 
    NEWMAT::ColumnVector&, NEWMAT::Matrix&); 
};

} // namespace OPTPP

#endif

