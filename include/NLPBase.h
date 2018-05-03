#ifndef NLPBase_h
#define NLPBase_h

/*----------------------------------------------------------------------
 Copyright (c) 2001, Sandia Corporation.   Under the terms of Contract 
 DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 work by or on behalf of the U.S. Government.

 P.J. Williams, Sandia National Laboratories, pwillia@sandia.gov
 ----------------------------------------------------------------------*/

#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#include <iostream>
#include <fstream>
#ifdef HAVE_STD
#include <cstring>
#include <cfloat>
#include <cmath>
#else
#include <string.h>
#include <float.h>
#include <math.h>
#endif

using std::ostream;

#include "globals.h"
#include "precisio.h"
#include "OptppArray.h"
#include "OptppFatalError.h"
#include "OptppExceptions.h"


namespace OPTPP {

/**
 * NLPBase is the Base Class for NonLinear Programming Problem
 *
 * @author P.J. Williams, Sandia National Laboratories, pwillia@sandia.gov
 * @date   Last modified 03/2007
 */

class NLPBase{

public:
// Destructor
  virtual ~NLPBase() {;}

// set and Accessor Methods
  virtual void setX(const int i, const real& x)  = 0;
  virtual void setX(const NEWMAT::ColumnVector& x)       = 0;

  virtual void setF(const real& fx)              = 0;

  virtual void setIsExpensive(const int e)       = 0;

  virtual void setFcnAccrcy(const int i, const real& accrcy)  = 0;
  virtual void setFcnAccrcy(const NEWMAT::ColumnVector& accrcy)       = 0;

  virtual int  getDim()         	const = 0;
  virtual int  getFevals()      	const = 0;
  virtual int  getIsExpensive() 	const = 0;
  virtual real getF()           	const = 0;
  virtual NEWMAT::ColumnVector getFcnAccrcy()   const = 0;
  virtual NEWMAT::ColumnVector getXc()  	const = 0;
  virtual real getFcnTime()     	const = 0;

// Constraint Accessor Methods
  virtual int  getNumOfCons()   	const = 0;
  virtual int  getNumOfNLCons()   	const = 0;
  virtual bool hasConstraints()	              = 0;
  virtual void printConstraints()             = 0; 

// Reset values to allow multiple instantiations 
  virtual void reset() = 0;

// Debugging tools
  virtual void setDebug()       = 0;
  virtual bool getDebug() const = 0;

// Function Evaluation Methods
  virtual void initFcn()  = 0;
  virtual real evalF()  = 0;
  virtual real evalF(const NEWMAT::ColumnVector &x)  = 0;
  virtual NEWMAT::ColumnVector evalG()  = 0;
  virtual NEWMAT::ColumnVector evalG(const NEWMAT::ColumnVector &x)  = 0;
  virtual NEWMAT::SymmetricMatrix evalH()  = 0;
  virtual NEWMAT::SymmetricMatrix evalH(NEWMAT::ColumnVector &x)  = 0;
  virtual void eval()  = 0;

// Constraint Evaluation Methods
  virtual NEWMAT::ColumnVector evalCF(const NEWMAT::ColumnVector &x)  = 0;
  virtual NEWMAT::Matrix evalCG(const NEWMAT::ColumnVector &x)  = 0;
  virtual NEWMAT::SymmetricMatrix evalCH(NEWMAT::ColumnVector &x)  = 0;
  virtual OptppArray<NEWMAT::SymmetricMatrix> evalCH(NEWMAT::ColumnVector &x, int darg)  = 0;
  virtual void evalC(const NEWMAT::ColumnVector &x)  = 0;

// Print Methods
  virtual void printState(char *s) = 0; 
  virtual void fPrintState(ostream *nlpout, char *s) = 0; 

};

} // namespace OPTPP
#endif
