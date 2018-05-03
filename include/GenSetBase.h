#ifndef GenSetBase_h
#define GenSetBase_h
//------------------------------------------------------------------------
// Generating Set Class - for use with OptGSS
//------------------------------------------------------------------------

/*------------------------------------------------------------------------
 Copyright (c) 2003,
 Ricardo Oliva (raoliva@lbl.gov)
 Lawrence Berkeley National Laboratory
 ------------------------------------------------------------------------*/

#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#define WANT_MATH

#include <string>
#include <iostream>
#ifdef HAVE_STD
#include <cfloat>
#include <cstdlib>
#else
#include <float.h>
#include <stdlib.h>
#endif

#include "newmatap.h"

using std::cerr;
using std::string;

namespace OPTPP {

class GenSetBase {  // Generating Set Base Class

 protected:
  int  Vdim;

  int  Size;
  int  nAct;

  NEWMAT::ColumnVector ActiveIDs;
  NEWMAT::ColumnVector InactiveIDs;

 public:
  virtual string classnm() { return "GenSetBase";}; 

  // default constructor;
  GenSetBase() : Vdim(0),  Size(0), nAct(0) {};

  // Constructor with specific vector-size
  GenSetBase(int n) : Vdim(n), Size(0), nAct(0) {};

  /// Destructor
  virtual ~GenSetBase() {;}  

  // Basic init method -- call after default constr.
  void init(int vd) { Vdim = vd; }

  // Basic set/get methods
  void setSize(int s) { Size = s; }
  void setVdim(int n) { Vdim = n; }
  int size() { return Size; }
  int vdim() { return Vdim; }

  //--
  // Generating Methods
  //--

  // -- wrt ALL Directions --

  //  virtual NEWMAT::ColumnVector generate(int i);
  ///< Returns  d_i = ith element of D

  //  virtual void generate(int i, NEWMAT::ColumnVector& y);
  ///< Stores d_i in y

  //  virtual NEWMAT::ColumnVector generate(int i, double a, NEWMAT::ColumnVector& x);
  ///< Returns the vector y_i =  x + a*d_i
  
  virtual void generate(int i, double a, NEWMAT::ColumnVector &x, 
    NEWMAT::ColumnVector &y) = 0;  
  ///< Stores in y the vector  x + a*d_i

  // -- wrt ACTIVE Directions --
  /*
  virtual 
    NEWMAT::ColumnVector generateActive(int i)
    { return generate(activeID(i)); }

  virtual 
    void generateActive(int i, NEWMAT::ColumnVector& y) 
    { generate(activeID(i), y); }

  virtual 
    NEWMAT::ColumnVector generateActive(int i, double s, NEWMAT::ColumnVector& x)
    { return generate(activeID(i), s, x); }
  */
  virtual void generateActive(int i, double s, NEWMAT::ColumnVector &x, 
    NEWMAT::ColumnVector &y)
    { generate(activeID(i), s, x, y); }

  // -- wrt INACTIVE Directions --
  /*
  virtual 
    NEWMAT::ColumnVector generateInactive(int i)
    ///< Returns b_i, the ith element of the INACTIVE subset of D
    { 
      NEWMAT::ColumnVector v; 
      v = generate(inactiveID(i)); 
      return v; 
    }

  virtual 
    void generateInactive(int i, NEWMAT::ColumnVector& y)
    ///< Stores b_i in y 
    { generate(inactiveID(i), y); }

  virtual 
    NEWMAT::ColumnVector generateInactive(int i, double s, NEWMAT::ColumnVector& x)
    ///< Returns the vector x + s*b_i,
    { return generate(inactiveID(i), s, x); }
  */
  virtual 
    void generateInactive(int i, double s, NEWMAT::ColumnVector &x, NEWMAT::ColumnVector &y)
    ///< Stores in y the vector x + s*b_i,
    { generate(inactiveID(i), s, x, y); }

  //--
  // Pruning methods
  //--
  virtual void initActive() {  
    // call this in constructor of derived class 
    // after size of derived class has been set
    if (Size==0) { 
      //cerr << "!!! ERROR: GenSetBase::initActive() called when size==0\n";
      return;
    }
    nAct = Size;
    ActiveIDs.ReSize(Size);
    for (int i=1; i<=Size; i++) ActiveIDs(i) = i; 
    InactiveIDs.ReSize(Size); 
    InactiveIDs = 0;
  }

  virtual int nActive() { return nAct; }
  virtual int nInactive() { return (Size - nAct); }
  virtual int activeID(int j) { return static_cast<int>(ActiveIDs(j)); }
  virtual int inactiveID(int j) { return static_cast<int>(InactiveIDs(j)); }

  virtual int init(){ return 0;}    ///< Computes initial generating set D
  virtual int init(NEWMAT::ColumnVector& pV){ return 0;}    

  virtual int update(){ return 0;}    ///< Updates D on each iteration
  virtual int update(NEWMAT::ColumnVector& pV){return 0;}        

  virtual bool prunes(){return false;} 
  ///< switch to true if implementing pruning in derived class

  bool   generateAll(NEWMAT::Matrix& M, NEWMAT::ColumnVector& X, double Delta=1.0);
  NEWMAT::Matrix generateAll(NEWMAT::ColumnVector& X, double D=1.0) {
    NEWMAT::Matrix M(Vdim,Size);
    generateAll(M,X,D);
    return M;
  }
  NEWMAT::Matrix generateAll(double Delta=1.0) {
    NEWMAT::ColumnVector X(Vdim);
    X = 0;
    return generateAll(X,Delta); 
  }

  bool   generateAllActive(NEWMAT::Matrix& M, NEWMAT::ColumnVector& X, double Delta=1.0);
  NEWMAT::Matrix generateAllActive(NEWMAT::ColumnVector& X, double D=1.0) {
    int n = nActive();
    int m = Vdim;
    NEWMAT::Matrix M(m,n);
    generateAllActive(M,X,D);
    return M;
  }
  NEWMAT::Matrix generateAllActive(double Delta=1.0) {
    NEWMAT::ColumnVector X(Vdim);
    X = 0;
    return generateAllActive(X,Delta); 
  }

  NEWMAT::Matrix pllMesh(int P, NEWMAT::ColumnVector& xc, NEWMAT::ColumnVector& xn, double d=0.0);

}; // end of GenSetBase class

} // namespace OPTPP

#endif
