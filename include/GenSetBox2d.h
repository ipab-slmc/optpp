#ifndef GenSetBox2d_h
#define GenSetBox2d_h
//------------------------------------------------------------------------
// Generating Set Classes - for use with OptGSS
//------------------------------------------------------------------------

/*------------------------------------------------------------------------
 Copyright (c) 2003,
 Ricardo Oliva (raoliva@lbl.gov)
 Lawrence Berkeley National Laboratory
 ------------------------------------------------------------------------*/

#include "GenSetBase.h"

using std::string;

namespace OPTPP {

//------------------------------------------------------------------------
// Standard Basis with 2*n vectors
//------------------------------------------------------------------------
class GenSetBox2d : public GenSetBase {  // Standard with 2*n vectors

 public:
  virtual string classnm() { return "GenSetBox2d";}

  // default constructor;
  GenSetBox2d(){};

  // Vsize specific constructor
  GenSetBox2d(int n) : GenSetBase(n) { 
    setSize(2*n+4);
    initActive();
  };

  void generate(int i, double a, NEWMAT::ColumnVector &x, NEWMAT::ColumnVector &y);
  ///< Stores the search direction in the vector y

  // overloaded pruning methods - virtual in base
  int init(){ return 0;}    ///< Computes initial generating set D
  int init(NEWMAT::ColumnVector& pV);  
  int update(){ return 0;}    ///< Updates D on each iteration
  int update(NEWMAT::ColumnVector& pV);   
  bool prunes() { return true; }

 
  /// Destructor
  virtual ~GenSetBox2d() {;}  
};

} // namespace OPTPP
#endif
