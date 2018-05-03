#ifndef OPT_DIRECT
#define OPT_DIRECT

//------------------------------------------------------------------------
// Direct -- A Derived Class from Optimize 
//------------------------------------------------------------------------

#include "Opt.h"

namespace OPTPP{

class OptDirect: public OptimizeClass {

 public:
  OptDirect(){}
  OptDirect(int n): OptimizeClass(n){}
  OptDirect(int n, TOLS t): OptimizeClass(n,t){}
  virtual ~OptDirect(){}

  /* these are virtual in Optimizeclass

  virtual void acceptStep(int, int) = 0;
  virtual int checkConvg() {return 0;}
  virtual void optimize() {}
  virtual void readOptInput() {}
  virtual void printStatus(char *) {}
  virtual void updateModel(int, int, NEWMAT::ColumnVector) = 0;
  */
};

} // namespace OPTPP
#endif
