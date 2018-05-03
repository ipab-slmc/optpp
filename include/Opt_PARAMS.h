// --
// parameter structures for GSS algorithms
//--

#ifndef OPT_PARAMS
#define OPT_PARAMS

namespace OPTPP {

class OptGSS_params {
 public:
  double Delta;
  double Delta_tol;
  double Phi;
  double Theta;
  int    Iter_max;
  bool   SearchAll;
  bool   printCOPYRIGHT;
  bool   printXiter;
  bool   printGiter;

  // constructor -- sets default parameters
  OptGSS_params() : Delta(0.0),  
    Delta_tol(1e-12), 
    Phi(1.3),  
    Theta(0.5), 
    Iter_max(10000), 
    SearchAll(true), 
    printCOPYRIGHT(false), 
    printXiter(false), 
    printGiter(false)
    {;}
};

} // namespace OPTPP

#endif
