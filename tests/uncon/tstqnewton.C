//
// Test program for Quasi-Newton optimization objects
//
// 1. Quasi Newton with trust regions on an NLF1
// 2. Quasi Newton with More-Thuente Line Search on an NLF1
// 3. Quasi Newton with Backtracking Line Search on an NLF1
//

#include <fstream>

#include "OptQNewton.h"
#include "NLF.h"
#include "tstfcn.h"

using NEWMAT::ColumnVector;
using namespace OPTPP;

#define true  1
#define false 0

void update_model(int, int, ColumnVector) {}

int main ()
{
  int n = 2;
  
  static char *status_file = {"tstqnewton.out"};

//----------------------------------------------------------------------------
// 1. Quasi-Newton with trust regions
//----------------------------------------------------------------------------

  //  Create a Nonlinear problem object

  NLF1 nlp(n,rosen,init_rosen);
  
  nlp.setIsExpensive(true);
  
  //  Build a Quasi-Newton object and optimize 

  OptQNewton objfcn(&nlp,update_model);   
  if (!objfcn.setOutputFile(status_file, 0))
    cerr << "main: output file open failed" << endl;
  objfcn.setTRSize(100.);
  objfcn.optimize();

#ifdef REG_TEST
  ColumnVector x_sol = nlp.getXc();
  double f_sol = nlp.getF();
  ostream* optout = objfcn.getOutputFile();
  if ((1.0 - x_sol(1) <= 1.e-2) && (1.0 - x_sol(2) <= 1.e-2) && (f_sol
								 <=
								 1.e-2))
    *optout << "Quasi-Newton 1 PASSED" << endl;
  else
    *optout << "Quasi-Newton 1 FAILED" << endl;
#endif

  objfcn.printStatus("Solution from quasi-newton: trust regions");
  objfcn.cleanup();	 

//----------------------------------------------------------------------------
// 2. Quasi-Newton with More and Thuente's line search
//----------------------------------------------------------------------------

  NLF1 nlp2(n,rosen,init_rosen);
  
  OptQNewton objfcn2(&nlp2,update_model);   
  objfcn2.setSearchStrategy(LineSearch);
  objfcn2.setOutputFile(status_file, 1);
  objfcn2.optimize();
  objfcn2.printStatus("Solution from quasi-newton: More-Thuente linesearch");

#ifdef REG_TEST
  x_sol = nlp2.getXc();
  f_sol = nlp2.getF();
  optout = objfcn2.getOutputFile();
  if ((1.0 - x_sol(1) <= 1.e-2) && (1.0 - x_sol(2) <= 1.e-2) && (f_sol
								 <=
								 1.e-2))
    *optout << "Quasi-Newton 2 PASSED" << endl;
  else
    *optout << "Quasi-Newton 2 FAILED" << endl;
#endif

  objfcn2.cleanup();	 
    
//----------------------------------------------------------------------------
// 3. Quasi-Newton with backtrack line search
//----------------------------------------------------------------------------

  NLF1 nlp3(n,rosen,init_rosen);
  nlp3.setIsExpensive(true);
  
  OptQNewton objfcn3(&nlp3, update_model);   
  objfcn3.setSearchStrategy(LineSearch);
  objfcn3.setOutputFile(status_file, 1);
  objfcn3.optimize();
  objfcn3.printStatus("Solution from quasi-newton: backtrack linesearch");

#ifdef REG_TEST
  x_sol = nlp3.getXc();
  f_sol = nlp3.getF();
  optout = objfcn3.getOutputFile();
  if ((1.0 - x_sol(1) <= 1.e-2) && (1.0 - x_sol(2) <= 1.e-2) && (f_sol
								 <=
								 1.e-2))
    *optout << "Quasi-Newton 3 PASSED" << endl;
  else
    *optout << "Quasi-Newton 3 FAILED" << endl;
#endif

  objfcn3.cleanup();	 
}
