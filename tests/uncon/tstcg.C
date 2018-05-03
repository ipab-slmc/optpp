/** \example tstcg.C
 * Test program for CG optimization objects
 *
 * 1. Nonlinear CG with More-Thuente Line Search on an NLF1
 *
 * 2. Nonlinear CG with backtracking Line Search on an NLF1
 *
 */


#include <fstream>

#include "OptCG.h"
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
  
  static char *status_file = {"tstcg.out"};

  //  Create a Nonlinear problem object

  NLF1 nlp(n,rosen,init_rosen);
  
  //  Build a CG object and optimize 

  OptCG objfcn(&nlp);   
  objfcn.setUpdateModel(update_model);
  if (!objfcn.setOutputFile(status_file, 0))
    cerr << "main: output file open failed" << endl;
  objfcn.setGradTol(1.e-6);
  objfcn.setMaxBacktrackIter(10);
  objfcn.optimize();
    
  objfcn.printStatus("Solution from CG: Fcn not Expensive");

#ifdef REG_TEST
  ColumnVector x_sol = nlp.getXc();
  double f_sol = nlp.getF();
  ostream* optout = objfcn.getOutputFile();
  if ((1.0 - x_sol(1) <= 1.e-2) && (1.0 - x_sol(2) <= 1.e-2) && (f_sol
								 <=
								 1.e-2))
    *optout << "CG 1 PASSED" << endl;
  else
    *optout << "CG 1 FAILED" << endl;
#endif

  objfcn.cleanup();

  //  Build a CG object and optimize 

  NLF1 nlp2(n,rosen,init_rosen);
  nlp2.setIsExpensive(true);

  OptCG objfcn2(&nlp2);   
  objfcn2.setUpdateModel(update_model);
  objfcn2.setOutputFile(status_file, 1);
  objfcn2.setGradTol(1.e-6);
  objfcn2.setLineSearchTol(1.e-4);
  objfcn2.setMaxBacktrackIter(10);
  objfcn2.optimize();
    
  objfcn2.printStatus("Solution from CG: Fcn Expensive");

#ifdef REG_TEST
  x_sol = nlp2.getXc();
  f_sol = nlp2.getF();
  optout = objfcn2.getOutputFile();
  if ((1.0 - x_sol(1) <= 1.e-2) && (1.0 - x_sol(2) <= 1.e-2) && (f_sol
								 <=
								 1.e-2))
    *optout << "CG 2 PASSED" << endl;
  else
    *optout << "CG 2 FAILED" << endl;
#endif

  objfcn2.cleanup();
}
