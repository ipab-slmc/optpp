/** \example tstLBFGS.C
 *
 * Test program for LBFGS optimization objects
 * 
 * 1. Limited Memory BFGS method with More-Thuente Search on an NLF1
 *
 */

#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#include <fstream>
#ifdef HAVE_STD
#include <cstdio>
#else
#include <stdio.h>
#endif

#include "OptLBFGS.h"
#include "NLF.h"

#include "tstfcn.h"

using NEWMAT::ColumnVector;

using namespace OPTPP;
void update_model(int, int, ColumnVector) {}

int main ()
{
  int n = 2;
  
  static char *status_file = {"tstLBFGS.out"};

  //  Create a Nonlinear problem object

  NLF1 nlp(n,rosen,init_rosen);
  
  //  Build a LBFGS object and optimize 

  OptLBFGS objfcn(&nlp);   
  objfcn.setUpdateModel(update_model);
  if (!objfcn.setOutputFile(status_file, 0))
    cerr << "main: output file open failed" << endl;
  objfcn.setGradTol(1.e-6);
  objfcn.setMaxBacktrackIter(10);
  objfcn.setPrintFinalX(true);
  objfcn.optimize();
    
  objfcn.printStatus("Solution from LBFGS: More and Thuente's linesearch");

#ifdef REG_TEST
  ColumnVector x_sol = nlp.getXc();
  double f_sol = nlp.getF();
  ostream* optout = objfcn.getOutputFile();
  if ((1.0 - x_sol(1) <= 1.e-2) && (1.0 - x_sol(2) <= 1.e-2) && (f_sol
								 <=
								 1.e-2))
    *optout << "LBFGS 1 PASSED" << endl;
  else
    *optout << "LBFGS 1 FAILED" << endl;
#endif

  objfcn.cleanup();

}
