/** 
 * Test 1: Newton's method with trust regions on an LSQNLF, no derivatives
 * Test 2: Newton's method with trust regions on an LSQNLF, analytic derivatives
 */

#include <fstream>

#include "LSQNLF.h"
#include "OptNewton.h"
#include "tstfcn.h"

using NEWMAT::ColumnVector;
using namespace OPTPP;

void update_model(int, int, ColumnVector) {}

int main ()
{
  int n = 2;
  
  static char *status_file = {"tstgnewton.out"};
//----------------------------------------------------------------------------
// Test Problem 1  - Finite Difference Gradients
//----------------------------------------------------------------------------

  //  Create a Nonlinear problem object

  LSQNLF nlp(n,n,rosen0_least_squares,init_rosen);
  
  //  Build a Newton object and optimize 

  //nlp.setIsExpensive(true);
  OptNewton objfcn(&nlp, update_model);
  objfcn.setOutputFile(status_file, 0);
  objfcn.setTRSize(1.0e3);
  nlp.setIsExpensive(false);
  //nlp.setDebug();
  objfcn.optimize();
  objfcn.printStatus("Solution from newton");

#ifdef REG_TEST
  ColumnVector x_sol = nlp.getXc();
  double f_sol = nlp.getF();
  ostream* optout = objfcn.getOutputFile();
  if ((1.0 - x_sol(1) <= 1.e-2) && (1.0 - x_sol(2) <= 1.e-2) && (f_sol
								 <=
								 1.e-2))
    *optout << "Gauss-Newton 1 PASSED" << endl;
  else
    *optout << "Gauss-Newton 1 FAILED" << endl;
#endif

  objfcn.cleanup();	 


//----------------------------------------------------------------------------
// Test Problem 2  - Analytic Gradients
//----------------------------------------------------------------------------

  //  Create a Nonlinear problem object

  LSQNLF nlp2(n,n,rosen_least_squares,init_rosen);
  
  //  Build a Newton object and optimize 

  OptNewton objfcn2(&nlp2, update_model);
  objfcn2.setOutputFile(status_file, 1);
  objfcn2.setTRSize(1.0e3);
  nlp2.setIsExpensive(false);
  //nlp2.setDebug();
  objfcn2.optimize();
  objfcn2.printStatus("Solution from newton");

#ifdef REG_TEST
  x_sol  = nlp2.getXc();
  f_sol  = nlp2.getF();
  optout = objfcn2.getOutputFile();
  if ((1.0 - x_sol(1) <= 1.e-2) && (1.0 - x_sol(2) <= 1.e-2) && (f_sol
								 <=
								 1.e-2))
    *optout << "Gauss-Newton 2 PASSED" << endl;
  else
    *optout << "Gauss-Newton 2 FAILED" << endl;
#endif

  objfcn2.cleanup();	 

}
