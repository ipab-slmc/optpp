/** 
 *  Hock-Schittkowski Test Problem 78.
 *
 *  Example of using NIPS algorithm on a NLF2.
 */

#include <iostream>
#include <fstream>

#include "NLF.h"
#include "OptNIPS.h"

#include "hockfcns.h"

using NEWMAT::ColumnVector;

using namespace OPTPP;

void update_model(int, int, ColumnVector) {}

int main ()
{
  int n = 5;

  static char *status_file = {"tsthock78.out"};

  //  Create a Constrained Nonlinear problem object
  NLF2 nips(n,hs78_2,init_hs78,create_constraint_hs78_2);

  //  Build a NIPS object and optimize
  OptNIPS objfcn(&nips, update_model);
  objfcn.setOutputFile(status_file, 0);

  //  Set function tolerance 
  objfcn.setFcnTol(1.0e-06);

  //  Set maximum allowable iterations for NIPS algorithm 
  objfcn.setMaxIter(150);

  //  Use a backtracking linesearch method to determine acceptable step 
  objfcn.setSearchStrategy(LineSearch);

  // Use the Argaez-Tapia merit function as a globalization strategy
  objfcn.setMeritFcn(ArgaezTapia);

  objfcn.optimize();
  objfcn.printStatus("Solution from nips");
  objfcn.cleanup();

#ifdef REG_TEST
  ColumnVector x_sol = nips.getXc();
  double f_sol = nips.getF();
  ostream* optout = objfcn.getOutputFile();
  if ((-1.0 - x_sol(1) <= 1.e-2) && (1.2683e-06 - x_sol(2) <= 1.e-2) && 
      (3.0 - x_sol(3) <= 1.e-2) && (-1.1096e-03 - x_sol(4) <= 1.e-2) && 
      (-1.1096e-03 - x_sol(5) <= 1.e-2) && (0 - f_sol <= 1.e-2))
    *optout << "Hock  78 PASSED" << endl;
  else
    *optout << "Hock  78 FAILED" << endl;
#endif
}

