/** 
 *  Hock-Schittkowski Test Problem 1.
 *   
 *  Example of using finite-difference NIPS algorithm on a NLF1.
 */

#include <iostream>
#include <fstream>

using namespace std;

#include "NLF.h"
#include "OptFDNIPS.h"

#include "hockfcns.h"

using NEWMAT::ColumnVector;
using namespace OPTPP;

void update_model(int, int, ColumnVector) {}

int main ()
{
  int n = 2;

  static char *status_file = {"tsthock1.out"};

  //  Create a Constrained Nonlinear problem object
  NLF1 nips(n,hs1,init_hs1,create_constraint_hs1);

  //  Build a finite-difference NIPS object and optimize
  OptFDNIPS objfcn(&nips, update_model);

  objfcn.setOutputFile(status_file, 0);
  // Set convergence tolerance 
  objfcn.setFcnTol(1.0e-06);

  //  Set maximum allowable iterations for FDNIPS algorithm 
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
  if ((1.0 - x_sol(1) <= 1.e-2) && (1.0 - x_sol(2) <= 1.e-2) && (0.0 -f_sol <=
								 1.e-2))
    *optout << "Hock  1 PASSED" << endl;
  else
    *optout << "Hock  1 FAILED" << endl;
#endif
}

