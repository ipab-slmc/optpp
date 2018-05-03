
#include <iostream>
#include <fstream>

#include "NLF.h"
#include "OptNIPS.h"

#include "hockfcns.h"

using NEWMAT::ColumnVector;

void update_model(int, int, ColumnVector) {}

int main ()
{
  int n = 3;

  static char *status_file = {"tsthock65.out"};

  //  Create a Constrained Nonlinear problem object
  NLF2 nips(n,hs65_2,init_hs65,create_constraint_hs65_2);

  //  Build a NIPS object and optimize
  OptNIPS objfcn(&nips, update_model);
  objfcn.setOutputFile(status_file, 0);
  objfcn.setFcnTol(1.0e-06);
  objfcn.setMaxIter(150);
  objfcn.setSearchStrategy(LineSearch);
  objfcn.setMeritFcn(ArgaezTapia);
  objfcn.optimize();
  objfcn.printStatus("Solution from nips");
  objfcn.cleanup();

#ifdef REG_TEST
  ColumnVector x_sol = nips.getXc();
  double f_sol = nips.getF();
  ostream* optout = objfcn.getOutputFile();
  if ((3.6505 - x_sol(1) <= 1.e-2) && (3.6505 - x_sol(2) <= 1.e-2) && 
      (4.6204 - x_sol(3) <= 1.e-2) && (9.5353e-01 - f_sol <= 1.e-2))
    *optout << "Hock  65 PASSED" << endl;
  else
    *optout << "Hock  65 FAILED" << endl;
#endif
}

