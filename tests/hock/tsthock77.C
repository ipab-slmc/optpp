
#include <iostream>
#include <fstream>

#include "NLF.h"
#include "OptFDNIPS.h"

#include "hockfcns.h"

using NEWMAT::ColumnVector;

using namespace OPTPP;

void update_model(int, int, ColumnVector) {}

int main ()
{
  int n = 5;

  static char *status_file = {"tsthock77.out"};

  //  Create a Constrained Nonlinear problem object
  NLF1 nips(n,hs77,init_hs77,create_constraint_hs77);

  //  Build a finite-difference NIPS object and optimize
  OptFDNIPS objfcn(&nips, update_model);
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
  if ((1.1662 - x_sol(1) <= 1.e-2) && (1.1821 - x_sol(2) <= 1.e-2) && 
      (1.3803 - x_sol(3) <= 1.e-2) && (1.5060 - x_sol(4) <= 1.e-2) && 
      (6.1092e-01 - x_sol(5) <= 1.e-2) && (2.4151e-01 - f_sol <= 1.e-2))
    *optout << "Hock  77 PASSED" << endl;
  else
    *optout << "Hock  77 FAILED" << endl;
#endif
}

