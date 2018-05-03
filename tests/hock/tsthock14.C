
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
  int n = 2;

  static char *status_file = {"tsthock14.out"};

  //  Create a Constrained Nonlinear problem object
  NLF1 nips(n,hs14,init_hs14,create_constraint_hs14);

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
  if ((8.2288e-01 - x_sol(1) <= 1.e-2) && (9.1144e-01 - x_sol(2) <= 1.e-2) && 
	  (1.3935 - f_sol <= 1.e-2))
    *optout << "Hock  14 PASSED" << endl;
  else
    *optout << "Hock  14 FAILED" << endl;
#endif
}

