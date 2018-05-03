
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
  int n = 3;

  static char *status_file = {"tsthock26.out"};

  //  Create a Constrained Nonlinear problem object
  NLF1 nips(n,hs26,init_hs26,create_constraint_hs26);

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
  if ((9.9885e-01 - x_sol(1) <= 1.e-2) && (9.9885e-01 - x_sol(2) <= 1.e-2) && 
	  (1.0011 - x_sol(3) <= 1.e-2) && (0 - f_sol <= 1.e-2))
    *optout << "Hock  26 PASSED" << endl;
  else
    *optout << "Hock  26 FAILED" << endl;
#endif
}

