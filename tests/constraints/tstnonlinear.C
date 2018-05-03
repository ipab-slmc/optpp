
#include <iostream>
#include <fstream>

#include "NLF.h"
#include "NLP.h"
#include "NonLinearEquation.h"
#include "NonLinearInequality.h"

#include "tstfcn.h"

using NEWMAT::ColumnVector;
using std::cout;
using namespace OPTPP;

int main ()
{
  int n = 2, ncnln = 1;
  ColumnVector x(n), b(ncnln), resid;

  //  Create a nonlinear problem object
  NLP trig_prob  = new NLF0(n,trig,init_trig);

  //  Initialize and evaluate the function at x
  trig_prob.initFcn();
  trig_prob.evalF();
  trig_prob.printState("Trigometric Function");

  // Create a nonlinear equation 
  b = 4;
  NLP* trig_constraint = new NLP( new NLF0(n,1,trig_as_a_constraint,init_trig));
  NonLinearEquation eqn(trig_constraint, b);

  x     = trig_prob.getXc();
  resid = eqn.evalResidual(x);

  cout << "========= Nonlinear Inequality  ================= \n";
  cout << "rhs" << "\t" << "residual" << "\n"; 
  cout << b(1)  << "\t" << resid(1)   << "\n";


}
