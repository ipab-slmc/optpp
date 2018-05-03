
#include <iostream>
#include <fstream>

#include "NLF.h"
#include "NLP.h"

#include "tstfcn.h"

using namespace OPTPP;

int main ()
{
  int n = 2;

  NLP trig_prob = new  NLF0(n, trig, init_trig);
  trig_prob.initFcn();
  trig_prob.eval();
  trig_prob.printState("Trigonometric Problem");

  NLP rosen_prob = new NLF1(n, rosen, init_rosen);
  rosen_prob.initFcn();
  rosen_prob.eval();
  rosen_prob.printState("RosenBrock Problem");

  NLP rosen2_prob = new NLF2(n, rosen2, init_rosen);
  rosen2_prob.initFcn();
  rosen2_prob.eval();
  rosen2_prob.printState("RosenBrock Problem with analytic Hessian ");
}
