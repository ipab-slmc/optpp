//------------------------------------------------------------------------
// Copyright (C) 1993, 1994: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//------------------------------------------------------------------------

#include <string>
#include <iostream>
#include <string>

#include "NLF.h"

#include "tstfcn.h"

using std::cout;
using std::string;

using namespace OPTPP;

void SetupTestProblem(string test_id, USERFCN0 *test_problem, 
		      INITFCN *init_problem)
{
  //
  // Given a character string test_id, set the function pointers
  // test_problem and init_problem 
  // 

  string tstfcn[21];
  USERFCN0 test_problems[21];
  INITFCN  init_problems[21];

  bool found = false;

  tstfcn       [0] = "erosen";
  test_problems[0] = erosen;
  init_problems[0] = init_erosen;

  tstfcn       [1] = "epowell";
  test_problems[1] = epowell;
  init_problems[1] = init_epowell;

  tstfcn       [2] = "trig";
  test_problems[2] = trig;
  init_problems[2] = init_trig;

  tstfcn       [3] = "penalty1";
  test_problems[3] = penalty1;
  init_problems[3] = init_penalty1;

  tstfcn       [4] = "penalty2";
  test_problems[4] = penalty2;
  init_problems[4] = init_penalty2;

  tstfcn       [5] = "vardim";
  test_problems[5] = vardim;
  init_problems[5] = init_vardim;

  tstfcn       [6] = "gen_brown";
  test_problems[6] = gen_brown;
  init_problems[6] = init_gen_brown;

  tstfcn       [7] = "broyden_tridiag";
  test_problems[7] = broyden_tridiag;
  init_problems[7] = init_broyden;

  tstfcn       [8] = "chain_singular";
  test_problems[8] = chain_singular;
  init_problems[8] = init_chain_singular;

  tstfcn       [9] = "gen_wood";
  test_problems[9] = gen_wood;
  init_problems[9] = init_gen_wood;

  tstfcn       [10] = "chain_wood";
  test_problems[10] = chain_wood;
  init_problems[10] = init_chain_wood;

  tstfcn       [11] = "broyden1a";
  test_problems[11] = broyden1a;
  init_problems[11] = init_broyden;

  tstfcn       [12] = "broyden1b";
  test_problems[12] = broyden1b;
  init_problems[12] = init_broyden;

  tstfcn       [13] = "broyden2a";
  test_problems[13] = broyden2a;
  init_problems[13] = init_broyden;

  tstfcn       [14] = "broyden2b";
  test_problems[14] = broyden2b;
  init_problems[14] = init_broyden;

  tstfcn       [15] = "tointbroy";
  test_problems[15] = tointbroy;
  init_problems[15] = init_broyden;

  tstfcn       [16] = "cragg_levy";
  test_problems[16] = gen_cragg_levy;
  init_problems[16] = init_gen_cragg_levy;

  tstfcn       [17] = "toint_trig";
  test_problems[17] = toint_trig;
  init_problems[17] = init_toint_trig;

  tstfcn       [18] = "chebyquad";
  test_problems[18] = chebyquad;
  init_problems[18] = init_chebyquad;

  tstfcn       [19] = "nelder";
  test_problems[19] = nelder;
  init_problems[19] = init_nelder;

  int maxprob = 20;

  //
  // Loop over possible test problems
  //
  int i = 0;
  while (!found & i < maxprob) {
    if (test_id == tstfcn[i]) {
      test_id      = tstfcn[i];
      *test_problem = test_problems[i];
      *init_problem = init_problems[i];
      found = true;
    }
    i++;
  }
  if (!found) {
    cout << "SetupTestProblem: No problem specified. Using rosen as default\n";
    *test_problem = rosen0;
    *init_problem = init_rosen;
  }
}

