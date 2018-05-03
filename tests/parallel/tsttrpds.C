/** 
 * Test program for "Bound Constrained" Quasi-Newton optimization objects.
 *
 * Example of using TRPDS algorithm on a NLF1 problem.  
 */

#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#include <string>
#include <iostream>
#include <fstream>
#ifdef HAVE_STD
#include <cstdio>
#else
#include <stdio.h>
#endif

#ifdef WITH_MPI
#include "mpi.h"
#endif

#include "NLF.h"
#include "NLP2.h"
#include "BoundConstraint.h"
#include "CompoundConstraint.h"
#include "OptConstrQNewton.h"
#include "ioformat.h"

#include "tstfcn.h"

using NEWMAT::ColumnVector;
using std::cerr;

using namespace OPTPP;

void SetupTestProblem(string test_id, USERFCN0 *test_problem, 
		      INITFCN *init_problem);
void update_model(int, int, ColumnVector) {}


int main (int argc, char* argv[])
{

  int ndim = 2;
  double time0, opt_time;

  // USERFCN0 test_problem;
  // INITFCN  init_problem;

  // string test_id;

  //
  // Setup the test problem
  // test_problem is a pointer to the function (fcn) to optimize
  // init_problem is a pointer to the function that initializes fcn
  // test_id is a character string identifying the test problem
  // ndim is the dimension of the problem

#ifdef WITH_MPI
  int me;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
#endif

  //if (argc != 3) {
  //  cout << "Usage: tstnewtpds problem_name ndim\n";
  //  exit(1);
  //}

  //test_id = argv[1];
  //ndim    = atoi(argv[2]);

  //SetupTestProblem(test_id, &test_problem, &init_problem);

  //
  // Generate file name based on the problem name
  //
  char status_file[80];
  //strcpy(status_file,test_id.c_str());
  strcpy(status_file, "tsttrpds");
#ifdef WITH_MPI
  sprintf(status_file,"%s.out.%d", status_file, me);
#else
  strcat(status_file,".out");
#endif

  //  Build a CompoundConstraint object
  ColumnVector bound(ndim );
  bound = -500.0;
  Constraint bc = new BoundConstraint(ndim,bound);

  //  Create an OptppArray of Constraints 
  OptppArray<Constraint> arrayOfConstraints;
  arrayOfConstraints.append(bc);

  //  Create a compound constraint 
  CompoundConstraint constraints(arrayOfConstraints);  
    
  //  Create a constrained Nonlinear problem object 

  //FDNLF1 nlp(ndim,test_problem, init_problem, &constraints);
  NLF1 nlp(ndim, erosen1, init_erosen);
  //  Initialize and evaluate the function at x

  nlp.setSpecOption(NoSpec);
  nlp.setIsExpensive(true);
  
  //  Build a Quasi-Newton object and optimize 

//----------------------------------------------------------------------------
// 1. Quasi-Newton with trust region PDS
//----------------------------------------------------------------------------


  OptConstrQNewton objfcn(&nlp,update_model);   
  if (!objfcn.setOutputFile(status_file, 0))
//    cerr << "main: output file open failed" << endl;
  ostream* optout = objfcn.getOutputFile();
  //*optout << "Test problem: " << test_id << endl;
  *optout << "Test problem: " << "erosen1" << endl;
  *optout << "Dimension   : " << ndim    << endl;
  objfcn.setSearchStrategy(TrustPDS);
  objfcn.setMaxFeval(10000);
  //  objfcn.setTRSize(100.);

  time0 = get_wall_clock_time();
  objfcn.optimize();
  opt_time = get_wall_clock_time() - time0;
  *optout << "wall clock time =" << e(opt_time,12,4) << endl;

  objfcn.printStatus("Solution from quasi-newton");
  objfcn.cleanup();

#ifdef REG_TEST
  ColumnVector x_sol = nlp.getXc();
  double f_sol = nlp.getF();
  if ((1.0 - x_sol(1) <= 1.e-2) && (1.0 - x_sol(2) <= 1.e-2) && (f_sol
								 <=
								 1.e-2))
    *optout << "TRPDS 1 PASSED" << endl;
  else
    *optout << "TRPDS 1 FAILED" << endl;
#endif

#ifdef WITH_MPI
  MPI_Finalize();
#endif

}
