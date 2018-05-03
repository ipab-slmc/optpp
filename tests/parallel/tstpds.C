/** 
 * Example of using a PDS optimize class on an NLF0 problem class.
 *
 * This class is used whenever the objective function does not
 * have analytic gradients.
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

#include "OptPDS.h"
#include "NLF.h"
#include "CompoundConstraint.h"
#include "BoundConstraint.h"
#include "OptppArray.h"
#include "cblas.h"
#include "ioformat.h"

#include "tstfcn.h"

using NEWMAT::ColumnVector;
using NEWMAT::Matrix;

using namespace OPTPP;

void SetupTestProblem(string test_id, USERFCN0 *test_problem, 
		      INITFCN *init_problem);
void update_model(int, int, ColumnVector) {}

int main (int argc, char* argv[])
{
  int i, j;
  int ndim = 2;
  double perturb;

  static char *schemefilename = {"myscheme"};

  // USERFCN0 test_problem;
  // INITFCN  init_problem;

  // string test_id;

#ifdef WITH_MPI
  int me;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
#endif

  //if (argc != 3) {
  //  cout << "Usage: tstpds problem_name ndim\n";
  //  exit(1);
  //}

  //test_id = argv[1];
  //ndim    = atoi(argv[2]);

  ColumnVector x(ndim);
  ColumnVector vscale(ndim);
  Matrix init_simplex(ndim,ndim+1);

  //SetupTestProblem(test_id, &test_problem, &init_problem);

  //
  // Generate file name based on the problem name
  //
  char status_file[80];
  //strcpy(status_file,test_id.c_str());
  strcpy(status_file, "tstpds");
#ifdef WITH_MPI
  sprintf(status_file,"%s.out.%d", status_file, me);
#else
  strcat(status_file,".out");
#endif

#ifdef TestBounds
  ColumnVector lbounds(ndim), ubounds(ndim);
  lbounds = -2.0;
  ubounds = 2.0;
  Constraint bc = new BoundConstraint(ndim, lbounds, ubounds);
#endif

  //  Create an OptppArray of Constraints 
  OptppArray<Constraint> arrayOfConstraints;

#ifdef TestBounds
  arrayOfConstraints.append(bc);
#endif

  //  Create a compound constraint 
  CompoundConstraint constraints(arrayOfConstraints);  
  
  //  Create a constrained Nonlinear problem object 
  //NLF0 nlp(ndim,test_problem, init_problem, &constraints);         
  NLF0 nlp(ndim,erosen, init_erosen, &constraints);         

  //  Build a PDS object and optimize 
  
  OptPDS objfcn(&nlp);
  objfcn.setOutputFile(status_file, 0);
  ostream* optout = objfcn.getOutputFile();
  //*optout << "Test problem: " << test_id << endl;
  *optout << "Test problem: " << "erosen" << endl;
  *optout << "Dimension   : " << ndim    << endl;
  objfcn.setFcnTol(1.49012e-8);
  objfcn.setMaxIter(500);
  objfcn.setMaxFeval(10000);
  objfcn.setSSS(256);

  vscale = 1.0;
  objfcn.setScale(vscale);

  x = nlp.getXc();
  for (i=1; i <= ndim; i++) {
    for (j=1; j <= ndim+1; j++) {
      init_simplex(i,j) = x(i);
    }
  }

  for (i=1; i<= ndim; i++) {
    perturb = x(i)*.01;
    init_simplex(i,i+1) = x(i) + perturb;
  }

  objfcn.setSimplexType(2);
  objfcn.setSimplex(init_simplex);

  objfcn.setCreateFlag();

  objfcn.setSchemeFileName(schemefilename);

  objfcn.optimize();
  
  objfcn.printStatus("Solution from PDS");

  objfcn.cleanup();

#ifdef REG_TEST
  ColumnVector x_sol = nlp.getXc();
  double f_sol = nlp.getF();
  if ((1.0 - x_sol(1) <= 1.e-2) && (1.0 - x_sol(2) <= 1.e-2) && (f_sol
								 <=
								 1.e-2))
    *optout << "PDS 1 PASSED" << endl;
  else
    *optout << "PDS 1 FAILED" << endl;
#endif

#ifdef WITH_MPI
  MPI_Finalize();
#endif    

}
