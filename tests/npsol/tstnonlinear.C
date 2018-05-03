//
// Test program for NPSOL objects
//

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

#include "OptNPSOL.h"

#include "tstfcn.h"

using NEWMAT::ColumnVector;

using namespace OPTPP;

double const LOW_BOUND =  -4.5;
double const UP_BOUND  =   4.5;

int main ()
{
  int          i, n=3;
  ColumnVector x(n);
  
  char *status_file = {"tstnonlinear.out"};
  ofstream opt_fp;
  opt_fp.open(status_file);

  int ndim  = n;
  int nclin = 0;
  int ncnln = 1;
  double ftol = 1.0e-06;
  ColumnVector lower(n+ncnln), upper(n+ncnln);

  //  Create boundary constraints 

  for(i=1; i<=ndim; i++) {lower(i) = LOW_BOUND;  upper(i) = UP_BOUND;}
  lower(n)       = -5.0;     upper(n)       = 5.0;
  lower(n+ncnln) =  -1.0E32; upper(n+ncnln) = 48;

  //  Build a optimization object and optimize 

  OptNPSOL objfcn(ndim,nclin,ncnln,hs65,init_hs65,ineq_hs65); 
  objfcn.setX(x);
  objfcn.setLower(lower);
  objfcn.setUpper(upper);
  objfcn.setFcnAccrcy(ftol);
  objfcn.setOutputFile(opt_fp);
  objfcn.optimize();
  objfcn.printStatus("Solution from NPSOL ");
  opt_fp.close();

}
