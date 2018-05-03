//
// Test program for NPSOL objects
//

#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#include <iostream.h>
#include <fstream.h>
#ifdef HAVE_STD
#include <cstdio>
#else
#include <stdio.h>
#endif

#include "OptNPSOL.h"
#include "ioformat.h"

#include "tstfcn.h"

using NEWMAT::ColumnVector;
using NEWMAT::Matrix;

using namespace OPTPP;

double const LOW_BOUND =  -2.0E0;
double const UP_BOUND  =   2.0E0;

int main ()
{
  int          i, n=2;
  ColumnVector x(n);
  
  char *status_file = {"tstlinear.out"};
  ofstream opt_fp;
  opt_fp.open(status_file);

  int ndim  = n;
  int nclin = 1;
  int ncnln = 0;
  double ftol = 1.0e-10;
  ColumnVector lower(n+nclin), upper(n+nclin);
  Matrix A(nclin, ndim);

  //  Create boundary constraints 

  for(i=1; i<=ndim; i++) {lower(i) = LOW_BOUND;  upper(i) = UP_BOUND;}
  lower(ndim+nclin) = 2*LOW_BOUND;
  upper(ndim+nclin) = 2*UP_BOUND;

  //  Create linear constraints 
  A = 1.0;

  //  Build a optimization object and optimize 

  OptNPSOL objfcn(ndim,nclin,ncnln,rosen,init_rosen);   
  objfcn.setX(x);
  objfcn.setLower(lower);
  objfcn.setUpper(upper);
  objfcn.setMatrix(A);
  objfcn.setFcnAccrcy(ftol);
  objfcn.setOutputFile(opt_fp);
  objfcn.optimize();
  objfcn.printStatus("Solution from NPSOL ");
  opt_fp.close();

}
