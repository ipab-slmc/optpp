
#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#ifdef HAVE_STD
#include <cstdio>
#include <cmath>
#else
#include <stdio.h>
#include <math.h>
#endif

#define WANT_MATH

#include "newmatap.h"

#include "NLP.h"
#include "NLF.h"

using NEWMAT::ColumnVector;
using NEWMAT::SymmetricMatrix;
using std::cout;

using namespace OPTPP;

static int    m, n, fcn_count=0;
static double **A, **patch, **lamp;
static SymmetricMatrix HH(1);

#if !(defined(__GNUC__) && __GNUC__ >=3)
  #ifndef SUN4
     const double M_PI = 3.14159265358979323846;
  #endif
#endif

/* A Test Example - The illumination example in Boyd */

void init_illum (int ndim, ColumnVector& x)
{
  int    i, j, k;
  double h,rij2,midx,midy,dx,dy,slope1,slope2,theta1,theta2,theta,scale;
  double dtmp, Ides;

  // allocate storage

  m = ndim;
  n = 11;
  Ides = 2.0;
  patch = new double*[n+1];
  for (i=0; i<=n; i++) patch[i] = new double[2];
  lamp = new double*[m+1];
  for (i=0; i<=m; i++) lamp[i] = new double[2];
  A = new double*[n+1];
  for (i=0; i<=n; i++) A[i] = new double[m+1];

  // initializing the patches and lamps

  patch[0][0] = 0.0;     patch[0][1] = 0.0;
  patch[1][0] = 0.0909;  patch[1][1] = 0.1;
  patch[2][0] = 0.1818;  patch[2][1] = 0.2;
  patch[3][0] = 0.2727;  patch[3][1] = 0.2;
  patch[4][0] = 0.3636;  patch[4][1] = 0.1;
  patch[5][0] = 0.4545;  patch[5][1] = 0.2;
  patch[6][0] = 0.5455;  patch[6][1] = 0.3;
  patch[7][0] = 0.6364;  patch[7][1] = 0.2;
  patch[8][0] = 0.7273;  patch[8][1] = 0.0;
  patch[9][0] = 0.8182;  patch[9][1] = 0.0;
  patch[10][0] = 0.9091; patch[10][1] = 0.2;
  patch[11][0] = 1.0;    patch[11][1] = 0.1;

  lamp[1][0] = 0.1;  lamp[1][1] = 1.0;
  lamp[2][0] = 0.3;  lamp[2][1] = 1.1;
  lamp[3][0] = 0.4;  lamp[3][1] = 0.6;
  lamp[4][0] = 0.6;  lamp[4][1] = 0.9;
  lamp[5][0] = 0.8;  lamp[5][1] = 0.9;
  lamp[6][0] = 0.9;  lamp[6][1] = 1.2;
  lamp[7][0] = 0.95; lamp[7][1] = 1.0;

  // initialize the A matrix

  for (i=1; i<=n; i++) {
    midx  = patch[i][0] - patch[i-1][0];
    midy  = patch[i][1] - patch[i-1][1];
    if (midx == 0.0)  {
      cout << "Error : patches cannot overlap each other. \n";
      exit(-1);
    }
    slope1 = midy / midx;
    theta1 = atan(slope1);
    midx  = 0.5 * midx + patch[i-1][0];
    midy  = 0.5 * midy + patch[i-1][1];
    for (j=1; j<=m; j++) {
      dx = lamp[j][0] - midx;
      dy = lamp[j][1] - midy;
      rij2 = dx * dx + dy * dy;
      if (dx == 0.0) theta2 = dy / fabs(dy) * M_PI * 0.5;
      else {
        slope2 = dy / dx;
        theta2 = atan(slope2);
      }
      if (theta2 < 0.0) theta2 = theta2 + M_PI;
      theta = M_PI * 0.5 - (theta2 - theta1);
      scale = cos(theta);
      if (scale < 0.0) A[i][j] = 0.0;
      else A[i][j] = scale / (rij2 * Ides);
    }
  }

  // initialize x 

  for (i=1; i<=m; i++) x(i) = 0.5;

  // initialize the Hessian

  HH.ReSize(m);
  for (i=1; i<=m; i++) {
    for (j=1; j<=m; j++) {
      dtmp = 0.0;
      for (k=1; k<=n; k++) dtmp += A[k][i] * A[k][j];
      HH(i,j) = dtmp;
    }
  }
}

void illum2(int mode, int nn, const ColumnVector& x, double& fx, 
	ColumnVector& g, SymmetricMatrix& H, int &result)
{ 
  fcn_count++;

  int    i, j;
  double dtmp;

  // output debug message

//    cout << "\n illum2: mode = " << mode
//       << " count = " << fcn_count << "\n";
//    for(i=1; i<=m; i++) 
//        cout << "x(" << i << ") = " << x(i) << "\n";

  // compute function value 

  if(mode & NLPFunction){
    fx = 0.0;
    for (i=1; i<=n; i++) {
      dtmp = 0.0;
      for (j=1; j<=m; j++)  dtmp += A[i][j] * x(j);
      dtmp = (1.0 - dtmp) * (1.0 - dtmp);
      fx = fx + dtmp;
    }
    fx = sqrt(fx);
    result = NLPFunction;
  }

  if(mode & NLPGradient){
  // compute gradient information 
    for (i=1; i<=m; i++) {
      dtmp = 0.0;
      for (j=1; j<=m; j++) dtmp += HH(i,j) * x(j);
      for (j=1; j<=n; j++) dtmp -= A[j][i];
      g(i) = 2.0 * dtmp;
    }
    result = NLPGradient;
  }


  if(mode & NLPHessian){
    H = HH;
    result = NLPHessian;
  }
}

void illum(int mode, int nn, const ColumnVector& x, double& fx, 
	ColumnVector& g, int &result)
{ 
  fcn_count++;

  int    i, j;
  double dtmp;


  // compute function value 

  if(mode & NLPFunction){
    fx = 0.0;
    for (i=1; i<=n; i++) {
      dtmp = 0.0;
      for (j=1; j<=m; j++) dtmp += A[i][j] * x(j);
        dtmp = (1.0 - dtmp) * (1.0 - dtmp);
      fx = fx + dtmp;
    }
    fx = sqrt(fx);
    result = NLPFunction;
  }

  // compute gradient information 
  
  if (nn > 0){
    if(mode & NLPGradient){
       for (i=1; i<=m; i++) {
          dtmp = 0.0;
          for (j=1; j<=m; j++) dtmp += HH(i,j) * x(j);
          for (j=1; j<=n; j++) dtmp -= A[j][i];
          g(i) = 2.0 * dtmp;
       }
       result = NLPGradient;
    }
  }
}

