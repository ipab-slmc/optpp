//------------------------------------------------------------------------
// Generating Set Class - for use with OptGSS
//------------------------------------------------------------------------

/*------------------------------------------------------------------------
 Copyright (c) 2003,
 Ricardo Oliva (raoliva@lbl.gov)
 Lawrence Berkeley National Laboratory
 ------------------------------------------------------------------------*/

#include "GenSetBase.h"

using NEWMAT::ColumnVector;
using NEWMAT::Matrix;

using std::cerr;
using std::endl;

namespace OPTPP {
//
// Methods for Generating Directions
//


bool GenSetBase::generateAll(Matrix& M, ColumnVector& X, double D){ 
  if (Size<=0 || Vdim<=0) {
//    cerr << "***ERROR: GenSetBase::generateAll(Matrix,...) "
//	 << "called with size=" << Size << ", vdim=" << Vdim << endl;
    return false;
  }
  if (M.Ncols() != Size || M.Nrows() != Vdim) {
//    cerr << "***ERROR: GenSetBase::generateAll(Matrix,...) "
//	 << "dimesion of M expected to be "
//	 << Vdim << "-by-" << Size
//	 << " but is " << M.Nrows() << "-by-" << M.Ncols()
//	 << endl;
    return false;
  }
  ColumnVector xi(Vdim);
  for (int i=1; i<=Size; i++) {
    generate(i, D, X, xi);
    M.Column(i) = xi;
  }
  return true;
}

bool GenSetBase::generateAllActive(Matrix& M, ColumnVector& X, double D){ 
  if (Size<=0 || Vdim<=0 || nActive()<=0) {
//    cerr << "***ERROR: GenSetBase::generateAllActive(Matrix,...) "
//	 << "called with size=" << Size << ", vdim=" << Vdim
//	 << " nActive = " << nActive()
//	 << endl;
    return false;
  }
  if (M.Ncols() != nActive() || M.Nrows() != Vdim ) {
//    cerr << "***ERROR: GenSetBase::generateAllActive(Matrix,...) "
//	 << "dimesion of M expected to be "
//	 << Vdim << "-by-" << nActive()
//	 << " but is " << M.Nrows() << "-by-" << M.Ncols()
//	 << endl;
    return false;
  }
  ColumnVector xi(Vdim);
  for (int i=1; i<=nActive(); i++) {
    generateActive(i, D, X, xi);
    M.Column(i) = xi;
  }
  return true;
}
/*
ColumnVector GenSetBase::generate(int i) { 
  ///< returns d_i, the ith basis element
  ColumnVector y(Vdim);
  y = 0;
  generate(i, 0, y, y); 
  return y;
}

void GenSetBase::generate(int i, ColumnVector& y) {
  ///< Stores d_i in y
  y = 0;
  generate(i, 0, y, y);
}

ColumnVector GenSetBase::generate(int i, double a, ColumnVector& x) {
  ///< Returns  x + a * d_i, with d_i = ith basis element
  ColumnVector y(Vdim);
  generate(i, a, x, y);
  return y;
}
*/
// -- purely virtual--  to be defined in  each derived class:
//
// generate(int i, double a, ColumnVector &x, ColumnVector &y)
// 
// ///< Set  y = x + a * d_i; should allow y and x to be same vector.
//


Matrix GenSetBase::pllMesh(int P, ColumnVector& xc, ColumnVector& xn, double r) 
  // P : num points we want ( ~ num of processors)
  // xc: the current point; 
  // xn: the "newton" point, tells the dir in which the mesh is grown
  //  r: ||xc-xn||, if known.
  //
  // return matrix M with genset generated at y_k = xc+k(xn-xc)
  // with radius r_k = k^n*r_0, r_0=||xc-xn||/2
  //
  // *** current implementation not optimized for efficiency ***
  //
{

    int k = 0;
    ColumnVector xk; 
    double       rk; 
    Matrix M, A, B, C;

    ColumnVector ns = xn - xc;  // newton step

    int m = Vdim;
    int n = Size;

    // return at least the base point
    M = xn;
    int nump = P-1; 

    //--
    // main loop:
    //--
    if (r <= 0) r = Norm2(ns);
    double r0 = 0.25*r;
    double pert = r0*.05;
    while (nump > 0) {

      // generate points xc + k*newton_step + r^k*d(i), i=1..size()
      ++k;
      xk = xc + k*ns;
      rk = k*r0;
      A = generateAll(xk,rk);

      // perturbation to avoid potential overlaps
      int RAND_HALF = RAND_MAX / 2;
      for (int i=1; i<=n; i++)
	for (int j=1; j<=m; j++) {
	  int sig = (rand() > RAND_HALF)? 1 : -1;
	  double amp = rand() / RAND_MAX * pert;	
	  A(j,i) = A(j,i) + sig * amp ;
	}

      // after the first iteration want to add xk to M
      if (k>1) {
	B = M | xk;
	M = B;
	--nump;
      }

      // addd to M as many columns as we can afford
      if (nump > n) {
	B = M | A;
      }
      else if (nump>0) { 
	C = A.SubMatrix(1,m,1,nump);
	B = M | C;
      }
      
      M = B;
      nump = nump - n;

    }
    return M;
}

} // namespace OPTPP
