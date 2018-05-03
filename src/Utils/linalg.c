
#include <math.h>

#include "cblas.h"

void daxpy(int *ndim, double *alpha, double *dx, int *inc1, 
	   double *dy, int *inc2)
     
{
  /* Local variables */
  static int i, m, ix, iy;
  int n = *ndim, incx = *inc1, incy = *inc2;
  double da = *alpha;


  /* Function Body */
  
  /*     constant times a vector plus a vector. */
  /*     uses unrolled loops for increments equal to one. */
  /*     jack dongarra, linpack, 3/11/78. */
  

  if (n <= 0) return;
  if (da == 0.) return;
  if (incx == 1 && incy == 1) {/* code for both increments equal to 1 */

    
    /* clean-up loop */

    m = n % 4;
    if (m != 0) {
      for (i = 0; i < m; ++i) dy[i] += da * dx[i];
      if (n < 4)  return;
    }

    for (i = m; i < n; i += 4) {
      dy[i] += da * dx[i];
      dy[i + 1] += da * dx[i + 1];
      dy[i + 2] += da * dx[i + 2];
      dy[i + 3] += da * dx[i + 3];
    }

    return;
  }

  /*        code for unequal increments or equal increments */
  /*        not equal to 1 */

  else {
    ix = 1;
    iy = 1;
    if (incx < 0) ix = (-n + 1) * incx + 1;
    if (incy < 0) iy = (-n + 1) * incy + 1;
    for (i = 0; i < n; ++i) {
      dy[iy] += da * dx[ix];
      ix += incx;
      iy += incy;
    }
  }
}


void dswap(int *ndim, double *dx, int *inc1, double *dy, int *inc2)

{
  /* Local variables */
  static int i, m;
  static double dtemp;
  static int ix, iy;
  int n = *ndim, incx = *inc1, incy = *inc2;
  
  /* Function Body */
  /*     interchanges two vectors. */
  /*     uses unrolled loops for increments equal one. */
  /*     jack dongarra, linpack, 3/11/78. */
  

  if (n <= 0) return;
  if (incx == 1 && incy == 1) {/* clean-up loop */
    m = n % 3;
    if (m != 0) {
      for (i = 0; i < m; ++i) {
	dtemp = dx[i];
	dx[i] = dy[i];
	dy[i] = dtemp;
      }
      if (n < 3) return;
    }
    for (i = m; i < n; i += 3) {
	dtemp     = dx[i];
	dx[i]     = dy[i];
	dy[i]     = dtemp;
	dtemp     = dx[i + 1];
	dx[i + 1] = dy[i + 1];
	dy[i + 1] = dtemp;
	dtemp     = dx[i + 2];
	dx[i + 2] = dy[i + 2];
	dy[i + 2] = dtemp;
      }
  }
  else {/* code for unequal increments or equal increments not equal to 1 */
    
    ix = 1;
    iy = 1;
    if (incx < 0) ix = (-n + 1) * incx + 1;
    if (incy < 0) iy = (-n + 1) * incy + 1;
    for (i = 0; i < n; ++i) {
      dtemp = dx[ix];
      dx[ix] = dy[iy];
      dy[iy] = dtemp;
      ix += incx;
      iy += incy;
    }
  }
}

/* 
 * Blas 1 routines written in C
 * Implemented: ddot, dnrm2, dscal
 * Blas 1: dasum, daxpy, dcopy, ddot, dnrm2, drot, drotg, dscal, dswap
 *
 * J.C. Meza
 * Sandia National Laboratories
 * meza@california.sandia.gov
 */
double ddot(int *ndim, double *dx, int *inc1, double *dy, int *inc2)

{
  /* Local variables */
  int i, m;
  double dtemp = 0.;
  int ix, iy;
  int n = *ndim, incx = *inc1, incy = *inc2;
  
  /* Function Body */
  /*     forms the dot product of two vectors. */
  /*     uses unrolled loops for increments equal to one. */
  /*     jack dongarra, linpack, 3/11/78. */
  

  dtemp = 0.;
  if (n <= 0) return 0.;

  if (incx == 1 && incy == 1) {/* code for both increments equal to 1 */
    
    /* clean-up loop */
    
    m = n % 5;
    if (m != 0) {
      for (i=0; i<m; ++i) dtemp += dx[i]*dy[i];
      if (n < 5) return dtemp;
    }
    
    for (i=m; i<n; i+=5) {
      dtemp = dtemp + dx[i]*dy[i] + dx[i+1]*dy[i+1] + dx[i+2]*dy[i+2] + 
	dx[i+3]*dy[i+3] + dx[i+4]*dy[i+4];
    }
  }
  
  else { /* code for unequal increments or equal increments not equal to 1 */
  
    ix = 1;    iy = 1;
    if (incx < 0) ix = (-n + 1) * incx + 1;
    if (incy < 0) iy = (-n + 1) * incy + 1;

    for (i=0; i<n; ++i) {
      dtemp += dx[ix] * dy[iy];
      ix += incx;
      iy += incy;
    }
  }
  return dtemp;
  
}

void dscal(int *ndim, double *alpha, double *dx, int *inc1)

{

  /* Local variables */
  static int i, m, ix;
  int n = *ndim, incx = *inc1;
  double da = *alpha;

  /* Function Body */
  /*     scales a vector by a constant. */
  /*     uses unrolled loops for increment equal to one. */
  /*     jack dongarra, linpack, 3/11/78. */
  /*     modified to correct problem with negative increment, 8/21/90. */
  

  if (n <= 0) return;
  if (incx == 1) {/* code for increment equal to 1 */

    /* clean-up loop */
    m = n % 5;
    if (m != 0) {
      for (i = 0; i < m; ++i) dx[i] = da * dx[i];
      if (n < 5) return;
    }

    for (i = m; i < n; i += 5) {
      dx[i] = da * dx[i];
      dx[i + 1] = da * dx[i + 1];
      dx[i + 2] = da * dx[i + 2];
      dx[i + 3] = da * dx[i + 3];
      dx[i + 4] = da * dx[i + 4];
    }
  }
  
  else {/* code for increment not equal to 1 */
      
    ix = 1;
    if (incx < 0) ix = (-n + + 1) * incx + 1;
    for (i = 0; i < n; ++i) {
      dx[ix] = da * dx[ix];
      ix += incx;
    }
  }
}

double dnrm2(int *ndim, double *dx, int *inc1)
{
  int i, ix;
  double sum;
  int n = *ndim, incx = *inc1;

  sum = 0.0;
  if (incx == 1) {
    for (i=0; i<n; ++i) {
      sum += dx[i]*dx[i];
    }
  }
  else {
    ix = 1;
    for (i=0; i<n; ++i) {
      sum += dx[ix];
      ix = ix + incx;
    }
  }
  
  return sqrt(sum);

}

void dcopy(int *ndim, double *dx, int *inc1, double *dy, int *inc2)

{
    /* Local variables */
    static int i, m, ix, iy;
    int n = *ndim, incx = *inc1, incy = *inc2;

    /* Function Body */

    /*     copies a vector, x, to a vector, y. */
    /*     uses unrolled loops for increments equal to one. */
    /*     jack dongarra, linpack, 3/11/78. */


    if (n <= 0) {
	return;
    }
    if (incx == 1 && incy == 1) { /* code for both increments equal to 1 */

      /* clean-up loop */

      m = n % 7;
      if (m != 0) {
	for (i = 0; i < m; ++i) dy[i] = dx[i];
	if (n < 7) return;
      }

      for (i = m; i < n; i += 7) {
	dy[i] = dx[i];
	dy[i + 1] = dx[i + 1];
	dy[i + 2] = dx[i + 2];
	dy[i + 3] = dx[i + 3];
	dy[i + 4] = dx[i + 4];
	dy[i + 5] = dx[i + 5];
	dy[i + 6] = dx[i + 6];
      }
      return;
    }

    else{ /* code for unequal increments or equal increments */
          /* not equal to 1 */

      ix = 1;
      iy = 1;
      if (incx < 0) ix = (-n + 1) * incx + 1;
      if (incy < 0) iy = (-n + 1) * incy + 1;
      for (i = 0; i < n; ++i) {
	dy[iy] = dx[ix];
	ix += incx;
	iy += incy;
      }
      return;
    }
  }

