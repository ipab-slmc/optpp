/*
//--------------------------------------------------------------------
// Copyright (C) 1993,1994: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//--------------------------------------------------------------------
*/

#include <math.h>

#include "pds.h"
#include "cblas.h"

double pdslen(int ndim, int type, double *s, double scale,
	      double *work)
{
  /*******************************************************************
   *
   * THIS IS A SERVICE FUNCTION USED TO DETERMINE THE LENGTH OF THE
   * LONGEST EDGE IN THE INITIAL SIMPLEX.  NOTE THAT THE JOB IS MADE
   * CONSIDERABLY EASIER IF THE USER SPECIFIED EITHER A RIGHT-ANGLED
   * OR REGULAR SIMPLEX.  WRITTEN BY VIRGINIA TORCZON.  LAST
   * MODIFICATION: MAY 18, 1992.
   *
   * INPUT 
   *    N             DIMENSION OF THE PROBLEM TO BE SOLVED
   *
   *    TYPE          TYPE OF SIMPLEX SPECIFIED BY THE USER CONTAINS
   *                  THE N+1 VERTICES IN THE INITIAL SIMPLEX (STORED
   *                  BY COLUMN)
   *
   *    SCALE         (RELATIVE) LENGTH OF EDGES IN THE INITIAL
   *                  SIMPLEX IF THE USER REQUESTED THAT THE SIMPLEX
   *                  BE CONSTRUCTED (RATHER THAN ENTERED DIRECTLY)
   *
   *    WORK          ONE-DIMENSIONAL WORK ARRAY OF LENGTH N USED TO
   *                  COMPUTE THE DIFFERENCE BETWEEN VERTICES 
   *
   *******************************************************************/

  /* System generated locals */

  double ret_val;

  /* Local variables */

  static double temp;
  static int i, j, k;
  int incx = 1;

  if (type == 1) {

    /* THE USER CHOSE A SIMPLE RIGHT-ANGLED SIMPLEX SO THE MAXIMUM
     * LENGTH IS EASY TO COMPUTE. */

    ret_val = fabs(scale) * sqrt(2.);
  } 
  else if (type == 2) {

    /* THE USER CHOSE A REGULAR SIMPLEX; I.E., THE LENGTH OF EVERY
     * EDGE IN THE SIMPLEX IS SPECIFIED BY THE CHOICE OF SCALE. */

    ret_val = fabs(scale);
  } 
  else {

    /* THE USER EITHER ENTERED THE SIMPLEX OR CHOSE A SCALED
     * RIGHT-ANGLED SIMPLEX.  COMPUTE THE LENGTHS OF ALL THE EDGES. */
    
    ret_val = 0.;

    for (i = 0; i <= ndim; i++) {

      for (j = i + 1; j <= ndim ; j++) {

	for (k = 0; k < ndim; k++) {
	  work[k] = s[k + j * ndim] - s[k + i * ndim];
	}

	temp = dnrm2(&ndim, work, &incx);

	if (temp > ret_val) {
	  ret_val = temp;
	}
      }
    }
  }

  return ret_val;
}

