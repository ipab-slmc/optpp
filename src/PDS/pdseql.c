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

int pdseql(int ndim, double scale, double *simplex)
{
  /*******************************************************************
   *
   * THIS IS A SERVICE SUBROUTINE USED TO CONSTRUCT AN INITIAL SIMPLEX
   * WHEN ONE IS NOT PROVIDED BY THE USER.  THIS ROUTINE CONSTRUCTS A
   * REGULAR SIMPLEX (I.E., ONE IN WHICH ALL OF THE EDGES ARE OF EQUAL
   * LENGTH) FOLLOWING AN ALGORITHM GIVEN BY JACOBY, KOWALIK, AND
   * PIZZO IN "ITERATIVE METHODS FOR NONLINEAR OPTIMIZATION PROBLEMS,"
   * PRENTICE-HALL (1972).  THIS ALGORITHM ALSO APPEARS IN SPENDLEY,
   * HEXT, AND HIMSWORTH, "SEQUENTIAL APPLICATION OF SIMPLEX DESIGNS
   * IN OPTIMISATION AND EVOLUTIONARY OPERATION," TECHNOMETRICS,
   * VOL. 4, NO. 4, NOVEMBER 1962, PAGES 441--461.  NOTE THAT THE USER
   * MUST PROVIDE AN INITIAL GUESS AT THE SOLUTION, WHICH IS USED AS
   * THE BASE VERTEX V0, AS WELL AS THE VARIABLE `SCALE', WHICH
   * SPECIFIES THE LENGTH OF THE EDGES IN THE INITIAL SIMPLEX.
   *
   * WRITTEN BY VIRGINIA TORCZON. 
   *
   * LAST MODIFICATION:  MAY 4, 1992.
   * 
   *    INPUT 
   *
   *       N             DIMENSION OF THE PROBLEM TO BE SOLVED
   *
   *       SCALE         DUAL PURPOSE; IF 
   *                        0    THEN USER SUPPLIED THE SIMPLEX TO BE
   *                             USED TO START THE SEARCH AND THIS
   *                             SUBROUTINE WILL NOT BE CALLED
   *                        ELSE THIS SUBROUTINE CONSTRUCTS A REGULAR
   *                             SIMPLEX WITH EDGES OF LENGTH `SCALE'
   *
   *       SIMPLEX       TWO-DIMENSIONAL ARRAY CONTAINING THE VERTEX 
   *                     SPECIFIED BY THE USER TO BEGIN THE SEARCH
   *
   *    OUTPUT
   *
   *       SIMPLEX       TWO-DIMENSIONAL ARRAY CONTAINING THE N+1
   *                     VERTICES OF THE SIMPLEX TO BE USED TO START
   *                     THE SEARCH.  THIS SIMPLEX WILL BE A REGULAR
   *                     SIMPLEX WITH EDGES OF LENGTH `SCALE'. 
   *    REMARKS
   *
   *       THE WORKSPACE FOR THE SIMPLEX IS DECLARED IN THE MAIN
   *       PROGRAM TO BE A ONE-DIMENSIONAL VECTOR `SIMPLEX'; HOWEVER,
   *       ALL SUBROUTINE CALLS REDEFINE THE ONE-DIMENSIONAL VECTOR
   *       `SIMPLEX' TO BE A TWO-DIMENSIONAL ARRAY `S' OF DIMENSION N
   *       BY N+1.  NOTE ALSO THAT THE N+1 VERTICES OF THE SIMPLEX ARE
   *       STORED BY COLUMN.
   *
   * CALCULATE THE AUXILIARY VALUES 'P' AND 'Q' AS FOLLOWS:
   *   
   *             SQRT(N + 1) - 1 + N 
   *      P  =  --------------------- * SCALE 
   *                 N * SQRT(2) 
   *
   *               SQRT(N + 1) - 1 
   *      Q  =    ----------------- * SCALE 
   *                 N * SQRT(2) 
   *
   *******************************************************************/

  static double temp;
  static int i, j;
  static double p, q;

  temp = ndim + 1.;
  q = ((sqrt(temp) - 1.) / (ndim * sqrt(2.))) * scale;
  p = q + 1. / sqrt(2.) * scale;

  /* NOW CALCULATE THE COORDINATES OF THE REMAINING VERTICES OF THE 
   * SIMPLEX AS FOLLOWS: 
   *
   * S(*,1)   = (S(1,1),     S(2,1),     S(3,1),     ..., S(N,1)) 
   * S(*,2)   = (S(1,1) + P, S(2,1) + Q, S(3,1) + Q, ..., S(N,1) + Q) 
   *                .
   *                . 
   *                . 
   *                . 
   * S(*,N+1) = (S(1,1) + Q, S(2,1) + Q, S(3,1) + Q, ..., S(N,1) + P) */

  for (j = 1; j <= ndim; j++) {

      for (i = 0; i < j - 1; i++) 
	simplex[i + j*ndim] = simplex[i] + q;

      i = j - 1;
      simplex[i + j*ndim] = simplex[i] + p;

      for (i = j; i < ndim; i++) 
	simplex[i + j*ndim] = simplex[i] + q;
    }
  
  return 0;
}

