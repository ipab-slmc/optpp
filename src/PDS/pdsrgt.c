/*
//--------------------------------------------------------------------
// Copyright (C) 1993,1994: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//--------------------------------------------------------------------
*/

#include "pds.h"

int pdsrgt(int ndim, double scale, double *simplex)
{
  /*******************************************************************
   *
   * THIS IS A SERVICE SUBROUTINE USED TO CONSTRUCT AN INITIAL SIMPLEX
   * WHEN ONE IS NOT PROVIDED BY THE USER.  THIS ROUTINE CONSTRUCTS A
   * RIGHT-ANGLE SIMPLEX (I.E., ONE IN WHICH ALL OF THE VERTICES ARE
   * OF EQUAL DISPLACEMENT ALONG THE COORDINATE AXES.  NOTE THAT THE
   * USER MUST PROVIDE THE VARIABLE `SCALE', WHICH SPECIFIES THE
   * LENGTH OF THE DISPLACEMENT, AS WELL AS THE INITIAL VERTEX USED TO
   * START THE SEARCH.
   *
   * WRITTEN BY VIRGINIA TORCZON. 
   *
   * LAST MODIFICATION:  MARCH 11, 1992. 
   *
   * INPUT
   *
   *    N             DIMENSION OF THE PROBLEM TO BE SOLVED
   *
   *    SCALE         DUAL PURPOSE; IF 
   *                     0    THEN USER SUPPLIED THE SIMPLEX TO BE
   *                          USED TO START THE SEARCH AND THIS
   *                          SUBROUTINE WILL NOT BE CALLED
   *                     ELSE THIS SUBROUTINE CONSTRUCTS A RIGHT-ANGLE 
   *                          SIMPLEX WITH EDGES OF LENGTH `SCALE'
   *
   *    SIMPLEX       TWO-DIMENSIONAL ARRAY CONTAINING THE VERTEX 
   *                  SPECIFIED BY THE USER TO BEGIN THE SEARCH 
   * OUTPUT 
   *
   *    SIMPLEX       TWO-DIMENSIONAL ARRAY CONTAINING THE N+1 VERTICES 
   *                  OF THE SIMPLEX TO BE USED TO START THE SEARCH.
   *                  THIS SIMPLEX WILL BE A RIGHT-ANGLE SIMPLEX WITH
   *                  DISPLACEMENTS ALONG THE COORDINATE AXES OF LENGTH 
   *                  (AND ORIENTATION) OF `SCALE'.
   *
   * REMARKS
   *
   *    THE WORKSPACE FOR THE SIMPLEX IS DECLARED IN THE MAIN PROGRAM
   *    TO BE A ONE-DIMENSIONAL VECTOR `SIMPLEX'; HOWEVER, ALL
   *    SUBROUTINE CALLS REDEFINE THE ONE-DIMENSIONAL VECTOR `SIMPLEX'
   *    TO BE A TWO-DIMENSIONAL ARRAY `S' OF DIMENSION N BY N+1.  NOTE
   *    ALSO THAT THE N+1 VERTICES OF THE SIMPLEX ARE STORED BY
   *    COLUMN.
   *
   *******************************************************************/

  static int i, j;

  for (j = 1; j <= ndim; j++) {

    for (i = 0; i < ndim; i++) {
      simplex[i + j * ndim] = simplex[i];
    }

    simplex[j - 1 + j * ndim] += scale;
  }

  return 0;
}

