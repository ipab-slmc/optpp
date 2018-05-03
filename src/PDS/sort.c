/*
//--------------------------------------------------------------------
// Copyright (C) 1993,1994: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//--------------------------------------------------------------------
*/

#include "pds.h"

int sort(int ndim, int *total, int *array, int *scheme, int *error)
{
  /*******************************************************************
   *
   * THIS IS A FORTRAN IMPLEMENTATION OF THE NON-RECURSIVE VERSION FOR
   * SORTING AN ARRAY OF NUMBERS. THE ALGORITHM IS THAT WHICH APPEARS
   * AS PROGRAM 2.11 IN WIRTH'S BOOK "ALGORITHMS + DATA STRUCTURES =
   * PASCAL PROGRAMS", PRENTICE-HALL, 1976. THE MODIFICATION TO
   * MINIMIZE THE STACK SIZE DESCRIBED IN (2.16) OF THIS REFERENCE IS
   * ALSO USED.
   *
   * -----ON INPUT-----
   *
   *       ARRAY  -- INDEX ARRAY OF POINTERS FOR THE N-TUPLES IN
   *                 `SCHEME'
   *
   *       SCHEME -- THE INTEGER N-TUPLES FOR THE SEARCH STRATEGY NOTE
   *                 THAT WE ARE ACTUALLY SORTING THE N-TUPLES IN
   *                 `SCHEME', BUT WE USE `ARRAY' TO STORE THIS
   *                 ORDERING
   *
   *       N      -- THE DIMENSION OF THE PROBLEM TO BE SOLVED
   *
   *       TOTAL  -- THE NUMBER OF "RAW" POINTS GENERATED FOR `SCHEME'
   *                 NOTE THAT THE NUMBER OF ITEMS TO BE SORTED IS
   *                 THUS N + 1 + TOTAL
   *
   * -----ON OUTPUT----
   *
   *       ARRAY  -- THE INDEX DENOTING A SORTING FOR THE N-TUPLES IN
   *                 `SCHEME'.  NOTE THAT THE ORIGINAL ARRAY OF
   *                 NUMBERS IS WRITTEN OVER AND LOST
   *
   *       ERROR  -- FLAG TO INDICATE THAT THE ARRAY TO BE SORTED IS
   *                 TOO LONG FOR THE INTERNAL STACK (SEE COMMENTS
   *                 BELOW)
   *
   * WRITTEN BY R.M.LEWIS JULY,1986
   *
   * LAST MODIFICATION BY VIRGINIA TORCZON FEBRUARY 27, 1992
   *
   * WIRTH'S REMEDY TO KEEP THE LENGTH OF THE STACK STACK (`M') LESS
   * THAN THE LENGTH OF THE LIST TO BE SORTED (`N') LIES IN STACKING
   * THE SORT REQUEST FOR THE LONGER PARTITION AND CONTINUING DIRECTLY
   * WITH THE FURTHER PARTITIONING OF THE SMALLER SECTIONS.  IF THIS
   * IS ENFORCED THEN THE LENGTH OF THE STACK CAN BE LIMITED TO M =
   * LOG N.
   *    2
   *
   * IN THIS IMPLEMENTATION, WE HAVE ALLOCATED THE STACK LOCALLY TO BE
   * OF LENGTH M, WHICH ALLOWS US TO SORT ARRAYS OF LENGTH UP TO 2**M.
   * TO MAKE SURE THAT WE DO NOT EXCEED THIS UPPER BOUND, THE
   * FOLLOWING TEST IS MADE BEFORE BEGINNING THE SEARCH.
   *
   *******************************************************************/

  /* System generated locals */

  int array_offset, scheme_dim1, scheme_offset;

  /* Local variables */

  static int temp, i, j, l, r, s, x, stack[64]	/* was [32][2] */;
  extern int order(int, int *, int *);
  static int power, remain, length;

  /* Parameter adjustments */

  array_offset = -(ndim);
  array -= array_offset;
  scheme_dim1 = ndim + 2;
  scheme_offset = scheme_dim1 * (-(ndim)) - 1;
  scheme -= scheme_offset;

  /* Function Body */

  power = 0;
  remain = 0;
  length = *total + ndim + 1;
L50:

  if (length > 1) {
    ++power;
    remain += length % 2;
    length /= 2;
    goto L50;
  }

  if (remain > 0) {
    ++power;
  }

  if (power <= 32) {

    /* INITIALIZE THE STACK : WE BEGIN THE PARTITIONING AT THE LEFT
     * AND RIGHT ENDS OF THE ARRAY */

    s = 1;
    stack[s - 1] = -(ndim);
    stack[s + 31] = *total;

    /* BEGIN THE LOOP TO SIMULATE RECURSION */

L100:

    /* TAKE THE TOP REQUEST OFF THE STACK */

    l = stack[s - 1];
    r = stack[s + 31];
    --s;

    /* SPLIT { ARRAY(L) , .... , ARRAY(R) } */

L200:
    i = l;
    j = r;
    x = array[(l + r) / 2];

L300:

    /* THE FOLLOWING IS EQUIVALENT TO: DO WHILE 
     *                                    (SCHEME(*,I) < SCHEME(*,X))
     *                                    I = I + 1
     *                                 END DO */

L330:
    if (order(ndim, &scheme[array[i] * scheme_dim1 + 1], &scheme[x * 
				       scheme_dim1 + 1]) < 0) {
      ++i;
      goto L330;
    }

    /* THE FOLLOWING IS EQUIVALENT TO: DO WHILE
     *                                    (SCHEME(*,J) > SCHEME(*,X))
     *                                    J = J + 1
     *                                 END DO */

L350:
    if (order(ndim, &scheme[array[j] * scheme_dim1 + 1], &scheme[x * 
				       scheme_dim1 + 1]) > 0) {
      --j;
      goto L350;
    }

    /* PERFORM A SWAP, AND THEN MOVE ON IN THE SCAN */

    if (i <= j) {
      temp = array[i];
      array[i] = array[j];
      array[j] = temp;

      ++i;
      --j;
    }

    if (i <= j) {
      goto L300;
    }

    /* PLACE ON THE STACK THE REQUEST TO SORT THE LONGER PARTITION AND
     * BEGIN FURTHER PARTIONING OF THE SMALLER PARTITIONS */

    if (j - l < r - i) {

      if (i < r) {
	++s;
	stack[s - 1] = i;
	stack[s + 31] = r;
      }

      r = j;
    } else {

      if (l < j) {
	++s;
	stack[s - 1] = l;
	stack[s + 31] = j;
      }

      l = i;
    }

    /* THE FOLLOWING TWO IF-THEN BLOCKS MARK THE ENDS OF TWO
     * REPEAT-UNTIL BLOCKS. CAN YOU SEE WHAT THEY ARE ??? */

    if (l < r) {
      goto L200;
    }

    if (s > 0) {
      goto L100;
    } else {
      return 0;
    }
  }

  *error = 1;
  
  return 0;
}

int quick(int ndim, int *array, int *error)
{
  /*******************************************************************
   *
   * THIS IS A FORTRAN IMPLEMENTATION OF THE NON-RECURSIVE VERSION FOR
   * SORTING AN ARRAY OF NUMBERS. THE ALGORITHM IS THAT WHICH APPEARS
   * AS PROGRAM 2.11 IN WIRTH'S BOOK "ALGORITHMS + DATA STRUCTURES =
   * PASCAL PROGRAMS", PRENTICE-HALL, 1976. THE MODIFICATION TO
   * MINIMIZE THE STACK SIZE DESCRIBED IN (2.16) OF THIS REFERENCE IS
   * ALSO USED.
   *
   * -----ON INPUT-----
   *
   *       ARRAY  -- THE ARRAY OF NUMBERS TO BE SORTED : THIS ROUTINE
   *                 ASSUMES THAT THESE NUMBERS ARE INTEGER.
   *
   *       N      -- THE NUMBER OF ITEMS TO BE SORTED
   *
   * -----ON OUTPUT----
   *
   *       ARRAY  -- THE NUMBERS SORTED IN ASCENDING ORDER : NOTE THAT
   *                 THE ORIGINAL ARRAY OF NUMBERS IS WRITTEN OVER AND
   *                 LOST
   *
   *       ERROR  -- FLAG TO INDICATE THAT THE ARRAY TO BE SORTED IS
   *                 TOO LONG FOR THE INTERNAL STACK (SEE COMMENTS
   *                 BELOW)
   *
   *     WRITTEN BY R.M.LEWIS JULY,1986
   *
   *     LAST MODIFICATION BY VIRGINIA TORCZON FEBRUARY 27, 1992
   *
   * WIRTH'S REMEDY TO KEEP THE LENGTH OF THE STACK STACK (`M') LESS
   * THAN THE LENGTH OF THE LIST TO BE SORTED (`N') LIES IN STACKING
   * THE SORT REQUEST FOR THE LONGER PARTITION AND CONTINUING DIRECTLY
   * WITH THE FURTHER PARTITIONING OF THE SMALLER SECTIONS.  IF THIS
   * IS ENFORCED THEN THE LENGTH OF THE STACK CAN BE LIMITED TO M =
   * LOG N.
   *    2
   *
   * IN THIS IMPLEMENTATION, WE HAVE ALLOCATED THE STACK LOCALLY TO BE
   * OF LENGTH M, WHICH ALLOWS US TO SORT ARRAYS OF LENGTH UP TO 2**M.
   * TO MAKE SURE THAT WE DO NOT EXCEED THIS UPPER BOUND, THE
   * FOLLOWING TEST IS MADE BEFORE BEGINNING THE SEARCH.
   *
   *******************************************************************/

  static int temp, i, j, l, r, s, x, stack[64]	/* was [32][2] */, 
    power, remain, length;

  /* Parameter adjustments */

  --array;

  /* Function Body */

  power = 0;
  remain = 0;
  length = ndim;

L50:
  if (length > 1) {
    ++power;
    remain += length % 2;
    length /= 2;
    goto L50;
  }

  if (remain > 0) {
    ++power;
  }

  if (power <= 32) {

    /* INITIALIZE THE STACK : WE BEGIN THE PARTITIONING AT THE LEFT
     * AND RIGHT ENDS OF THE ARRAY */

    s = 1;
    stack[s - 1] = 1;
    stack[s + 31] = ndim;


    /* BEGIN THE LOOP TO SIMULATE RECURSION */

L100:

    /* TAKE THE TOP REQUEST OFF THE STACK */

    l = stack[s - 1];
    r = stack[s + 31];
    --s;

    /* SPLIT { ARRAY(L) , .... , ARRAY(R) } */

L200:
    i = l;
    j = r;
    x = array[(l + r) / 2];

L300:


    /* THE FOLLOWING IS EQUIVALENT TO :   DO WHILE( ARRAY(I) < X )
     *                                       I = I + 1
     *                                    END DO */

L330:
    if (array[i] < x) {
      ++i;
      goto L330;
    }

    /* THE FOLLOWING IS EQUIVALENT TO :   DO WHILE( ARRAY(J) > X )
     *                                        J = J + 1
     *                                     END DO */

L350:
    if (array[j] > x) {
      --j;
      goto L350;
    }

    /* PERFORM A SWAP, AND THEN MOVE ON IN THE SCAN */

    if (i <= j) {
      temp = array[i];
      array[i] = array[j];
      array[j] = temp;

      ++i;
      --j;
    }

    if (i <= j) {
      goto L300;
    }


    /* PLACE ON THE STACK THE REQUEST TO SORT THE LONGER PARTITION AND
     * BEGIN FURTHER PARTIONING OF THE SMALLER PARTITIONS */

    if (j - l < r - i) {

      if (i < r) {
	++s;
	stack[s - 1] = i;
	stack[s + 31] = r;
      }

      r = j;
    } else {

      if (l < j) {
	++s;
	stack[s - 1] = l;
	stack[s + 31] = j;
      }

      l = i;
    }

    /* THE FOLLOWING TWO IF-THEN BLOCKS MARK THE ENDS OF TWO
     * REPEAT-UNTIL BLOCKS. CAN YOU SEE WHAT THEY ARE ??? */

    if (l < r) {
      goto L200;
    }

    if (s > 0) {
      goto L100;
    } else {
      return 0;
    }
  }

  *error = 1;

  return 0;
}

