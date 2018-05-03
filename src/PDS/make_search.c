/*
//--------------------------------------------------------------------
// Copyright (C) 1993,1994: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//--------------------------------------------------------------------
*/

#include "pds.h"

int make_search(int ndim, int out, int *max, int *scheme, int *index,
		int *list, int *unique, int *factor, int *error)
{
  /*******************************************************************
   *
   * THIS IS A CONTROLLING ROUTINE TO CONSTRUCT AN ACTUAL SEARCH
   * STRATEGY.  FOR A FURTHER DESCRIPTION OF THIS ALGORITHM, AS IT
   * RELATES TO THE PARALLEL DIRECT SEARCH METHODS, SEE J. E. DENNIS,
   * JR. AND VIRGINIA TORCZON, "DIRECT SEARCH METHODS ON PARALLEL
   * MACHINES," SIAM J. OPTIMIZATION, VOL. 1, NO. 4, PP. 448-474,
   * NOVEMBER 1991.
   *
   * THE CONSTRUCTION OF THIS SEARCH STRATEGY IS DISCUSSED ON
   * PP. 459-462 AND PP. 463-466; THE ALGORITHM IMPLEMENTED BELOW IS
   * FORMALLY STATED IN TABLE 3, P. 465.
   *
   * WRITTEN BY VIRGINIA TORCZON
   *
   * LAST MODIFICATION:  MARCH 11, 1992.
   *
   * PARAMETERS
   *
   *    NOTE THAT THESE RATHER MESSY DIMENSION STATEMENTS ARE TO ALLOW
   *    THE MAXIMUM POSSIBLE NUMBER OF COLUMNS IN `SCHEME' GIVEN `N',
   *    THE DIMENSION OF THE PROBLEM, AND `MAX', THE SIZE OF `SCHEME'
   *    DECLARED IN THE MAIN PROGRAM.
   *
   *    INPUT
   *
   *       N              THE DIMENSION OF THE PROBLEM TO BE SOLVED
   *
   *       MAX            THE DECLARED DIMENSION OF THE VECTOR `SCHEME'
   *
   *       OUT            UNIT NUMBER FOR THE FILE TO WHICH THE POINTS
   *                      IN THE SEARCH SCHEME ARE TO BE WRITTEN
   *    WORK
   *
   *       SCHEME        A WORK VECTOR PASSED FROM THE MAIN PROGRAM
   *                     THAT IS REPARTITIONED INTO A MATRIX WITH ROWS
   *                     FROM -1 TO N AND AS MANY COLUMNS AS THE
   *                     ORIGINAL DIMENSION OF THE VECTOR WILL ALLOW.
   *                     THIS ARRAY IS USED TO HOLD THE N-TUPLES (IN
   *                     ROWS 1 THROUGH N) FOR EVERY POINT GENERATED
   *                     AS A POSSIBLE POINT IN THE SEARCH STRATEGY.
   *                     IN ROWS -1 AND 0 ARE THE SCALAR `A' AND THE
   *                     POINT `BEST' NECESSARY TO RECONSTRUCT THE
   *                     SIMPLEX ASSOCIATED WITH THAT N-TUPLE.
   *
   *       INDEX         ARRAY TO KEEP TRACK OF EACH COLUMN OF
   *                     `SCHEME'
   *
   *       LIST          ARRAY TO KEEP TRACK OF A LIST OF DISTINCT
   *                     N-TUPLES IN SCHEME (IF AN N-TUPLE OCCURS MORE
   *                     THAN ONCE IN `SCHEME', `LIST' RECORDS ONLY
   *                     THE FIRST OCCURRENCE)
   *    OUTPUT
   *
   *       UNIQUE        THE FINAL LENGTH OF `LIST' (I.E., THE TOTAL
   *                     NUMBER OF DISTINCT POINTS GENERATED)
   *
   *       FACTOR        THE SCALING FACTOR USED TO GENERATE THE LIST;
   *                     `FACTOR' IS CHOSEN TO ALLOW ALL THE WORK TO
   *                     BE DONE IN INTEGER ARITHMETIC BOTH TO SAVE
   *                     STORAGE AND TO ELIMINATE ANY AMBIGUITIES
   *                     ABOUT THE DISTINCTIVENESS OF
   *                     N-TUPLE---AMBIGUITIES THAT COULD POSSIBLE
   *                     OCCUR WHEN USING FLOATING POINT ARITHMETIC
   *
   *       ERROR         ERROR FLAG TO SIGNAL WHETHER OR NOT PREMATURE
   *                     TERMINATION OCCURRED.  IF `ERROR' IS SET TO
   *                     ZERO, THEN NO ERRORS HAVE BEEN FLAGGED.
   *                     CURRENTLY, THE ONLY ERROR FLAGGED IS WHEN THE
   *                     QUICKSORT ROUTINES RETURN BECAUSE THE LENGTH
   *                     OF THE ARRAY PASSED TO SORT EXCEEDS THE
   *                     CAPACITY OF THE INTERNAL STACK.  FOR FURTHER
   *                     INFORMATION ON THIS ERROR, SEE THE COMMENTS
   *                     IN THE ROUTINES `SORT' AND `QUICK'.
   *
   * LOCAL CONSTANTS
   *
   *    THE PARAMETERS ALPHA, BETA, AND GAMMA SET THE SIZE OF THE *
   *    REFLECTION, CONTRACTION, AND EXPANSION STEPS, RESPECTIVELY.  *
   *    THESE CONSTANTS ARE CURRENTLY SET TO THE DEFAULT VALUES *
   *    RECOMMENDED BY NELDER AND MEAD (NOTE THAT WE DIVIDE BY BETA!).
   *    * THE ALGORITHM BELOW ASSUMES THAT EACH CONSTANT IS A POWER OF
   *    TWO * SO THAT THE WORK CAN BE DONE IN INTEGER ARITHMETIC.  THE
   *    CRITICAL * "PIECES" OF THE SEARCH SCHEME RESIDE IN THE DATA
   *    STRUCTURE * `SCHEME'.  THE MATRIX `SCHEME' CAN BE THOUGHT OF
   *    AS TWO VECTORS, * `A', WHICH IS STORED IN THE -1 ROW OF
   *    `SCHEME', AND `BEST', WHICH * IS STORED IN THE 0 ROW OF
   *    `SCHEME', AND ONE MATRIX `NTUPLES', * WHICH IS STORED IN THE
   *    REMAINING ROWS OF `SCHEME'.  THE PURPOSE * OF THESE PIECES ARE
   *    AS FOLLOWS:
   *
   *       A        AN INTEGER ARRAY WHICH CONTAINS THE SCALAR NECESSARY
   *                TO RECONSTRUCT THE SIMPLEX ASSOCIATED WITH A GIVEN
   *                POINT I
   *
   *       BEST     AN INTEGER ARRAY WHICH CONTAINS THE POINTER TO THE
   *                VERTEX IN THE ORIGINAL SIMPLEX ASSOCIATED WITH A
   *                GIVEN POINT I
   *
   *       NTUPLES  AN INTEGER MATRIX WHICH CONTAINS THE SCALARS
   *                NECESSARY TO CONSTRUCT POINT I FROM THE N EDGES
   *                ADJACENT TO THE BEST VERTEX IN THE ORIGINAL SIMPLEX;
   *                THE N-TUPLES ARE STORED BY COLUMN
   *
   *    NOTE THAT FOR ALL THREE STRUCTURES WE ALSO KEEP INFORMATION
   *    ABOUT * THE N+1 VERTICES IN THE ORIGINAL SIMPLEX IN
   *    POSITIONS/COLUMNS -N * THROUGH 0 SO THAT WE HAVE ALL THE
   *    INFORMATION WE NEED BOTH TO * GENERATE THE NEW POINTS AND TO
   *    DETECT DUPLICATION.
   *
   * INITIALIZE THE OPTIMIZATION PROCEDURE BY DETERMINING THE POINTS
   * TO BE COMPUTED AT EACH ITERATION ON THE HYPERCUBE.
   *
   * THERE IS AN UPPER LIMIT ON THE NUMBER OF POINTS WE CAN CONSIDER
   * WHICH IS BASED ON THE DIMENSION OF THE PROBLEM `N' AND THE SIZE
   * OF `SCHEME'.  CALCULATE THAT LIMIT.
   *
   * WE WILL ALSO RESCALE BY THE "GREATEST COMMON DENOMINATOR"
   * (DETERMINED BY `BETA').  IN ESSENCE, WE ARE TRYING TO GUESS HOW
   * MANY "LOOK-AHEADS" WE MUST DO TO FILL OUT THE "RAW" LIST OF
   * POINTS WE WILL GENERATE.  THE SIZE OF THE RAW LIST IS DETERMINED
   * BY `LIMIT' (I.E., THE MAXIMUM SIZE OF `SCHEME').  ONCE WE RESCALE
   * ACCORDINGLY, WE CAN WORK STRICTLY IN INTEGER ARITHMETIC.
   *
   *******************************************************************/

  /* System generated locals */

  int scheme_dim1, scheme_offset, index_offset, i_1, i_2;

  /* Local variables */

  int beta;

  static int leaf;
  static int c, i, j, k;

  static int limit, total;

  /* Parameter adjustments */

  scheme_dim1 = ndim + 2;
  scheme_offset = scheme_dim1 * (-(ndim)) - 1;
  scheme -= scheme_offset;
  index_offset = -(ndim);
  index -= index_offset;
  --list;

  /* Function Body */

  limit = (*max - ndim * ndim - ndim * 3 - 2) / (ndim + 2);
  *factor = depth(ndim, 2, limit);

  /* INITIALIZE THE ROOT OF THE TREE, WHICH IS THE CURRENT BEST
   * VERTEX. */

  leaf = 0;
  j = 0;
  index[j] = 0;
  scheme[j * scheme_dim1 - 1] = *factor;
  scheme[j * scheme_dim1] = 0;
  i_1 = ndim;

  for (k = 1; k <= i_1; ++k) {
    scheme[k + leaf * scheme_dim1] = 0;
  }

  /* INITIALIZE THE REMAINING VERTICES IN THE SIMPLEX SO THAT THEY ARE
   * NOT LATER DUPLICATED. */

  i_1 = ndim;

  for (i = 1; i <= i_1; ++i) {
    index[-i] = -i;
    scheme[-i * scheme_dim1 - 1] = *factor;
    scheme[-i * scheme_dim1] = 0;
    i_2 = ndim;

    for (k = 1; k <= i_2; ++k) {
      scheme[k + -i * scheme_dim1] = 0;
    }

    scheme[i + -i * scheme_dim1] = *factor;
  }

  /* CHECK TO PREVENT OVERFLOW (I.E., SO THAT `LIMIT' IS NOT
   * EXCEEDED). */

L4000:

  if (j <= limit - (ndim * 3 + 1)) {
    c = scheme[leaf * scheme_dim1];

    /* FIRST COMPUTE ALL N REFLECTION POINTS. */

    i_1 = c - 1;

    for (i = 0; i <= i_1; ++i) {
      ++j;
      index[j] = j;
      scheme[j * scheme_dim1 - 1] = -scheme[leaf * scheme_dim1 - 1];
      scheme[j * scheme_dim1] = i;
      i_2 = ndim;

      for (k = 1; k <= i_2; ++k) {
	scheme[k + j * scheme_dim1] = scheme[k + leaf * scheme_dim1];
      }

      if (i != 0) {
	scheme[i + j * scheme_dim1] += scheme[j * scheme_dim1 - 1];
      }

      scheme[c + j * scheme_dim1] -= scheme[j * scheme_dim1 - 1];
    }

    i_1 = ndim;

    for (i = c + 1; i <= i_1; ++i) {
      ++j;
      index[j] = j;
      scheme[j * scheme_dim1 - 1] = -scheme[leaf * scheme_dim1 - 1];
      scheme[j * scheme_dim1] = i;
      i_2 = ndim;

      for (k = 1; k <= i_2; ++k) {
	scheme[k + j * scheme_dim1] = scheme[k + leaf * scheme_dim1];
      }

      scheme[i + j * scheme_dim1] += scheme[j * scheme_dim1 - 1];

      if (c != 0) {
	scheme[c + j * scheme_dim1] -= scheme[j * scheme_dim1 - 1];
      }

    }

    /* NEXT COMPUTE ALL N + 1 CONTRACTION POINTS.  WE INCLUDE THE
     * POSSIBILITY THAT THE BEST IS NOT REPLACED. */

    i_1 = ndim;

    for (i = 0; i <= i_1; ++i) {
      ++j;
      index[j] = j;
      scheme[j * scheme_dim1 - 1] = scheme[leaf * scheme_dim1 - 1] / 2;
      scheme[j * scheme_dim1] = i;
      i_2 = ndim;

      for (k = 1; k <= i_2; ++k) {
	scheme[k + j * scheme_dim1] = scheme[k + leaf * scheme_dim1];
      }

      if (i != 0) {
	scheme[i + j * scheme_dim1] += scheme[j * scheme_dim1 - 1];
      }

      if (c != 0) {
	scheme[c + j * scheme_dim1] -= scheme[j * scheme_dim1 - 1];
      }

    }

    /* FINALLY, COMPUTE ALL N EXPANSION POINTS. */

    i_1 = c - 1;

    for (i = 0; i <= i_1; ++i) {
      ++j;
      index[j] = j;
      scheme[j * scheme_dim1 - 1] = scheme[leaf * scheme_dim1 - 1] * -2;

      scheme[j * scheme_dim1] = i;
      i_2 = ndim;

      for (k = 1; k <= i_2; ++k) {
	scheme[k + j * scheme_dim1] = scheme[k + leaf * scheme_dim1];
      }

      if (i != 0) {
	scheme[i + j * scheme_dim1] += scheme[j * scheme_dim1 - 1];
      }

      scheme[c + j * scheme_dim1] -= scheme[j * scheme_dim1 - 1];
    }

    i_1 = ndim;

    for (i = c + 1; i <= i_1; ++i) {
      ++j;
      index[j] = j;
      scheme[j * scheme_dim1 - 1] = scheme[leaf * scheme_dim1 - 1] * -2;

      scheme[j * scheme_dim1] = i;
      i_2 = ndim;

      for (k = 1; k <= i_2; ++k) {
	scheme[k + j * scheme_dim1] = scheme[k + leaf * scheme_dim1];
      }

      scheme[i + j * scheme_dim1] += scheme[j * scheme_dim1 - 1];

      if (c != 0) {
	scheme[c + j * scheme_dim1] -= scheme[j * scheme_dim1 - 1];
      }

    }

    ++leaf;
    goto L4000;
  }

  total = j;

  /* NOW, ELIMINATE ALL THE DUPLICATE POINTS.  WHEN THERE ARE
   * DUPLICATES, KEEP THE POINT THAT WAS GENERATED FIRST. */

  *error = 0;

  /* FIRST SORT THE LIST OF "RAW" POINTS; THE SORT IS KEYED ON THE
   * N-TUPLES IN `SCHEME'. */

  sort(ndim, &total, &index[index_offset], &scheme[scheme_offset], error);

  if (*error == 0) {

    /* NOW GO THROUGH THIS SORTED LIST AND BUILD A NEW LIST OF UNIQUE
     * N-TUPLES (THE INDEX FOR THIS LIST IS STORED IN THE ARRAY
     * `LIST'). */

    *unique = 0;
    i = -(ndim);
    i_1 = total;

    for (j = -(ndim) + 1; j <= i_1; ++j) {

      if (order(ndim, &scheme[index[i] * scheme_dim1 + 1],
		&scheme[index[j] * scheme_dim1 + 1]) != 0) {

	/* YOU HAVE ANOTHER (UNIQUE) POINT TO ADD TO THE LIST.  (MAKE
	 * SURE YOU ARE NOT ADDING ONE OF THE ORIGINAL VERTICES!) */

	if (index[i] > 0) {
	  ++(*unique);
	  list[*unique] = index[i];
	}

	i = j;

      } else {

	/* YOU HAVE A DUPLICATE POINT. */

	if (index[i] > index[j]) {

	  /* KEEP THE POINT WHICH WAS GENERATED "FIRST." */

	  i = j;
	}
      }
    }

    /* ADD THE LAST POINT TO THE LIST. */

    if (index[i] > 0) {

      /* AGAIN, MAKE SURE YOU ARE NOT ADDING ONE OF THE ORIGINAL
       * VERTICES! */

      ++(*unique);
      list[*unique] = index[i];
    }

    /* NOW RESORT THE LIST OF UNIQUE POINTS, KEYED ON `LIST', SO THAT
     * THE POINTS CAN BE LOADED IN THE CORRECT ORDER. */

    quick(*unique, &list[1], error);

    if (*error == 0) {

      /* WRITE ALL THIS INFORMATION OUT TO A FILE FOR LATER USE. */

      beta = (int) 2.0;
      *error = writes(out, ndim, *unique, *factor, beta, 
		     &scheme[scheme_offset], &list[1]);
    }
  }

  return *error;
}

int depth(int ndim, int beta, int max)
{
  /*******************************************************************
   *
   * THIS A SERVICE FUNCTION USED TO FIGURE OUT HOW MANY LOOK-AHEADS
   * WILL BE NEEDED TO GENERATE A LIST WITH A TOTAL OF `MAX' POINTS.
   * USING THE CONTRACTION FACTOR, WE CALCULATE A RESCALING FACTOR SO
   * THAT ALL WORK TO GENERATE THIS LIST OF POINTS CAN BE DONE IN
   * INTEGER ARITHMETIC.  THE RETURN VALUE DEPENDS ON THE CONTRACTION
   * FACTOR AND THE MAXIMUM NUMBER OF POINTS THAT CAN BE GENERATED
   * GIVEN THE AMOUNT OF WORKSPACE ALLOCATED IN THE MAIN PROGRAM.
   *
   * WRITTEN BY VIRGINIA TORCZON
   *
   * LAST MODIFICATION:  MARCH 11, 1992.
   *
   * INPUT
   *
   *    N             THE DIMENSION OF THE PROBLEM TO BE SOLVED
   *
   *    BETA          THE CONTRACTION FACTOR USED TO REDUCE THE
   *                  SIZE OF THE SIMPLEX
   *
   *    MAX           UPPER LIMIT ON THE NUMBER OF POINTS THAT CAN
   *                  BE CONSIDERED FOR THE TEMPLATE
   *
   *******************************************************************/

  /* System generated locals */

  int ret_val;

  /* Local variables */

  static int temp, factor, growth, sum;

  factor = beta;
  growth = ndim * 3 + 1;
  temp = growth;
  sum = temp;

L1000:
  if (sum < max) {
    temp *= growth;
    sum += temp;
    factor *= beta;
    goto L1000;
  }

  ret_val = factor;

  return ret_val;
}

int order(int ndim, int *x, int *y)
{
  /*******************************************************************
   *
   * THIS IS A SERVICE FUNCTION USED TO COMPARE TWO INTEGER N-TUPLES
   * TO ALLOW FOR A PARTIAL ORDERING.  THE RETURN VALUES ARE:
   *      -1   IF X < Y
   *       0   IF X = Y
   *       1   IF X > Y
   * NOTE THAT THE PARTIAL ORDERING OF THE N-TUPLES IS DETERMINED BY
   * SCANNING THE TUPLES FROM THE FIRST ELEMENT TO THE NTH ELEMENT.
   *
   * WRITTEN BY VIRGINIA TORCZON
   *
   * LAST MODIFICATION:  MARCH 11, 1992.
   *
   * INPUT
   *
   *    N             DIMENSION OF THE PROBLEM TO BE SOLVED
   *
   *    X             FIRST N-TUPLE FOR COMPARISON
   *
   *    Y             SECOND N-TUPLE FOR COMPARISON
   *
   *******************************************************************/

  /* System generated locals */

  int ret_val;

  /* Local variables */

  static int same, i;

  /* Parameter adjustments */

  --x;
  --y;

  same = 0;
  i = 1;

L1000:
  if (same == 0) {

    if (x[i] < y[i]) {
      same = -1;
    } else if (x[i] > y[i]) {
      same = 1;
    }

    if (i < ndim) {
      ++i;
      goto L1000;
    }

  }

  ret_val = same;

  return ret_val;
}

