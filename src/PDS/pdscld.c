
#include "pds.h"

int pdscld(int ndim, double scale, double *simplex)
{
  /*******************************************************************
   *
   * This is a service subroutine used to construct an initial
   * simplex when one is not provided by the user.  This routine
   * constructs a scaled right-angle simplex (i.e., one in which the
   * edges are scaled by the value of the original coordinates to give
   * the displacement along the coordinate axes).  Note that the user
   * must provide the variable SCALE, which specifies the original
   * length (and orientation) of the displacement, as well as the
   * initial vertex used to start the search (and scale the actual
   * displacements).
   *
   * Written by Virginia Torczon. 
   * Last modification:  January 6, 1994. 
   *
   * Parameters
   * 
   *    Input 
   *
   *       N             dimension of the problem to be solved
   *
   *       SCALE         dual purpose; if
   *                        0    then user supplied the simplex to be
   *                             used to start the search and this
   *                             subroutine will not be called
   *                        else this subroutine constructs a scaled
   *                             right-angle simplex with edges of
   *                             length
   *                             |SCALE*V0(I)|(sign(SCALE*V0(I))
   *                             determines orientation)
   *
   *       S             two-dimensional array containing the vertex 
   *                     specified by the user to begin the search
   * 
   *    Output
   *
   *       S             two-dimensional array containing the N+1
   *                     vertices of the simplex to be used to start
   *                     the search.  This simplex will be a scaled
   *                     right-angle simplex with displacements along
   *                     the coordinate axes with base length (and
   *                     orientation) of SCALE, further scaled (and
   *                     reoriented) by the coordinates of V0.
   *
   * Remarks
   *
   *    The workspace for the simplex is declared in the main program
   *    to be a one-dimensional vector SIMPLEX; however, all
   *    subroutine calls redefine the one-dimensional vector SIMPLEX
   *    to be a two-dimensional array S of dimension N by N+1.  Note
   *    also that the N+1 vertices of the simplex are stored by
   *    column.
   *
   *******************************************************************/

  static int i, j;

  for (j = 1; j <= ndim; j++) {

    for (i = 0; i < ndim; i++) {
      simplex[i + j * ndim] = simplex[i];
    }

    /* Make sure that you do not scale by zero! */

    if (simplex[j - 1] + 1. != 1.) {
      simplex[j - 1 + j * ndim] += scale * simplex[j - 1];
    } else {
      simplex[j - 1 + j * ndim] += scale;
    }
  }

  return 0;
}

