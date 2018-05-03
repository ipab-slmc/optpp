/*********************************************************************
*   File:  dogleg.c
*
*   The source code in this file computes a dogleg step using
*   Powell's method.
*
*   Code development:
*      12 Oct 94 - Originated by T. Plantenga, Sandia National Labs.
*       4 Nov 94 - Converted to C++ by J.C. Meza, Sandia National Labs.
**********************************************************************/
#include "Opt.h"
#include "ioformat.h"

using NEWMAT::ColumnVector;
using NEWMAT::SymmetricMatrix;
using NEWMAT::DiagonalMatrix;

namespace OPTPP {

int dogleg(NLP1* nlp, ostream *fout, 
	   SymmetricMatrix& Hessian, ColumnVector& grad, 
	   ColumnVector& dogleg_step, ColumnVector& sx, 
	   real& dnorm, real& TR_size, real step_max)       
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *   This routine solves the trust region subproblem by Powell's
 *   dogleg method.  It constructs the Cauchy and Newton steps
 *      d_cp = -alpha grad
 *      d_nw = -H grad
 *   using L-BFGS Hessian approximations.  The dogleg path runs from
 *   dogleg_step = 0 to d_cp to d_nw, and its intersection with the
 *   spherical trust region constraint is the dogleg solution point.
 *   The Newton step should always be further away than the Cauchy
 *   point, but it may lie inside the trust region.
 *
 *   Parameters
 *     nlp          -->  pointer to nonlinear problem object 
 *
 *     Hessian      -->  Hessian Matrix
 *
 *     grad         -->  Vector of length n which specifies the 
 *                       gradient on input
 *     sx           -->  Scaling vector for x
 *                       gradient on input
 *     dogleg_step <-->  Vector of length n which specifies the 
 *                       newton direction on input. On output it will
 *                       contain the dogleg step 
 *     dnorm       <--   Norm of the dogleg step
 *     TR_size      -->  Trust region size
 *     step_max     -->  Maximum step size allowed
 *     debug        -->  Debug flag
 *     outfile      -->  File pointer for output messages
 *
*********************************************************************/
{
// Local variables
  int   i;
  real  gBg, alpha;
  real  norm_cp, norm_nw, gnorm;
  real  a, b, c, tau;
  int n = nlp->getDim();
  ColumnVector tmp_vec(n);
  ColumnVector scaled_grad(n);
  real CP_length;
  DiagonalMatrix Dx(n);
  int trace = 0;

// Return the Newton step if it's inside the trust region. 

  norm_nw = Norm2(dogleg_step);

  if (norm_nw <= TR_size) {
    dnorm = norm_nw;
    return(Newton_Step);
  }

  if (trace) *fout << "\nNewton step outside of trust region\n";
//
// Compute Scaled gradient
//
  for (i=1; i<=n; i++) Dx(i,i) = sx(i);
  scaled_grad = Dx*grad;

// Compute the Cauchy point in tmp_vec.  If it lies outside
// the trust region, then use take a step to the boundary

  tmp_vec = Dx*Hessian*Dx*scaled_grad;
  gBg = Dot(scaled_grad, tmp_vec);
  gnorm = Norm2(scaled_grad);
  CP_length = (gnorm * gnorm * gnorm) / gBg;

// If TR_size has not been set
// set it equal to the minimum of CP length and max step allowed

  if (TR_size == 0.0) TR_size = min(CP_length,step_max);

// Check to make sure that the CP is inside the trust region
// If it's outside, then just take the scaled gradient to the TR boundary

  if (trace) 
    *fout << "Compute Cauchy Point\n" 
         << "gBg       = " << e(gBg,12,4) << "\n"
         << "gnorm     = " << e(gnorm,12,4) << "\n"
         << "CP_length = " << e(CP_length,12,4) << "\n";

  if (CP_length >= TR_size) {
    alpha = -TR_size / gnorm;
    dogleg_step = scaled_grad*alpha;
    dnorm = Norm2(dogleg_step);
    if (trace) 
      *fout  <<"Cauchy point outside trust region\n"
            <<"alpha = " << e(alpha,12,4) << "\n"
            <<"dnorm = " << e(dnorm,12,4) << "\n";
    return(Cauchy_Step);
  } 

// O.K.,  the CP is inside the trust region
//
  alpha = -(gnorm * gnorm / gBg);
  tmp_vec = scaled_grad * alpha;
  norm_cp = Norm2(tmp_vec);

// The final case is to return the intersection of the segment
// connecting the Cauchy point to the Newton point with the
// trust region boundary; i.e., find tau such that
//     || d_cp + tau*(d_nw - d_cp) ||_2 = TR_size.
//  This requires solving a single quadratic equation.

  dogleg_step = dogleg_step - tmp_vec;
  a = Dot(dogleg_step, dogleg_step);
  b = 2.0 * Dot(dogleg_step, tmp_vec);
  c = (TR_size * TR_size) - (norm_cp * norm_cp);
  tau = (-b + sqrt(b*b + 4.0*a*c)) / (2.0*a);
  dogleg_step = dogleg_step*tau;
  dogleg_step = tmp_vec + dogleg_step;
  dnorm = Norm2(dogleg_step);

  if (trace) *fout << "Taking a dogleg step\n"
                  << "dnorm = " << e(dnorm,12,4) << "\n";
  return(Dogleg_Step);
}

} // namespace OPTPP

