
#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#ifdef HAVE_STD
#include <cmath>
#include <cerrno>
#include <cstring>
#else
#include <math.h>
#include <errno.h>
#include <string.h>
#endif

#include "OptPDS.h"
#include "pds.h"
#include "newmat.h"
#include "common.h"
#include "cblas.h"

using NEWMAT::ColumnVector;

/* Structures for constraints and parallel configuration. */

extern struct pdscon pdscon;
extern struct conbcmni conbcmni;

namespace OPTPP {

int pdschk(NLP0* nlp, int ndim, double *xc, double *xt, double tr_size,
	   double *dist, int trpds, double feas_tol)
{
  /* Check trust region constraint for TRPDS. */

  int i ;

  *dist = 0.0;
  
// PJW
  if(nlp->hasConstraints()){
    CompoundConstraint* constraints = nlp->getConstraints();

    ColumnVector xtrial(ndim);
    for (i = 0; i < ndim; i++) {
      xtrial(i+1) = xt[i]; 
    }

    bool feasible = constraints->amIFeasible(xtrial, feas_tol);
    if(!feasible){
      //       printf("pdschk: Current point violates the constraints. \n");
       return feasible;
    }
  }
  if(trpds){
    ColumnVector diff(ndim);
  
// PJW

    for (i = 0; i < ndim; i++)
      diff(i+1) = xc[i] - xt[i];

    *dist = Norm2(diff);

    if(*dist < 0.0) {
      printf("pdschk: Distance is negative: %e\n", *dist);
    }

    if (*dist <= tr_size)
      return 1;
    else{
      return 0;
    }
  }
  /* The problem is unconstrained AND the solution method is not TRPDS */
  else
    return 1;
}

} // namespace OPTPP
