
//--------------------------------------------------------------------
// Copyright (C) 1993,1994:
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//--------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#include <iostream>
#include <fstream>
#ifdef HAVE_STD
#include <cmath>
#else
#include <math.h>
#endif

#include "newmat.h"

#include "pds.h"

using namespace std;
using NEWMAT::ColumnVector;

extern "C" {
int bin_open(char *filename, int *fd);
int bin_close(int fd);
}

namespace OPTPP {
int create_scheme (ostream *fout, int ndim, int scheme_limit, char
		   *scheme_name, int *scheme, int debug)
{
  /*******************************************************************
   *
   * create a search strategy for the parallel direct search method
   *
   * since the size of the workspace plays a critical role in the
   * total number of points generated for the search scheme, a brief
   * explanation of how this workspace is used is included here for
   * those who wish to optimize the use of memory.
   *
   * upon entry into the subroutine `search' the vector `scheme' is
   * partitioned into a matrix of the following form
   *                integer         scheme(-1:n,-n:?) 
   * where the total number of columns depends on the amount of space
   * allocated in the calling program.
   *
   * the vector `index' is used as a permutation array to keep track
   * of every column of `scheme'.  in `search' it is an array of the
   * form
   *           integer         index(-n:?) 
   *
   * the vector `list' is used to keep track of each unique n-tuple in
   * `scheme'.  and thus is declared to be an array of the form
   *           integer         list(?) 
   *
   * thus, `index' and `list' really only need to be large enough to
   * track the total number of columns in scheme but are declared to
   * be as large as `scheme' to prevent any possible overflow.  for
   * the most efficient use of space---which may become an issue when
   * generating very large search schemes on a processor with a
   * limited amount of memory, the constant `dim' can be set equal to
   * the dimension of the problem(s) for which the search scheme is
   * being generated and the constant `max' can be set equal to the
   * number of columns of workspace to be allowed in `scheme' so that
   * the actual amount of space for `scheme' becomes
   *           integer         scheme(-1:dim,-dim:max) 
   *
   * note that the constant `limit' automatically takes care of this
   * in the driver.  the workspace for index and list can then be
   * redefined---in the driver--as
   *           integer         index(1+dim+max), list(max) 
   * without any danger of overflow.  thus all space created in the
   * calling program will be used in the subroutine `make_search'.
   *  
   * Original version due to Virginia Torczon
   * This version hacked by J.C. Meza
   *
   *******************************************************************/

  int error;
  int factor, unique;
  int fpscheme;
  ColumnVector list(scheme_limit), index(scheme_limit);

  (*fout) << "Creating SCHEME file: " << scheme_name << "\n";

  error = bin_open(scheme_name, &fpscheme);
  if (error != 0) {
//    cerr << "create_scheme: error opening scheme file for writing.   \n"
//         << "The TMP environment variable may need to be set to a    \n"
//         << "valid temporary file system.  Otherwise, PDS and TRPDS  \n"
//         << "will not run correctly.  Please set the TMP environment \n"
//         << "variable and re-run the problem. \n" << endl;
    return error;
  }

  make_search(ndim, fpscheme, &scheme_limit, scheme, (int *)index.Store(),
	      (int *)list.Store(), &unique, &factor, &error);
  
  if (error == 0) {
    if (debug) {
      (*fout) << "Successfully completed a search strategy.\n";
      (*fout) << "Dimension of the problem = " << ndim << "\n";
      (*fout) << "Number of unique points  = " << unique << "\n";
      (*fout) << "Restoration factor       = " << factor << "\n";
      (*fout) << "Initialization phase finished.\n\n";
    }
  }
  else {
    (*fout) << "Returned without a completed search strategy. \n";
    (*fout) << "Internal stack overflow in quicksort routines.\n";
    (*fout) << "Check the documentation for further details.\n" << endl;
    return error;
  }

  error = bin_close(fpscheme);
  return error;  
}

} // namespace OPTPP

