//------------------------------------------------------------------------
// Generating Set Standard 2*n Basis Class D = [I -I]
//--------------------------------------------------
// Basis vectors represented by integers:
//------------------------------------------------------------------------

/*------------------------------------------------------------------------
 Copyright (c) 2003,
 Ricardo Oliva (raoliva@lbl.gov)
 Lawrence Berkeley National Laboratory
 ------------------------------------------------------------------------*/

#include "GenSetStd.h" 

using NEWMAT::ColumnVector;

namespace OPTPP {
///> Stores the search direction in the vector y
void GenSetStd::generate(int i, double a, ColumnVector &x, ColumnVector &y)
{
  //  sets y = x + a * d[i] 
  if (i<1 || i>Size) {
//    cerr << classnm() << "of size " << Size << ". Basis index out of range: " << i << "\n";
    return;
  }

  y = x;

  if (i<=Vdim)
    y(i) += a;  
  else
    y(i-Vdim) -= a;
}

//--
// the pruning methods
//--

int GenSetStd::init(ColumnVector& gX)  {

  ActiveIDs.ReSize(Size);
  for (int i=1; i<=Size; i++) ActiveIDs(i) = i; 

  return update(gX);
}

int GenSetStd::update(ColumnVector& gX)  {
  if (Size<1) {
//    cerr << "GenSetStd1 Error: update() called on an empty GenSet\n";
    return -1;
  }

  //--
  // Update == Pruning
  //--
  // determine which search directions are descending
  // and sets Active accordingly
  //--
  int nIna = 0;
  nAct = 0; 
  ActiveIDs = 0; 
  InactiveIDs = 0;
  double gradangle = 0.0;

  // {I} ==> gX*d = gX(i)
  for (int i=1; i<=Vdim; i++) {
    if (gX(i) <= gradangle) {
      ActiveIDs(++nAct) = i;
    } 
    else {
      InactiveIDs(++nIna) = i;
    }
  }

  // { -I }
  for (int i=Vdim+1; i<=Size; i++) {
    if (gX(i-Vdim) >= gradangle) {
      ActiveIDs (++nAct) = i;
    }
    else {
      InactiveIDs(++nIna) = i;
    }
  }

  return 0;
}

} // namespace OPTPP
