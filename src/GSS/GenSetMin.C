//------------------------------------------------------------------------
// Generating Set Minimum Size Basis Class D = [I -1]
//--------------------------------------------------
// Basis vectors represented by integers:
//------------------------------------------------------------------------

/*------------------------------------------------------------------------
 Copyright (c) 2003,
 Ricardo Oliva (raoliva@lbl.gov)
 Lawrence Berkeley National Laboratory
 ------------------------------------------------------------------------*/

#include "GenSetMin.h" 

using NEWMAT::ColumnVector;

namespace OPTPP {

///> Stores the search direction in the vector y
void GenSetMin::generate(int i, double a, ColumnVector &x, ColumnVector &y)
{

  //  sets y = x + a * d[i] 
  if (i<1 || i>Size) {
//    cerr << "GenSetMin: size=" << Size << "; Basis index out of range: " << i << "\n";
    return;
  }

  y << x;
  if (i<Size)
    y(i) += a; 
  else
    y -= 1.0;
   
}

//--
// the pruning methods
//--

int GenSetMin::init(ColumnVector& gX)  {

  return update(gX);
}

int GenSetMin::update(ColumnVector& gX)  {
  if (Size<1) {
//    cerr << "GenSetMin Error: update() called on an empty set\n";
    return -1;
  }

  //--
  // Update == Pruning
  //--
  // determine which search directions are descending
  // and sets Active accordingly
  //--
  int nIna = 0; 
  nAct = 0; // all inactive;
  ActiveIDs = 0; 
  InactiveIDs = 0;
  double gradangle = 0.0;

  // d_i = I(:,i) ==> -gX dot d_i = -gX(i), 
  for (int i=1; i<=Vdim; i++) {
    if (gX(i) <= gradangle) {
      ActiveIDs(++nAct) = i;
    } 
    else 
      InactiveIDs(++nIna) = i;
  }

  // { -1 } ==> -gX dot d_i = sum(gX)
  double dot = gX.Sum();
  if (dot >= gradangle) {
    ActiveIDs(++nAct) = Size;
  }
  else 
    InactiveIDs(++nIna) = Size;

  return 0;
}

} // namespace OPTPP
