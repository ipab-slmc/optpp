//------------------------------------------------------------------------
// Generating Set Basis Class D = [I -I (1;1) (1;-1) (-1;1) (-1;-1)]
// for 2-d only!!!!
//-----------------------------------------------------------------------
// Basis vectors represented by integers
//------------------------------------------------------------------------

/*------------------------------------------------------------------------
 Copyright (c) 2003,
 Ricardo Oliva (raoliva@lbl.gov)
 Lawrence Berkeley National Laboratory
 ------------------------------------------------------------------------*/

#include "GenSetBox2d.h" 

using NEWMAT::ColumnVector;
using std::cerr;

namespace OPTPP {

///> Stores the search direction in the vector y
void GenSetBox2d::generate(int i, double a, ColumnVector &x, ColumnVector &y)
{
  //  sets y = x + a * d[i] 
  if (i<1 || i>Size) {
//    cerr << "Gen_Set_Box2d: Basis index out of range: " << i << "\n";
    return;
  }

  y << x;

  if (i<=Vdim)
    y(i) += a;  
  else if (i<=2*Vdim)
    y(i-Vdim) -= a;
  else {
    double w = a / sqrt(2.0);

    switch (i-2*Vdim) {
    case 1:
      y(1) += w;   y(2) += w;  break;
    case 2:
      y(1) += w;   y(2) -= w;  break;
    case 3:
      y(1) -= w;   y(2) += w;  break;
    case 4:
      y(1) -= w;   y(2) -= w;  break;
    }
  }
}

//--
// the pruning methods
//--

int GenSetBox2d::init(ColumnVector& gX)  {

  ActiveIDs.ReSize(Size);
  for (int i=1; i<=Size; i++) ActiveIDs(i) = i; 

  return update(gX);
}

int GenSetBox2d::update(ColumnVector& gX)  {
  if (Size<1) {
//    cerr << "GenSetBox2d Error: update() called on an empty GenSet\n";
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
  for (int i=Vdim+1; i<=2*Vdim; i++) {
    if (gX(i-Vdim) >= gradangle) {
      ActiveIDs (++nAct) = i;
    }
    else {
      InactiveIDs(++nIna) = i;
    }
  }

  // { Corner Directions }
  for (int i=2*Vdim+1; i<=Size; i++) {
    int s = i-2*Vdim;
    double dot;
    switch (s) { 
    case 1:  dot =  gX(1) + gX(2);  break;
    case 2:  dot =  gX(1) - gX(2);  break;
    case 3:  dot = -gX(1) + gX(2);  break;
    case 4:  dot = -gX(1) - gX(2);  break;
    default: dot = 0.0;
    }
    if (gradangle != 0.0) dot /= sqrt(2.0);

    if (dot < gradangle) {
      ActiveIDs (++nAct) = i;
    }
    else 
      InactiveIDs(++nIna) = i;

  } // for

  return 0;
}

} // namespace OPTPP
