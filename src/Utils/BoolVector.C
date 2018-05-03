//------------------------------------------------------------------------
// P.J. Williams
// Sandia National Laboratories
// pwillia@sandia.gov
// Last modified 11/05/1999
//------------------------------------------------------------------------

#include "BoolVector.h"

using namespace std;

namespace OPTPP {

// Constructors
BoolVector::BoolVector(int sz):
    size(sz), p(0)
    {p = new bool[size = sz];}

BoolVector::BoolVector(int sz, const bool& val) {
     p = new bool[size = sz];
     for (int i=0; i<size; i++) p[i] = val;
  }

BoolVector::BoolVector(int sz, const BoolVector& val) {
     p = new bool[size = sz];
     for (int i=0; i<size; i++) p[i] = val.p[i];
  }

BoolVector::BoolVector(const BoolVector & val) {  
     p = new bool[size=val.size];
     for (int i=0; i<size; i++) p[i] = val.p[i];
  }

BoolVector& BoolVector::operator=(const BoolVector& val) 
{
    if (this != &val) {
      p = new bool[size=val.size];
      for (int i=0; i< val.size; i++) p[i] = val.p[i];
    }
    return *this;
}
//PJW


bool& BoolVector::operator()( int index ) 
{
    if (!(1 <= index && index <= size)) {
//      cerr << "BoolVector: out of bounds\n";
      exit(1);
    }
    return p[index-1];
}

const bool& BoolVector::operator()( int index ) const 
{
    if (!(1 <= index && index <= size)) {
//      cerr << "BoolVector: out of bounds\n";
      exit(1);
    }
    return p[index-1];
}

} // namespace OPTPP
