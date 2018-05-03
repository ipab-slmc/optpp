
#if (defined(__sgi) || defined(__xlc__) || defined(__xlC__))
#define WANT_MATH
#else
#define WANT_STREAM
#define WANT_MATH
#endif

#include <string>
#include <iostream>
#include <fstream>

using namespace std;

#include "ioformat.h"
#include "newmat.h"

using NEWMAT::LowerTriangularMatrix;
using NEWMAT::UpperTriangularMatrix;
using NEWMAT::DiagonalMatrix;
using NEWMAT::SymmetricMatrix;
using NEWMAT::Matrix;
using NEWMAT::Real;

namespace OPTPP {

class PrintCounter
{
   int count;
   char* s;
public:
   ~PrintCounter();
   PrintCounter(char * sx) : count(0), s(sx) {}
   void operator++() { count++; }
};

PrintCounter PCZ("Number of non-zero matrices (should be 1) = ");
PrintCounter PCN("Number of matrices tested                 = ");

PrintCounter::~PrintCounter()
//jcm{ cout << s << count << "\n"; }
{}


void Print(const Matrix& X)
{
   ++PCN;
   cout << "\nPrint::Matrix type: " << X.Type().Value() << " (";
   cout << X.Nrows() << ", ";
   cout << X.Ncols() << ")\n\n";
   if (X.IsZero()) { cout << "All elements are zero\n" << flush; return; }
   int nr=X.Nrows(); int nc=X.Ncols();
   for (int i=1; i<=nr; i++)
   {
      for (int j=1; j<=nc; j++)  cout << e(X(i,j),14,6) << "\t";
      cout << "\n";
   }
   cout << flush; ++PCZ;
}

void Print(const UpperTriangularMatrix& X)
{
   ++PCN;
   cout << "\nMatrix type: " << X.Type().Value() << " (";
   cout << X.Nrows() << ", ";
   cout << X.Ncols() << ")\n\n";
   if (X.IsZero()) { cout << "All elements are zero\n" << flush; return; }
   int nr=X.Nrows(); int nc=X.Ncols();
   int i, j;
   for (i=1; i<=nr; i++)
   {
      for (j=1; j<i; j++) cout << "\t";
      for (j=i; j<=nc; j++)  cout << e(X(i,j),14,6)  << "\t";
      cout << "\n";
   }
   cout << flush; ++PCZ;
}

void Print(const DiagonalMatrix& X)
{
   ++PCN;
   cout << "\nMatrix type: " << X.Type().Value() << " (";
   cout << X.Nrows() << ", ";
   cout << X.Ncols() << ")\n\n";
   if (X.IsZero()) { cout << "All elements are zero\n" << flush; return; }
   int nr=X.Nrows(); int nc=X.Ncols();
   for (int i=1; i<=nr; i++)
   {
      for (int j=1; j<i; j++) cout << "\t";
      if (i<=nc) cout << e(X(i,i),14,6) << "\t";
      cout << "\n";
   }
   cout << flush; ++PCZ;
}

void Print(const SymmetricMatrix& X)
{
   ++PCN;
   cout << "\nMatrix type: " << X.Type().Value() << " (";
   cout << X.Nrows() << ", ";
   cout << X.Ncols() << ")\n\n";
   if (X.IsZero()) { cout << "All elements are zero\n" << flush; return; }
   int nr=X.Nrows(); int nc=X.Ncols();
   int i, j;
   for (i=1; i<=nr; i++)
   {
      for (j=1; j<i; j++) cout << e(X(j,i),14,6) << "\t";
      for (j=i; j<=nc; j++)  cout << e(X(i,j),14,6) << "\t";
      cout << "\n";
   }
   cout << flush; ++PCZ;
}

void Print(const LowerTriangularMatrix& X)
{
   ++PCN;
   cout << "\nMatrix type: " << X.Type().Value() << " (";
   cout << X.Nrows() << ", ";
   cout << X.Ncols() << ")\n\n";
   if (X.IsZero()) { cout << "All elements are zero\n" << flush; return; }
   int nr=X.Nrows();
   for (int i=1; i<=nr; i++)
   {
      for (int j=1; j<=i; j++) cout << e(X(i,j),14,6) << "\t";
      cout << "\n";
   }
   cout << flush; ++PCZ;
}


void Clean(Matrix& A, Real c)
{
   int nr = A.Nrows(); int nc = A.Ncols();
   for (int i=1; i<=nr; i++)
   {
      for ( int j=1; j<=nc; j++)
      { Real a = A(i,j); if ((a < c) && (a > -c)) A(i,j) = 0.0; }
   }
}

void Clean(DiagonalMatrix& A, Real c)
{
   int nr = A.Nrows();
   for (int i=1; i<=nr; i++)
   { Real a = A(i,i); if ((a < c) && (a > -c)) A(i,i) = 0.0; }
}

void FPrint(ostream *fout, const Matrix& X)
{
   (*fout) << "\nFPrint::Matrix type: " << X.Type().Value() << " (";
   (*fout) << X.Nrows() << ", ";
   (*fout) << X.Ncols() << ")\n\n";
   if (X.IsZero()) { (*fout) << "All elements are zero\n" << flush; return; }
   int nr=X.Nrows(); int nc=X.Ncols();
   for (int i=1; i<=nr; i++)
   {
      for (int j=1; j<=nc; j++)  (*fout) << e(X(i,j),14,6) << "\t";
      (*fout) << "\n";
   }
   (*fout) << flush; ++PCZ;
}

void FPrint(ostream *fout, const SymmetricMatrix& X)
{
   ++PCN;
   (*fout) << "\nFPrint::Matrix type: " << X.Type().Value() << " (";
   (*fout) << X.Nrows() << ", ";
   (*fout) << X.Ncols() << ")\n\n";
   if (X.IsZero()) { (*fout) << "All elements are zero\n" << flush; return; }
   int nr=X.Nrows(); int nc=X.Ncols();
   int i, j;
   for (i=1; i<=nr; i++)
   {
      for (j=1; j<i; j++) (*fout) << e(X(j,i),14,6) << "\t";
      for (j=i; j<=nc; j++)  (*fout) << e(X(i,j),14,6) << "\t";
      (*fout) << "\n";
   }
   (*fout) << flush; ++PCZ;
}
void FPrint(ostream *fout, const DiagonalMatrix& X)
{
   ++PCN;
   (*fout) << "\nMatrix type: " << X.Type().Value() << " (";
   (*fout) << X.Nrows() << ", ";
   (*fout) << X.Ncols() << ")\n\n";
   if (X.IsZero()) { (*fout) << "All elements are zero\n" << flush; return; }
   int nr=X.Nrows(); int nc=X.Ncols();
   for (int i=1; i<=nr; i++)
   {
      for (int j=1; j<i; j++) (*fout) << "\t";
      if (i<=nc) (*fout) << e(X(i,i),14,6) << "\t";
      (*fout) << "\n";
   }
   (*fout) << flush; ++PCZ;
}

} // namespace OPTPP
