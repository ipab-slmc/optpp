// Test program for Linear Equation objects
//
// Modified by P.J. Williams
// October 28, 1999

#include <iostream>
#include <fstream>

#include "LinearInequality.h"
#include "LinearEquation.h"

using NEWMAT::ColumnVector;
using NEWMAT::Matrix;
using NEWMAT::Real;

using std::cout;

using namespace OPTPP;

void PrintConstr(LinearInequality& lineq, ColumnVector& x);

int main ()
{
  int num_constr = 3, num_var    = 2;
  ColumnVector xc(num_var), b(num_constr);
  Matrix       A(num_constr, num_var); 

//----------------------------------------------------------------------------
//
//  Case Ax >= 0
//
//----------------------------------------------------------------------------
  Real a[] = {1.0, 0.0, 0.0, 1.0, 1.0, 1.0};
  A        << a;
  for(int i = 1; i<= num_var; i++)
     xc(i)  = 1.0*i;

  b = 0.0;

  //  Declare the object
  LinearInequality first_lineq(A,b);        
  cout << "*****  Ax >= b *****\n";
  PrintConstr(first_lineq,xc);

//----------------------------------------------------------------------------
//
//  Case Ax <= b
//
//----------------------------------------------------------------------------
    b(1)    = 1.0;
    b(2)    = 1.0;
    b(3)    = 2.0;

  //  Declare the object
  
  bool flag = false;
  LinearInequality second_lineq(A,b,flag);        
  cout << "*****  Ax <= b  *****\n";
  PrintConstr(second_lineq,xc);
  
}

//----------------------------------------------------------------------------
void PrintConstr(LinearInequality& lineq, ColumnVector& x) {
 
  double eps = 1.0e-08;
// Retrieve info
  int num_constr        = lineq.getNumOfCons();
  int num_var           = lineq.getNumOfVars();
  ColumnVector lower    = lineq.getLower();
  ColumnVector upper    = lineq.getUpper();
  ColumnVector residual = lineq.evalResidual(x);
  ColumnVector Ax       = lineq.evalAx(x);
  bool         feasible = lineq.amIFeasible(x,eps);

  cout << "The current point x : \n";
  for (int i = 1; i <= num_var; i++) {
    cout << i <<  "\t" <<x(i) << "\n";
  }

  cout << "Index \t  Ax \t Lower \t Upper \t Residual \n";
  for (int i = 1; i <= num_constr; i++) {
     cout << i <<  "\t" << Ax(i) << "\t" 
	  << lower(i) << "\t" << upper(i) << "\t" 
	  << residual(i) << "\n";
  }

  if(feasible)
    cout << "Constraints are feasible." <<  "\n";
  else
    cout << "Constraints are infeasible." <<  "\n";

  cout << "\n";
}

