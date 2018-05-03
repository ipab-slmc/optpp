// Test program for Simple Bound objects
//
// Modified by P.J. Williams
// October 28, 1999

#include <iostream>
#include <fstream>

#include "BoundConstraint.h"

using NEWMAT::ColumnVector;
using std::cout;
using namespace OPTPP;

void PrintConstr(BoundConstraint& bc, ColumnVector& x);

int main ()
{
  int i, num_constr = 2;
  ColumnVector xc(num_constr), lower(num_constr), upper(num_constr);
  BoolVector   std_form(num_constr), hard(num_constr);

//----------------------------------------------------------------------------
//
//  Case x <= b 
//
//----------------------------------------------------------------------------

  lower = 20;
  upper = 50;
  for(i = 1; i<= num_constr; i++){
    std_form(i) = false;
    hard(i) = false;
  }

  xc(1) = 48; xc(2) = 25;

//  Declare the bound constraint
  BoundConstraint bc(num_constr,upper,std_form);
  cout << "***********  Upper Bounds           *********************\n";
  PrintConstr(bc,xc);

  BoundConstraint both(num_constr,lower,upper);
  cout << "***********  Lower and Upper Bounds *********************\n";
  PrintConstr(both,xc);

//----------------------------------------------------------------------------
//
//  Case x >= a 
//
//----------------------------------------------------------------------------
 xc    = 2.0;
 lower = 3.0;

//  Declare the bound constraint
  BoundConstraint bc2(num_constr,lower);
  cout << "***********  Lower Bounds           *********************\n";
  PrintConstr(bc2,xc);

}
//----------------------------------------------------------------------------
//
//----------------------------------------------------------------------------
void PrintConstr(BoundConstraint& bc, ColumnVector& x) {

// Retrieve info
  bool feasible;
  double eps;
  int i, index, num_cons;
  ColumnVector low, up, residual;
  OptppArray<int> indices;

  eps       = 1.0e-08;
  num_cons  = bc.getNumOfCons();
  low       = bc.getLower();
  up        = bc.getUpper();
  residual  = bc.evalResidual(x);
  feasible  = bc.amIFeasible(x, eps);
  indices   = bc.getConstraintMappingIndices();

  cout << "X \t Lower \t Upper \t Residual\n";
  for (i=1; i<=num_cons; i++) {
       index = indices[i-1];
       cout << x(index) << "\t" << low(index) << "\t" 
            << up(index) << "\t" << residual(i) <<"\n" ;
  }

  if(feasible)
    cout << "Constraints are feasible." <<  "\n";
  else
    cout << "Constraints are infeasible." <<  "\n";

  cout << "\n";

}
