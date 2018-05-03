// Test program for Compound Constraints objects
//
// Modified by P.J. Williams
// November 11, 1999

#include <iostream>
#include <fstream>

#include "OptppArray.h"
#include "Constraint.h"
#include "BoundConstraint.h"
#include "LinearEquation.h"
#include "LinearInequality.h"
#include "CompoundConstraint.h"

using NEWMAT::ColumnVector;
using NEWMAT::Matrix;

using std::cout;

using namespace OPTPP;


void PrintConstr(Constraint& constraint_, ColumnVector& x);
void PrintCompConstr(CompoundConstraint& constraint_, ColumnVector& x);

int main ()
{

  int bc_num_constr  = 2;
  ColumnVector xc(bc_num_constr), lower(bc_num_constr);

  int num_constr = 3;
  int num_var    = 2;
  Matrix       A(num_constr, num_var), B(num_constr, num_var);
  ColumnVector rhs(num_constr), b(num_constr);
  ColumnVector Ax(num_constr);


//----------------------------------------------------------------------------
//
//  Case x >= 3
//
//----------------------------------------------------------------------------

 for(int i = 1; i <= num_var; i++){
     xc(i)  = 1.0*i;
 }
 lower = 3;

 Constraint bc =  new BoundConstraint(bc_num_constr,lower);        
 cout << "*****  Bound Constraint  *****\n";
 PrintConstr(bc,xc);

//----------------------------------------------------------------------------
//
//  Case Ax = 0
//
//----------------------------------------------------------------------------

 for(int i = 1; i <= num_var; i++){
     A(i,i) = 1.0;
 }
 A(1,2) = 0.0;
 A(2,1) = 0.0;
 A(3,1) = 1.0;
 A(3,2) = 1.0;

//  Declare the object

 b = 0.0;
 Constraint leqn = new LinearEquation(A,b);
 cout << "*****  Linear Equation, b = 0  *****\n";
 PrintConstr(leqn,xc);

//----------------------------------------------------------------------------
//
//  Case Ax <= b
//
//----------------------------------------------------------------------------

 b  = 2.0;

 Constraint ineq = new LinearInequality(A,b);
 cout << "*****  Linear Inequality,  *****\n";
 PrintConstr(ineq,xc);
 
//----------------------------------------------------------------------------
//
//  Compound Constraints
//
//---------------------------------------------------------------------------
 OptppArray<Constraint> ac(0);
 ac.append(bc);
 ac.append(leqn);
 ac.append(ineq);
 cout << "*****  Compound Constraints,  *****\n";
 CompoundConstraint constraints(ac);
 PrintCompConstr(constraints,xc);

}

void PrintConstr(Constraint& eqn, ColumnVector& x) {

// Retrieve info
double eps;
bool feasibility;
int i, index, num_cons;
ColumnVector lower, upper, residual;
OptppArray<int> mapping;

eps         = 1.0e-08;
num_cons    = eqn.getNumOfCons();
residual    = eqn.evalResidual(x);
lower       = eqn.getLower();
upper       = eqn.getUpper();
feasibility = eqn.amIFeasible(x,eps);
mapping     = eqn.getConstraintMappingIndices();	 

cout << "Index  \t Resid. \t Lower \t Upper \n";
for (i=1; i<=num_cons; i++){ 
  index = mapping[i-1];
  cout << index <<  "\t" << residual(i) << "\t" 
       << lower(index) <<"\t" <<  upper(index) <<"\n";
}

  if(feasibility)
    cout << "Constraints are feasible." <<  "\n";
  else
    cout << "Constraints are infeasible." <<  "\n";
}

void PrintCompConstr(CompoundConstraint& eqn, ColumnVector& x) {

// Retrieve info
double eps = 1.0e-08;
int num_sets             = eqn.getNumOfSets();
int num_constr           = eqn.getNumOfCons();
ColumnVector type        = eqn.getConstraintType();
ColumnVector residual    = eqn.evalResidual(x);
ColumnVector lower       = eqn.getLower();
ColumnVector upper       = eqn.getUpper();
bool         feasibility = eqn.amIFeasible(x,eps);

cout << "\nThere are " << num_sets << " sets of equations. \n";
cout << "Constraints are re-ordered so that equations "
     << "appear before inequalities. \n\n";
cout << "Index \t Resid. \t Lower \t Upper \n";
for (int i=1; i<= num_constr; i++) 
  cout << i << "\t" << residual(i) << "\t" 
       << lower(i) <<"\t" <<  upper(i) <<"\n";

  if(feasibility)
    cout << "Constraints are feasible." <<  "\n";
  else
    cout << "Constraints are infeasible." <<  "\n";

}

