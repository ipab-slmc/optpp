//------------------------------------------------------------------------
// Copyright (C) 1993: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//------------------------------------------------------------------------

#include "NLP2.h"
#include "TOLS.h"
#include "cblas.h"
#include "ioformat.h"

using namespace std;

namespace OPTPP {

//-------------------------------------------------------------------------
// Output Routines
//-------------------------------------------------------------------------

void NLP2::printState(char * s) 
{ // Print out current state: x current, gradient and Function value
  cout << "\n\n=========  " << s << "  ===========\n\n";
  cout << "\n    i\t    xc \t\t grad  \t\t fcn_accrcy \n";
  for (int i=1; i<=dim; i++) 
    cout << d(i,6) << e(mem_xc(i),12,4)<< "\t" << e(mem_grad(i),12,4) 
         << "\t" << e(mem_fcn_accrcy(i),12,4) << "\n";
  cout <<"Function Value     = " << e(fvalue,12,4) << "\n";
  double gnorm = Norm2(mem_grad);
  cout <<"Norm of gradient   = " << e(gnorm,12,4) << "\n";
  cout <<"\n\n==============================================\n\n";
}

void NLP2::fPrintState(ostream *nlpout, char * s) 
{ // Print out current state: x current, gradient and Function value
  (*nlpout) << "\n\n=========  " << s << "  ===========\n\n";
  (*nlpout) << "\n    i\t    xc \t\t grad \t\t fcn_accrcy \n";
  for (int i=1; i<=dim; i++) 
    (*nlpout)<< d(i,6) << e(mem_xc(i),12,4) << "\t" << e(mem_grad(i),12,4) 
             << "\t" << e(mem_fcn_accrcy(i),12,4) << "\n";
  (*nlpout) <<"Function Value     = " << e(fvalue,12,4) << "\n";
  double gnorm = Norm2(mem_grad);
  (*nlpout) <<"Norm of gradient   = " << e(gnorm,12,4) << "\n";
 // (*nlpout) <<"Function Accuracy  = " << e(mem_fcn_accrcy,12,4) << "\n";
  (*nlpout) <<"\n\n==============================================\n\n";
}

} // namespace OPTPP
