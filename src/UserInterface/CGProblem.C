
using namespace std;

#include "xercesc/util/XMLString.hpp"

#include"CGProblem.h"

namespace OPTPP {

DOMElement* CGProblem::FindParameterXML()
{
  DOMElement* solverXML = GetSolverXML();
  DOMElement* parameterXML;
  DOMNodeList* nodeList;
  DOMNode* tmpNode;

  if((nodeList =
      solverXML->getElementsByTagName(XMLString::transcode("CGvar")))->getLength() > 0)
    {
      tmpNode = nodeList->item(0);
      parameterXML = (DOMElement *) tmpNode;
    }

  return parameterXML;
}

OptError CGProblem::CreateFunctionOptimizer(OptimizeClass* &objfcn,
					    NLP0* &func)
{
  OptError error;

  int numVar = GetNumVar();
  int derivOrder = GetDerivOrder();
  VariableList* variables = GetVariables();

  if ((variables->upperExists()) || (variables->lowerExists()) ||
      (variables->linearExists()) || (variables->nonlinearExists()))
  {
    cout << "CG does not support constraints." << endl;
    cout << "Constraints being ignored." << endl;
  }

  // construct a function with 2nd derivative info

  if(derivOrder == 2)
    {
      USERFCN2 userFcn = GetUserFunction2();
      INITFCN initFcn = GetInitFunction();
      if(userFcn == NULL || initFcn == NULL)	
	{
	  error.value = -4;
	  error.msg = "Error loading function";
	  return error;
	}

      NLF2* myFunc = new NLF2(numVar, userFcn, initFcn);
      myFunc->setIsExpensive(true);
      myFunc->setX(variables->GetInitialValues());
      func = myFunc;
    }
  else if(derivOrder == 1)
    {
      // construct a function with 1st derivative info

      USERFCN1 userFcn = GetUserFunction1();
      INITFCN initFcn = GetInitFunction();
      if(userFcn == NULL || initFcn == NULL)	
	{
	  error.value = -4;
	  error.msg = "Error loading function";
	  return error;
	}

      NLF1* myFunc = new NLF1(numVar, userFcn, initFcn);
      myFunc->setIsExpensive(true);
      myFunc->setX(variables->GetInitialValues());
      func = myFunc;
    }
  else if(derivOrder == 0)
    {
      USERFCN0 userFcn = GetUserFunction0();
      INITFCN initFcn = GetInitFunction();
      if(userFcn == NULL || initFcn == NULL)	
	{
	  error.value = -4;
	  error.msg = "Error loading function";
	  return error;
	}

      // construct a function with no derivative info 
      // but use finite differences to approximate a first derivative

      FDNLF1* myFunc = new FDNLF1(numVar, userFcn, initFcn);
      myFunc->setIsExpensive(true);
      myFunc->setX(variables->GetInitialValues());
      func = myFunc;
    }

  objfcn = new OptCG((NLP1 *)func);
  error.value = 0;
  error.msg = "Success";
  return error;
}

OptError CGProblem::CreateApplicationOptimizer(OptimizeClass* &objfcn,
					       NLP0* &func)
{
  OptError error;

  int numVar = GetNumVar();
  VariableList* variables = GetVariables();
  AppLauncher* launcher = GetAppLauncher();
  USERFCN0APP userFcn = &(launcher->run_app);
  INITFCNAPP initFcn = &(launcher->init_app);

  if ((variables->upperExists()) || (variables->lowerExists()) ||
      (variables->linearExists()) || (variables->nonlinearExists()))
  {
    cout << "CG does not support constraints." << endl;
    cout << "Constraints being ignored" << endl;
  }

  FDNLF1APP * myFunc = new FDNLF1APP(numVar, userFcn, initFcn,
				     launcher);
  myFunc->setIsExpensive(true);
  func = myFunc;

  objfcn = new OptCG(myFunc);

  error.value = 0;
  error.msg = "Success";
  return error;
}	


void CGProblem::SetParameters(OptimizeClass* objfcn)
{
  Problem::SetParameters(objfcn);
}

} // namespace OPTPP
