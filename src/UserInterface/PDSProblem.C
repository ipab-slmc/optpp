
using namespace std;

#include "xercesc/util/XMLString.hpp"

#include "PDSProblem.h"

namespace OPTPP {

DOMElement* PDSProblem::FindParameterXML()
{
  DOMElement* solverXML = GetSolverXML();

  DOMNodeList* nodeList;
  DOMNode* tmpNode;
  DOMElement* parameterXML;

  if((nodeList =
      solverXML->getElementsByTagName(XMLString::transcode("PDSvar")))->getLength() > 0)
    {
      tmpNode = nodeList->item(0);
      parameterXML = (DOMElement *) tmpNode;
    }

  return parameterXML;
}

OptError PDSProblem::CreateFunctionOptimizer(OptimizeClass * &objfcn,
					     NLP0* &func)
{
  OptError error;
  int numVar = GetNumVar();
  int derivOrder = GetDerivOrder();
  VariableList* variables = GetVariables();

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

      NLF2 * myFunc = new NLF2(numVar, userFcn, initFcn, GetConstraints());
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

      NLF1 * myFunc = new NLF1(numVar, userFcn, initFcn, GetConstraints());
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
      NLF0 * myFunc = new NLF0(numVar, userFcn, initFcn, GetConstraints());

      myFunc->setIsExpensive(true);
      myFunc->setX(variables->GetInitialValues());

      func = myFunc;
    }
  objfcn = new OptPDS(func);
  error.value = 0; 
  error.msg = "Success";
  return error;
}

OptError PDSProblem::CreateApplicationOptimizer(OptimizeClass * &objfcn,
						NLP0* &func)
{
  OptError error;
  int numVar = GetNumVar();
  AppLauncher * launcher = GetAppLauncher();
  USERFCN0APP userFcn = &(launcher->run_app);
  INITFCNAPP initFcn = &(launcher->init_app);

  NLF0APP * myFunc = new NLF0APP(numVar, userFcn, initFcn, launcher,
				 GetConstraints());

  myFunc->setIsExpensive(true);
  func = myFunc;

  objfcn = new OptPDS(myFunc);
  error.value = 0;
  error.msg = "Success";
  return error;
}	

void PDSProblem::SetParameters(OptimizeClass* objfcn)
{
  Problem::SetParameters(objfcn);

  DOMElement* pdsXML = GetParameterXML();

  string searchSize, simplexType;

  OptPDS * objfcnPDS = (OptPDS *) objfcn;

  searchSize =
    XMLString::transcode(pdsXML->getAttribute(XMLString::transcode("searchSize")));
  if (searchSize != "" )
    objfcnPDS->setSSS(atoi(searchSize.c_str()));

  simplexType = 
    XMLString::transcode(pdsXML->getAttribute(XMLString::transcode("simpType")));
  if (simplexType != "")
  {
    objfcnPDS->setSimplexType(atoi(simplexType.c_str()));
  }
}

} // namespace OPTPP
