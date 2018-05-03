
using namespace std;

#include "xercesc/util/XMLString.hpp"

#include "NewtonProblem.h"

namespace OPTPP {

const int NewtonProblem::lineSearch = 1;
const int NewtonProblem::trustRegion = 2;
const int NewtonProblem::trustPDS = 3;

DOMElement* NewtonProblem::FindParameterXML()
{
  DOMElement* solverXML = GetSolverXML();

  DOMElement* parameterXML;
  DOMNodeList* nodeList;
  DOMNode* tmpNode;

  string trName = "trustRegion";
  string lsName = "lineSearch";
  string pdsName = "trustPDS";

  // now check to see which search method is used.

  if((nodeList =
      solverXML->getElementsByTagName(XMLString::transcode(trName.c_str())))->getLength() > 0)
    {
      searchType_ = trustRegion;
      tmpNode = nodeList->item(0);
      parameterXML = (DOMElement *) tmpNode;
    }
  else if((nodeList =
	   solverXML->getElementsByTagName(XMLString::transcode(lsName.c_str())))->getLength() > 0)
    {
      searchType_ = lineSearch;
      tmpNode = nodeList->item(0);
      parameterXML = (DOMElement *) tmpNode;
    }
  else if((nodeList =
	   solverXML->getElementsByTagName(XMLString::transcode(pdsName.c_str())))->getLength() > 0)
  {
    searchType_ = trustPDS;
    tmpNode = nodeList->item(0);
    parameterXML = (DOMElement *) tmpNode;
  }
  else
  {
//    cerr << "CANNOT determine search strategy type" << endl;
    exit(1);
  }
  return parameterXML;
}

OptError NewtonProblem::CreateFunctionOptimizer(OptimizeClass *
						&objfcn, NLP0 * &func)
{
  OptError error;

  int numVar = GetNumVar();
  int derivOrder = GetDerivOrder();
  VariableList* variables = GetVariables();
  CompoundConstraint* constraints = 0;

  if ((variables->linearExists()) || (variables->nonlinearExists()))
  {
    cout << "Newton does not support linear or nonlinear constraints." << endl;
    cout << "Linear and nonlinear constraints being ignored." << endl;
  }
  if ((variables->upperExists()) || (variables->lowerExists()))
  {
    constraints = new CompoundConstraint(Constraint(variables->GetBoundConstraints()));
  }

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
		
      NLF2 * myFunc = new NLF2(numVar, userFcn, initFcn, constraints);
      myFunc->setIsExpensive(true);
      myFunc->setX(variables->GetInitialValues());

      if ((variables->upperExists()) || (variables->lowerExists()))
      {
	objfcn = new OptBCNewton(myFunc);
      }
      else
      {
	objfcn = new OptNewton(myFunc);
      }

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

      NLF1 * myFunc = new NLF1(numVar, userFcn, initFcn, constraints);
      myFunc->setIsExpensive(true);
      myFunc->setX(variables->GetInitialValues());

      if ((variables->upperExists()) || (variables->lowerExists()))
      {
	objfcn = new OptBCQNewton(myFunc);
      }
      else
      {
	objfcn = new OptQNewton(myFunc);
      }

      func = myFunc;
    }
  else if(derivOrder == 0)
    {
      INITFCN initFcn = GetInitFunction();
      USERFCN0 userFcn = GetUserFunction0();
      if(userFcn == NULL)
	{
	  error.value = -4;
	  error.msg = "Error loading  user function";
	  return error;
	}
      else if(initFcn == NULL)	
	{
	  error.value = -4;
	  error.msg = "Error loading init function";
	  return error;
	}

      // construct a function with no derivative info 
      // but use finite differences to approximate a first derivative
		
      FDNLF1 * myFunc = new FDNLF1(numVar, userFcn, initFcn, constraints);
      myFunc->setIsExpensive(true);
      myFunc->setX(variables->GetInitialValues());
      if(searchType_ == trustPDS)
	myFunc->setSpecOption(NoSpec);

      if ((variables->upperExists()) || (variables->lowerExists()))
      {
	objfcn = new OptBCQNewton(myFunc);
      }
      else
      {
	objfcn = new OptQNewton(myFunc);
      }

      func = myFunc;
    }

  error.value = 0; 
  error.msg = "Success";
  return error;
}

OptError NewtonProblem::CreateApplicationOptimizer(OptimizeClass * &objfcn,
						   NLP0 * &func)
{
  OptError error;
  int numVar = GetNumVar();
  VariableList* variables = GetVariables();
  CompoundConstraint* constraints = 0;

  AppLauncher * launcher = GetAppLauncher();
  USERFCN0APP userFcn = &(launcher->run_app);
  INITFCNAPP initFcn = &(launcher->init_app);

  if(userFcn == NULL)
  {
    error.value = -4;
    error.msg = "Error loading  user function";
    return error;
  }
  else if(initFcn == NULL)	
  {
    error.value = -4;
    error.msg = "Error loading init function";
    return error;
  }

  if ((variables->linearExists()) || (variables->nonlinearExists()))
  {
    cout << "Newton does not support linear or nonlinear constraints." << endl;
    cout << "Linear and nonlinear constraints being ignored." << endl;
  }
  if ((variables->upperExists()) || (variables->lowerExists()))
  {
    constraints = new CompoundConstraint(Constraint(variables->GetBoundConstraints()));
  }

  FDNLF1APP * myFunc = new FDNLF1APP(numVar, userFcn, initFcn,
				     launcher, constraints);
  myFunc->setIsExpensive(true);
  if(searchType_ == trustPDS)
    myFunc->setSpecOption(NoSpec);
  func = myFunc;

  if ((variables->upperExists()) || (variables->lowerExists()))
  {
    objfcn = new OptBCQNewton(myFunc);
  }
  else
  {
    objfcn = new OptQNewton(myFunc);
  }

  error.value = 0;
  error.msg = "Success";
  return error;
}	

void NewtonProblem::SetParameters(OptimizeClass* objfcn)
{
  Problem::SetParameters(objfcn);

  DOMElement* searchXML = GetParameterXML();
  VariableList* variables = GetVariables();

  string gradMult, searchSize, maxBack;

  OptNewtonLike * objfcnNewt = (OptNewtonLike *) objfcn;

  if(searchType_ == trustPDS)
  {
    objfcnNewt->setSearchStrategy(TrustPDS);

    gradMult =
      XMLString::transcode(searchXML->getAttribute(XMLString::transcode("gradMult")));
    if (gradMult != "")
      objfcnNewt->setGradMult(atof(gradMult.c_str()));

    searchSize =
      XMLString::transcode(searchXML->getAttribute(XMLString::transcode("searchSize")));
    if (searchSize != "")
      objfcnNewt->setSearchSize(atoi(searchSize.c_str()));
  }
  else if (searchType_ == lineSearch)
  {
    objfcnNewt->setSearchStrategy(LineSearch);

    maxBack =
      XMLString::transcode(searchXML->getAttribute(XMLString::transcode("maxBTIter")));
    if (maxBack != "")
      objfcnNewt->setMaxBacktrackIter(atoi(maxBack.c_str()));
  }
  else if (searchType_ == trustRegion)
  {
    if ((variables->upperExists()) || (variables->lowerExists()))
    {
      cout << "Newton with Trust Region does not support bounds." << endl;
      cout << "Using Line Search instead." << endl;
      objfcnNewt->setSearchStrategy(LineSearch);
    }
    else
    {
      objfcnNewt->setSearchStrategy(TrustRegion);
      gradMult =
	XMLString::transcode(searchXML->getAttribute(XMLString::transcode("gradMult")));
      if (gradMult != "")
	objfcnNewt->setGradMult(atof(gradMult.c_str()));
    }
  }
  else
  {
//    cerr << "Unrecognized search strategy type" << endl;
    exit(1);
  }
}

}
