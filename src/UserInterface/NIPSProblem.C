
using namespace std;

#include "xercesc/util/XMLString.hpp"

#include "NIPSProblem.h"


namespace OPTPP {

DOMElement* NIPSProblem::FindParameterXML()
{	
  DOMElement* solverXML = GetSolverXML();

  DOMNodeList* nodeList;
  DOMNode* tmpNode;
  DOMElement* parameterXML;

  string atName = "argaezTapia";
  string nfmName = "NormFMu";
  string vsName = "VanShanno";

  if((nodeList =
      solverXML->getElementsByTagName(XMLString::transcode(atName.c_str())))->getLength() > 0)
    {
      meritFunction_ = ArgaezTapia;
      tmpNode = nodeList->item(0);
      parameterXML = (DOMElement *) tmpNode;
    }
  else if((nodeList =
	   solverXML->getElementsByTagName(XMLString::transcode(nfmName.c_str())))->getLength() > 0)
    {
      meritFunction_ = NormFmu;
      tmpNode = nodeList->item(0);
      parameterXML = (DOMElement *) tmpNode;
    }
  else if((nodeList = solverXML->getElementsByTagName(XMLString::transcode(vsName.c_str())))->getLength() > 0)
    {
      meritFunction_ = VanShanno;
      tmpNode = nodeList->item(0);
      parameterXML = (DOMElement *) tmpNode;
    }
  else{
//    cerr << "CANNOT determine merit function type" << endl;
    exit(1);
  }

  return parameterXML;
}

OptError NIPSProblem::CreateFunctionOptimizer(OptimizeClass * &objfcn,
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

      objfcn = new OptNIPS(myFunc);

      func = myFunc;
    }
  else if(derivOrder == 1)
    {
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
      objfcn = new OptFDNIPS(myFunc);
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

      FDNLF1 * myFunc = new FDNLF1(numVar, userFcn, initFcn, GetConstraints());

      myFunc->setIsExpensive(true);
      myFunc->setX(variables->GetInitialValues());

      objfcn = new OptFDNIPS(myFunc);

      func = myFunc;
    }	
  else
    {
      error.value = -3;
      error.msg = "Invalid Derivative Order has been given"; 
      return error;
    }

  error.value = 0;
  error.msg = "Success";
  return error;
}

OptError NIPSProblem::CreateApplicationOptimizer(OptimizeClass * &objfcn,
						 NLP0* &func)
{
  OptError error;
  int numVar = GetNumVar();
  AppLauncher * launcher = GetAppLauncher();
  USERFCN0APP userFcn = &(launcher->run_app);
  INITFCNAPP initFcn = &(launcher->init_app);

  FDNLF1APP * myFunc = new FDNLF1APP(numVar, userFcn, initFcn,
				     launcher, GetConstraints());

  myFunc->setIsExpensive(true);
  func = myFunc;

  objfcn = new OptFDNIPS(myFunc);
  error.value = 0;
  error.msg = "Success";
  return error;
}	

void NIPSProblem::SetParameters(OptimizeClass* objfcn)
{
  Problem::SetParameters(objfcn);

  DOMElement* meritXML = GetParameterXML();

  string centParm, lenBound, maxBack;

  OptNIPSLike * objfcnNIPS = (OptNIPSLike *) objfcn;
  objfcnNIPS->setMeritFcn(meritFunction_);

  centParm =
    XMLString::transcode(meritXML->getAttribute(XMLString::transcode("centParm")));
  if (centParm != "")
    objfcnNIPS->setCenteringParameter(atof(centParm.c_str()));

  lenBound =
    XMLString::transcode(meritXML->getAttribute(XMLString::transcode("lenBound")));
  if (lenBound != "")
    objfcnNIPS->setStepLengthToBdry(atof(lenBound.c_str()));

  maxBack =
    XMLString::transcode(meritXML->getAttribute(XMLString::transcode("maxBack")));
  if (maxBack != "")
    objfcnNIPS->setMaxBacktrackIter(atoi(maxBack.c_str()));
  objfcnNIPS->setSearchStrategy(LineSearch);
}

} // namespace OPTPP
