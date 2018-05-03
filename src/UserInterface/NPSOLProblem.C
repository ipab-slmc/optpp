
using namespace std;

#include "xercesc/util/XMLString.hpp"

#include "NPSOLProblem.h"

namespace OPTPP {

DOMElement* NPSOLProblem::FindParameterXML()
{
  DOMElement* solverXML = GetSolverXML();

  DOMElement* parameterXML;
  DOMNodeList* nodeList;
  DOMNode* tmpNode;

  if((nodeList =
      solverXML->getElementsByTagName(XMLString::transcode("NPSOLvar")))->getLength() > 0)
    {
      tmpNode = nodeList->item(0);
      parameterXML = (DOMElement *) tmpNode;
    }

  return parameterXML;
}

OptError NPSOLProblem::CreateFunctionOptimizer(OptimizeClass *
					       &objfcn, NLP0* &func)
{
  OptError error;
  error.value = -5;
  error.msg = "This method is not implemented"; 
  return error;
}

OptError NPSOLProblem::CreateApplicationOptimizer(OptimizeClass * &objfcn,
						  NLP0* &func)
{
  OptError error;
  error.value = -5;
  error.msg = "This method is not implemented"; 
  return error;
}	

void NPSOLProblem::SetParameters(OptimizeClass* objfcn)
{
  Problem::SetParameters(objfcn);
  //  DOMElement* parameterXML = GetParameterXML();
}

} // namespace OPTPP
