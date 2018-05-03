
#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

using namespace std;

#ifdef WITH_MPI
#include "mpi.h"
#endif

#include "xercesc/util/XMLString.hpp"

#include "Problem.h"

namespace OPTPP {

void * Problem::GetFunction(string libName, string funcName)
{
  char * error;
  void * handle;
  void* function;

  // Open up the specified library, or look in the global table if 
  // no library has been specified 
  if(libName == "")
    {
      dlopen(NULL, RTLD_LAZY);
    }
  else
    {
      handle = dlopen(libName.c_str(), RTLD_LAZY);
		
      // add the library handle to a vector 
      // so it will be closed at program completion
      loadedLibs_.push_back(handle);
    }
	
  // if something went wrong, print out the error
  if(!handle)
    {
//      cerr << dlerror() << endl;
      return NULL;
    }

  // Get the specified Function
  function = dlsym(handle, funcName.c_str());

  // if something went wrong, print out the error
  if((error = dlerror()) != NULL)
    {
//      cerr << error << endl;
//      cerr << "Possible causes are: Incorrect library" << endl;
//      cerr << "                     Incorrect function name" << endl;
//      cerr << "                     Function not delcared with \"C\" linkage"
//	   << endl;
      return NULL;
    }

  return function;
}

void Problem::CloseLoadedLibs()
{
  // loop throgh all the loaded libraries and close them 
  for(int i = 0; i < loadedLibs_.size(); i++)
    {
      dlclose(loadedLibs_[i]);
    }

  loadedLibs_.clear();
}

int Problem::GetNumVar()
{
  // first check if we have already read in this value yet.
  if(!hasNumVar_)
    {
      DOMNodeList* nodeList;
      DOMNode* tmpNode;
      DOMElement* variableXML;
      DOMElement* tmpXML;

      // if not, count the number of variables given in the
      // Variable Array
      if((nodeList =
	  applicationXML_->getElementsByTagName(XMLString::transcode("VariableClass")))->getLength() > 0)
	{
	  tmpNode = nodeList->item(0);
	  tmpXML = (DOMElement *)tmpNode;
	}
      if((nodeList =
	  applicationXML_->getElementsByTagName(XMLString::transcode("Array")))->getLength() > 0)
	{
	  tmpNode = nodeList->item(0);
	  variableXML = (DOMElement *)tmpNode;
	}

      numVar_ =
	variableXML->getElementsByTagName(XMLString::transcode("Variables"))->getLength();
      hasNumVar_ = true;
    }
  return numVar_;
}

int Problem::GetDerivOrder()
{
  // first check to see if we already know the derivative order
  if(!hasDerivOrder_)
    {
      // if not, look in the XML and see if we can find out.
      DOMElement* subroutineXML = GetSubroutineXML();
      string second =
	XMLString::transcode(subroutineXML->getAttribute(XMLString::transcode("Second")));
      string first =
	XMLString::transcode(subroutineXML->getAttribute(XMLString::transcode("First")));
      if(second == "true")
	{
	  derivOrder_ = 2;
	}else if(first == "true")
	  {
	    derivOrder_ = 1;
	  }else{
	    derivOrder_ = 0;
	  }
      hasDerivOrder_ = true;
    }

  return derivOrder_;
}

DOMElement* Problem::GetSubroutineXML()
{
  // if we already have found the Subroutine XML, just return it
  if(hasSubroutineXML_)
    return subroutineXML_;

  DOMNodeList* nodeList;
  DOMNode* tmpNode;

  // if we have not found it before, go get it and save the reference
  if((nodeList =
      applicationXML_->getElementsByTagName(XMLString::transcode("Library")))->getLength() > 0)
    {
      tmpNode = nodeList->item(0);
      subroutineXML_ = (DOMElement *) tmpNode;
      hasSubroutineXML_ = true;
    }

  // then return the xml
  return subroutineXML_;
}

DOMElement* Problem::GetSolverXML()
{
  return solverXML_;
}

DOMElement* Problem::GetParameterXML()
{
  // if we have not found the parameter XML yet, go and find it.
  if(!hasParameterXML_)
    {
      parameterXML_ = FindParameterXML();
      hasParameterXML_ = true;
    }
  // return the parameter XML
  return parameterXML_;
}


USERFCN0 Problem::GetUserFunction0()
{
  // load in the library function and cast appropriately before returning
  DOMElement* subroutineXML = GetSubroutineXML();
  string fEval =
    XMLString::transcode(subroutineXML->getAttribute(XMLString::transcode("FEval")));
  string libName =
    XMLString::transcode(subroutineXML->getAttribute(XMLString::transcode("LibName")));
  void* function = GetFunction(libName, fEval);
  return (USERFCN0) function;
}

USERFCN1 Problem::GetUserFunction1()
{
  // load in the library function and cast appropriately before returning
  DOMElement* subroutineXML = GetSubroutineXML();
  string fEval =
    XMLString::transcode(subroutineXML->getAttribute(XMLString::transcode("FEval")));
  string libName =
    XMLString::transcode(subroutineXML->getAttribute(XMLString::transcode("LibName")));
  void* function = GetFunction(libName, fEval);
  return (USERFCN1) function;
}

USERFCN2 Problem::GetUserFunction2()
{
  // load in the library function and cast appropriately before returning
  DOMElement* subroutineXML = GetSubroutineXML();
  string fEval =
    XMLString::transcode(subroutineXML->getAttribute(XMLString::transcode("FEval")));
  string libName =
    XMLString::transcode(subroutineXML->getAttribute(XMLString::transcode("LibName")));
  void* function = GetFunction(libName, fEval);
  return (USERFCN2) function;
}

INITFCN Problem::GetInitFunction()
{
  // load in the library function and cast appropriately before returning
  DOMElement* subroutineXML = GetSubroutineXML();
  string init =
    XMLString::transcode(subroutineXML->getAttribute(XMLString::transcode("Init")));
  string libName =
    XMLString::transcode(subroutineXML->getAttribute(XMLString::transcode("LibName")));
  void* function = GetFunction(libName, init);
  return (INITFCN) function;
}

AppLauncher * Problem::GetAppLauncher()
{
  return &launcher_;
}

void Problem::SetParameters(OptimizeClass* objfcn)
{
  // try to find the parameters that are common to all algorithms
  DOMElement* parameterXML = GetSolverXML();

  DOMNodeList* nodeList;
  DOMNode* tmpNode;
  DOMElement* optionXML;

  string outFile, maxIter, maxFeval, debug;
  string fcnTol, stepTol, gradTol, conTol, minStep, maxStep;

  if ((nodeList = parameterXML->getElementsByTagName(XMLString::transcode("BasicOptions")))->getLength() >0)
  {
    tmpNode = nodeList->item(0);
    optionXML = (DOMElement *) tmpNode;

    outFile =
      XMLString::transcode(optionXML->getAttribute(XMLString::transcode("outFile")));
    maxIter =
      XMLString::transcode(optionXML->getAttribute(XMLString::transcode("maxIter")));
    maxFeval =
      XMLString::transcode(optionXML->getAttribute(XMLString::transcode("maxFeval")));
    debug =
      XMLString::transcode(optionXML->getAttribute(XMLString::transcode("Debug")));
  }

  if ((nodeList = parameterXML->getElementsByTagName(XMLString::transcode("AdvancedOptions")))->getLength() >0)
  {
    tmpNode = nodeList->item(0);
    optionXML = (DOMElement *) tmpNode;

    fcnTol =
      XMLString::transcode(optionXML->getAttribute(XMLString::transcode("fcnTol")));
    stepTol =
      XMLString::transcode(optionXML->getAttribute(XMLString::transcode("stepTol")));
    minStep =
      XMLString::transcode(optionXML->getAttribute(XMLString::transcode("minStep")));
    maxStep =
      XMLString::transcode(optionXML->getAttribute(XMLString::transcode("maxStep")));
    conTol =
      XMLString::transcode(optionXML->getAttribute(XMLString::transcode("conTol")));
    gradTol =
      XMLString::transcode(optionXML->getAttribute(XMLString::transcode("gradTol")));
  }

#ifdef WITH_MPI
  int me;
#endif

  // if a parameter was present, then set it
  if(outFile != "")
#ifdef WITH_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    char status_file[80];
    sprintf(status_file,"%s.%d",outFile.c_str(), me);
    objfcn->setOutputFile(status_file, 0);
#else
    objfcn->setOutputFile(outFile.c_str(), 0);
#endif

  if(maxIter != "")
    objfcn->setMaxIter(atoi(maxIter.c_str()));

  if(maxFeval != "")
    objfcn->setMaxFeval(atoi(maxFeval.c_str()));

  if(fcnTol != "")
    objfcn->setFcnTol(atof(fcnTol.c_str()));

  if(stepTol != "")
    objfcn->setStepTol(atof(stepTol.c_str()));

  if(minStep != "")
    objfcn->setMinStep(atof(minStep.c_str()));

  if(maxStep != "")
    objfcn->setMaxStep(atof(maxStep.c_str()));

  if(debug == "true")
    objfcn->setDebug();

  if(conTol != "")
    objfcn->setConTol(atof(conTol.c_str()));

  if(gradTol != "")
    objfcn->setGradTol(atof(gradTol.c_str()));
}

VariableList * Problem::GetVariables()
{
  return &variables_;
}

CompoundConstraint * Problem::GetConstraints()
{
  return variables_.GetGeneralConstraints();
}


OptError Problem::optimize()
{	
  OptError error;

  // first create the appropriate optimizer.

  if(useApplicationOptimizer_)
    {
      error = CreateApplicationOptimizer(objfcn_, func_);
    }else{
      error = CreateFunctionOptimizer(objfcn_, func_);
    }

  if(error.value < 0)
    return error;

  //then set the algorithm parameters.
  SetParameters(objfcn_);

  // optimize the solution	
  objfcn_->optimize();
  // print out the results
  objfcn_->printStatus("Solution " );

  // clean up all of the junk
  objfcn_->cleanup();

  delete objfcn_;
  delete func_;
  if(!useApplicationOptimizer_)
    {
      CloseLoadedLibs();
    }

  return error;
}

void Problem::SetProblemXML(DOMElement* applicationXML)
{
  // save the problem XML
  applicationXML_ = applicationXML;
  variables_ = VariableList(applicationXML_);
  DOMNodeList* nodeList;
  DOMNode* tmpNode;
  DOMElement* appXML;

  // if the function is from an external application,
  // create a launcher and set the problem to use an Application
  // Optimizer
  if((nodeList =
      applicationXML_->getElementsByTagName(XMLString::transcode("Application")))->getLength() > 0)
    {
      tmpNode = nodeList->item(0);
      appXML = (DOMElement *)tmpNode;
      launcher_ = AppLauncher(appXML, variables_, true);
      useApplicationOptimizer_ = true;
    }
}

} // namespace OPTPP
