
#include <dlfcn.h>

using namespace std;

#include "xercesc/util/XMLString.hpp"

#include "VariableList.h"
#include "AppLauncher.h"
#include "NLFAPP.h"

using NEWMAT::ColumnVector;
using NEWMAT::Matrix;
using NEWMAT::SymmetricMatrix;

namespace OPTPP {

VariableList::VariableList(DOMElement* problemXML)
{	
  compound_ = NULL;
  bound_ = NULL;
  linEquality_ = NULL;
  linInEquality_ = NULL;
	
  DOMNodeList* nodeList;
  DOMNode* tmpNode;
  DOMElement* variableXML;
  DOMElement* tmpXML;
  lboundsExist_ = false;
  uboundsExist_ = false;
  lconstraintsExist_ = false;
  nlconstraintsExist_ = false;

  // Get the Variable List XML
  if((nodeList =
      problemXML->getElementsByTagName(XMLString::transcode("VariableClass")))->getLength() > 0)
    {
      tmpNode = nodeList->item(0);
      tmpXML = (DOMElement *)tmpNode;
    }
  else
    {
//      cerr << "Error in XML file" << endl;
      exit(2);
    }
  if((nodeList =
      tmpXML->getElementsByTagName(XMLString::transcode("Array")))->getLength() > 0)
    {
      tmpNode = nodeList->item(0);
      variableListXML_ = (DOMElement *)tmpNode;
    }
  else
    {
//      cerr << "Error in XML file" << endl;
      exit(2);
    }

  // Read in the Variables.
  nodeList =
    variableListXML_->getElementsByTagName(XMLString::transcode("Variables"));

  numVars_ = nodeList->getLength();

  for(int i = 0; i < numVars_; i++)
    {
      tmpNode = nodeList->item(i);
      variableXML = (DOMElement *) tmpNode;

      inputVariable newVar;

      newVar.varName =
	XMLString::transcode(variableXML->getAttribute(XMLString::transcode("theName")));
      string initialValue =
	XMLString::transcode(variableXML->getAttribute(XMLString::transcode("initVal")));
      string upperBound =
	XMLString::transcode(variableXML->getAttribute(XMLString::transcode("upper")));
      string lowerBound =
	XMLString::transcode(variableXML->getAttribute(XMLString::transcode("lower")));

      newVar.initialValue = atof(initialValue.c_str());

      if(upperBound != "")
	{
	  newVar.hasUpperBound = true;
	  newVar.upperBound = atof(upperBound.c_str());
	  uboundsExist_ = true;
	}else{
	  newVar.hasUpperBound = false;
	}

      if(lowerBound != "")
	{
	  newVar.hasLowerBound = true;
	  newVar.lowerBound = atof(lowerBound.c_str());
	  lboundsExist_ = true;
	}else{
	  newVar.hasLowerBound = false;
	}

      variables_.push_back(newVar);	
    }

  // Get the Linear Constraint List XML
  if((nodeList =
      problemXML->getElementsByTagName(XMLString::transcode("LConstraintClass")))->getLength() > 0)
    {
      tmpNode = nodeList->item(0);
      tmpXML = (DOMElement *)tmpNode;
    }
  else
    {
//      cerr << "Error in XML file" << endl;
      exit(2);
    }
  if((nodeList =
      tmpXML->getElementsByTagName(XMLString::transcode("Array")))->getLength() > 0)
    {
      tmpNode = nodeList->item(0);
      linConstraintListXML_ = (DOMElement *)tmpNode;
      if (linConstraintListXML_->getElementsByTagName(XMLString::transcode("LinearConstraint"))->getLength() > 0)
	lconstraintsExist_ = true;
    }
  else
    {
//      cerr << "Error in XML file" << endl;
      exit(2);
    }
	
  // Get the Non-Linear Constraint List XML
  if((nodeList =
      problemXML->getElementsByTagName(XMLString::transcode("NLConstraintClass")))->getLength() > 0)
    {
      tmpNode = nodeList->item(0);
      tmpXML = (DOMElement *)tmpNode;
    }
  else
    {
//      cerr << "Error in XML file" << endl;
      exit(2);
    }
  if((nodeList =
      tmpXML->getElementsByTagName(XMLString::transcode("Array")))->getLength() > 0)
    {
      tmpNode = nodeList->item(0);
      nlConstraintListXML_ = (DOMElement *)tmpNode;
      if ((nlConstraintListXML_->getElementsByTagName(XMLString::transcode("Library"))->getLength() > 0) || (nlConstraintListXML_->getElementsByTagName(XMLString::transcode("Application"))->getLength() > 0))
	nlconstraintsExist_ = true;
    }
  else
    {
//      cerr << "Error in XML file" << endl;
      exit(2);
    }
}

VariableList::~VariableList()
{
  for(int i = 0; i < loadedLibs_.size(); i++)
    {
      dlclose(loadedLibs_[i]);
    }

  loadedLibs_.clear();

  for(int i = 0; i < constraintFunctions_.size(); i++)
    {
      delete constraintFunctions_[i];
    }

  for(int i = 0; i < constraintBase_.size(); i++)
    {
      delete constraintBase_[i];
    }

  for(int i = 0; i < nlEquations_.size(); i++)
    {
      delete nlEquations_[i];
    }

  for(int i = 0; i < launchers_.size(); i++)
    {
      delete launchers_[i];
    }

  if(compound_ != NULL)
    delete compound_;

  if(bound_ != NULL)
    delete bound_;

  if(linEquality_ != NULL)
    delete linEquality_;

  if(linInEquality_ != NULL)
    delete linInEquality_;
}

ColumnVector VariableList::GetInitialValues()
{
  ColumnVector init(variables_.size());
  for(int i = 0; i < variables_.size(); i++)
    {
      init(i+1) = variables_[i].initialValue;
    }
  return init;
}

ColumnVector VariableList::GetUpperBounds()
{
  ColumnVector upper(variables_.size());
  for(int i = 0; i < variables_.size(); i++)
    {
      if (variables_[i].hasUpperBound)
	upper(i+1) = variables_[i].upperBound;
      else
	upper(i+1) = MAX_BND;
    }
  return upper;
}

ColumnVector VariableList::GetLowerBounds()
{
  ColumnVector lower(variables_.size());
  for(int i = 0; i < variables_.size(); i++)
    {
      if (variables_[i].hasLowerBound)
	lower(i+1) = variables_[i].lowerBound;
      else
	lower(i+1) = MIN_BND;
    }
  return lower;
}

string VariableList::GetVariableName(int i)
{
  return variables_[i].varName;
}

BoundConstraint * VariableList::GetBoundConstraints()
{
  if ((lboundsExist_) && (uboundsExist_))
  {
    bound_ = new BoundConstraint(numVars_, 
				 GetLowerBounds(),
				 GetUpperBounds());
  }
  else if (lboundsExist_)
  {
    bound_ = new BoundConstraint(numVars_, GetLowerBounds());
  }
  else if (uboundsExist_)
  {
    bound_ = new BoundConstraint(numVars_, GetUpperBounds(), false);
  }

  return bound_;
}

CompoundConstraint * VariableList::GetGeneralConstraints()
{
  int numConstraints = 0;

  BoundConstraint* bounds = GetBoundConstraints();
  if (bounds != NULL)
  {
    numConstraints++;
  }

  DOMNodeList* nlNodeListLib =
    nlConstraintListXML_->getElementsByTagName(XMLString::transcode("Library"));

  DOMNodeList* nlNodeListApp =
    nlConstraintListXML_->getElementsByTagName(XMLString::transcode("Application"));

  DOMNodeList* linNodeList =
    linConstraintListXML_->getElementsByTagName(XMLString::transcode("LinearConstraint"));

  DOMNodeList* nodeList;
  DOMNode* tmpNode;
  DOMElement* tmpXML;

  // generate the constraints array.

  // check for existence of linear constraints
	
  int numLinear = linNodeList->getLength();
  int numEquality = 0;
  int numInequality = 0;
  int equalityCount = 0;
  int inequalityCount = 0;

  for(int i = 0; i < numLinear; i++)
  {
    tmpNode = linNodeList->item(i);
    DOMElement* constraintXML = (DOMElement *)tmpNode;

    string comparison =
      XMLString::transcode(constraintXML->getAttribute(XMLString::transcode("operator")));

    // choose the appropriate array (equality/inequality)
    if(comparison == "=")
    {
      numEquality++;
    }else{
      numInequality++;
    }
  }

  Matrix equalityMat(numEquality, numVars_);
  Matrix inequalityMat(numInequality, numVars_);
  ColumnVector equalityRHS(numEquality);
  ColumnVector inequalityRHS(numInequality);

  // Read in the linear Constraints
	
  // loop through each element in the array,
  for(int i = 0; i < numLinear; i++)
  {
    tmpNode = linNodeList->item(i);
    DOMElement* constraintXML = (DOMElement *)tmpNode;

    string comparison =
      XMLString::transcode(constraintXML->getAttribute(XMLString::transcode("operator")));

    Matrix * A;
    ColumnVector * b;

    int index1;

    // choose the appropriate array (equality/inequality)
    if(comparison == "=")
    {
      equalityCount++;
      index1 = equalityCount;

      A = &equalityMat;
      b = &equalityRHS;
    }else{
      inequalityCount++;
      index1 = inequalityCount;

      A = &inequalityMat;
      b = &inequalityRHS;
    }

    /*      if(index1 > numVars_)
	    {
	    cerr << "The system is over constrained, ignoring further linear constraints." << endl;
	    continue;
	    } */

    string rhs =
      XMLString::transcode(constraintXML->getAttribute(XMLString::transcode("rhs")));

    // get the rhs and add to the rhs vector
    (*b)(index1) = atof(rhs.c_str());

    if((nodeList =
	constraintXML->getElementsByTagName(XMLString::transcode("Table")))->getLength() > 0)
    {
      tmpNode = nodeList->item(0);
      tmpXML = (DOMElement *)tmpNode;
    }

    DOMNodeList* tableNodeList =
      tmpXML->getElementsByTagName(XMLString::transcode("Entry"));

    DOMElement* coefficientXML;

    // loop through the table, 
    for(int j  = 0; j < tableNodeList->getLength(); j++)
    {
      tmpNode = tableNodeList->item(j);
      coefficientXML = (DOMElement *) tmpNode;
      string variable =
	XMLString::transcode(coefficientXML->getAttribute(XMLString::transcode("Variable")));
      string coefficient =
	XMLString::transcode(coefficientXML->getAttribute(XMLString::transcode("Coefficient")));

      // find the appropriate variable index for each name given
      int index2;
      bool foundK = false;
      for(int k = 0; k < numVars_; k++)
      {
	if(variable == GetVariableName(k))
	{
	  foundK = true;
	  index2 = k+1; // increment to get index of Matrix correct
	  break;
	}
      }

      // add the coefficient in the proper array spot
      if(foundK)
      {
	(*A)(index1, index2) = atof(coefficient.c_str());
      }
    }
  }

  if(numInequality > 0)
  {
    linInEquality_ = new LinearInequality(inequalityMat, inequalityRHS);
    numConstraints++;
  }
	
  if(numEquality > 0)
  {
    linEquality_ = new LinearEquation(equalityMat, equalityRHS);
    numConstraints++;
  }
	
  // check for existence of non-linear constraints.

  int numNonLinearLib = nlNodeListLib->getLength();
  int numNonLinearApp = nlNodeListApp->getLength();
  int numNLEquality = 0;
  int numNLInequality = 0;

  NonLinearEquation * nonLinearE;
  NonLinearInequality * nonLinearI;
  ColumnVector nlRHS(1);

  numConstraints += numNonLinearLib;
  numConstraints += numNonLinearApp;
	
  OptppArray<Constraint> constraints(numConstraints);

  // Read in the Non-Linear Constraints
	
  int constraintIndex;
	
  // Loop through the function call constraints 
  for(constraintIndex = 0; constraintIndex < numNonLinearLib;
      constraintIndex++)
  {
    tmpNode = nlNodeListLib->item(constraintIndex);
    DOMElement* constraintXML = (DOMElement *)tmpNode;

    string comparison =
      XMLString::transcode(constraintXML->getAttribute(XMLString::transcode("operator")));

    string rhs =
      XMLString::transcode(constraintXML->getAttribute(XMLString::transcode("rhs")));

    // get the rhs and add to the rhs vector
    nlRHS(1) = atof(rhs.c_str());
		
    NLP * myFunc;
		
    // check the derivative order and instatiate the appropriate
    // function type
    // i.e. NLF0, NLF1, NLF2
    // and store the function, and NLFBase for future deletion
    int derivOrder = GetDerivOrder(constraintXML);

    if(derivOrder == 2)
    {
      USERNLNCON2 userFcn = GetConstraintFunction2(constraintXML);
      INITFCN initFcn = GetConstraintInitFunction(constraintXML);
      if(userFcn == NULL || initFcn == NULL)	
      {
	continue;
      }
		
      NLF2 * func = new NLF2(numVars_, 1, userFcn, initFcn);
      constraintBase_.push_back(func);
      func->initFcn();

      myFunc = new NLP(func);
      constraintFunctions_.push_back(myFunc);
    }
    else if(derivOrder == 1)
    {
      USERNLNCON1 userFcn = GetConstraintFunction1(constraintXML);
      INITFCN initFcn = GetConstraintInitFunction(constraintXML);
      if(userFcn == NULL || initFcn == NULL)	
      {
	continue;
      }
		
      NLF1 * func = new NLF1(numVars_, 1, userFcn, initFcn);
      constraintBase_.push_back(func);
      func->initFcn();

      myFunc = new NLP(func);
      constraintFunctions_.push_back(myFunc);
    }
    else 
    {
      USERNLNCON0 userFcn = GetConstraintFunction(constraintXML);
      INITFCN initFcn = GetConstraintInitFunction(constraintXML);
      if(userFcn == NULL || initFcn == NULL)	
      {
	continue;
      }
		
      NLF0 * func = new NLF0(numVars_, 1, userFcn, initFcn);
      constraintBase_.push_back(func);
      func->initFcn();

      myFunc = new NLP(func);
      constraintFunctions_.push_back(myFunc);
    }

    // store the equation for later deletion
    if(comparison == "=")
    {
      numNLEquality++;
      // create a nonlinear equation from the function
      nonLinearE = new NonLinearEquation(myFunc, nlRHS);
      nlEquations_.push_back(nonLinearE);
      // add the constriat to the CompoundConstraints.
      constraints.put(constraintIndex,Constraint(nonLinearE));
    }
    else
    {
      numNLInequality++;
      // create a nonlinear equation from the function
      nonLinearI = new NonLinearInequality(myFunc, nlRHS);
      nlInequalities_.push_back(nonLinearI);
      // add the constriat to the CompoundConstraints.
      constraints.put(constraintIndex,Constraint(nonLinearI));
    }
  }

  // what about App called non-linear constraints? (Hard problem)
  for(int i = 0; i < numNonLinearApp; i++)
  {
    tmpNode = nlNodeListApp->item(i);
    DOMElement* constraintXML = (DOMElement *)tmpNode;

    string comparison =
      XMLString::transcode(constraintXML->getAttribute(XMLString::transcode("operator")));

    string rhs =
      XMLString::transcode(constraintXML->getAttribute(XMLString::transcode("rhs")));

    // get the rhs and add to the rhs vector
    nlRHS(1) = atof(rhs.c_str());

    // Generate a new Applauncher for the constraint
    AppLauncher * launcher = new AppLauncher(constraintXML, *this, false);

    // store the launcher for future deletion.
    launchers_.push_back(launcher);

    // get the appropriate function calls
    USERNLNCON0APP userFcn = &(launcher->run_app_nlncon);
    INITFCNAPP initFcn = &(launcher->init_app);
		
    // cereate the functions and store them for future deletion
    FDNLF1APP * func = new FDNLF1APP(numVars_, 1, userFcn,
				     initFcn, launcher);
    constraintBase_.push_back(func);
    func->initFcn();

    NLP * myFunc = new NLP(func);
    constraintFunctions_.push_back(myFunc);

    // store the equation for later deletion
    if(comparison == "=")
    {
      numNLEquality++;
      // create a nonlinear equation from the function
      nonLinearE = new NonLinearEquation(myFunc, nlRHS);
      nlEquations_.push_back(nonLinearE);
      // add the constriat to the CompoundConstraints.
      constraints.put(constraintIndex,Constraint(nonLinearE));
    }
    else
    {
      numNLInequality++;
      // create a nonlinear equation from the function
      nonLinearI = new NonLinearInequality(myFunc, nlRHS);
      nlInequalities_.push_back(nonLinearI);
      // add the constriat to the CompoundConstraints.
      constraints.put(constraintIndex,Constraint(nonLinearI));
    }
    constraintIndex++;
  }

  // create the compound constraints from bound/linear/non-linear constraints.
	
  if(numInequality > 0){
    constraints.put(constraintIndex, Constraint(linInEquality_));
    constraintIndex++;
  }
	
  if(numEquality > 0){
    constraints.put(constraintIndex, Constraint(linEquality_));
    constraintIndex++;
  }

  if (bounds != NULL)
  {
    constraints.put(constraintIndex, Constraint(bounds));
  }

  if (numConstraints > 0)
  {
    compound_ = new CompoundConstraint(constraints);
  }
  else
  {
    compound_ = 0;
  }

  return compound_;
}


USERNLNCON0 VariableList::GetConstraintFunction(DOMElement* constraintXML)
{
  string fEval =
    XMLString::transcode(constraintXML->getAttribute(XMLString::transcode("FEval")));
  string libName =
    XMLString::transcode(constraintXML->getAttribute(XMLString::transcode("LibName")));
  void* function = GetFunction(libName, fEval);
  return (USERNLNCON0) function;
}

USERNLNCON1 VariableList::GetConstraintFunction1(DOMElement* constraintXML)
{
  string fEval =
    XMLString::transcode(constraintXML->getAttribute(XMLString::transcode("FEval")));
  string libName =
    XMLString::transcode(constraintXML->getAttribute(XMLString::transcode("LibName")));
  void* function = GetFunction(libName, fEval);
  return (USERNLNCON1) function;
}

USERNLNCON2 VariableList::GetConstraintFunction2(DOMElement* constraintXML)
{
  string fEval =
    XMLString::transcode(constraintXML->getAttribute(XMLString::transcode("FEval")));
  string libName =
    XMLString::transcode(constraintXML->getAttribute(XMLString::transcode("LibName")));
  void* function = GetFunction(libName, fEval);
  return (USERNLNCON2) function;
}

INITFCN VariableList::GetConstraintInitFunction(DOMElement* constraintXML)
{
  string init =
    XMLString::transcode(constraintXML->getAttribute(XMLString::transcode("Init")));
  string libName =
    XMLString::transcode(constraintXML->getAttribute(XMLString::transcode("LibName")));
  void* function = GetFunction(libName, init);
  return (INITFCN) function;
}

void * VariableList::GetFunction(string libName, string funcName)
{
  char * error;
  void * handle;
  void* function;

  if(libName == "")
    {
      dlopen(NULL, RTLD_LAZY);
    }
  else
    {
      handle = dlopen(libName.c_str(), RTLD_LAZY);
      loadedLibs_.push_back(handle);
    }

  if(!handle)
    {
//      cerr << dlerror() << endl;
      return NULL;
    }

  function = dlsym(handle, funcName.c_str());
  if((error = dlerror()) != NULL)
    {
//      cerr << error << endl;
//      cerr << "Possible causes are: Incorrect library" << endl;
//      cerr << "                     Incorrect function name" << endl;
//      cerr << "                     Function not delcared with \"C\" linkage" << endl;
      return NULL;
    }

  return function;
}

int VariableList::GetDerivOrder(DOMElement* subroutineXML)
{
  int derivOrder;

  string second =
    XMLString::transcode(subroutineXML->getAttribute(XMLString::transcode("Second")));
  string first =
    XMLString::transcode(subroutineXML->getAttribute(XMLString::transcode("First")));
  if(second == "true")
    {
      derivOrder = 2;
    }else if(first == "true")
      {
	derivOrder = 1;
      }else{
	derivOrder = 0;
      }

  return derivOrder;
}

bool VariableList::lowerExists()
{
  return lboundsExist_;
}

bool VariableList::upperExists()
{
  return uboundsExist_;
}

bool VariableList::linearExists()
{
  return lconstraintsExist_;
}

bool VariableList::nonlinearExists()
{
  return nlconstraintsExist_;
}

} // namespace OPTPP
