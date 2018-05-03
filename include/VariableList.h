#ifndef VARIABLELIST_H
#define VARIABLELIST_H

#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#include <iostream>
#include <string>
#include <vector>
#ifdef HAVE_STD
#include <cstdlib>
#else
#include <stdlib.h>
#endif

#include "xercesc/dom/DOM.hpp"
#include "xercesc/util/PlatformUtils.hpp"

#include "newmat.h"

#include "Opt.h"
#include "BoundConstraint.h"
#include "NonLinearEquation.h"
#include "NonLinearInequality.h"
#include "CompoundConstraint.h"
#include "LinearInequality.h"
#include "LinearEquation.h"

using std::string;

XERCES_CPP_NAMESPACE_USE

namespace OPTPP {

typedef struct{
	string varName;
	double initialValue;
	bool hasUpperBound;
	double upperBound;
	bool hasLowerBound;
	double lowerBound;
}inputVariable;

class AppLauncher;

class VariableList
{
	private:

		vector<inputVariable> variables_;
		DOMElement* variableListXML_;
		DOMElement* nlConstraintListXML_;
		DOMElement* linConstraintListXML_;
		int numVars_;
		vector<void *> loadedLibs_;
		vector<NLP *> constraintFunctions_;
		vector<NLPBase *> constraintBase_;
		vector<NonLinearEquation *> nlEquations_;
		vector<NonLinearInequality *> nlInequalities_;
		vector<AppLauncher *> launchers_;
		CompoundConstraint * compound_;
		BoundConstraint * bound_;
		LinearEquation * linEquality_;
		LinearInequality * linInEquality_;
		bool lboundsExist_, uboundsExist_;
		bool lconstraintsExist_, nlconstraintsExist_;
		

	public:
		/** no-arg constructor */
		VariableList(){	compound_ = NULL;
			bound_ = NULL;
			linEquality_ = NULL;
			linInEquality_ = NULL;
			lboundsExist_ = false;
			uboundsExist_ = false;
			lconstraintsExist_ = false;
			nlconstraintsExist_ = false;
		}
		/** main constructor */
		VariableList(DOMElement* problemXML);

		/** descructor */
		~VariableList();
	
		/** gets the initial Values */
		NEWMAT::ColumnVector GetInitialValues();

		/** Gets the upper bounds */
		NEWMAT::ColumnVector GetUpperBounds();

		/** Gets the lower bounds */
		NEWMAT::ColumnVector GetLowerBounds();
		
		/** Gets the name of the specified variable */
		string GetVariableName(int i);
		
		/** returns the number of variables */
		int size(){return numVars_;}
		
		/** Retuns the Boundary Constraints */
		BoundConstraint * GetBoundConstraints();

		/** Returns all of the general constraints */
		CompoundConstraint * GetGeneralConstraints();
		
		/** Returns the non-differentiable constraint function */
		USERNLNCON0 GetConstraintFunction(DOMElement* constraintXML);
		
		/** Returns the first order differentiable constraint function */
		USERNLNCON1 GetConstraintFunction1(DOMElement* constraintXML);
		
		/** Returns the second order differentiable constraint function */
		USERNLNCON2 GetConstraintFunction2(DOMElement* constraintXML);
		
		/** Returns the constraint initialization function */
		INITFCN GetConstraintInitFunction(DOMElement* constraintXML);

		/** Dynamically loads a function from the given library */
		void * GetFunction(string libName, string funcName);

		/** Gets the order of differentiation */
		int GetDerivOrder(DOMElement* subroutineXML);

		bool lowerExists();
		bool upperExists();
		bool linearExists();
		bool nonlinearExists();
};

} // namespace OPTPP
#endif
