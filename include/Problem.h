#ifndef PROBLEM_H
#define PROBLEM_H

#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#include <iostream>
#include <string>
#include <vector>
#include <dlfcn.h>
#include <sys/wait.h>
#include <unistd.h>
#ifdef HAVE_STD
#include <cstdlib>
#else
#include <stdlib.h>
#endif

#include "xercesc/dom/DOM.hpp"
#include "xercesc/util/PlatformUtils.hpp"

#include "Opt.h"
#include "AppLauncher.h"
#include "VariableList.h"
#include "NLFAPP.h"
#include "BoundConstraint.h"
#include "NonLinearEquation.h"
#include "CompoundConstraint.h"

using std::string;

namespace OPTPP {

typedef struct{
	int value;
	string msg;
} OptError;

/**
 * Class Problem is the parent for all the different types of solvers.
 */
class Problem
{
	private:
		DOMElement* solverXML_;
		bool hasSubroutineXML_;
		bool hasParameterXML_;
		bool useApplicationOptimizer_;
		bool hasNumVar_;
		bool hasDerivOrder_;
		OptimizeClass * objfcn_;
		NLP0 * func_;
		int numVar_;
		int derivOrder_;
		DOMElement* subroutineXML_;
		DOMElement* applicationXML_;
		VariableList variables_;
		
		

		AppLauncher launcher_;

		vector<void *> loadedLibs_;

		/**
		 * Dymanically loads the given Library function
		 * and returns int as a void pointer
		 */
		void * GetFunction(string libName, string funcName);
	
		/**
		 * Closes all of the Libraries that have been dynamically Loaded 
		 */
		void CloseLoadedLibs();

	protected:
		// this needs to be set by each individual subclass
		DOMElement* parameterXML_;
		
		/** Returns the number of Variables */
		int GetNumVar();

		/** Reterns the order of Differentiation */
		int GetDerivOrder();

		/**
		 * Returns the XML Element containing Subroutine and library info
		 */
		DOMElement* GetSubroutineXML();

		/**
		 * Returns the XML Element Which contains the solver properites
		 */
		DOMElement* GetSolverXML();

		/**
		 * Returns the XML Element containing the Algorithm parameters
		 */
		DOMElement* GetParameterXML();
			
		/**
		 * Initially Retrieves the XML Element containing Algorithm parameters
		 * from the Solver XML
		 */
		virtual DOMElement* FindParameterXML() = 0;
		
		/** Returns the non-differentiable function to be optimized */
		USERFCN0 GetUserFunction0();

		/** Returns the first order differentiable function to be optimized */
		USERFCN1 GetUserFunction1();

		/** Returns the second order differentiable function to be optimized */
		USERFCN2 GetUserFunction2();

		/** Returns the initialization funciton for the Optimization */
		INITFCN GetInitFunction();

		/** 
		 * Returns the Application Launcher for Optimizing an external 
		 * Function Evaluation
		 */
		AppLauncher * GetAppLauncher();
		
		/**
		 * Sets the Algorithm Parameters for the Optimization run
		 */
		virtual void SetParameters(OptimizeClass* objfcn);

		/**
		 * Gets a VariableList object containing the
		 * optimization  variables and related information.
		 */
		VariableList * GetVariables();

		/**
		 * Gets a CompoundConstraint object based upon whether the 
		 * optimization is unconstrained, bounded, or general
		 */
		CompoundConstraint * GetConstraints();

		/**
		 * Creates an Optimizer that will optimize a dynamically loaded
		 * function from an external library
		 */
		virtual OptError CreateFunctionOptimizer(OptimizeClass * &objfcn, NLP0* &func) = 0;

		/**
		 * Creates an Optimizer that will optimize a function run by an
		 * external Application
		 */
		virtual OptError CreateApplicationOptimizer(OptimizeClass * &objfcn, NLP0* &func) = 0;
		
	public:
		/** constructor */
		Problem(DOMElement* solverXML):
			solverXML_(solverXML),
			hasSubroutineXML_(false),
			hasParameterXML_(false),
			useApplicationOptimizer_(false),
			hasNumVar_(false),
			hasDerivOrder_(false)
		{;}

		/* Destructor */

		virtual ~Problem(){}
		
		/** Optimize the Function */
		OptError optimize();
		
		/** Sets the XML that defines the Problem to be solved */
		void SetProblemXML(DOMElement* applicationXML);
};

} // namespace OPTPP

#endif
