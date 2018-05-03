#ifndef APPLAUNCHER_H
#define APPLAUNCHER_H

#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#include <iostream>
#include <string>

#include <dlfcn.h>
#include <unistd.h>
#include <sys/wait.h>
#ifdef HAVE_STD
#include <cstdio>
#include <cstring>
#else
#include <stdio.h>
#include <string.h>
#endif

#include "xercesc/dom/DOM.hpp"
#include "xercesc/util/PlatformUtils.hpp"

#include "BoundConstraint.h"
#include "NonLinearEquation.h"
#include "CompoundConstraint.h"
#include "VariableList.h"

using std::string;

namespace OPTPP {

class AppLauncher
{
	private: 

		string appName_;
		 string appInput_;
		 string appOutput_;
		 static string appDir_;
		 VariableList * variables_;

	public:
		/** no-arg Constructor*/
		AppLauncher(){;}

		/** no-variable Constructor */
		AppLauncher(DOMElement* appXML, bool createDir);

		/** main constructor */
		AppLauncher(DOMElement* appXML,
			    VariableList& variables, bool createDir);
		
		
		/** Returns the initial Values in x */
		void init_app(int ndim, NEWMAT::ColumnVector& x);

		/** Static method for init_app which can be used by OPT++ */
		static void init_app(int ndim, NEWMAT::ColumnVector& x, AppLauncher * launcher)
		{
			launcher->init_app(ndim, x);
		}

		/** Evaluates the Application function*/
		void run_app(int ndim, const NEWMAT::ColumnVector& x, double& fx, int& result);
		
		/** static method for run_app which can be used by OPT++ */
		static void run_app(int ndim, const NEWMAT::ColumnVector& x, double& fx, 
					int& result, AppLauncher * launcher)
		{
			launcher->run_app(ndim, x, fx, result);
		}

		/** Evaluates an Application NonLinearConstraint */
		void run_app_nlncon(int ndim, int nlncons, const NEWMAT::ColumnVector& x, NEWMAT::ColumnVector& fx, int& result);

		/** Static method for run_app_nlncon which can be used by OPT++*/
		static void run_app_nlncon(int ndim, int nlncons, const NEWMAT::ColumnVector& x, NEWMAT::ColumnVector& fx, 
					int& result, AppLauncher * launcher)
		{
			launcher->run_app_nlncon(ndim, nlncons, x, fx, result);
		}

		/** Actually Launches the Application and gathers results */
		void RunFunctionEvaluation(int ndim, const NEWMAT::ColumnVector & x);
		int setupin(int ndim, const NEWMAT::ColumnVector& x,
			    const char *fileName);
		int substitute_value(char *newLine, char *line,
				     const char *pattern, double value);
};	

} // namespace OPTPP
#endif
