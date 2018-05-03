#ifndef NewtonPROBLEM_H
#define NewtonPROBLEM_H

#include "xercesc/dom/DOM.hpp"
#include "xercesc/util/PlatformUtils.hpp"

#include "OptNewton.h"
#include "OptQNewton.h"
#include "OptBCNewton.h"
#include "OptBCQNewton.h"
#include "Problem.h"

namespace OPTPP {

class NewtonProblem:public Problem
{
	private:
		int searchType_;
		static const int lineSearch;
		static const int trustRegion;
		static const int trustPDS;

	virtual DOMElement* FindParameterXML();

		virtual OptError CreateFunctionOptimizer(OptimizeClass * &objfcn, NLP0* &func);
		virtual OptError CreateApplicationOptimizer(OptimizeClass * &objfcn, NLP0* &func);

		virtual void SetParameters(OptimizeClass* objfcn);

	public:
		NewtonProblem(DOMElement* solverXML):Problem(solverXML){;}
};

} // namespace OPTPP
#endif

