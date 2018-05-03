#ifndef CGPROBLEM_H
#define CGPROBLEM_H

#include "xercesc/dom/DOM.hpp"
#include "xercesc/util/PlatformUtils.hpp"

#include "OptCG.h"
#include "Problem.h"

namespace OPTPP {

class CGProblem:public Problem
{
	private:

		virtual DOMElement* FindParameterXML();

		virtual OptError CreateFunctionOptimizer(OptimizeClass * &objfcn, NLP0* &func);
		virtual OptError CreateApplicationOptimizer(OptimizeClass * &objfcn, NLP0* &func);

		virtual void SetParameters(OptimizeClass* objfcn);

	public:
		CGProblem(DOMElement* solverXML):Problem(solverXML){;}
};

} // namespace OPTPP
#endif

