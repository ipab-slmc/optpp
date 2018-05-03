#ifndef PDSPROBLEM_H
#define PDSPROBLEM_H

#include "xercesc/dom/DOM.hpp"
#include "xercesc/util/PlatformUtils.hpp"

#include "OptPDS.h"
#include "Problem.h"

namespace OPTPP {

class PDSProblem:public Problem
{
	private:
		virtual DOMElement* FindParameterXML();

		virtual OptError CreateFunctionOptimizer(OptimizeClass * &objfcn, NLP0* &func);
		virtual OptError CreateApplicationOptimizer(OptimizeClass * &objfcn, NLP0* &func);

		virtual void SetParameters(OptimizeClass* objfcn);

	public:
		PDSProblem(DOMElement* solverXML):Problem(solverXML){;}
};
} // namespace OPTPP
#endif

