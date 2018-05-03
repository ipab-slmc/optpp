#ifndef NPSOLPROBLEM_H
#define NPSOLPROBLEM_H

#include "xercesc/dom/DOM.hpp"
#include "xercesc/util/PlatformUtils.hpp"

#include "OptNPSOL.h"
#include "Problem.h"

namespace OPTPP {

class NPSOLProblem:public Problem
{
	private:

		virtual DOMElement* FindParameterXML();

		virtual OptError CreateFunctionOptimizer(OptimizeClass * &objfcn, NLP0* &func);
		virtual OptError CreateApplicationOptimizer(OptimizeClass * &objfcn, NLP0* &func);

		virtual void SetParameters(OptimizeClass* objfcn);

	public:
		NPSOLProblem(DOMElement* solverXML):Problem(solverXML){;}
};
} // namespace OPTPP
#endif

