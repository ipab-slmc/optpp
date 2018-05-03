#ifndef NIPSPROBLEM_H
#define NIPSPROBLEM_H

#include "xercesc/dom/DOM.hpp"
#include "xercesc/util/PlatformUtils.hpp"

#include "OptNIPS.h"
#include "OptQNIPS.h"
#include "OptFDNIPS.h"
#include "Problem.h"

namespace OPTPP {

class NIPSProblem:public Problem
{
	private:
		MeritFcn meritFunction_;

	virtual DOMElement* FindParameterXML();

		virtual OptError CreateFunctionOptimizer(OptimizeClass * &objfcn, NLP0* &func);
		virtual OptError CreateApplicationOptimizer(OptimizeClass * &objfcn, NLP0* &func);

		virtual void SetParameters(OptimizeClass* objfcn);

	public:
		NIPSProblem(DOMElement* solverXML):Problem(solverXML){;}
};

} // namespace OPTPP
#endif

