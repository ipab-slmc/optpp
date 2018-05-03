
#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#include <iostream>
#ifdef HAVE_STD
#include <cstdlib>
#else
#include <stdlib.h>
#endif

#include "OptppFatalError.h"

using namespace std;

namespace OPTPP{

// user-level error-invocation functions. These will throw an exception
// if exceptions are supported, otherwise they will print a mesg and
// exit. 

void OptppfatalError(const char* mesg)
{
#ifdef USEEXCEPTIONS
	throw OptppExceptions(mesg);
#else
//	cerr << "fatal error: " << mesg << endl;
	exit(1);
#endif
}

void OptppmemoryError(const char* mesg)
{
#ifdef USEEXCEPTIONS
	throw OptppMemoryError(mesg);
#else
//	cerr << "memory error: " << mesg << endl;
	exit(1);
#endif
}

void OptpprangeError(const char* mesg, int i, int low, int high)
{
#ifdef USEEXCEPTIONS
	throw OptppRangeError(mesg, i, low, high);
#else
//	cerr << "range error: " << mesg << " index=" << i << " bounds:["
//			 << low << ", " << high << "]" << endl;
	exit(1);
#endif
}

void OptppmathError(const char* mesg)
{
#ifdef USEEXCEPTIONS
	throw OptppMathError(mesg);
#else
//	cerr << "math error: " << mesg << endl;
	exit(1);
#endif
}

void OptppdomainError(const char* mesg, const double& badValue)
{
#ifdef USEEXCEPTIONS
	throw OptppDomainError(mesg, badValue);
#else
//	cerr << "domain error: " << mesg << " bad value = " << badValue << endl;
	exit(1);
#endif
}

void OptppzeroDivide(const char* mesg)
{
#ifdef USEEXCEPTIONS
	throw OptppZeroDivide(mesg);
#else
//	cerr << "division by zero error: " << mesg << endl;
	exit(1);
#endif
}

} // namespace OPTPP











