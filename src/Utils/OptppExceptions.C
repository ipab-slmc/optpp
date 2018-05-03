
#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#include <iostream>
#ifdef HAVE_STD
#include <cstring>
#include <cstdlib>
#else
#include <string.h>
#include <stdlib.h>
#endif

#include "OptppExceptions.h"

using namespace std;

namespace OPTPP {

// general bailout function

void bailout(const OptppExceptions& e)
{
	e.print();
//	cerr << "bailing out... " << endl;
	exit(1);
}

// base OptppExceptions class

OptppExceptions::OptppExceptions(const char *mesg)
{
	strcpy(mesg_, mesg);
}

void OptppExceptions::print() const 
{
//	cerr << "Unspecified exception detected: " << mesg_ << endl;
}


// memory errors

OptppMemoryError::OptppMemoryError(const char* mesg)
	: OptppExceptions(mesg)
{;}

void OptppMemoryError::print() const 
{
//	cerr << "Memory exception detected: " << mesg_ << endl;
}	


// bounds error

OptppRangeError::OptppRangeError(const char* mesg, int i, int low, int high)
	: OptppExceptions(mesg), i_(i), low_(low), high_(high)
{;}

void OptppRangeError::print() const 
{
//	cerr << "range exception: " << mesg_ << " index=" << i_ << " bounds:["
//			 << low_ << ", " << high_ << "]" << endl;
}

// math errors

// general math error

OptppMathError::OptppMathError(const char* mesg)
	: OptppExceptions(mesg)
{;}

void OptppMathError::print() const 
{
//	cerr << "Math exception detected: " << mesg_ << endl;
}	

// domain error (e.g., sqrt(-1), log(0))

OptppDomainError::OptppDomainError(const char* mesg, const double& badValue)
	: OptppMathError(mesg), badValue_(badValue)
{;}

void OptppDomainError::print() const 
{
//	cerr << "Math domain error detected: " << mesg_
//			 << " bad value = " << badValue_ << endl;
}	

// division by zero

OptppZeroDivide::OptppZeroDivide(const char* mesg)
	: OptppMathError(mesg)
{;}

void OptppZeroDivide::print() const 
{
//	cerr << "Division by zero exception detected: " << mesg_ << endl;
}	


// recoverable exception

RecoverableOptppExceptions::RecoverableOptppExceptions(const OptppExceptions& e)
	: OptppExceptions(), e_(e)
{
	;
}

void RecoverableOptppExceptions::print() const 
{
//	cerr << "recoverable exception detected: " << endl;
	e_.print();
}	


} // namespace OPTPP
