#ifndef OPTPPEXCEPTIONS_H
#define OPTPPEXCEPTIONS_H

namespace OPTPP {

/**
 * OptppExceptions is the base class for OptppMemoryError, OptppRangeError,
 * OptppMathError, OptppDomainError, OptppZeroDivide.
 *
 * @note  Modified by P.J. Williams (Nomenclature Changes),
 * Sandia National Laboratories, pwillia@sandia.gov
 * @date  02/2006
 */ 


// CPJW
// Changed the class name from Exception to Exceptions to prevent clash with 
// class Exception previously defined in newmat08/myexcepts.h 
// Changed the class name from Exceptions to OptppExceptions to prevent clash 
// with class CSRMException previously defined in the DDace/CPPUtilities.h 
// 

class OptppExceptions
{
 public:
	OptppExceptions() {;}
	OptppExceptions(const char* mesg);

	virtual void print() const ;

	virtual const char* message() const {return mesg_;}

 protected:
	char mesg_[1000];
};


class OptppMemoryError : public OptppExceptions
{
 public:
	OptppMemoryError(const char* mesg);

	virtual void print() const ;
	
 private:
};

class OptppRangeError : public OptppExceptions
{
 public:
	OptppRangeError(const char* mesg, int i, int low, int high);

	virtual void print() const ;
	
 private:
	int i_;
	int low_;
	int high_;
};


class OptppMathError : public OptppExceptions
{
 public:
	OptppMathError(const char* mesg);

	virtual void print() const ;
	
 private:
};

class OptppDomainError : public OptppMathError
{
 public:
	OptppDomainError(const char* mesg, const double& badValue);

	virtual void print() const ;
	
 private:
	double badValue_;
};

class OptppZeroDivide : public OptppMathError 
{
 public:
	OptppZeroDivide(const char* mesg);

	virtual void print() const ;
	
 private:
};

class RecoverableOptppExceptions : public OptppExceptions
{
 public:
	RecoverableOptppExceptions(const OptppExceptions& e);

	virtual void print() const ;
 private:
	OptppExceptions e_;
};

void bailout(const OptppExceptions& e);

} // namespace OPTPP

#endif
