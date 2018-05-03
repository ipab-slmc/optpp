#ifndef OPTPPFATALERROR_H
#define OPTPPFATALERROR_H

#include "OptppExceptions.h"
/**
 * Portable error routines: under the hood, they throw an exception
 * on machines where exceptions are supported, or simply print a mesg and
 * then terminate on machines w/o exception support. 
 *
 * to switch compilation of exceptions on/off, toggle the macro USEEXCEPTIONS
 */

namespace OPTPP {

void OptppfatalError(const char* mesg);
void OptppmemoryError(const char* mesg);
void OptpprangeError(const char* mesg, int i, int low, int high);
void OptppmathError(const char* mesg);
void OptppdomainError(const char* mesg, const double& badValue);
void OptppzeroDivide(const char* mesg);

} // namespace OPTPP

#endif
