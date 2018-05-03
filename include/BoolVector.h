#ifndef bcv_h
#define bcv_h

#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#include <iostream>
#ifdef HAVE_STD
#include <cstdlib>
#else
#include <stdlib.h>
#endif

namespace OPTPP {

/**
 * BoolVector defines a vector of bools.
 * @author Charles Tong
 * @note Copy constructor and assignment operator added by P.J. Williams
 * @date Last modified 02/2006
 */
class BoolVector
{
  /// Length of vector
  int  size;	
  /// Pointer to a bool 
  bool *p;	

 public:
 /**
  * @param sz an integer argument
  * @see BoolVector(int sz, cost bool& val)
  * @see BoolVector(int sz, const BoolVector& val)
  */
  BoolVector(int sz);
 /**
  * @param sz an integer argument
  * @param val a bool 
  * @see BoolVector(int sz)
  * @see BoolVector(int sz, const BoolVector& val)
  */
  BoolVector(int sz, const bool & val); 
 /**
  * @param sz an integer argument
  * @param val a BoolVector 
  * @see BoolVector(int sz)
  * @see BoolVector(int sz, cost bool& val)
  */
  BoolVector(int sz, const BoolVector & val);
 /**
  * Copy Constructor
  */
  BoolVector(const BoolVector & val);

 /**
  * Default Constructor
  * @note Creates a bool vector of length 100
  */
  BoolVector(void  ) {p = new bool[size = 100];}
 /**
  * Destructor
  */
 ~BoolVector(void  ) {delete []p;}

  /** 
   * @return The length of the vector
   */
  int Size(void) {return size;}

  /**
   * @return Reference to the BoolVector element indexed by index
   */
  bool &operator()( int index ); 
  /**
   * @return Const reference to the BoolVector element indexed by index
   */
  const bool &operator()( int index ) const ; 

 /**
  * Assignment operator
  */
  BoolVector& operator=(const BoolVector& val); 
};

} // namespace OPTPP
#endif
