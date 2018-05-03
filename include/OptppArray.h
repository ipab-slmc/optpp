#ifndef OPTPPARRAY_H
#define OPTPPARRAY_H

#include <iostream>

#include "OptppFatalError.h"

using std::ostream;

namespace OPTPP {

/**
 * Simple array class. 
 *
 * Bounds checking is ON by default. To turn it off (for optimal performance)
 * do -DNOBOUNDSCHECK on the compilation command line.
 *
 *
 * Parallel support for primitive types relies on template specialization. 
 * If your compiler supports template specialization, define
 * TEMPLATE_SPECIALIZATION and use the code here. 
 * Otherwise, use the workaround routines in 
 * BadCompilerHacks.[cpp,h].
 * IRIX CC 7.2 and egcs both support specialization.
 * IRIX CC 7.1 appears not to (though the documentation says otherwise).
 */

template<class T>
class OptppArray {
public:
       /**
        * Default Constructor
        * @see OptppArray(int n)
        * @see OptppArray(int n, const T* cOptppArray)
        * @see OptppArray(int n, const T& t)
        */
	OptppArray();
        /**
         * @param n an integer argument
         */
	OptppArray(int n); 

        /**
         * @param n an integer argument
         * @param cOptppArray a pointer to class T
         */
	OptppArray(int n, const T* cOptppArray);

        /**
         * @param n an integer argument
         * @param t a reference to class T
         */
	OptppArray(int n, const T& t); 

        /// explicit copy destructor needed to prevent memory corruption
	~OptppArray();
        /// explicit copy constructor needed to prevent memory corruption
	OptppArray(const OptppArray<T>& other);  
        /// explicit assignment needed to prevent memory corruption
	const OptppArray<T>& operator=(const OptppArray<T>& other); 
  
	/// resize the OptppArray 
	void resize(int newN); 
	/// reserve n slots in OptppArray 
	void reserve(int n); 
	/// return numbered of reserved slots in OptppArray 
	int  reserve() const ;
	
	/// add a new entry. 
	OptppArray<T>& append(const T& rhs);

  /**	
   * simple accessors.
   * If NOBOUNDSCHECK is not set and a bounds error occurs, crash.
   * In cases where error handling is to be used, and exceptions are
   * not supported, use get and put instead.
   */
	T& operator[](int i);
	const T& operator[](int i) const;

  /**	
   * accessors that return error codes upon bounds violations (if
   * NOBOUNDSCHECK is not set):
   * @return 1 OK
   * @return 0 error
   */
	bool get(int i, T& value) const ;
	bool put(int i, const T& value) ;
	
   /**
    * @return Length of array
    */
	int length() const; 
  /**	
   * accessors that return error codes if there are problems
   * with parallel processing 
   * @return 1 OK
   * @return 0 error
   */
	bool bcast(int sender);
	bool send(int tag, int dest);
	bool recv(int tag, int src);
 private:
  T* data_;   	///< T class pointer to the data    
  int len_;   	///< Length of array   
  int reserve_;	///< Amount of reserved space 
  
  void indexCheckCrash(int i) const;
  bool indexCheckNoCrash(int i) const;
};

//template<class T> ostream& operator<<(ostream& os,const OptppArray<T>& array);
// print in form (), (1), or (1,2)
template<class T> ostream& operator<<(ostream& os, const OptppArray<T>& item)
{
  const OptppArray<T>& myItem = item;
  T myElem;
 

  os << "(";
  int ultimate = myItem.length()-1;        // the largest valid index
 
  for (int i = 0; i < ultimate; ++i) {
    myElem = myItem[i];
    os << myElem;
    os << ", ";
  }
  if (ultimate >= 0)
    os << myItem[ultimate];
  os << ")";
  return os;
}

inline int bump(int start, int to) {
  if (start == 0)
    start = 1;
  while (start < to)
    start *= 2;
  return start;
}

template<class T>
OptppArray<T> sliceOptppArray(const OptppArray<OptppArray<T> >& array, int index)
{
	OptppArray<T> rtn(array.length());

	for (int i=0; i<array.length(); i++)
		{
			rtn[i] = array[i][index];
		}
	return rtn;
}


// utilities for building small arrays

template<class T>
OptppArray<T> tuple(const T& a)
{
	OptppArray<T> rtn(1, a);
	return rtn;
}

template<class T>
OptppArray<T> tuple(const T& a, const T& b)
{
	OptppArray<T> rtn(2);
	rtn[0] = a;
	rtn[1] = b;
	return rtn;
}

template<class T>
OptppArray<T> tuple(const T& a, const T& b, const T& c)
{
	OptppArray<T> rtn(3);
	rtn[0] = a;
	rtn[1] = b;
	rtn[2] = c;
	return rtn;
}

template<class T>
OptppArray<T> tuple(const T& a, const T& b, const T& c, const T& d)
{
	OptppArray<T> rtn(4);
	rtn[0] = a;
	rtn[1] = b;
	rtn[2] = c;
	rtn[3] = d;
	return rtn;
}

template<class T>
OptppArray<T> tuple(const T& a, const T& b, const T& c, const T& d, const T& e)
{
	OptppArray<T> rtn(5);
	rtn[0] = a;
	rtn[1] = b;
	rtn[2] = c;
	rtn[3] = d;
	rtn[4] = e;
	return rtn;
}

template<class T>
OptppArray<T> tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
							 const T& f)
{
	OptppArray<T> rtn(6);
	rtn[0] = a;
	rtn[1] = b;
	rtn[2] = c;
	rtn[3] = d;
	rtn[4] = e;
	rtn[5] = f;
	return rtn;
}

template<class T> inline OptppArray<T>::OptppArray()
  : data_(0),
    len_(0), 
    reserve_(0)
{
}

template<class T> inline  OptppArray<T>::OptppArray(int n)
  : data_(0),
    len_(n), 
    reserve_(n)
{
  if (len_ < 0)
    OptppfatalError("Negative length passed to OptppArray<T>::OptppArray(int n)");
  if (len_ > 0) 
    {
      data_ = new T [len_];
      if (!data_)
	OptppmemoryError("OptppArray constructor out of memory");
    }
}

template<class T>  inline OptppArray<T>::OptppArray(int n, const T* cOptppArray)
  : data_(0),
    len_(n),
    reserve_(n)
{
  if (len_ < 0)
    OptppfatalError("Negative lenght passed to OptppArray::OptppArray(n, cOptppArray)");
  if (len_ > 0)
    {
      data_ = new T [len_];
      if (!data_)
	OptppmemoryError("OptppArray constructor out of memory");
    }
  for (int i = 0; i < len_; i++) data_[i] = cOptppArray[i];
}

template<class T>  inline OptppArray<T>::OptppArray(int n, const T& t)
  : data_(0),
    len_(n), 
    reserve_(n)
{
  if (len_ < 0)
    OptppfatalError("Negative length passed to OptppArray<T>::OptppArray(int n)");
  if (len_ > 0) 
    {
      data_ = new T [len_];
      if (!data_)
	OptppmemoryError("OptppArray constructor out of memory");
    }
  for (int i = 0; i < len_; ++ i) 
    data_[i] = t;
}


template<class T>  inline OptppArray<T>::OptppArray(const OptppArray<T>& arr)
  : data_(0),
    len_(arr.len_), 
    reserve_(arr.len_)
{
  if (len_ > 0) 
    {
      data_ = new T [reserve_];
      if (!data_)
	OptppmemoryError("OptppArray constructor out of memory");
      for (int i = 0; i < len_; ++i)
	data_[i] = arr.data_[i];
    }
}

template<class T> inline OptppArray<T>::~OptppArray() 
{
  delete [] data_;
}

template<class T> 
inline OptppArray<T>& OptppArray<T>::append(const T& rhs) 
{
  resize(len_+1);
  data_[len_-1] = rhs;
  return *this;
}

template<class T>  
inline const OptppArray<T>& OptppArray<T>::operator=(const OptppArray<T>& arr) 
{
  if (this != &arr) 
    { // don't bother to assign if they're already identical
      if (reserve_ < arr.len_) 
	{ //If the reserved space is too small to hold arr
	  delete [] data_;
	  data_ = 0;
	  reserve_ = arr.len_;
	  if (reserve_ > 0) {
	    data_ = new T[reserve_];
	    if (!data_)
	      OptppmemoryError("OptppArray constructor out of memory");
	  }
	}
      len_ = arr.len_;
      for (int i = 0; i < len_; ++i)
	data_[i] = arr[i];
    }
  return *this;
}

template<class T> 
inline void OptppArray<T>::resize(int newN) {
  if (len_ != newN) { // do not do anything if new size is not different
    if (newN < 0)
      OptppfatalError("Negative length passed to OptppArray<T>::resize(int newN)");
    if(newN > reserve_)
      reserve(bump(reserve_, newN));
    len_ = newN;
  }
}

template<class T>
inline int OptppArray<T>::length() const {
  return len_;
}

template<class T> 
inline void OptppArray<T>::reserve(int N){
  if(reserve_ != N){
    if(N < 0){
      OptppfatalError("Negative length passed to OptppArray<T>::reserve(int N)");
    }
    if(N < len_){ len_ = N;}
    reserve_ = N;
    T* oldData = data_;
    data_ = 0;
    data_ = new T [reserve_];
    if (!data_)
      OptppmemoryError("OptppArray<T>::reserve(int N) out of memory");
    for (int i = 0; i < len_; i++)
      data_[i] = oldData[i];
    delete [] oldData;
  }
}

template<class T> 
inline int OptppArray<T>::reserve() const{
  return reserve_;
}

template<class T> 
inline T& OptppArray<T>::operator[](int i) {
#ifndef NOBOUNDSCHECK
  indexCheckCrash(i);
#endif
  return data_[i];
}

template<class T> 
inline const T& OptppArray<T>::operator[](int i) const {
#ifndef NOBOUNDSCHECK
  indexCheckCrash(i);
#endif
  return data_[i];
}

template<class T> 
inline void OptppArray<T>::indexCheckCrash(int i) const 
{
  if (i<0 || i>=len_)
    OptpprangeError("OptppArray<T>", i, 0, len_-1);
}

template<class T> 
inline bool OptppArray<T>::indexCheckNoCrash(int i) const 
{
  if (i<0 || i>=len_) return false;
  return true;
}

template<class T>
inline bool OptppArray<T>::get(int i, T& value) const
{
#ifndef NOBOUNDSCHECK
  if (!indexCheckNoCrash(i)) return false;
#endif
  value = data_[i];
  return true;
}

template<class T>
inline bool OptppArray<T>::put(int i, const T& value)
{
#ifndef NOBOUNDSCHECK
  if (!indexCheckNoCrash(i)) return false;
#endif
  data_[i]= value;
  return true;
}

} // namespace OPTPP

#endif
