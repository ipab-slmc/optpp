#ifndef SMARTPTR_H
#define SMARTPTR_H

#include "OptppFatalError.h"

namespace OPTPP {

/**
 *	Howard Hinnant's reference counting handle class.
 *	
 *	This is to be used as a pointer to class T.  
 *	This will feel and smell just like a built-in pointer except:
 *	1.  There is no need to call delete on the pointer.
 *	2.  The default copy constructor and assignment implement ref-counting.
 *	3.  The user may call isNonUnique to determine if this pointer is
 *	    the only pointer to the data.  This can be used to hide the
 *	    ref-counting behavior of a class.
 *	4.  Checks for dereference of a null pointer.
 */


template <class T>
class SmartPtr {
public:
	inline SmartPtr(T* ptr = 0);
	inline SmartPtr(const SmartPtr<T>& other);
	inline ~SmartPtr();
	const SmartPtr<T>& operator =(const SmartPtr<T>& rhs);
	inline T* operator ->();
	inline bool isNull() const;
	inline const T* operator ->() const;
	inline T& operator *();
	inline const T& operator *() const;
	inline operator const T* () const;
	inline bool isNonUnique() const;

	inline int refCount() const;

protected:
	T* ptr_;
	int* refCount_;
};

template <class T> inline
SmartPtr<T>::SmartPtr(T* ptr)
	: ptr_(ptr),
	refCount_(0)
{
	if (ptr_) {
		refCount_ = new int;
		if (refCount_ == 0)
			OptppmemoryError("SmartPtr::SmartPtr(T* ptr) out of memory");
		*refCount_ = 1;
	}
}

template <class T> inline
SmartPtr<T>::SmartPtr(const SmartPtr<T>& other)
	: ptr_(other.ptr_),
	refCount_(other.refCount_)
{
	if (refCount_ != 0)
		++(*refCount_);
}

template <class T> 
inline SmartPtr<T>::~SmartPtr() {
	if (refCount_ != 0 && --(*refCount_) == 0) {
		delete ptr_;
		ptr_ = 0;
		delete refCount_;
		refCount_ = 0;
	}
}

template <class T> inline
bool SmartPtr<T>::isNull() const {
  return (ptr_ == 0);
}

template <class T> inline
T* SmartPtr<T>::operator ->() {
	if (ptr_ == 0)
		OptppfatalError("SmartPtr<T>::operator ->() on null pointer");
	return ptr_;
}

template <class T> inline
const T* SmartPtr<T>::operator ->() const {
	if (ptr_ == 0)
		OptppfatalError("SmartPtr<T>::operator ->() on null pointer");
	return ptr_;
}

template <class T> inline
T& SmartPtr<T>::operator *() {
	if (ptr_ == 0)
		OptppfatalError("SmartPtr<T>::operator *() on null pointer");
	return *ptr_;
}

template <class T> inline
const T& SmartPtr<T>::operator *() const {
	if (ptr_ == 0)
		OptppfatalError("SmartPtr<T>::operator *() on null pointer");
	return *ptr_;
}

template <class T> inline
SmartPtr<T>::operator const T* () const {
	return ptr_;
}

template <class T> inline
bool SmartPtr<T>::isNonUnique() const {
	return refCount_ == 0 ? false : *refCount_ != 1;
}

template <class T> inline
int SmartPtr<T>::refCount() const {
	return refCount_ == 0 ? 0 : *refCount_;
}

template <class T> 
inline const SmartPtr<T>& SmartPtr<T>::operator =(const SmartPtr<T>& rhs) {
	if (ptr_ != rhs.ptr_) {
		if (refCount_ != 0 && --(*refCount_) == 0) {
			delete ptr_;
			delete refCount_;
		}
		ptr_ = rhs.ptr_;
		refCount_ = rhs.refCount_;
		if (refCount_ != 0)
			++(*refCount_);
	}
	return *this;
}

} // namespace OPTPP
#endif
