/**************************************************************************

	TABLE.H - Types for data tables (1D, 2D and 3D) 

	The contents of this file are subject to the Academic Development
	and Distribution License Version 1.0 (the "License"); you may not
	use this file except in compliance with the License.  You should
	have received a copy of the License with this file.  If you did not
	then please contact whoever distributed this file too you, since
	they may be in violation of the License, and this may affect your
	rights under the License.

	Software distributed under the License is distributed on an "AS IS"
	basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See
	the License for the specific language governing rights and
	limitations under the License.

	The Initial Developer of the Original Software is Jeremy Lea.

	Portions Copyright (C) 2016 OpenPave.org.

	Contributor(s): Jeremy Lea <reg@openpave.org>.

	Purpose:
		This header implements classes for storing data tables.  Each axis
		of the table is defined by an axis class that provides the sheet,
		column or row index.

	Design:
		These classes make use of many new C++11 features, including
		variadic templates and SFINAE.  This makes them possible, but also
		complex.
		
		The basic class is table.  Because the issues with variadic
		templates with multiple lists of classes in multiple axes is not
		easily addressed, there are three verions with duplicate code.  If
		you need more dimensions then you need to clone and add that...
		But remember that the axes can have many layers of keys, which is
		a much cleaner solution in many cases.
		
	History:
		2016/03/05 - Created a basic implementation.

**************************************************************************/

#ifndef __TABLE_H
#define __TABLE_H

#include <cstring>
#include <stdexcept>
#include <type_traits>
#include "axis.h"

/*
 * Class table
 *
 */
template<typename A, typename V>
class table : protected listener
{
public:
	table(A & ax, const unsigned b = DFLT_BLK)
	  : x(ax), len(0), buflen(0), blklen(b), buffer(nullptr) {
		unsigned s = ax.length();

		listen(ax,message<axis_message,unsigned>(
			[this](axis_message e, unsigned p){
				switch (e) {
				case axis_message::add:
					this->add(p);
					break;
				case axis_message::remove:
					this->remove(p);
					break;
				case axis_message::deleting:
					throw std::runtime_error("Out of order destructors!");
				default:
					break;
				}
		}));
		allocate(s); // Creates enough space
		copy(s,nullptr);   // Constructs the elements and increases size
	}
	~table() {
		deallocate();
	}
	// Behave like an array. Zero indexed.
	V & operator[] (const unsigned p) const {
#if defined(DEBUG)
		if (~inbounds(p))
			throw std::out_of_range("Index out of range!");
#endif
		return buffer[p];
	}

private:
	A & x;                     // The axis for this table
	unsigned len;              // The number of elements...
	unsigned buflen;           // The allocated buffer size...
	unsigned blklen;           // The minimum block size.
	struct _V {                // Placement new wrapper
		V _v;
		explicit _V() : _v() {}
		explicit _V(const V & v) : _v(v) {}
		operator V () const {
			return _v;
		}
		void * operator new(size_t, void * p) {
			return p;
		}
		void operator delete(void * , void *) {
		}
	};
	V * buffer;                // The actual storage.

	// Check if an index is within the set...
	inline bool inbounds(const unsigned p) const {
		return (p < len ? true : false);
	}
	// Calculate the buffer size.
	inline unsigned bufsize(unsigned s) {
		while (s > 8*blklen)
			blklen *= 8;
		return blklen*(s/blklen+(s%blklen?1:0));
	}
	void allocate(const unsigned s) {
		unsigned b = bufsize(s);
		if (b == buflen)
			return;
		if (b == 0) {
			free(buffer);
			buffer = nullptr;
			buflen = 0;
			return;
		}
		V * temp = static_cast<V *>(realloc(buffer,b*sizeof(V)));
		if (temp == nullptr)
			throw std::bad_alloc();
		buffer = temp;
		buflen = b;
	}
	void deallocate() {
		if (buffer != nullptr) {
			for (unsigned i = 0; i < len; i++)
				buffer[i].~V();
			free(buffer);
			buffer = nullptr;
		}
		len = 0;
		buflen = 0;
	}
	// POD constructor
	template<typename T>
	typename std::enable_if<std::is_pod<T>::value>::type
	init(const unsigned i, const T * v) {
		if (v)
			buffer[i] = *v;
	}
	// not POD
	template<typename T>
	typename std::enable_if<!std::is_pod<T>::value>::type
	init(const unsigned i, const T * v) {
		if (v)
			new(&buffer[i]) _V(*v);
		else
			new(&buffer[i]) _V();
	}
	void copy(const unsigned s, const V * v) {
		for (unsigned i = 0; i < s; i++)
			init(len++,(v ? &v[i] : nullptr));
	}
	// Add an element at position p.
	void add(unsigned p, const V * v = nullptr) {
		allocate(len+1);
		if (p < len)
			std::memmove(&buffer[p+1],&buffer[p],(len-p)*sizeof(V));
		len++;
		init(p,v);
	}
	// Remove element position p.
	void remove(unsigned p) {
		buffer[p].~V();
		if (p+1 < len)
			std::memmove(&buffer[p],&buffer[p+1],(len-p-1)*sizeof(V));
		len--;
		allocate(len);
	}
};
 
#endif // TABLE_H
