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

	Portions Copyright (C) 2016-2022 OpenPave.org.

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

#pragma once
#ifndef __TABLE_H
#define __TABLE_H

#if !defined(DFLT_BLK)
#define DFLT_BLK    64
#endif

#include <cstring>
#include <stdexcept>
#include <type_traits>
#include "axis.h"

#if defined(DEBUG)
#define NOEXCEPT
#else
#define NOEXCEPT noexcept
#endif

namespace OP {

/*
 * Class table
 *
 */
template<typename V, typename...As>
class table : protected listener
{
public:
	table(As &...as)
	  : axes(as...), buflen(0), blklen(DFLT_BLK), buffer(nullptr) {
		const unsigned s = init(as...);
		allocate(s); // Creates enough space
		//copy(s,nullptr);   // Constructs the elements and increases size
	}
	~table() {
		allocate(0);
	}
	// Behave like an array. Zero indexed.
	V & operator [] (unsigned p) const NOEXCEPT {
#if defined(DEBUG)
		if (!inbounds(p))
			throw std::out_of_range("Index out of range!");
#endif
		return buffer[p];
	}
	// Behave like an array. Zero indexed.
	const V & operator () (const typename As::index_t...ps) const {
		return buffer[make_index(0,ps...)];
	}
	V & operator () (const typename As::index_t...ps) noexcept {
		return buffer[make_index(0,ps...)];
	}
	// Behave like an array using the axis keys
	const V & operator () (const typename As::key_t &...ks) const {
		return buffer[make_fromkey(0,ks...)];
	}
	V & operator () (const typename As::key_t &...ks) {
		return buffer[make_fromkey(0,ks...)];
	}
	// Return the size of array below dimension d.
	unsigned length(unsigned d = sizeof...(As)) const noexcept {
		unsigned s = 1;
		for (unsigned i = d; i > 0; i--)
			s *= sizes[i-1];
		return s;
	}

private:
	// The axes are stored in a tuple for ease of processing.
	std::tuple<const As &...> axes;
	unsigned sizes[sizeof...(As)];
	// 1D storage for the data
	unsigned buflen;           // The allocated buffer size...
	unsigned blklen;           // The minimum block size.
	struct _V {                // Placement new wrapper
		V _v;
		_V() : _v() {}
		explicit _V(const V & v) : _v(v) {}
		operator V () const {
			return _v;
		}
		void * operator new (size_t, void * p) noexcept {
			return p;
		}
		void operator delete (void * , void *) noexcept {
		}
	};
	V * buffer;                // The actual storage.

	template<typename...Ks>
	unsigned init(Ks &...) noexcept {
		return 1;
	}
	template<typename K, typename...Ks>
	unsigned init(K & ax, Ks &...ks) {
		constexpr const unsigned d = sizeof...(Ks);
		const unsigned s = ax.length();
		listen(ax,message<axis_message,unsigned>(
			[=](axis_message e, unsigned p){
				switch (e) {
				case axis_message::add:
					this->add(d,p);
					break;
				case axis_message::remove:
					this->remove(d,p);
					break;
				case axis_message::empty:
					this->empty();
					break;
				case axis_message::deleting:
					throw std::runtime_error("Out of order destructors!");
				}
		}));
		sizes[d] = s;
		return s*init(ks...);
	}
	template<typename...Is>
	unsigned make_index(unsigned p, Is...) const noexcept {
		return p;
	}
	template<typename I, typename...Is>
	unsigned make_index(unsigned p, I i, Is...is) const noexcept {
		return make_index(p*sizes[sizeof...(Is)]+i,is...);
	}
	template<typename...Ks>
	unsigned make_fromkey(unsigned p, Ks...) const noexcept {
		return p;
	}
	template<typename K, typename...Ks>
	unsigned make_fromkey(unsigned p, K k, Ks...ks) const {
		const auto & a = std::get<sizeof...(As)-sizeof...(Ks)-1>(axes);
		const unsigned i = a.getorderof(k);
		return make_fromkey(p*sizes[sizeof...(Ks)]+i,ks...);
	}
	// Return the step size for dimension d.
	unsigned step(unsigned d) const noexcept {
		unsigned s = 1;
		for (unsigned i = d+1; i < sizeof...(As); i++)
			s *= sizes[i];
		return s;
	}
	// Check if an index is within the set...
	bool inbounds(unsigned p) const noexcept {
		return (p < length() ? true : false);
	}
	// Calculate the buffer size.
	unsigned bufsize(unsigned s) noexcept {
		while (s > 8*blklen)
			blklen *= 8;
		return blklen*(s/blklen+((s%blklen)?1:0));
	}
	void allocate(unsigned s) {
		const unsigned b = bufsize(s);
		if (b == buflen)
			return;
		if (b == 0) {
			for (unsigned i = 0; i < length(); i++)
				buffer[i].~V();
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
	void empty() {
		allocate(0);
		for (unsigned i = 0; i < sizeof...(As); i++)
			sizes[i] = 0;
	}
	// POD constructor
	template<typename T>
	typename std::enable_if<std::is_pod<T>::value>::type
	initelem(unsigned i, const T * v) noexcept {
		if (v)
			buffer[i] = *v;
	}
	// not POD
	template<typename T>
	typename std::enable_if<!std::is_pod<T>::value>::type
	initelem(unsigned i, const T * v) {
		if (v)
			new(&buffer[i]) _V(*v);
		else
			new(&buffer[i]) _V();
	}
	//void copy(unsigned s, const V * v) {
	//	for (unsigned i = 0; i < s; i++)
	//		initelem(len++,(v ? &v[i] : nullptr));
	//}
	// Add an element at position p.
	void add(unsigned d, unsigned p, const V * v = nullptr) {
		const unsigned bs = length(d);
		const unsigned ss = step(d);
		const unsigned os = sizes[d];
		sizes[d]++;
		allocate(length());
		unsigned nl = length()/sizes[d];
		for (unsigned i = ss; i > 0; i--) {
			const unsigned ob = ((i-1)*os+p)*bs;
			const unsigned nb = ((i-1)*(os+1)+p+1)*bs;
			const unsigned ln = (i == ss ? os-p : os);
			std::memmove(&buffer[nb],&buffer[ob],ln*bs*sizeof(V));
			for (unsigned j = 0; j < bs; j++)
				initelem(nb-bs+j,v ? &v[--nl] : nullptr);
		}
	}
	// Remove element position p.
	void remove(unsigned d, unsigned p) {
		const unsigned bs = length(d);
		const unsigned ss = step(d);
		const unsigned os = sizes[d];
		sizes[d]--;
		for (unsigned i = 0; i < ss; i++) {
			const unsigned ob = (i*os+p+1)*bs;
			const unsigned nb = (i*(os-1)+p)*bs;
			const unsigned ln = (i == ss-1 ? os-p-1 : os);
			for (unsigned j = 0; j < bs; j++)
				buffer[nb+j].~V();
			std::memmove(&buffer[nb],&buffer[ob],ln*bs*sizeof(V));
		}
		allocate(length());
	}
};

} // namespace OP

#endif // TABLE_H
