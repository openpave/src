/**************************************************************************

	SET.H - Simple set template classes

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

	Portions Copyright (C) 2006-2022 OpenPave.org.

	Contributor(s): Jeremy Lea <reg@openpave.org>.

	Purpose:
		This header implements a number of templated C++ classes for
		various types of sets. The code is fairly compact, but not
		optimised for big sets (particularly sorting).

	Design:
		There are two types of sets. The first are normal sets with
		one variable. Then there are associative sets, which store
		a key and a value.

		The hierarchy goes from a fixed size set, to a resizeable set,
		to a sorted resizeable set.

		All of the templated sets take a minimum block size.  Major
		performance enhancements will be seen on resizeable sets if this
		value is bigger than the maximum size of the set, because the
		memory is reused.

	Status:
		The swaps in the sorts should be replaced by memmove.

	History:
		1994       - Created by Jeremy Lea <reg@openpave.org>
		2002/01/23 - Modifications for use on a modern C++ compiler
		2007/09/11 - Modified to use realloc and placement new.
		2008/02/11 - Removed the idea of a default element.

**************************************************************************/

#ifndef __SET_H
#define __SET_H

#if !defined(DFLT_BLK)
#define DFLT_BLK    64
#endif

#include <algorithm>
#include <climits>
#include <cstring>
#include <new>
#include <stdexcept>

namespace OP {

/*
 * class set - Set base class
 *
 * This class provides the basic working for the rest of the set classes,
 * including the calculation of the buffer size.
 */
class set {
public:
	// The length. Nice for lots of things...
	unsigned length() const noexcept {
		return size;
	}
	// Check if an index is within the set...
	bool inbounds(unsigned p) const noexcept {
		return (p < size ? true : false);
	}

protected:
	unsigned size;             // The size of the set...
	unsigned buffer;           // The allocated buffer size...

	// Simple constructor...
	explicit set(unsigned b = DFLT_BLK) noexcept
	  : size(0), buffer(0), block(b > 1 ? b : 1) {
	}
	set(const set & s) noexcept
	  : size(0), buffer(0), block(s.block) {
	}
	set & operator = (const set & s) noexcept {
		size = 0;
		buffer = 0;
		block = s.block;
		return *this;
	}
	~set() {
	}
	// Calculate the buffer size.
	unsigned bufsize(unsigned s) noexcept {
		while (s > 8*block)
			block *= 8;
		//while (64*s < block)
		//	block /= 8;
		return block*(s/block+(s%block?1:0));
	}

private:
	unsigned block;            // The minimum block size.
};

/*
 * class fset - Fixed size set of type V.
 *
 * This class acts like a fixed size array.
 */
template <class V>
class fset : public set {
public:
	// Nice simple constructor...
	fset(unsigned s, unsigned b)
	  : set(b), value(nullptr) {
		allocate(s); // Creates enough space
		copy(s,nullptr);   // Constructs the elements and increases size
	}
	explicit fset(unsigned s, const V * v = nullptr, unsigned b = DFLT_BLK)
	  : set(b), value(nullptr) {
		allocate(s);
		copy(s,v);
	}
	// Copy constructor.
	fset(const fset<V> & v)
	  : set(v), value(nullptr) {
		allocate(v.size);
		copy(v.size,v.value);
	}
	// Wow, a destructor...
	~fset() {
		deallocate();
	}

	// Assignment operator.
	fset<V> & operator = (const fset<V> & v) {
		deallocate();
		allocate(v.size);
		copy(v.size,v.value);
		return *this;
	}
	// Allow the size to be changed, even for 'fixed' sets.
	void resize(unsigned s) {
		deallocate();
		allocate(s);
		copy(s,nullptr);
	}
	// Behave like an array. Zero indexed.
	V & operator [] (unsigned p) const {
#if defined(DEBUG)
		if (!inbounds(p))
			throw std::out_of_range("Index out of range!");
#endif
		return value[p];
	}
	// Flatten the set into a traditional array.  Caller is responsible
	// for correct size.
	void copyout(V * v) const noexcept {
		for (unsigned i = 0; i < size; i++)
			v[i] = value[i];
	}

protected:
	V * value;                 // The buffer.
	struct _V {                // Placement new wrapper
		V _v;
		_V() noexcept : _v() {}
		explicit _V(const V & v) : _v(v) {}
		void * operator new (size_t, void * p) noexcept {
			return p;
		}
		void operator delete (void * , void *) noexcept {
		}
	};

	// Since this is a fixed size set, hide the allocator...
	void allocate(unsigned s) {
		unsigned b = bufsize(s);
		if (b == buffer)
			return;
		if (b == 0) {
			free(value);
			value = nullptr;
			buffer = 0;
			return;
		}
		V * temp = static_cast<V *>(realloc(value,b*sizeof(V)));
		if (temp == nullptr)
			throw std::bad_alloc();
		value = temp;
		buffer = b;
	}
	void deallocate() noexcept {
		if (value) {
			for (unsigned i = 0; i < size; i++)
				value[i].~V();
			free(value);
			value = nullptr;
		}
		size = 0;
		buffer = 0;
	}
	void init(unsigned i, const V * v) {
		if (v != nullptr)
			new(&value[i]) _V(*v);
		else
			new(&value[i]) _V();
	}
	void copy(unsigned s, const V * v) {
		for (unsigned i = 0; i < s; i++)
			init(size++,(v ? &v[i] : nullptr));
	}
	// Also hide the null constructor.
	fset()
	  : set(), value(nullptr) {
		allocate(size);
	}
};

/*
 * Special case these to avoid constructors
 */
template<>
inline void
fset<int>::init(unsigned i, const int * v) {
	if (v)
		value[i] = *v;
}
template<>
inline void
fset<unsigned>::init(unsigned i, const unsigned * v) {
	if (v)
		value[i] = *v;
}
template<>
inline void
fset<float>::init(unsigned i, const float * v) {
	if (v)
		value[i] = *v;
}
template<>
inline void
fset<double>::init(unsigned i, const double * v) {
	if (v)
		value[i] = *v;
}

/*
 * class sset - Sizeable set of type V.
 *
 * This class builds on fset to make a set which can be resized
 * at will. Lots of nice add() and remove() options are available.
 */
template <class V>
class sset : public fset<V> {
public:
	// Provide a null constructor for empty sets.
	sset() noexcept
	  : fset<V>() {
	}
	sset(unsigned s, unsigned b)
	  : fset<V>(s,b) {
	}
	// Simple constructor.
	explicit sset(unsigned s, const V * v = nullptr, unsigned b = DFLT_BLK)
	  : fset<V>(s,v,b) {
	}
	// Up conversion constructor.
	explicit sset(const fset<V> & v)
	  : fset<V>(v) {
	}
	// Copy constructor.
	sset(const sset<V> & v)
	  : fset<V>(v) {
	}
	// We let someone else clean up...
	~sset() {
	}

	// Assignment operator.
	sset<V> & operator = (const sset<V> & v) {
		fset<V>::operator=(v);
		return *this;
	}

	// Add one element, at the end.
	void add(const V & v) {
		add(this->size,&v,1);
	}
	// Add a whole set, at the end.
	void add(const fset<V> & v) {
		add(this->size,&(v[0]),v.length());
	}
	// Add an array, at the end.
	void add(const V * v, unsigned s = 1) {
		add(this->size,v,s);
	}
	// Insert one element at position p.
	void add(unsigned p, const V & v) {
		add(p,&v,1);
	}
	// Insert a set at position p.
	void add(unsigned p, const fset<V> & v) {
		add(p,&(v[0]),v.length());
	}
	// Add an array at position p. (Actually do the work too).
	void add(unsigned p, const V * v, unsigned s = 1) {
		unsigned i;
		if (s == 0)
			throw std::invalid_argument("Cannot add zero elements into a set!");
		if (p >= this->size) {
			this->allocate(p+s);
			for (i = this->size; p > 0 && i < p-1; i++)
				this->init(i,nullptr);
			this->size = p+s;
		} else {
			this->allocate(this->size+s);
			memmove(&this->value[p+s],&this->value[p],
					(this->size-p)*sizeof(V));
			this->size += s;
		}
		for (i = 0; i < s; i++)
			this->init(p+i,v ? &v[i] : nullptr);
	}
	// Remove the last element.
	void remove() {
		remove(this->size-1,1);
	}
	// Remove s elements, starting at position p.
	void remove(unsigned p, unsigned s = 1) {
		if (!this->inbounds(p))
			throw std::out_of_range("Cannot remove out of range element from set!");
		if (s == 0)
			throw std::invalid_argument("Cannot remove zero elements from set!");
		s = (p+s > this->size ? this->size-p : s);
		for (unsigned i = 0; i < s; i++)
			this->value[p+i].~V();
		if (this->size-p > s)
			memmove(&this->value[p],&this->value[p+s],
				(this->size-p-s)*sizeof(V));
		this->size -= s;
		this->allocate(this->size);
	}
	void empty() noexcept {
		this->deallocate();
		this->allocate(0);
	}
};

/*
 * class oset - An (lazy) ordered sizeable set.
 *
 * This is a sizeable set which can be sorted. Because we don't limit
 * the sets to read only, this has to be lazy, so you need to call sort
 * before you expect the set to be sorted.
 *
 * If you want to store your own classes, then define operator's >,
 * <=, and =.
 */
template <class V>
class oset : public sset<V> {
public:
	// Null constructor.
	oset()
	  : sset<V>() {
	}
	oset(unsigned s, unsigned b)
	  : sset<V>(s,b) {
	}
	// Simple constructor.
	explicit oset(unsigned s, const V * v = nullptr, unsigned b = DFLT_BLK)
	  : sset<V>(s,v,b) {
		if (v)
			sort();
	}
	// Copy constructor.
	explicit oset(const fset<V> & v)
	  : sset<V>(v) {
		sort();
	}
	// Copy constructor.
	oset(const oset<V> & v)
	  : sset<V>(v) {
		sort();
	}
	// Destructor.
	~oset() {
	}

	// Guess what?
	void sort() {
		qsort(0,this->size);
	}
	// Do a lookup, and return -1 if the value is not found.
	// You must sort the set first!
	unsigned findvalue(const V & v) const noexcept {
		unsigned l = 0, r = this->size;
		while (l < r) {
			unsigned i = l + (r-l)/2;
			if (this->value[i] < v)
				l = i+1;
			else
				r = i;
		}
		if (l < this->size && this->value[l] == v)
			return l;
		else
			return UINT_MAX;
	}

protected:
	// A highly optimised quick sort. Don't touch...
	void qsort(unsigned l, unsigned r) {
		if (r <= l)
			return;
		unsigned i, j, k, p = l+(r-1-l)/2;
		if (this->value[l] > this->value[p])
			std::swap(this->value[l],this->value[p]);
		if (this->value[l] > this->value[r-1])
			std::swap(this->value[l],this->value[r-1]);
		if (this->value[p] > this->value[r-1])
			std::swap(this->value[p],this->value[r-1]);
		if (r-1-l <= 2)
			return;
		for (i = l, j = r-1, k = 0; ;
				p = (p==i?j++:(p==j?i--:p)), k++) {
			while (++i < p && !(this->value[i] > this->value[p])) {};
			while (p < --j && !(this->value[p] > this->value[j])) {};
			if (i >= j)
				break;
			std::swap(this->value[i],this->value[j]);
		}
		if (k == 0) {
			isort(l,p);
			isort(p+1,r);
		} else {
			qsort(l,p);
			qsort(p+1,r);
		}
	}
	// Insertion sort for if the set looks sorted already.
	void isort(unsigned l, unsigned r) noexcept {
		for (unsigned i = l; l < r && i < r-1; i++) {
			for (unsigned j = i+1; j > l
					&& this->value[j-1] > this->value[j]; j--)
				std::swap(this->value[j-1],this->value[j]);
		}
	}
};

/*
 * class cset - An compact ordered sizeable set.
 *
 * Build on the ordered set to make a set which take unique
 * values only. This is also lazy.
 *
 * Custom classes should define operator!=.
 */
template <class V>
class cset : public oset<V> {
public:
	// Getting the hang of this yet?
	cset()
	  : oset<V>() {
	}
	cset(unsigned s, unsigned b)
	  : oset<V>(s,b) {
	}
	explicit cset(unsigned s, const V * v = nullptr, unsigned b = DFLT_BLK)
	  : oset<V>(s,v,b) {
		if (v)
			compact();
	}
	explicit cset(const fset<V> & v)
	  : oset<V>(v) {
		compact();
	}
	cset(const cset<V> & v)
	  : oset<V>(v) {
		compact();
	}
	~cset() {
	}

	// Sort then compact the set.
	void sort() {
		this->qsort(0,this->size);
		compact();
	}

protected:
	void compact() {
		unsigned i, j, s;
		for (i = 1, j = 0; this->size > 1 && i < this->size; i++) {
			if (this->value[j] != this->value[i])
				j++;
			if (j != i) {
				do
					this->value[i].~V();
				while (++i < this->size
				 && this->value[j] == this->value[i]);
				s = i;
				while (s < this->size-1
				 && this->value[s] != this->value[s+1])
					s++;
				if (i < this->size) {
					memmove(&this->value[j+1],&this->value[i],
						(s-i+1)*sizeof(V));
					j = s-(i-(j+1)); i = s;
				}
			}
		}
		this->size = j+1;
		this->allocate(this->size);
	}
};

/*
 * class iset - Indexed set
 *
 * This is like a sset, expect that sorting only sorts an index.
 * As a result hasvalue() is much faster, but at the expense of having
 * the set be read only...
 */
template <class V>
class iset : public set {
public:
	// Make one...
	iset()
	  : set(), idx(nullptr), value(nullptr) {
		allocate(size);
		copy(size,nullptr,true);
	}
	iset(unsigned s, unsigned b)
	  : set(b), idx(nullptr), value(nullptr) {
		allocate(s);
		copy(s,nullptr,true);
	}
	// Basic constructor
	explicit iset(unsigned s, const V * v = nullptr, unsigned b = DFLT_BLK)
	  : set(b), idx(nullptr), value(nullptr) {
		allocate(s);
		copy(s,v,true);
		allocate(size);
	}
	// Copy constructor.
	iset(const iset<V> & v)
	  : set(v), idx(nullptr), value(nullptr) {
		allocate(v.size);
		copy(v.size,v.value,false);
	}
	// Copy from an fset.
	explicit iset(const fset<V> & v)
	  : set(v), idx(nullptr), value(nullptr) {
		allocate(v.size);
		copy(v.size,v.value,true);
		allocate(size);
	}
	// Clean up.
	~iset() {
		deallocate();
	}

	// Do a lookup, and return -1 if the value is not found.
	unsigned hasvalue(const V & v) const {
		unsigned p = findvalue(v);
		if (p < size && value[idx[p]] == v)
			return idx[p];
		else
			return UINT_MAX;
	}
	// Assignment operator.
	iset<V> & operator = (const iset<V> & v) {
		deallocate();
		allocate(v.size);
		copy(v.size,v.value,false);
		return *this;
	}
	iset<V> & operator = (const fset<V> & v) {
		deallocate();
		allocate(v.size);
		copy(v.size,v.value,true);
		allocate(size);
		return *this;
	}
	// Only integer keys make sense.
	const V & operator [] (unsigned p) const {
#if defined(DEBUG)
		if (!inbounds(p))
			throw std::out_of_range("Index out of range!");
#endif
		return value[p];
	}
	// Allow sorted access.
	const V & getindex(unsigned i) const {
#if defined(DEBUG)
		if (!inbounds(i))
			throw std::out_of_range("Index out of range!");
#endif
		return value[idx[i]];
	}
	// Get the position of an element in the sort.
	unsigned getorder(unsigned i) const {
#if defined(DEBUG)
		if (!inbounds(i))
			throw std::out_of_range("Index out of range!");
#endif
		return idx[i];
	}
	// Flatten the set into a traditional array.  Caller is responsible
	// for correct size.
	void copyout(V * v) const {
		for (unsigned i = 0; i < size; i++)
			v[i] = value[i];
	}


	// Add one value at the end.
	void add(const V & v) {
		return add(&v,1);
	}
	// Add a whole set, at the end.
	void add(const iset<V> & v) {
		return add(v.value,v.size);
	}
	void add(const fset<V> & v) {
		return add(v.value,v.size);
	}
	// Add an array of values... There is no point in
	// a position based addition.
	void add(const V * v, unsigned s = 1) {
		if (s == 0)
			throw std::invalid_argument("Cannot add zero elements into a set!");
		if (v == nullptr)
			throw std::invalid_argument("Cannot add a null to a set!");
		allocate(size+s);
		copy(s,v,true);
		allocate(size);
	}
	// Now start removing them...
	void remove(const V & v) {
		unsigned p = hasvalue(v), q = UINT_MAX;
		if (p == UINT_MAX)
			throw std::invalid_argument("Cannot remove value not in set!");
		value[p].~V();
		if (p < --size)
			memmove(&value[p],&value[p+1],(size-p)*sizeof(V));
		for (unsigned i = 0; i <= size; i++) {
			if (idx[i] == p)
				q = i;
			else if (idx[i] > p)
				idx[i]--;
		}
		if (q < size)
			memmove(&idx[q],&idx[q+1],(size-q)*sizeof(unsigned));
		allocate(size);
	}
	void empty() noexcept {
		deallocate();
		allocate(0);
	}

protected:
	unsigned * idx;            // The index.
	V * value;                 // Take a guess...
	struct _V {                // Placement new wrapper
		V _v;
		explicit _V(const V & v) : _v(v) {}
		void * operator new (size_t, void * p) {
			return p;
		}
		void operator delete (void * , void *) {
		}
	};

	// Make some space...
	void allocate(unsigned s) {
		unsigned b = bufsize(s);
		if (b == buffer)
			return;
		if (b == 0) {
			free(idx);
			idx = nullptr;
			free(value);
			value = nullptr;
			buffer = 0;
			return;
		}
		unsigned * itemp = static_cast<unsigned *>(
				realloc(idx,b*sizeof(unsigned)));
		V * vtemp = static_cast<V *>(realloc(value,b*sizeof(V)));
		if (itemp == nullptr || vtemp == nullptr)
			throw std::bad_alloc();
		idx = itemp;
		value = vtemp;
		buffer = b;
	}
	void deallocate() noexcept {
		if (idx) {
			free(idx);
			idx = nullptr;
		}
		if (value) {
			for (unsigned i = 0; i < size; i++)
				value[i].~V();
			free(value);
			value = nullptr;
		}
		size = 0;
		buffer = 0;
	}
	// Find the position which is either equal or greater...
	unsigned findvalue(const V & v) const {
		unsigned l = 0, r = size;
		while (l < r) {
			unsigned i = l + (r-l)/2;
			if (value[idx[i]] < v)
				l = i+1;
			else
				r = i;
		}
		return l;
	}
	void init(const V * v) {
		if (v == nullptr)
			throw std::invalid_argument("Cannot insert null element into set!");
		new(&value[size]) _V(*v);
		unsigned p = findvalue(*v);
		if (p < size)
			memmove(&idx[p+1],&idx[p],(size-p)*sizeof(unsigned));
		idx[p] = size++;
	}
	void copy(unsigned s, const V * v, bool checkdups) {
		for (unsigned i = 0, p; v && i < s; i++) {
			if (checkdups && (p = hasvalue(v[i])) != UINT_MAX)
				value[p] = v[i];
			else
				init(&v[i]);
		}
	}
};

/*
 * class kfset - Keyed fixed set.
 *
 * This implements a set where the values contain a base
 * class which acts as a key. The values and keys are the same
 * data, just treated as values or keys depending on the role.
 * The key must be unique. This is intended for basic database
 * like functionality.
 */
template <class K, class V>
class kfset : public set {
public:
	// Simple constructor.
	kfset(unsigned s, unsigned b)
	  : set(b), value(nullptr) {
		allocate(s);
		copy(s,nullptr,true);
	}
	explicit kfset(unsigned s, const V * v = nullptr, unsigned b = DFLT_BLK)
	  : set(b), value(nullptr) {
		allocate(s);
		copy(s,v,true);
		allocate(size);
	}
	// Copy constructor.
	kfset(const kfset<K,V> & v)
	  : set(v), value(nullptr) {
		allocate(v.size);
		copy(v.size,v.value,false);
	}
	// Clean up.
	~kfset() {
		deallocate();
	}

	// Do a key lookup, and return UINT_MAX if the key is not found.
	unsigned haskey(const K & k) const noexcept {
		for (unsigned i = 0; i < size; i++) {
			if (static_cast<K &>(value[i]) == k)
				return i;
		}
		return UINT_MAX;
	}
	// Assignment operator.
	kfset<K,V> & operator = (const kfset<K,V> & v) {
		deallocate();
		allocate(v.size);
		copy(v.size,v.value,false);
		return *this;
	}
	// Return data based on a key lookup.
	V & operator [] (const K & k) const {
		unsigned p = haskey(k);
#if defined(DEBUG)
		if (!inbounds(p))
			throw std::out_of_range("Key out of range!");
#endif
		return value[p];
	}
	// Allow integer keys.
	V & operator [] (unsigned p) const {
#if defined(DEBUG)
		if (!inbounds(p))
			throw std::out_of_range("Index out of range!");
#endif
		return value[p];
	}
	// Flatten the set into a traditional array.  Caller is responsible
	// for correct size.
	void copyout(V * v) const {
		for (unsigned i = 0; i < size; i++)
			v[i] = value[i];
	}


protected:
	V * value;                 // The data.
	struct _V {                // Placement new wrapper
		V _v;
		explicit _V(const V & v) : _v(v) {}
		void * operator new (size_t, void * p) noexcept {
			return p;
		}
		void operator delete (void * , void *) noexcept {
		}
	};

	// Hide the allocation function.
	void allocate(unsigned s) {
		unsigned b = bufsize(s);
		if (b == buffer)
			return;
		if (b == 0) {
			free(value);
			value = nullptr;
			buffer = 0;
			return;
		}
		V * temp = static_cast<V *>(realloc(value,b*sizeof(V)));
		if (temp == nullptr)
			throw std::bad_alloc();
		value = temp;
		buffer = b;
	}
	void deallocate() noexcept {
		if (value) {
			for (unsigned i = 0; i < size; i++)
				value[i].~V();
			free(value);
			value = nullptr;
		}
		size = 0;
		buffer = 0;
	}
	void init(const V * v) {
		if (v == nullptr)
			throw std::invalid_argument("Cannot insert null element into set!");
		new(&value[size++]) _V(*v);
	}
	void copy(unsigned s, const V * v, bool checkdups) {
		for (unsigned i = 0, p; v && i < s; i++) {
			if (checkdups && (p = haskey(v[i])) != UINT_MAX)
				value[p] = v[i];
			else
				init(&v[i]);
		}
	}
	// And the null constructor.
	kfset()
	  : set(), value(nullptr) {
		allocate(size);
	}
};

/*
 * class ksset - Keyed sizeable set
 *
 * Think combination of cset and kfset.
 */
template <class K, class V>
class ksset : public kfset<K,V> {
public:
	// C++ sucks...
	ksset()
	  : kfset<K,V>() {
	}
	ksset(unsigned s, unsigned b)
	  : kfset<K,V>(s,b) {
	}
	explicit ksset(unsigned s, const V * v = nullptr, unsigned b = DFLT_BLK)
	  : kfset<K,V>(s,v,b) {
	}
	explicit ksset(const kfset<K,V> & v)
	  : kfset<K,V>(v) {
	}
	ksset(const ksset<K,V> & v)
	  : kfset<K,V>(v) {
	}
	ksset<K,V> & operator = (const ksset<K,V> & v) {
		kfset<K,V>::operator=(v);
		return *this;
	}
	~ksset() {
	}

	// Add one value at the end.
	void add(const V & v) {
		add(&v,1);
	}
	// Add a whole set, at the end.
	void add(const kfset<K,V> & v) {
		add(v.value,v.size);
	}
	// Add an array of values... There is no point in
	// a position based addition.
	void add(const V * v, unsigned s = 1) {
		if (s == 0)
			throw std::invalid_argument("Cannot add zero elements to a set!");
		if (v == nullptr)
			throw std::invalid_argument("Cannot add a null to a set!");
		this->allocate(this->size+s);
		this->copy(s,v,true);
		this->allocate(this->size);
	}
	// Remove based on key.
	void remove(const K & k) {
		unsigned p = this->haskey(k);
		if (p == UINT_MAX)
			throw std::invalid_argument("Cannot remove key not in set!");
		this->value[p].~V();
		if (p < --this->size)
			memmove(&this->value[p],&this->value[p+1],
				(this->size-p)*sizeof(V));
		this->allocate(this->size);
	}
	// Replace key/value with another.
	void replace(const V & v) {
		unsigned p = this->haskey(v);
		if (p == UINT_MAX)
			throw std::invalid_argument("Cannot replace value not in set!");
		this->value[p] = v;
	}
	void empty() noexcept {
		this->deallocate();
		this->allocate(0);
	}
};

/*
 * class koset - Keyed ordered set
 *
 * Think oset and ksset.
 */
template <class K, class V>
class koset : public ksset<K,V> {
public:
	koset()
	  : ksset<K,V>() {
	}
	koset(unsigned s, unsigned b)
	  : ksset<K,V>(s,b) {
	}
	explicit koset(unsigned s, const V * v = nullptr, unsigned b = DFLT_BLK)
	  : ksset<K,V>(s,v,b) {
		if (v)
			sort();
	}
	explicit koset(const kfset<K,V> & v)
	  : ksset<K,V>(v) {
		sort();
	}
	~koset() {
	}

	void sort() {
		qsort(0,this->size);
	}
protected:
	void qsort(unsigned l, unsigned r) {
		if (r <= l)
			return;
		unsigned i, j, k, p = l+(r-1-l)/2;
		if (static_cast<K &>(this->value[l])
				> static_cast<K &>(this->value[p]))
			std::swap(this->value[l],this->value[p]);
		if (static_cast<K &>(this->value[l])
				> static_cast<K &>(this->value[r-1]))
			std::swap(this->value[l],this->value[r-1]);
		if (static_cast<K &>(this->value[p])
				> static_cast<K &>(this->value[r-1]))
			std::swap(this->value[p],this->value[r-1]);
		if (r-1-l <= 2)
			return;
		for (i = l, j = r-1, k = 0; ;
						p = (p==i?j++:(p==j?i--:p)), k++) {
			while (++i < p && !(static_cast<K &>(this->value[i])
					> static_cast<K &>(this->value[p]))) {};
			while (p < --j && !(static_cast<K &>(this->value[p])
					> static_cast<K &>(this->value[j]))) {};
			if (i >= j)
				break;
			std::swap(this->value[i],this->value[j]);
		}
		if (k == 0) {
			isort(l,p);
			isort(p+1,r);
		} else {
			qsort(l,p);
			qsort(p+1,r);
		}
	}
	void isort(unsigned l, unsigned r) {
		for (unsigned i = l; l < r && i < r-1; i++) {
			for (unsigned j = i+1; j > l
					 && static_cast<K &>(this->value[j-1])
					 > static_cast<K &>(this->value[j]); j--) {
				std::swap(this->value[j-1],this->value[j]);
			}
		}
	}
};

/*
 * class kiset - Keyed indexed set
 *
 * This is like a ksset, expect that sorting only sorts an index.
 * As a result haskey() is much faster, but at the expense of having
 * the set be read only...
 */
template <class K, class V>
class kiset : public set {
public:
	// Make one...
	kiset()
	  : set(), idx(nullptr), value(nullptr) {
		allocate(size);
	}
	kiset(unsigned s, unsigned b)
	  : set(b), idx(nullptr), value(nullptr) {
		allocate(s);
		copy(s,nullptr,true);
	}
	// Basic constructor
	explicit kiset(unsigned s, const V * v = nullptr, unsigned b = DFLT_BLK)
	  : set(b), idx(nullptr), value(nullptr) {
		allocate(s);
		copy(s,v,true);
		allocate(size);
	}
	// Copy constructor.
	kiset(const kiset<K,V> & v)
	  : set(v), idx(nullptr), value(nullptr) {
		allocate(v.size);
		copy(v.size,v.value,false);
	}
	// Copy from a kfset.
	explicit kiset(const kfset<K,V> & v)
	  : set(v), idx(nullptr), value(nullptr) {
		allocate(v.size);
		copy(v.size,v.value,true);
		allocate(size);
	}
	// Clean up.
	~kiset() {
		deallocate();
	}

	// Do a key lookup, and return -1 if the key is not found.
	unsigned haskey(const K & k) const noexcept {
		unsigned p = findkey(k);
		if (p < size && static_cast<K &>(value[idx[p]]) == k)
			return idx[p];
		else
			return UINT_MAX;
	}
	// Assignment operator.
	kiset<K,V> & operator = (const kiset<K,V> & v) {
		deallocate();
		allocate(v.size);
		copy(v.size,v.value,false);
		return *this;
	}
	kiset<K,V> & operator = (const kfset<K,V> & v) {
		deallocate();
		allocate(v.size);
		copy(v.size,v.value,true);
		allocate(size);
		return *this;
	}
	// Return data based on a key lookup.
	const V & operator [] (const K & k) const {
		unsigned p = haskey(k);
#if defined(DEBUG)
		if (!inbounds(p))
			throw std::out_of_range("Key out of range!");
#endif
		return value[p];
	}
	// Allow integer keys.
	const V & operator [] (unsigned p) const {
#if defined(DEBUG)
		if (!inbounds(p))
			throw std::out_of_range("Index out of range!");
#endif
		return value[p];
	}
	// Allow sorted access.
	const V & getindex(unsigned i) const {
#if defined(DEBUG)
		if (!inbounds(i))
			throw std::out_of_range("Index out of range!");
#endif
		return value[idx[i]];
	}
	// Get the position of an element in the sort.
	unsigned getorder(unsigned i) const {
#if defined(DEBUG)
		if (!inbounds(i))
			throw std::out_of_range("Index out of range!");
#endif
		return findkey(value[i]);
	}
	// Flatten the set into a traditional array.  Caller is responsible
	// for correct size.
	void copyout(V * v) const noexcept {
		for (unsigned i = 0; i < size; i++)
			v[i] = value[i];
	}


	// Add one value at the end.
	void add(const V & v) {
		add(&v,1);
	}
	// Add a whole set, at the end.
	void add(const kfset<K,V> & v) {
		add(v.value,v.size);
	}
	// Add an array of values... There is no point in
	// a position based addition.
	void add(const V * v, unsigned s = 1) {
		if (s == 0)
			throw std::invalid_argument("Cannot add zero elements into a set!");
		if (v == nullptr)
			throw std::invalid_argument("Cannot add a null to a set!");
		allocate(size+s);
		copy(s,v,true);
		allocate(size);
	}
	// Now start removing them...
	void remove(const K & k) {
		unsigned p = haskey(k), q = UINT_MAX;
		if (p == UINT_MAX)
			throw std::invalid_argument("Cannot remove key not in set!");
		value[p].~V();
		if (p < --size)
			memmove(&value[p],&value[p+1],(size-p)*sizeof(V));
		for (unsigned i = 0; i <= size; i++) {
			if (idx[i] == p)
				q = i;
			else if (idx[i] > p)
				idx[i]--;
		}
		if (q < size)
			memmove(&idx[q],&idx[q+1],(size-q)*sizeof(unsigned));
		allocate(size);
	}
	// Replace key/value with another.
	void replace(const V & v) {
		unsigned p = haskey(v);
		if (p == UINT_MAX)
			throw std::invalid_argument("Cannot replace key not in set!");
		value[p] = v;
	}
	void empty() noexcept {
		deallocate();
		allocate(0);
	}

protected:
	unsigned * idx;            // The index.
	V * value;                 // Take a guess...
	struct _V {                // Placement new wrapper
		V _v;
		explicit _V(const V & v) : _v(v) {}
		void * operator new (size_t, void * p) noexcept {
			return p;
		}
		void operator delete (void * , void *) noexcept {
		}
	};

	// Make some space...
	void allocate(unsigned s) {
		unsigned b = bufsize(s);
		if (b == buffer)
			return;
		if (b == 0) {
			free(idx);
			idx = nullptr;
			free(value);
			value = nullptr;
			buffer = 0;
			return;
		}
		unsigned * itemp = static_cast<unsigned *>(
				realloc(idx,b*sizeof(unsigned)));
		V * vtemp = static_cast<V *>(realloc(value,b*sizeof(V)));
		if (itemp == nullptr || vtemp == nullptr)
			throw std::bad_alloc();
		idx = itemp;
		value = vtemp;
		buffer = b;
	}
	void deallocate() noexcept {
		if (idx) {
			free(idx);
			idx = nullptr;
		}
		if (value) {
			for (unsigned i = 0; i < size; i++)
				value[i].~V();
			free(value);
			value = nullptr;
		}
		size = 0;
		buffer = 0;
	}
	// Find the position which is either equal or greater...
	unsigned findkey(const K & k) const noexcept {
		unsigned l = 0, r = size;
		while (l < r) {
			unsigned i = l + (r-l)/2;
			if (static_cast<K &>(value[idx[i]]) < k)
				l = i+1;
			else
				r = i;
		}
		return l;
	}
	void init(const V * v) {
		if (v == nullptr)
			throw std::invalid_argument("Cannot insert a null into a set");
		new(&value[size]) _V(*v);
		unsigned p = findkey(static_cast<const K &>(*v));
		if (p < size)
			memmove(&idx[p+1],&idx[p],(size-p)*sizeof(unsigned));
		idx[p] = size++;
	}
	void copy(unsigned s, const V * v, bool checkdups) {
		for (unsigned i = 0, p; v && i < s; i++) {
			if (checkdups && (p = haskey(v[i])) != UINT_MAX)
				value[p] = v[i];
			else
				init(&v[i]);
		}
	}
};

/*
 * class afset - Associative fixed set
 *
 * Instead of having the key as a base class of the values, these
 * classes have separate keys and values. Now we have to manage
 * two arrays of things...
 */
template <class K, class V>
class afset : public set {
public:
	// Make one...
	afset(unsigned s, unsigned b)
	  : set(b), key(nullptr), value(nullptr) {
		allocate(s);
		copy(s,nullptr,nullptr,true);
	}
	explicit afset(unsigned s, const K * k = nullptr,
			const V * v = nullptr, unsigned b = DFLT_BLK)
	  : set(b), key(nullptr), value(nullptr) {
		allocate(s);
		copy(s,k,v,true);
		allocate(size);
	}
	// Copy one...
	afset(const afset<K,V> & v)
	  : set(v), key(nullptr), value(nullptr) {
		allocate(v.size);
		copy(v.size,v.key,v.value,false);
	}
	// Kill one...
	~afset() {
		deallocate();
	}

	// Return an index for a key.
	unsigned haskey(const K & k) const noexcept {
		for (unsigned i = 0; i < size; i++) {
			if (key[i] == k)
				return i;
		}
		return UINT_MAX;
	}
	// Assignment operator.
	afset<K,V> & operator = (const afset<K,V> & v) {
		deallocate();
		allocate(v.size);
		copy(v.size,v.key,v.value,false);
		return *this;
	}
	// Behave like an indexed array...
	V & operator [] (const K & k) const {
		unsigned p = haskey(k);
#if defined(DEBUG)
		if (!inbounds(p))
			throw std::out_of_range("Key out of range!");
#endif
		return value[p];
	}
	// Linear access.
	K & getkey(unsigned p) const {
#if defined(DEBUG)
		if (!inbounds(p))
			throw std::out_of_range("Index out of range!");
#endif
		return key[p];
	}
	// More linear access.
	V & getvalue(unsigned p) const {
#if defined(DEBUG)
		if (!inbounds(p))
			throw std::out_of_range("Index out of range!");
#endif
		return value[p];
	}

protected:
	K * key;                   // The keys.
	V * value;                 // Take a guess...
	struct _K {                // Placement new wrapper
		K _k;
		explicit _K(const K & k) : _k(k) {}
		void * operator new (size_t, void * p) noexcept {
			return p;
		}
		void operator delete (void * , void *) noexcept {
		}
	};
	struct _V {                // Placement new wrapper
		V _v;
		_V() : _v() {}
		explicit _V(const V & v) : _v(v) {}
		void * operator new (size_t, void * p) {
			return p;
		}
		void operator delete (void * , void *) {
		}
	};

	// Make some space...
	void allocate(unsigned s) {
		unsigned b = bufsize(s);
		if (b == buffer)
			return;
		if (b == 0) {
			free(key);
			key = nullptr;
			free(value);
			value = nullptr;
			buffer = 0;
			return;
		}
		K * ktemp = static_cast<K *>(realloc(key,b*sizeof(K)));
		V * vtemp = static_cast<V *>(realloc(value,b*sizeof(V)));
		if (ktemp == nullptr || vtemp == nullptr)
			throw std::bad_alloc();
		key = ktemp;
		value = vtemp;
		buffer = b;
	}
	void deallocate() noexcept {
		if (key) {
			for (unsigned i = 0; i < size; i++)
				key[i].~K();
			free(key);
			key = nullptr;
		}
		if (value) {
			for (unsigned i = 0; i < size; i++)
				value[i].~V();
			free(value);
			value = nullptr;
		}
		size = 0;
		buffer = 0;
	}
	void init(const K * k, const V * v) {
		if (k == nullptr)
			throw std::invalid_argument("Cannot insert a null key into a set");
		if (v)
			new(&value[size]) _V(*v);
		else
			new(&value[size]) _V();
		new(&key[size++]) _K(*k);
	}
	void copy(unsigned s, const K * k, const V * v, bool checkdups) {
		for (unsigned i = 0, p; k && i < s; i++) {
			if (checkdups && (p = haskey(k[i])) != UINT_MAX)
				value[p] = (v ? v[i] : V());
			else
				init(&k[i],(v ? &v[i] : nullptr));
		}
	}
	// Don't allow yobos to make empty sets...
	afset()
	  : set(), key(nullptr), value(nullptr) {
		allocate(size);
	}
};

/*
 * class asset - Associative sizeable set
 */
template <class K, class V>
class asset : public afset<K,V> {
public:
	// Empty sizeable sets are OK.
	asset()
	  : afset<K,V>() {
	}
	asset(unsigned s, unsigned b)
	  : afset<K,V>(s,b) {
	}
	// Simple constructor.
	explicit asset(unsigned s, const K * k = nullptr,
			const V * v = nullptr, unsigned b = DFLT_BLK)
	  : afset<K,V>(s,k,v,b) {
	}
	// Copy constructor.
	explicit asset(const afset<K,V> & v)
	  : afset<K,V>(v) {
	}
	// Destructor.
	~asset() {
	}

	// Add one...
	void add(const K & k, const V & v) {
		add(&k,&v,1);
	}
	// Add a whole set, at the end.
	void add(const afset<K,V> & v) {
		add(v.key,v.value,v.size);
	}
	// Add a whole bunch...
	void add(const K * k, const V * v = nullptr,
			unsigned s = 1) {
		if (s == 0)
			throw std::invalid_argument("Cannot add zero elements into a set!");
		if (k == nullptr)
			throw std::invalid_argument("Cannot add a null key into a set!");
		this->allocate(this->size+s);
		this->copy(s,k,v,true);
		this->allocate(this->size);
	}
	// Now start removing them...
	void remove(const K & k) {
		unsigned p = this->haskey(k);
		if (p == UINT_MAX)
			throw std::invalid_argument("Cannot remove value not in set!");
		this->key[++p].~K();
		this->value[p].~V();
		if (p < --this->size) {
			memmove(&this->key[p],&this->key[p+1],
					(this->size-p)*sizeof(K));
			memmove(&this->value[p],&this->value[p+1],
					(this->size-p)*sizeof(V));
		}
		this->allocate(this->size);
	}
	// Or replacing them...
	void replace(const K & k, const V & v) {
		unsigned p = this->haskey(k);
		if (p == UINT_MAX)
			throw std::invalid_argument("Cannot replace key not in set!");
		this->value[p] = v;
	}
	void empty() noexcept {
		this->deallocate();
		this->allocate(0);
	}
};

/*
 * class aoset - Associative ordered set
 *
 * Ordered by key set.
 */
template <class K, class V>
class aoset : public asset<K,V> {
public:
	aoset()
	  : asset<K,V>() {
	}
	aoset(unsigned s, unsigned b)
	  : asset<K,V>(s,b) {
	}
	explicit aoset(unsigned s, const K * k = nullptr,
			const V * v = nullptr, unsigned b = DFLT_BLK)
	  : asset<K,V>(s,k,v,b) {
		if (k)
			sort();
	}
	explicit aoset(const afset<K,V> & v)
	  : asset<K,V>(v) {
		sort();
	}
	~aoset() {
	}

	void sort() {
		qsort(0,this->size);
	}
protected:
	void qsort(unsigned l, unsigned r) {
		if (r <= l)
			return;
		unsigned i, j, k, p = l+(r-1-l)/2;
		if (this->key[l] > this->key[p]) {
			std::swap(this->key[l],this->key[p]);
			std::swap(this->value[l],this->value[p]);
		}
		if (this->key[l] > this->key[r-1]) {
			std::swap(this->key[l],this->key[r-1]);
			std::swap(this->value[l],this->value[r-1]);
		}
		if (this->key[p] > this->key[r-1]) {
			std::swap(this->key[p],this->key[r-1]);
			std::swap(this->value[p],this->value[r-1]);
		}
		if (r-1-l <= 2)
			return;
		for (i = l, j = r-1, k = 0; ;
				p = (p==i?j++:(p==j?i--:p)), k++) {
			while (++i < p && !(this->key[i] > this->key[p])) {};
			while (p < --j && !(this->key[p] > this->key[j])) {};
			if (i >= j)
				break;
			std::swap(this->key[i],this->key[j]);
			std::swap(this->value[i],this->value[j]);
		}
		if (k == 0) {
			isort(l,p);
			isort(p+1,r);
		} else {
			qsort(l,p);
			qsort(p+1,r);
		}
	}
	void isort(unsigned l, unsigned r) {
		for (unsigned i = l; l < r && i < r-1; i++) {
			for (unsigned j = i+1; j > l
					&& this->key[j-1] > this->key[j]; j--) {
				std::swap(this->key[j-1],this->key[j]);
				std::swap(this->value[j-1],this->value[j]);
			}
		}
	}
};

/*
 * class avoset - Associative value ordered set
 */
template <class K, class V>
class avoset : public asset<K,V> {
public:
	avoset()
	  : asset<K,V>() {
	}
	avoset(unsigned s, unsigned b)
	  : asset<K,V>(s,b) {
	}
	explicit avoset(unsigned s, const K * k = nullptr,
			const V * v = nullptr, unsigned b = DFLT_BLK)
	  : asset<K,V>(s,k,v,b) {
		if (k)
			sort();
	}
	explicit avoset(const afset<K,V> & v)
	  : asset<K,V>(v) {
		sort();
	}
	~avoset() {
	}

	void sort() {
		qsort(0,this->size);
	}

protected:
	void qsort(unsigned l, unsigned r) {
		if (r <= l)
			return;
		unsigned i, j, k, p = l+(r-1-l)/2;
		if (this->value[l] > this->value[p]) {
			std::swap(this->key[l],this->key[p]);
			std::swap(this->value[l],this->value[p]);
		}
		if (this->value[l] > this->value[r-1]) {
			std::swap(this->key[l],this->key[r-1]);
			std::swap(this->value[l],this->value[r-1]);
		}
		if (this->value[p] > this->value[r-1]) {
			std::swap(this->key[p],this->key[r-1]);
			std::swap(this->value[p],this->value[r-1]);
		}
		if (r-1-l <= 2)
			return;
		for (i = l, j = r-1, k = 0; ;
				p = (p==i?j++:(p==j?i--:p)), k++) {
			while (++i < p && !(this->value[i] > this->value[p])) {};
			while (p < --j && !(this->value[p] > this->value[j])) {};
			if (i >= j)
				break;
			std::swap(this->key[i],this->key[j]);
			std::swap(this->value[i],this->value[j]);
		}
		if (k == 0) {
			isort(l,p);
			isort(p+1,r);
		} else {
			qsort(l,p);
			qsort(p+1,r);
		}
	}
	void isort(unsigned l, unsigned r) {
		for (unsigned i = l; l < r && i < r-1; i++) {
			for (unsigned j = i+1; j > l
					&& this->value[j-1] > this->value[j]; j--) {
				std::swap(this->key[j-1],this->key[j]);
				std::swap(this->value[j-1],this->value[j]);
			}
		}
	}
};

} // namespace OP

#endif // SET_H
