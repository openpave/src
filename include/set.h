/**************************************************************************

	SET.H - Simple set template classes

	$OpenPave$

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

	Portions Copyright (C) 2006-2008 OpenPave.org.

	Contributor(s): Jeremy Lea <reg@openpave.org>.

	Purpose:
		This header implements a number of templated C++ classes for
		various types of sets. The code is fairly compact, but not
		optimised for big sets (particularly sorting).

	Design:
		There are two types of sets. The first are normal sets with
		one variable. Then there are associative sets, which store
		a key and a value.

		The heirarchy goes from a fixed size set, to a resizeable set,
		to a sorted resizable set.

		All of the templated sets take a minimum block size.  Major
		performance enhancements will be seen on resizable sets if this
		value is bigger than the maximum size of the set, because the
		memory is reused.
		
	Status:
		The swaps in the sorts should be replaced by memmove.  It should
		also learn to throw expections...

	History:
		1994       - Created by Jeremy Lea <reg@openpave.org>
		2002/01/23 - Modifications for use on a modern C++ compiler
		2007/09/11 - Modified to use realloc and placement new.
		2008/02/11 - Removed the idea of a default element.

**************************************************************************/

#ifndef __SET_H
#define __SET_H

#include "event.h"
#include "mathplus.h"
#include <string.h>
#include <limits.h>
#include <assert.h>

#define DFLT_BLK	64

/*
 * class set - Set base class
 *
 * This class provides the basic working for the rest of the set classes,
 * including the calcualtion of the buffer size.
 */
class set {
public:
	// The length. Nice for lots of things...
	inline unsigned length() const {
		return size;
	}
	// Check if an index is within the set...
	inline bool inbounds(const unsigned p) const {
		return (p < size ? true : false);
	}

protected:
	unsigned size;             // The size of the set...
	unsigned buffer;           // The allocated buffer size...

	// Simple constructor...
	inline explicit set(const unsigned b = DFLT_BLK)
	  : size(0), buffer(0), block(b > 1 ? b : 1) {
	}
	inline explicit set(const set & s)
	  : size(0), buffer(0), block(s.block) {
	}
	inline ~set() {
	}
	// Calculate the buffer size.
	inline unsigned bufsize(unsigned s) {
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
 * This class acts like a fixed size array, except that it returns
 * the default member if the index is out of bounds.
 */
template <class V>
class fset : public set {
public:
	// Nice simple constructor...
	inline explicit fset(const unsigned s, const unsigned b)
	  : set(b), value(0) {
		allocate(s);
	}
	inline explicit fset(const unsigned s, const V * v = 0,
			const unsigned b = DFLT_BLK)
	  : set(b), value(0) {
		allocate(s);
		copy(s,v);
	}
	// Copy constructor.
	inline explicit fset(const fset<V> & v)
	  : set(v), value(0) {
		allocate(v.size);
		copy(v.size,v.value);
	}
	// Wow, a destructor...
	inline ~fset() {
		deallocate();
	}

	// Assignment operator.
	inline fset<V> & operator= (const fset<V> & v) {
		deallocate();
		allocate(v.size);
		copy(v.size,v.value);
		return *this;
	}
	// Allow the size to be changed, even for 'fixed' sets.
	void resize(const unsigned s) {
		deallocate();
		allocate(s);
		copy(s,0);
	}
	// Behave like an array. Zero indexed.
	V & operator[] (const unsigned p) const {
		assert(inbounds(p));
		return value[p];
	}

protected:
	V * value;                 // The buffer.
	struct _V {                // Placement new wrapper
		V _v;
		explicit _V() : _v() {}
		explicit _V(const V & v) : _v(v) {}
		void * operator new(size_t, void * p) {
			return p;
		} 
		void operator delete(void * , void *) {
		} 
	};

	// Since this is a fixed size set, hide the allocator...
	void allocate(const unsigned s) {
		unsigned b = bufsize(s);
		if (b == buffer)
			return;
		if (b == 0) {
			free(value);
			value = 0;
			buffer = 0;
			return;
		}
		V * temp = static_cast<V *>(realloc(value,b*sizeof(V)));
		if (temp == 0)
			throw std::bad_alloc();
		value = temp;
		buffer = b;
	}
	void deallocate() {
		if (value) {
			for (unsigned i = 0; i < size; i++)
				value[i].~V();
			free(value);
			value = 0;
		}
		size = 0;
		buffer = 0;
	}
	void init(const unsigned i, const V * v) {
		if (v != 0)
			new(&value[i]) _V(*v);
		else
			new(&value[i]) _V();
	}
	void copy(const unsigned s, const V * v) {
		for (unsigned i = 0; i < s; i++)
			init(size++,(v ? &v[i] : 0));
	}
	// Also hide the null constructor.
	inline explicit fset()
	  : set(), value(0) {
		allocate(size);
	}
};

/*
 * Special case these to avoid constructors
 */
template<>
inline void fset<int>::init(const unsigned i, const int * v) {
	if (v)
		value[i] = *v;
}
template<>
inline void fset<unsigned>::init(const unsigned i, const unsigned * v) {
	if (v)
		value[i] = *v;
}
template<>
inline void fset<float>::init(const unsigned i, const float * v) {
	if (v)
		value[i] = *v;
}
template<>
inline void fset<double>::init(const unsigned i, const double * v) {
	if (v)
		value[i] = *v;
}

/*
 * class sset - Sizable set of type V.
 *
 * This class builds on fset to make a set which can be resized
 * at will. Lots of nice add() and remove() options are available.
 */
template <class V>
class sset : public fset<V> {
public:
	// Provide a null constuctor for empty sets.
	inline explicit sset()
	  : fset<V>() {
	}
	inline explicit sset(const unsigned s, const unsigned b)
	  : fset<V>(s,b) {
	}
	// Simple constructor.
	inline explicit sset(const unsigned s, const V * v = 0,
			const unsigned b = DFLT_BLK)
	  : fset<V>(s,v,b) {
	}
	// Copy constructor.
	inline explicit sset(const fset<V> & v)
	  : fset<V>(v) {
	}
	// We let someone else clean up...
	inline ~sset() {
	}

	// Add one element, at the end.
	inline void add(const V & v) {
		add(this->size,&v,1);
	}
	// Add a whole set, at the end.
	inline void add(const fset<V> & v) {
		add(this->size,&(v[0]),v.length());
	}
	// Add an array, at the end.
	inline void add(const V * v, const unsigned s = 1) {
		add(this->size,v,s);
	}
	// Insert one element at position p.
	inline void add(const unsigned p, const V & v) {
		add(p,&v,1);
	}
	// Insert a set at position p.
	inline void add(const unsigned p, const fset<V> & v) {
		add(p,&(v[0]),v.length());
	}
	// Add an array at position p. (Actually do the work too).
	void add(unsigned p, const V * v, const unsigned s = 1) {
		unsigned i;
		assert(s > 0 && v != 0);
		if (p >= this->size) {
			this->allocate(p+s);
			for (i = this->size; p > 0 && i < p-1; i++)
				this->init(i,0);
			this->size = p+s;
		} else {
			this->allocate(this->size+s);
			memmove(&this->value[p+s],&this->value[p],
					(this->size-p)*sizeof(V));
			this->size += s;
		}
		for (i = 0; i < s; i++)
			this->init(p+i,&v[i]);
	}
	// Remove the last element.
	inline void remove() {
		remove(this->size-1,1);
	}
	// Remove s elements, starting at position p.
	void remove(unsigned p, unsigned s = 1) {
		assert(this->inbounds(p) && s > 0);
		s = (p+s > this->size ? this->size-p : s);
		for (unsigned i = 0; i < s; i++)
			this->value[p+i].~V();
		if (this->size-p > s)
			memmove(&this->value[p],&this->value[p+s],
				(this->size-p-s)*sizeof(V));
		this->size -= s;
		this->allocate(this->size);
	}
	inline void empty() {
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
	inline explicit oset()
	  : sset<V>() {
	}
	inline explicit oset(const unsigned s, const unsigned b)
	  : sset<V>(s,b) {
	}
	// Simple constructor.
	inline explicit oset(const unsigned s, const V * v = 0,
			const unsigned b = DFLT_BLK)
	  : sset<V>(s,v,b) {
		if (v)
			sort();
	}
	// Copy constructor.
	inline explicit oset(const fset<V> & v)
	  : sset<V>(v) {
		sort();
	}
	// Destructor.
	inline ~oset() {
	}

	// Guess what?
	inline void sort() {
		qsort(0,this->size);
	}
	// Do a lookup, and return -1 if the value is not found.
	// You must sort the set first!
	inline unsigned findvalue(const V & v) const {
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
	void qsort(const unsigned l, const unsigned r) {
		if (r <= l)
			return;
		unsigned i, j, k, p = l+(r-1-l)/2;
		if (this->value[l] > this->value[p])
			swap(this->value[l],this->value[p]);
		if (this->value[l] > this->value[r-1])
			swap(this->value[l],this->value[r-1]);
		if (this->value[p] > this->value[r-1])
			swap(this->value[p],this->value[r-1]);
		if (r-1-l <= 2)
			return;
		for (i = l, j = r-1, k = 0; ;
				p = (p==i?j++:(p==j?i--:p)), k++) {
			while (++i < p && !(this->value[i] > this->value[p])) {};
			while (p < --j && !(this->value[p] > this->value[j])) {};
			if (i >= j)
				break;
			swap(this->value[i],this->value[j]);
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
	void isort(const unsigned l, const unsigned r) {
		for (unsigned i = l; l < r && i < r-1; i++) {
			for (unsigned j = i+1; j > l
					&& this->value[j-1] > this->value[j]; j--)
				swap(this->value[j-1],this->value[j]);
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
	inline explicit cset()
	  : oset<V>() {
	}
	inline explicit cset(const unsigned s, const unsigned b)
	  : oset<V>(s,b) {
	}
	inline explicit cset(const unsigned s, const V * v = 0,
			const unsigned b = DFLT_BLK)
	  : oset<V>(s,v,b) {
		if (v)
			compact();
	}
	inline explicit cset(const fset<V> & v)
	  : oset<V>(v) {
		compact();
	}
	inline ~cset() {
	}

	// Sort then compact the set.
	inline void sort() {
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
	inline explicit iset()
	  : set(), idx(0), value(0) {
		allocate(size);
	}
	inline explicit iset(const unsigned s, const unsigned b)
	  : set(b), idx(0), value(0) {
		allocate(s);
	}
	// Basic constructor
	inline explicit iset(const unsigned s, const V * v = 0,
			const unsigned b = DFLT_BLK)
	  : set(b), idx(0), value(0) {
		allocate(s);
		copy(s,v,true);
		allocate(size);
	}
	// Copy constuctor.
	inline explicit iset(const iset<V> & v)
	  : set(v), idx(0), value(0) {
		allocate(v.size);
		copy(v.size,v.value,false);
	}
	// Copy from an fset.
	inline explicit iset(const fset<V> & v)
	  : set(v), idx(0), value(0) {
		allocate(v.size);
		copy(v.size,v.value,true);
		allocate(size);
	}
	// Clean up.
	inline ~iset() {
		deallocate();
	}

	// Do a lookup, and return -1 if the value is not found.
	inline unsigned hasvalue(const V & v) const {
		unsigned p = findvalue(v);
		if (p < size && value[idx[p]] == v)
			return idx[p];
		else
			return UINT_MAX;
	}
	// Assignment operator.
	inline iset<V> & operator= (const iset<V> & v) {
		deallocate();
		allocate(v.size);
		copy(v.size,v.value,false);
		return *this;
	}
	inline iset<V> & operator= (const fset<V> & v) {
		deallocate();
		allocate(v.size);
		copy(v.size,v.value,true);
		allocate(size);
		return *this;
	}
	// Only integer keys make sense.
	inline const V & operator[] (const unsigned p) const {
		assert(inbounds(p));
		return value[p];
	}
	// Allow sorted acess.
	inline const V & getindex(const unsigned i) const {
		assert(inbounds(i));
		return value[idx[i]];
	}
	// Get the position of an element in the sort.
	inline unsigned getorder(const unsigned i) const {
		assert(inbounds(i));
		return findvalue(value[i]);
	}

	// Add one value at the end.
	inline void add(const V & v) {
		return add(&v,1);
	}
	// Add a whole set, at the end.
	inline void add(const iset<V> & v) {
		return add(v.value,v.size);
	}
	inline void add(const fset<V> & v) {
		return add(v.value,v.size);
	}
	// Add an array of values... There is no point in
	// a position based addition. 
	void add(const V * v, const unsigned s = 1) {
		assert(s > 0 && v != 0);
		allocate(size+s);
		copy(s,v,true);
		allocate(size);
	}
	// Now start removing them...
	void remove(const V & v) {
		unsigned p = hasvalue(v), q = UINT_MAX;
		assert(p != UINT_MAX);
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
	inline void replace(const V & v) {
		unsigned p = hasvalue(v);
		assert (p != UINT_MAX);
		value[p] = v;
	}
	inline void empty() {
		deallocate();
		allocate(0);
	}

protected:
	unsigned * idx;            // The index.
	V * value;                 // Take a guess...
	struct _V {                // Placement new wrapper
		V _v;
		explicit _V(const V & v) : _v(v) {}
		void * operator new(size_t, void * p) {
			return p;
		} 
		void operator delete(void * , void *) {
		} 
	};

	// Make some space...
	void allocate(const unsigned s) {
		unsigned b = bufsize(s);
		if (b == buffer)
			return;
		if (b == 0) {
			free(idx);
			idx = 0;
			free(value);
			value = 0;
			buffer = 0;
			return;
		}
		unsigned * itemp = static_cast<unsigned *>(
				realloc(idx,b*sizeof(unsigned)));
		V * vtemp = static_cast<V *>(realloc(value,b*sizeof(V)));
		if (itemp == 0 || vtemp == 0)
			throw std::bad_alloc();
		idx = itemp;
		value = vtemp;
		buffer = b;
	}
	void deallocate() {
		if (idx) {
			free(idx);
			idx = 0;
		}
		if (value) {
			for (unsigned i = 0; i < size; i++)
				value[i].~V();
			free(value);
			value = 0;
		}
		size = 0;
		buffer = 0;
	}
	// Find the position which is either equal or greater...
	inline unsigned findvalue(const V & v) const {
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
		assert(v != 0);
		new(&value[size]) _V(*v);
		unsigned p = findvalue(*v);
		if (p < size)
			memmove(&idx[p+1],&idx[p],(size-p)*sizeof(unsigned));
		idx[p] = size++;
	}
	void copy(const unsigned s, const V * v, bool checkdups) {
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
	inline explicit kfset(const unsigned s, const unsigned b)
	  : set(b), value(0) {
		allocate(s);
	}
	inline explicit kfset(const unsigned s, const V * v = 0,
			const unsigned b = DFLT_BLK)
	  : set(b), value(0) {
		allocate(s);
		copy(s,v,true);
		allocate(size);
	}
	// Copy constuctor.
	inline explicit kfset(const kfset<K,V> & v)
	  : set(v), value(0) {
		allocate(v.size);
		copy(v.size,v.value,false);
	}
	// Clean up.
	inline ~kfset() {
		deallocate();
	}

	// Do a key lookup, and return UINT_MAX if the key is not found.
	inline unsigned haskey(const K & k) const {
		for (unsigned i = 0; i < size; i++) {
			if (static_cast<K &>(value[i]) == k)
				return i;
		}
		return UINT_MAX;
	}
	// Assignment operator.
	inline kfset<K,V> & operator= (const kfset<K,V> & v) {
		deallocate();
		allocate(v.size);
		copy(v.size,v.value,false);
		return *this;
	}
	// Return data based on a key lookup.
	inline V & operator[] (const K & k) const {
		unsigned p = haskey(k);
		assert(inbounds(p));
		return value[p];
	}
	// Allow integer keys.
	inline V & operator[] (const unsigned p) const {
		assert(inbounds(p));
		return value[p];
	}

protected:
	V * value;                 // The data.
	struct _V {                // Placement new wrapper
		V _v;
		explicit _V(const V & v) : _v(v) {}
		void * operator new(size_t, void * p) {
			return p;
		} 
		void operator delete(void * , void *) {
		} 
	};

	// Hide the allocation function.
	void allocate(const unsigned s) {
		unsigned b = bufsize(s);
		if (b == buffer)
			return;
		if (b == 0) {
			free(value);
			value = 0;
			buffer = 0;
			return;
		}
		V * temp = static_cast<V *>(realloc(value,b*sizeof(V)));
		if (temp == 0)
			throw std::bad_alloc();
		value = temp;
		buffer = b;
	}
	void deallocate() {
		if (value) {
			for (unsigned i = 0; i < size; i++)
				value[i].~V();
			free(value);
			value = 0;
		}
		size = 0;
		buffer = 0;
	}
	void init(const V * v) {
		assert(v != 0);
		new(&value[size++]) _V(*v);
	}
	void copy(const unsigned s, const V * v, bool checkdups) {
		for (unsigned i = 0, p; v && i < s; i++) {
			if (checkdups && (p = haskey(v[i])) != UINT_MAX)
				value[p] = v[i];
			else
				init(&v[i]);
		}
	}
	// And the null constructor.
	inline explicit kfset()
	  : set(), value(0) {
		allocate(size);
	}
};

/*
 * class ksset - Keyed sizable set
 *
 * Think combination of sset and kfset.
 */
template <class K, class V>
class ksset : public kfset<K,V> {
public:
	// C++ sucks...
	inline explicit ksset()
	  : kfset<K,V>() {
	}
	inline explicit ksset(const unsigned s, const unsigned b)
	  : kfset<K,V>(s,b) {
	}
	inline explicit ksset(const unsigned s, const V * v = 0,
			const unsigned b = DFLT_BLK)
	  : kfset<K,V>(s,v,b) {
	}
	inline explicit ksset(const kfset<K,V> & v)
	  : kfset<K,V>(v) {
	}
	inline ~ksset() {
	}

	// Add one value at the end.
	inline void add(const V & v) {
		add(&v,1);
	}
	// Add a whole set, at the end.
	inline void add(const kfset<K,V> & v) {
		add(v.value,v.size);
	}
	// Add an array of values... There is no point in
	// a position based addition. 
	void add(const V * v, const unsigned s = 1) {
		assert(s > 0 && v != 0);
		this->allocate(this->size+s);
		this->copy(s,v,true);
		this->allocate(this->size);
	}
	// Remove based on key.
	void remove(const K & k) {
		unsigned p = this->haskey(k);
		assert(p != UINT_MAX);
		this->value[p].~V();
		if (p < --this->size)
			memmove(&this->value[p],&this->value[p+1],
				(this->size-p)*sizeof(V));
		this->allocate(this->size);
	}
	// Replace key/value with another.
	inline void replace(const V & v) {
		unsigned p = this->haskey(v);
		assert(p != UINT_MAX);
		this->value[p] = v;
	}
	inline void empty() {
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
	inline explicit koset()
	  : ksset<K,V>() {
	}
	inline explicit koset(const unsigned s, const unsigned b)
	  : ksset<K,V>(s,b) {
	}
	inline explicit koset(const unsigned s, const V * v = 0,
			const unsigned b = DFLT_BLK)
	  : ksset<K,V>(s,v,b) {
		if (v)
			sort();
	}
	inline explicit koset(const kfset<K,V> & v)
	  : ksset<K,V>(v) {
		sort();
	}
	inline ~koset() {
	}

	inline void sort() {
		qsort(0,this->size);
	}
protected:
	void qsort(const unsigned l, const unsigned r) {
		if (r <= l)
			return;
		unsigned i, j, k, p = l+(r-1-l)/2;
		if (this->value[l] > this->value[p])
			swap(this->value[l],this->value[p]);
		if (this->value[l] > this->value[r-1])
			swap(this->value[l],this->value[r-1]);
		if (this->value[p] > this->value[r-1])
			swap(this->value[p],this->value[r-1]);
		if (r-1-l <= 2)
			return;
		for (i = l, j = r-1, k = 0; ;
						p = (p==i?j++:(p==j?i--:p)), k++) {
			while (++i < p && !(this->value[i] > this->value[p])) {};
			while (p < --j && !(this->value[p] > this->value[j])) {};
			if (i >= j)
				break;
			swap(this->value[i],this->value[j]);
		}
		if (k == 0) {
			isort(l,p);
			isort(p+1,r);
		} else {
			qsort(l,p);
			qsort(p+1,r);
		}
	}
	void isort(const unsigned l, const unsigned r) {
		for (unsigned i = l; l < r && i < r-1; i++) {
			for (unsigned j = i+1; j > l
					 && this->value[j-1] > this->value[j]; j--) {
				swap(this->value[j-1],this->value[j]);
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
	inline explicit kiset()
	  : set(), idx(0), value(0) {
		allocate(size);
	}
	inline explicit kiset(const unsigned s, const unsigned b)
	  : set(b), idx(0), value(0) {
		allocate(s);
	}
	// Basic constructor
	inline explicit kiset(const unsigned s, const V * v = 0,
			const unsigned b = DFLT_BLK)
	  : set(b), idx(0), value(0) {
		allocate(s);
		copy(s,v,true);
		allocate(size);
	}
	// Copy constuctor.
	inline explicit kiset(const kiset<K,V> & v)
	  : set(v), idx(0), value(0) {
		allocate(v.size);
		copy(v.size,v.value,false);
	}
	// Copy from a kfset.
	inline explicit kiset(const kfset<K,V> & v)
	  : set(v), idx(0), value(0) {
		allocate(v.size);
		copy(v.size,v.value,true);
		allocate(size);
	}
	// Clean up.
	inline ~kiset() {
		deallocate();
	}

	// Do a key lookup, and return -1 if the key is not found.
	inline unsigned haskey(const K & k) const {
		unsigned p = findkey(k);
		if (p < size && static_cast<K &>(value[idx[p]]) == k)
			return idx[p];
		else
			return UINT_MAX;
	}
	// Assignment operator.
	inline kiset<K,V> & operator= (const kiset<K,V> & v) {
		deallocate();
		allocate(v.size);
		copy(v.size,v.value,false);
		return *this;
	}
	inline kiset<K,V> & operator= (const kfset<K,V> & v) {
		deallocate();
		allocate(v.size);
		copy(v.size,v.value,true);
		allocate(size);
		return *this;
	}
	// Return data based on a key lookup.
	inline const V & operator[] (const K & k) const {
		unsigned p = haskey(k);
		assert(inbounds(p));
		return value[p];
	}
	// Allow integer keys.
	inline const V & operator[] (const unsigned p) const {
		assert(inbounds(p));
		return value[p];
	}
	// Allow sorted acess.
	inline const V & getindex(const unsigned i) const {
		assert(inbounds(i));
		return value[idx[i]];
	}
	// Get the position of an element in the sort.
	inline unsigned getorder(const unsigned i) const {
		assert(inbounds(i));
		return findkey(value[i]);
	}

	// Add one value at the end.
	inline void add(const V & v) {
		add(&v,1);
	}
	// Add a whole set, at the end.
	inline void add(const kfset<K,V> & v) {
		add(v.value,v.size);
	}
	// Add an array of values... There is no point in
	// a position based addition. 
	void add(const V * v, const unsigned s = 1) {
		assert(s > 0 && v != 0);
		allocate(size+s);
		copy(s,v,true);
		allocate(size);
	}
	// Now start removing them...
	void remove(const K & k) {
		unsigned p = haskey(k), q = UINT_MAX;
		assert(p != UINT_MAX);
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
	inline void replace(const V & v) {
		unsigned p = haskey(v);
		assert(p != UINT_MAX);
		value[p] = v;
	}
	inline void empty() {
		deallocate();
		allocate(0);
	}

protected:
	unsigned * idx;            // The index.
	V * value;                 // Take a guess...
	struct _V {                // Placement new wrapper
		V _v;
		explicit _V(const V & v) : _v(v) {}
		void * operator new(size_t, void * p) {
			return p;
		} 
		void operator delete(void * , void *) {
		} 
	};

	// Make some space...
	void allocate(const unsigned s) {
		unsigned b = bufsize(s);
		if (b == buffer)
			return;
		if (b == 0) {
			free(idx);
			idx = 0;
			free(value);
			value = 0;
			buffer = 0;
			return;
		}
		unsigned * itemp = static_cast<unsigned *>(
				realloc(idx,b*sizeof(unsigned)));
		V * vtemp = static_cast<V *>(realloc(value,b*sizeof(V)));
		if (itemp == 0 || vtemp == 0)
			throw std::bad_alloc();
		idx = itemp;
		value = vtemp;
		buffer = b;
	}
	void deallocate() {
		if (idx) {
			free(idx);
			idx = 0;
		}
		if (value) {
			for (unsigned i = 0; i < size; i++)
				value[i].~V();
			free(value);
			value = 0;
		}
		size = 0;
		buffer = 0;
	}
	// Find the position which is either equal or greater...
	inline unsigned findkey(const K & k) const {
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
		assert(v != 0);
		new(&value[size]) _V(*v);
		unsigned p = findkey(static_cast<const K &>(*v));
		if (p < size)
			memmove(&idx[p+1],&idx[p],(size-p)*sizeof(unsigned));
		idx[p] = size++;
	}
	void copy(const unsigned s, const V * v, bool checkdups) {
		for (unsigned i = 0, p; v && i < s; i++) {
			if (checkdups && (p = haskey(v[i])) != UINT_MAX)
				value[p] = v[i];
			else
				init(&v[i]);
		}
	}
};

/*
 * class afset - Assosiative fixed set
 *
 * Instead of having the key as a base class of the values, these
 * classes have seperate keys and values. Now we have to manage
 * two arrays of things...
 */
template <class K, class V>
class afset : public set {
public:
	// Make one...
	inline explicit afset(const unsigned s, const unsigned b)
	  : set(b), key(0), value(0) {
		allocate(s);
	}
	inline explicit afset(const unsigned s,	const K * k = 0,
			const V * v = 0, const unsigned b = DFLT_BLK)
	  : set(b), key(0), value(0) {
		allocate(s);
		copy(s,k,v,true);
		allocate(size);
	}
	// Copy one...
	inline explicit afset(const afset<K,V> & v)
	  : set(v), key(0), value(0) {
		allocate(v.size);
		copy(v.size,v.key,v.value,false);
	}
	// Kill one...
	inline ~afset() {
		deallocate();
	}

	// Return an index for a key.
	inline unsigned haskey(const K & k) const {
		for (unsigned i = 0; i < size; i++) {
			if (key[i] == k)
				return i;
		}
		return UINT_MAX;
	}
	// Assignement operator.
	inline afset<K,V> & operator= (const afset<K,V> & v) {
		deallocate();
		allocate(v.size);
		copy(v.size,v.key,v.value,false);
		return *this;
	}
	// Behave like an indexed array...
	inline V & operator[] (const K & k) const {
		unsigned p = haskey(k);
		assert(inbounds(p));
		return value[p];
	}
	// Linear access.
	inline K & getkey(const unsigned p) const {
		assert(inbounds(p));
		return key[p];
	}
	// More linear access.
	inline V & getvalue(const unsigned p) const {
		assert(inbounds(p));
		return value[p];
	}
protected:
	K * key;                   // The keys.
	V * value;                 // Take a guess...
	struct _K {                // Placement new wrapper
		K _k;
		explicit _K(const K & k) : _k(k) {}
		void * operator new(size_t, void * p) {
			return p;
		} 
		void operator delete(void * , void *) {
		} 
	};
	struct _V {                // Placement new wrapper
		V _v;
		explicit _V() : _v() {}
		explicit _V(const V & v) : _v(v) {}
		void * operator new(size_t, void * p) {
			return p;
		} 
		void operator delete(void * , void *) {
		} 
	};

	// Make some space...
	void allocate(const unsigned s) {
		unsigned b = bufsize(s);
		if (b == buffer)
			return;
		if (b == 0) {
			free(key);
			key = 0;
			free(value);
			value = 0;
			buffer = 0;
			return;
		}
		K * ktemp = static_cast<K *>(realloc(key,b*sizeof(K)));
		V * vtemp = static_cast<V *>(realloc(value,b*sizeof(V)));
		if (ktemp == 0 || vtemp == 0)
			throw std::bad_alloc();
		key = ktemp;
		value = vtemp;
		buffer = b;
	}
	void deallocate() {
		if (key) {
			for (unsigned i = 0; i < size; i++)
				key[i].~K();
			free(key);
			key = 0;
		}
		if (value) {
			for (unsigned i = 0; i < size; i++)
				value[i].~V();
			free(value);
			value = 0;
		}
		size = 0;
		buffer = 0;
	}
	void init(const K * k, const V * v) {
		assert(k != 0);
		if (v)
			new(&value[size]) _V(*v);
		else
			new(&value[size]) _V();
		new(&key[size++]) _K(*k);
	}
	void copy(const unsigned s, const K * k, const V * v, bool checkdups) {
		for (unsigned i = 0, p; k && i < s; i++) {
			if (checkdups && (p = haskey(k[i])) != UINT_MAX)
				value[p] = (v ? v[i] : V());
			else
				init(&k[i],(v ? &v[i] : 0));
		}
	}
	// Don't allow yobos to make empty sets...
	inline explicit afset()
	  : set(), key(0), value(0) {
		allocate(size);
	}
};

/*
 * class asset - Associative sizable set
 */
template <class K, class V>
class asset : public afset<K,V> {
public:
	// Empty sizeable sets are OK.
	inline explicit asset()
	  : afset<K,V>() {
	}
	inline explicit asset(const unsigned s, const unsigned b)
	  : afset<K,V>(s,b) {
	}
	// Simple constructor.
	inline explicit asset(const unsigned s, const K * k = 0,
			const V * v = 0, const unsigned b = DFLT_BLK)
	  : afset<K,V>(s,k,v,b) {
	}
	// Copy constructor.
	inline explicit asset(const afset<K,V> & v)
	  : afset<K,V>(v) {
	}
	// Destructor.
	inline ~asset() {
	}

	// Add one...
	inline void add(const K & k,const V & v) {
		add(&k,&v,1);
	}
	// Add a whole set, at the end.
	inline void add(const afset<K,V> & v) {
		add(v.key,v.value,v.size);
	}
	// Add a whole bunch...
	void add(const K * k, const V * v = 0,
			const unsigned s = 1) {
		assert(s > 0 && k != 0);
		this->allocate(this->size+s);
		this->copy(s,k,v,true);
		this->allocate(this->size);
	}
	// Now start removing them...
	void remove(const K & k) {
		unsigned p = this->haskey(k);
		assert(p != UINT_MAX);
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
	inline void replace(const K & k, const V & v) {
		unsigned p = this->haskey(k);
		assert(p != UINT_MAX);
		this->value[p] = v;
	}
	inline void empty() {
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
	inline explicit aoset()
	  : asset<K,V>() {
	}
	inline explicit aoset(const unsigned s, const unsigned b)
	  : asset<K,V>(s,b) {
	}
	inline explicit aoset(const unsigned s, const K * k = 0,
			const V * v = 0, const unsigned b = DFLT_BLK)
	  : asset<K,V>(s,k,v,b) {
		if (k)
			sort();
	}
	inline explicit aoset(const afset<K,V> & v)
	  : asset<K,V>(v) {
		sort();
	}
	inline ~aoset() {
	}

	inline void sort() {
		qsort(0,this->size);
	}
protected:
	void qsort(const unsigned l, const unsigned r) {
		if (r <= l)
			return;
		unsigned i, j, k, p = l+(r-1-l)/2;
		if (this->key[l] > this->key[p]) {
			swap(this->key[l],this->key[p]);
			swap(this->value[l],this->value[p]);
		}
		if (this->key[l] > this->key[r-1]) {
			swap(this->key[l],this->key[r-1]);
			swap(this->value[l],this->value[r-1]);
		}
		if (this->key[p] > this->key[r-1]) {
			swap(this->key[p],this->key[r-1]);
			swap(this->value[p],this->value[r-1]);
		}
		if (r-1-l <= 2)
			return;
		for (i = l, j = r-1, k = 0; ;
				p = (p==i?j++:(p==j?i--:p)), k++) {
			while (++i < p && !(this->key[i] > this->key[p])) {};
			while (p < --j && !(this->key[p] > this->key[j])) {};
			if (i >= j)
				break;
			swap(this->key[i],this->key[j]);
			swap(this->value[i],this->value[j]);
		}
		if (k == 0) {
			isort(l,p);
			isort(p+1,r);
		} else {
			qsort(l,p);
			qsort(p+1,r);
		}
	}
	void isort(const unsigned l, const unsigned r) {
		for (unsigned i = l; l < r && i < r-1; i++) {
			for (unsigned j = i+1; j > l
					&& this->key[j-1] > this->key[j]; j--) {
				swap(this->key[j-1],this->key[j]);
				swap(this->value[j-1],this->value[j]);
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
	inline explicit avoset()
	  : asset<K,V>() {
	}
	inline explicit avoset(const unsigned s, const unsigned b)
	  : asset<K,V>(s,b) {
	}
	inline explicit avoset(const unsigned s, const K * k = 0,
			const V * v = 0, const unsigned b = DFLT_BLK)
	  : asset<K,V>(s,k,v,b) {
		if (k)
			sort();
	}
	inline explicit avoset(const afset<K,V> & v)
	  : asset<K,V>(v) {
		sort();
	}
	inline ~avoset() {
	}

	inline void sort() {
		qsort(0,this->size);
	}
protected:
	void qsort(const unsigned l, const unsigned r) {
		if (r <= l)
			return;
		unsigned i, j, k, p = l+(r-1-l)/2;
		if (this->value[l] > this->value[p]) {
			swap(this->key[l],this->key[p]);
			swap(this->value[l],this->value[p]);
		}
		if (this->value[l] > this->value[r-1]) {
			swap(this->key[l],this->key[r-1]);
			swap(this->value[l],this->value[r-1]);
		}
		if (this->value[p] > this->value[r-1]) {
			swap(this->key[p],this->key[r-1]);
			swap(this->value[p],this->value[r-1]);
		}
		if (r-1-l <= 2)
			return;
		for (i = l, j = r-1, k = 0; ;
				p = (p==i?j++:(p==j?i--:p)), k++) {
			while (++i < p && !(this->value[i] > this->value[p])) {};
			while (p < --j && !(this->value[p] > this->value[j])) {};
			if (i >= j)
				break;
			swap(this->key[i],this->key[j]);
			swap(this->value[i],this->value[j]);
		}
		if (k == 0) {
			isort(l,p);
			isort(p+1,r);
		} else {
			qsort(l,p);
			qsort(p+1,r);
		}
	}
	void isort(const unsigned l, const unsigned r) {
		for (unsigned i = l; l < r && i < r-1; i++) {
			for (unsigned j = i+1; j > l
					&& this->value[j-1] > this->value[j]; j--) {
				swap(this->key[j-1],this->key[j]);
				swap(this->value[j-1],this->value[j]);
			}
		}
	}
};

#endif // SET_H
