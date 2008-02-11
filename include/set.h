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
	// The first element. Nice for loops...
	inline const int start() const {
		return 0;
	}
	// The last element. Also nice for loops...
	inline const int end() const {
		return size-1;
	}
	// The length. Nice for lots of things...
	inline const int length() const {
		return size;
	}
	// Check if an index is within the set...
	inline bool inbounds(const int p) const {
		return (p >= 0 && p < size ? true : false);
	}

protected:
	int size;                  // The size of the set...
	int buffer;                // The allocated buffer size...

	// Simple constructor...
	inline explicit set(const int b = DFLT_BLK)
	  : size(0), buffer(0), block(b > 1 ? b : 1) {
	}
	inline explicit set(const set & s)
	  : size(0), buffer(0), block(s.block) {
	}
	inline ~set() {
	}
	// Calculate the buffer size.
	inline int bufsize(int s) {
		assert(s >= 0);
		while (s > 8*block)
			block *= 8;
		//while (64*s < block)
		//	block /= 8;
		return block*(s/block+(s%block?1:0));
	}

private:
	int block;                 // The minimum block size.
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
	inline explicit fset(const int s, const int b = DFLT_BLK,
			const V * v = 0)
	  : set(b), value(0) {
		if (!allocate(s))
			return;
		copy(s,v);
	}
	// Copy constructor.
	inline explicit fset(const fset<V> & v)
	  : set(v), value(0) {
		if (!allocate(v.size))
			return;
		copy(v.size,v.value);
	}
	// Wow, a destructor...
	inline ~fset() {
		deallocate();
	}

	// Assignment operator.
	inline fset<V> & operator= (const fset<V> & v) {
		deallocate();
		if (!allocate(v.size))
			return *this;
		copy(v.size,v.value);
		return *this;
	}
	// Allow the size to be changed, even for 'fixed' sets.
	bool resize(const int s) {
		deallocate();
		if (!allocate(s))
			return false;
		copy(s,0);
		return true;
	}
	// Behave like an array. Zero indexed.
	V & operator[] (const int p) const {
		assert(inbounds(p));
		return value[p];
	}

protected:
	V * value;					// The buffer.
	struct _V {                 // Placement new wrapper
		V _v;
		explicit _V() : _v() {}
		explicit _V(const V & v) : _v(v) {}
		void *operator new(size_t, void * p) {
			return p;
		} 
	};

	// Since this is a fixed size set, hide the allocator...
	bool allocate(const int s) {
		int b = bufsize(s);
		if (b == buffer)
			return true;
		if (b == 0) {
			free(value);
			value = 0;
			buffer = 0;
			return true;
		}
		V * temp = static_cast<V *>(realloc(value,b*sizeof(V)));
		if (temp == 0) {
			event_msg(EVENT_ERROR,"Out of memory in fset::allocate()!");
			return false;
		}
		value = temp;
		buffer = b;
		return true;
	}
	void deallocate() {
		if (value) {
			for (int i = 0; i < size; i++)
				value[i].~V();
			free(value);
			value = 0;
		}
		size = 0;
		buffer = 0;
	}
	void init(const int i, const V * v) {
		if (v != 0)
			new(&value[i]) _V(*v);
		else
			new(&value[i]) _V();
	}
	void copy(const int s, const V * v) {
		for (int i = 0; i < s; i++)
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
inline void fset<int>::init(const int i, const int * v) {
	if (v)
		value[i] = *v;
}
template<>
inline void fset<double>::init(const int i, const double * v) {
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
	// Simple constructor.
	inline explicit sset(const int s, const int b = DFLT_BLK,
			const V * v = 0)
	  : fset<V>(s,b,v) {
	}
	// Copy constructor.
	inline explicit sset(const fset<V> & v)
	  : fset<V>(v) {
	}
	// We let someone else clean up...
	inline ~sset() {
	}

	// Add one element, at the end.
	inline bool add(const V & v) {
		return add(this->size,&v,1);
	}
	// Add a whole set, at the end.
	inline bool add(const fset<V> & v) {
		return add(this->size,v.value,v.size);
	}
	// Add an array, at the end.
	inline bool add(const V * v, const int s = 1) {
		return add(this->size,v,s);
	}
	// Insert one element at position p.
	inline bool add(const int p, const V & v) {
		return add(p,&v,1);
	}
	// Insert a set at position p.
	inline bool add(const int p, const fset<V> & v) {
		return add(p,v.value,v.size);
	}
	// Add an array at position p. (Actually do the work too).
	bool add(int p, const V * v, const int s = 1) {
		int i;
		if (p < 0 || s <= 0 || v == 0)
			return false;
		if (p >= this->size) {
			if (!this->allocate(p+s))
				return false;
			for (i = this->size; i < p-1; i++)
				this->init(i,0);
			this->size = p+s;
		} else {
			if (!allocate(this->size+s))
				return false;
			memmove(&this->value[p+s],&this->value[p],
					(this->size-p)*sizeof(V));
			this->size += s;
		}
		for (i = 0; i < s; i++)
			this->init(p+i,&v[i]);
		return true;
	}
	// Remove the last element.
	inline bool remove() {
		return remove(this->size-1,1);
	}
	// Remove s elements, starting at position p.
	bool remove(int p, int s = 1) {
		if (!this->inbounds(p) || s == 0)
			return false;
		s = (p+s > this->size ? this->size-p : s);
		for (int i = 0; i < s; i++)
			this->value[p+i].~V();
		if (this->size-p > s)
			memmove(&this->value[p],&this->value[p+s],
				(this->size-p-s)*sizeof(V));
		this->size -= s;
		return allocate(this->size);
	}
	inline bool empty() {
		this->deallocate();
		return this->allocate(0);
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
	// Simple constructor.
	inline explicit oset(const int s, const int b = DFLT_BLK,
			const V * v = 0)
	  : sset<V>(s,b,v) {
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
		qsort(0,this->size-1);
	}
	// Do a lookup, and return -1 if the value is not found.
	// You must sort the set first!
	inline int findvalue(const V & v) const {
		int l = 0, r = this->size;
		while (l < r) {
			int i = l + (r-l)/2;
			if (this->value[i] < v)
				l = i+1;
			else
				r = i;
		}
		if (l < this->size && this->value[l] == v)
			return l;
		else
			return -1;
	}

protected:
	// A highly optimised quick sort. Don't touch...
	void qsort(const int l, const int r) {
		if (r > l) {
			int i, j, k, p = l+(r-l)/2;
			if (this->value[l] > this->value[p])
				swap(this->value[l],this->value[p]);
			if (this->value[l] > this->value[r])
				swap(this->value[l],this->value[r]);
			if (this->value[p] > this->value[r])
				swap(this->value[p],this->value[r]);
			for (i = l, j = r, k = 0; ;
								p = (p==i?j++:(p==j?i--:p)), k++) {
				while (++i < p && !(this->value[i] > this->value[p]));
				while (p < --j && !(this->value[p] > this->value[j]));
				if (i >= j)
					break;
				swap(this->value[i],this->value[j]);
			}
			if (k == 0) {
				isort(l,p-1);
				isort(p+1,r);
			} else {
				qsort(l,p-1);
				qsort(p+1,r);
			}	
		}
	}
	// Insertion sort for if the set looks sorted already.
	void isort(const int l, const int r) {
		for (int i = l; l < r && i < r; i++) {
			for (int j = i; i >= l && j >= l
					&& this->value[j] > this->value[j+1]; j--)
				swap(this->value[j],this->value[j+1]);
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
	inline explicit cset(const int s, const int b = DFLT_BLK,
			const V * v = 0)
	  : oset<V>(s,b,v) {
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
		this->qsort(0,this->size-1);
		compact();
	}

protected:
	void compact() {
		int i, j, s;
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
		allocate(this->size);
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
	// Basic constructor
	inline explicit iset(const int s, const int b = DFLT_BLK,
			const V * v = 0)
	  : set(b), idx(0), value(0) {
		if (!allocate(s))
			return;
		copy(s,v,true);
		allocate(size);
	}
	// Copy constuctor.
	inline explicit iset(const iset<V> & v)
	  : set(v), idx(0), value(0) {
		if (!allocate(v.size))
			return;
		copy(v.size,v.value,false);
	}
	// Copy from an fset.
	inline explicit iset(const fset<V> & v)
	  : set(v), idx(0), value(0) {
		if (!allocate(v.size))
			return;
		copy(v.size,v.value,true);
		allocate(size);
	}
	// Clean up.
	inline ~iset() {
		deallocate();
	}

	// Do a lookup, and return -1 if the value is not found.
	inline int hasvalue(const V & v) const {
		int p = findvalue(v);
		if (p < size && value[idx[p]] == v)
			return idx[p];
		else
			return -1;
	}
	// Assignment operator.
	inline iset<V> & operator= (const iset<V> & v) {
		deallocate();
		if (!allocate(v.size))
			return *this;
		copy(v.size,v.value,false);
		return *this;
	}
	inline iset<V> & operator= (const fset<V> & v) {
		deallocate();
		if (!allocate(v.size))
			return *this;
		copy(v.size,v.value,true);
		allocate(size);
		return *this;
	}
	// Only integer keys make sense.
	inline const V & operator[] (const int p) const {
		assert(inbounds(p));
		return value[p];
	}
	// Allow sorted acess.
	inline const V & getindex(const int i) const {
		assert(inbounds(i));
		return value[idx[i]];
	}
	// Get the position of an element in the sort.
	inline int getorder(const int i) const {
		assert(inbounds(i));
		return findvalue(value[i]);
	}

	// Add one value at the end.
	inline bool add(const V & v) {
		return add(&v,1);
	}
	// Add a whole set, at the end.
	inline bool add(const iset<V> & v) {
		return add(v.value,v.size);
	}
	inline bool add(const fset<V> & v) {
		return add(v.value,v.size);
	}
	// Add an array of values... There is no point in
	// a position based addition. 
	bool add(const V * v, const int s = 1) {
		if (s <= 0 || v == 0)
			return false;
		if (!allocate(size+s))
			return false;
		copy(s,v,true);
		return allocate(size);
	}
	// Now start removing them...
	bool remove(const V & v) {
		int p = hasvalue(v), q = -1;
		if (p == -1)
			return false;
		value[p].~V();
		if (p < --size)
			memmove(&value[p],&value[p+1],(size-p)*sizeof(V));
		for (int i = 0; i <= size; i++) {
			if (idx[i] == p)
				q = i;
			else if (idx[i] > p)
				idx[i]--;
		}
		if (q < size)
			memmove(&idx[q],&idx[q+1],(size-q)*sizeof(int));
		return allocate(size);
	}
	// Replace key/value with another.
	inline bool replace(const V & v) {
		int p = hasvalue(v);
		if (p == -1)
			return false;
		value[p] = v;
		return true;
	}
	inline bool empty() {
		deallocate();
		return allocate(0);
	}

protected:
	int * idx;                 // The index.
	V * value;                 // Take a guess...
	struct _V {                // Placement new wrapper
		V _v;
		explicit _V() : _v() {}
		explicit _V(const V & v) : _v(v) {}
		void *operator new(size_t, void * p) {
			return p;
		} 
	};

	// Make some space...
	bool allocate(const int s) {
		int b = bufsize(s);
		if (b == buffer)
			return true;
		if (b == 0) {
			free(idx);
			idx = 0;
			free(value);
			value = 0;
			buffer = 0;
			return true;
		}
		int * itemp = static_cast<int *>(realloc(idx,b*sizeof(int)));
		V * vtemp = static_cast<V *>(realloc(value,b*sizeof(V)));
		if (itemp == 0 || vtemp == 0) {
			event_msg(EVENT_ERROR,"Out of memory in iset::allocate()!");
			return false;
		}
		idx = itemp;
		value = vtemp;
		buffer = b;
		return true;
	}
	void deallocate() {
		if (idx) {
			free(idx);
			idx = 0;
		}
		if (value) {
			for (int i = 0; i < size; i++)
				value[i].~V();
			free(value);
			value = 0;
		}
		size = 0;
		buffer = 0;
	}
	// Find the position which is either equal or greater...
	inline int findvalue(const V & v) const {
		int l = 0, r = size;
		while (l < r) {
			int i = l + (r-l)/2;
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
		int p = findvalue(*v);
		if (p < size)
			memmove(&idx[p+1],&idx[p],(size-p)*sizeof(int));
		idx[p] = size++;
	}
	void copy(const int s, const V * v, bool checkdups) {
		for (int i = 0, p; v && i < s; i++) {
			if (checkdups && (p = hasvalue(v[i])) != -1)
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
	inline explicit kfset(const int s, const int b = DFLT_BLK,
			const V * v = 0)
	  : set(b), value(0) {
		if (!allocate(s))
			return;
		copy(s,v,true);
		allocate(size);
	}
	// Copy constuctor.
	inline explicit kfset(const kfset<K,V> & v)
	  : set(v), value(0) {
		if (!allocate(v.size))
			return;
		copy(v.size,v.value,false);
	}
	// Clean up.
	inline ~kfset() {
		deallocate();
	}

	// Do a key lookup, and return -1 if the key is not found.
	inline int haskey(const K & k) const {
		for (int i = 0; i < size; i++) {
			if (static_cast<K &>(value[i]) == k)
				return i;
		}
		return -1;
	}
	// Assignment operator.
	inline kfset<K,V> & operator= (const kfset<K,V> & v) {
		deallocate();
		if (!allocate(v.size))
			return *this;
		copy(v.size,v.value,false);
		return *this;
	}
	// Return data based on a key lookup.
	inline V & operator[] (const K & k) const {
		int p = haskey(k);
		assert(inbounds(p));
		return value[p];
	}
	// Allow integer keys.
	inline V & operator[] (const int p) const {
		assert(inbounds(p));
		return value[p];
	}

protected:
	V * value;                 // The data.
	struct _V {                // Placement new wrapper
		V _v;
		explicit _V() : _v() {}
		explicit _V(const V & v) : _v(v) {}
		void *operator new(size_t, void * p) {
			return p;
		} 
	};

	// Hide the allocation function.
	bool allocate(const int s) {
		int b = bufsize(s);
		if (b == buffer)
			return true;
		if (b == 0) {
			free(value);
			value = 0;
			buffer = 0;
			return true;
		}
		V * temp = static_cast<V *>(realloc(value,b*sizeof(V)));
		if (temp == 0) {
			event_msg(EVENT_ERROR,"Out of memory in kfset::allocate()!");
			return false;
		}
		value = temp;
		buffer = b;
		return true;
	}
	void deallocate() {
		if (value) {
			for (int i = 0; i < size; i++)
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
	void copy(const int s, const V * v, bool checkdups) {
		for (int i = 0, p; v && i < s; i++) {
			if (checkdups && (p = haskey(v[i])) != -1)
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
	inline explicit ksset(const int s, const int b = DFLT_BLK,
			const V * v = 0)
	  : kfset<K,V>(s,b,v) {
	}
	inline explicit ksset(const kfset<K,V> & v)
	  : kfset<K,V>(v) {
	}
	inline ~ksset() {
	}

	// Add one value at the end.
	inline bool add(const V & v) {
		return add(&v,1);
	}
	// Add a whole set, at the end.
	inline bool add(const kfset<K,V> & v) {
		return add(v.value,v.size);
	}
	// Add an array of values... There is no point in
	// a position based addition. 
	bool add(const V * v, const int s = 1) {
		if (s <= 0 || v == 0)
			return false;
		if (!allocate(this->size+s))
			return false;
		copy(s,v,true);
		return allocate(this->size);
	}
	// Remove based on key.
	bool remove(const K & k) {
		int p = haskey(k);
		if (p == -1)
			return false;
		this->value[p].~V();
		if (p < --this->size)
			memmove(&this->value[p],&this->value[p+1],
				(this->size-p)*sizeof(V));
		return allocate(this->size);
	}
	// Replace key/value with another.
	inline bool replace(const V & v) {
		int p = haskey(v);
		if (p == -1)
			return false;
		this->value[p] = v;
		return true;
	}
	inline bool empty() {
		this->deallocate();
		return this->allocate(0);
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
	inline explicit koset(const int s, const int b = DFLT_BLK,
			const V * v = 0)
	  : ksset<K,V>(s,b,v) {
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
		qsort(0,this->size-1);
	}
protected:
	void qsort(const int l, const int r) {
		if (r > l) {
			int i, j, k, p = l+(r-l)/2;
			if (this->value[l] > this->value[p])
				swap(this->value[l],this->value[p]);
			if (this->value[l] > this->value[r])
				swap(this->value[l],this->value[r]);
			if (this->value[p] > this->value[r])
				swap(this->value[p],this->value[r]);
			for (i = l, j = r, k = 0; ;
							p = (p==i?j++:(p==j?i--:p)), k++) {
				while (++i < p && !(this->value[i] > this->value[p]));
				while (p < --j && !(this->value[p] > this->value[j]));
				if (i >= j)
					break;
				swap(this->value[i],this->value[j]);
			}
			if (k == 0) {
				isort(l,p-1);
				isort(p+1,r);
			} else {
				qsort(l,p-1);
				qsort(p+1,r);
			}
		}
	}
	void isort(const int l, const int r) {
		for (int i = l; l < r && i < r; i++) {
			for (int j = i; i >= l && j >= l
					 && this->value[j] > this->value[j+1]; j--) {
				swap(this->value[j],this->value[j+1]);
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
	// Basic constructor
	inline explicit kiset(const int s, const int b = DFLT_BLK,
			const V * v = 0)
	  : set(b), idx(0), value(0) {
		if (!allocate(s))
			return;
		copy(s,v,true);
		allocate(size);
	}
	// Copy constuctor.
	inline explicit kiset(const kiset<K,V> & v)
	  : set(v), idx(0), value(0) {
		if (!allocate(v.size))
			return;
		copy(v.size,v.value,false);
	}
	// Copy from a kfset.
	inline explicit kiset(const kfset<K,V> & v)
	  : set(v), idx(0), value(0) {
		if (!allocate(v.size))
			return;
		copy(v.size,v.value,true);
		allocate(size);
	}
	// Clean up.
	inline ~kiset() {
		deallocate();
	}

	// Do a key lookup, and return -1 if the key is not found.
	inline int haskey(const K & k) const {
		int p = findkey(k);
		if (p < size && static_cast<K &>(value[idx[p]]) == k)
			return idx[p];
		else
			return -1;
	}
	// Assignment operator.
	inline kiset<K,V> & operator= (const kiset<K,V> & v) {
		deallocate();
		if (!allocate(v.size))
			return *this;
		copy(v.size,v.value,false);
		return *this;
	}
	inline kiset<K,V> & operator= (const kfset<K,V> & v) {
		deallocate();
		if (!allocate(v.size))
			return *this;
		copy(v.size,v.value,true);
		allocate(size);
		return *this;
	}
	// Return data based on a key lookup.
	inline const V & operator[] (const K & k) const {
		int p = haskey(k);
		assert(inbounds(p));
		return value[p];
	}
	// Allow integer keys.
	inline const V & operator[] (const int p) const {
		assert(inbounds(p));
		return value[p];
	}
	// Allow sorted acess.
	inline const V & getindex(const int i) const {
		assert(inbounds(i));
		return value[idx[i]];
	}
	// Get the position of an element in the sort.
	inline int getorder(const int i) const {
		assert(inbounds(i));
		return findkey(value[i]);
	}

	// Add one value at the end.
	inline bool add(const V & v) {
		return add(&v,1);
	}
	// Add a whole set, at the end.
	inline bool add(const kfset<K,V> & v) {
		return add(v.value,v.size);
	}
	// Add an array of values... There is no point in
	// a position based addition. 
	bool add(const V * v, const int s = 1) {
		if (s <= 0 || v == 0)
			return false;
		if (!allocate(size+s))
			return false;
		copy(s,v,true);
		return allocate(size);
	}
	// Now start removing them...
	bool remove(const K & k) {
		int p = haskey(k), q = -1;
		if (p == -1)
			return false;
		value[p].~V();
		if (p < --size)
			memmove(&value[p],&value[p+1],(size-p)*sizeof(V));
		for (int i = 0; i <= size; i++) {
			if (idx[i] == p)
				q = i;
			else if (idx[i] > p)
				idx[i]--;
		}
		if (q < size)
			memmove(&idx[q],&idx[q+1],(size-q)*sizeof(int));
		return allocate(size);
	}
	// Replace key/value with another.
	inline bool replace(const V & v) {
		int p = haskey(v);
		if (p == -1)
			return false;
		value[p] = v;
		return true;
	}
	inline bool empty() {
		deallocate();
		return allocate(0);
	}

protected:
	int * idx;                 // The index.
	V * value;                 // Take a guess...
	struct _V {                // Placement new wrapper
		V _v;
		explicit _V() : _v() {}
		explicit _V(const V & v) : _v(v) {}
		void *operator new(size_t, void * p) {
			return p;
		} 
	};

	// Make some space...
	bool allocate(const int s) {
		int b = bufsize(s);
		if (b == buffer)
			return true;
		if (b == 0) {
			free(idx);
			idx = 0;
			free(value);
			value = 0;
			buffer = 0;
			return true;
		}
		int * itemp = static_cast<int *>(realloc(idx,b*sizeof(int)));
		V * vtemp = static_cast<V *>(realloc(value,b*sizeof(V)));
		if (itemp == 0 || vtemp == 0) {
			event_msg(EVENT_ERROR,"Out of memory in kiset::allocate()!");
			return false;
		}
		idx = itemp;
		value = vtemp;
		buffer = b;
		return true;
	}
	void deallocate() {
		if (idx) {
			free(idx);
			idx = 0;
		}
		if (value) {
			for (int i = 0; i <= size; i++)
				value[i].~V();
			free(value);
			value = 0;
		}
		size = 0;
		buffer = 0;
	}
	// Find the position which is either equal or greater...
	inline int findkey(const K & k) const {
		int l = 0, r = size;
		while (l < r) {
			int i = l + (r-l)/2;
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
		int p = findkey(static_cast<const K &>(*v));
		if (p < size)
			memmove(&idx[p+1],&idx[p],(size-p)*sizeof(int));
		idx[p] = size++;
	}
	void copy(const int s, const V * v, bool checkdups) {
		for (int i = 0, p; v && i < s; i++) {
			if (checkdups && (p = haskey(v[i])) != -1)
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
	inline explicit afset(const int s, const int b = DFLT_BLK,
			const K * k = 0, const V * v = 0)
	  : set(b), key(0), value(0) {
		if (!allocate(s))
			return;
		copy(s,k,v,true);
		allocate(size);
	}
	// Copy one...
	inline explicit afset(const afset<K,V> & v)
	  : set(v), key(0), value(0) {
		if (!allocate(v.size))
			return;
		copy(v.size,v.key,v.value,false);
	}
	// Kill one...
	inline ~afset() {
		deallocate();
	}

	// Return an index for a key.
	inline int haskey(const K & k) const {
		for (int i = 0; i < size; i++) {
			if (key[i] == k)
				return i;
		}
		return -1;
	}
	// Assignement operator.
	inline afset<K,V> & operator= (const afset<K,V> & v) {
		deallocate();
		if (!allocate(v.size))
			return *this;
		copy(v.size,v.key,v.value,false);
		return *this;
	}
	// Behave like an indexed array...
	inline V & operator[] (const K & k) const {
		int p = haskey(k);
		assert(inbounds(p));
		return value[p];
	}
	// Linear access.
	inline K & getkey(const int p) const {
		assert(inbounds(p));
		return key[p];
	}
	// More linear access.
	inline V & getvalue(const int p) const {
		assert(inbounds(p));
		return value[p];
	}
protected:
	K * key;                   // The keys.
	V * value;                 // Take a guess...
	struct _K {                // Placement new wrapper
		K _k;
		explicit _K() : _k() {}
		explicit _K(const K & k) : _k(k) {}
		void *operator new(size_t, void * p) {
			return p;
		} 
	};
	struct _V {                // Placement new wrapper
		V _v;
		explicit _V() : _v() {}
		explicit _V(const V & v) : _v(v) {}
		void *operator new(size_t, void * p) {
			return p;
		} 
	};

	// Make some space...
	bool allocate(const int s) {
		int b = bufsize(s);
		if (b == buffer)
			return true;
		if (b == 0) {
			free(key);
			key = 0;
			free(value);
			value = 0;
			buffer = 0;
			return true;
		}
		K * ktemp = static_cast<K *>(realloc(key,b*sizeof(K)));
		V * vtemp = static_cast<V *>(realloc(value,b*sizeof(V)));
		if (ktemp == 0 || vtemp == 0) {
			event_msg(EVENT_ERROR,"Out of memory in afset::allocate()!");
			return false;
		}
		key = ktemp;
		value = vtemp;
		buffer = b;
		return true;

	}
	void deallocate() {
		if (key) {
			for (int i = 0; i < size; i++)
				key[i].~K();
			free(key);
			key = 0;
		}
		if (value) {
			for (int i = 0; i < size; i++)
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
	void copy(const int s, const K * k, const V * v, bool checkdups) {
		for (int i = 0, p; k && i < s; i++) {
			if (checkdups && (p = haskey(k[i])) != -1)
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
	// Simple constructor.
	inline explicit asset(const int s, const int b = DFLT_BLK,
			const K * k = 0, const V * v = 0)
	  : afset<K,V>(s,b,k,v) {
	}
	// Copy constructor.
	inline explicit asset(const afset<K,V> & v)
	  : afset<K,V>(v) {
	}
	// Destructor.
	inline ~asset() {
	}

	// Add one...
	inline bool add(const K & k,const V & v) {
		return add(&k,&v,1);
	}
	// Add a whole set, at the end.
	inline bool add(const afset<K,V> & v) {
		return add(v.key,v.value,v.size);
	}
	// Add a whole bunch...
	bool add(const K * k, const V * v = 0,
			const int s = 1) {
		if (s <= 0 || k == 0)
			return false;
		if (!allocate(this->size+s))
			return false;
		copy(s,k,v,true);
		return allocate(this->size);
	}
	// Now start removing them...
	bool remove(const K & k) {
		int p = haskey(k) + 1;
		if (p == 0)
			return false;
		this->key[p].~K();
		this->value[p].~V();
		if (p < --this->size) {
			memmove(&this->key[p],&this->key[p+1],
					(this->size-p)*sizeof(K));
			memmove(&this->value[p],&this->value[p+1],
					(this->size-p)*sizeof(V));
		}
		return allocate(this->size);
	}
	// Or replacing them...
	inline bool replace(const K & k, const V & v) {
		int p = haskey(k);
		if (p == -1)
			return false;
		this->value[p] = v;
		return true;
	}
	inline bool empty() {
		this->deallocate();
		return this->allocate(0);
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
	inline explicit aoset(const int s, const int b = DFLT_BLK,
			const K * k = 0, const V * v = 0)
	  : asset<K,V>(s,b,k,v) {
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
		qsort(0,this->size-1);
	}
protected:
	void qsort(const int l, const int r) {
		if (r > l) {
			int i, j, k, p = l+(r-l)/2;
			if (this->key[l] > this->key[p]) {
				swap(this->key[l],this->key[p]);
				swap(this->value[l],this->value[p]);
			}
			if (this->key[l] > this->key[r]) {
				swap(this->key[l],this->key[r]);
				swap(this->value[l],this->value[r]);
			}
			if (this->key[p] > this->key[r]) {
				swap(this->key[p],this->key[r]);
				swap(this->value[p],this->value[r]);
			}
			for (i = l, j = r, k = 0; ;
								p = (p==i?j++:(p==j?i--:p)), k++) {
				while (++i < p && !(this->key[i] > this->key[p]));
				while (p < --j && !(this->key[p] > this->key[j]));
				if (i >= j)
					break;
				swap(this->key[i],this->key[j]);
				swap(this->value[i],this->value[j]);
			}
			if (k == 0) {
				isort(l,p-1);
				isort(p+1,r);
			} else {
				qsort(l,p-1);
				qsort(p+1,r);
			}
		}
	}
	void isort(const int l, const int r) {
		for (int i = l; l < r && i < r; i++) {
			for (int j = i; i >= l && j >= l
					&& this->key[j] > this->key[j+1]; j--) {
				swap(this->key[j],this->key[j+1]);
				swap(this->value[j],this->value[j+1]);
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
	inline explicit avoset(const int s, const int b = DFLT_BLK,
			const K * k = 0, const V * v = 0)
	  : asset<K,V>(s,b,k,v) {
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
		qsort(0,this->size-1);
	}
protected:
	void qsort(const int l, const int r) {
		if (r > l) {
			int i, j, k, p = l+(r-l)/2;
			if (this->value[l] > this->value[p]) {
				swap(this->key[l],this->key[p]);
				swap(this->value[l],this->value[p]);
			}
			if (this->value[l] > this->value[r]) {
				swap(this->key[l],this->key[r]);
				swap(this->value[l],this->value[r]);
			}
			if (this->value[p] > this->value[r]) {
				swap(this->key[p],this->key[r]);
				swap(this->value[p],this->value[r]);
			}
			for (i = l, j = r, k = 0; ;
								p = (p==i?j++:(p==j?i--:p)), k++) {
				while (++i < p && !(this->value[i] > this->value[p]));
				while (p < --j && !(this->value[p] > this->value[j]));
				if (i >= j)
					break;
				swap(this->key[i],this->key[j]);
				swap(this->value[i],this->value[j]);
			}
			if (k == 0) {
				isort(l,p-1);
				isort(p+1,r);
			} else {
				qsort(l,p-1);
				qsort(p+1,r);
			}
		}
	}
	void isort(const int l, const int r) {
		for (int i = l; l < r && i < r; i++) {
			for (int j = i; i >= l && j >= l
					&& this->value[j] > this->value[j+1]; j--) {
				swap(this->key[j],this->key[j+1]);
				swap(this->value[j],this->value[j+1]);
			}
		}
	}
};

#endif // SET_H
