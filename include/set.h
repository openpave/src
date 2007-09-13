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

	The Original Code is OpenPave.org Core Libraries.

	The Initial Developer of the Original Code is OpenPave.org.

	Portions Copyright (C) 2006 OpenPave.org.

	Contributor(s): Jeremy Lea <reg@openpave.org>.

	Purpose:
		This header implements a number of templated C++ classes for
		various types of sets. The code is fairly compact, but not
		optimised for big sets (particularly sorting).

	Design:
		The sets are all one based, not zero based like C arrays, which
		make them more usefull logically, and allows the zero element to
		be used as a default. However, the operator[] functions are
		zero based, so as not to be confusing to consumers.

		There are two types of sets. The first are normal sets with
		one variable. Then there are associative sets, which store
		a key and a value.

		The heirarchy goes from a fixed size set, to a resizeable set,
		to a sorted resizable set.

		All of the templated sets take a minimum block size, which is
		defaults to 100. Major performance enhancements will be seen
		on resizable sets if this value is bigger than the maximum size
		of the set, because the memory is reused.
		
	Status:
		The swaps in the sorts should be replaced by memmove.  It should
		also learn to throw expections...

	History:
		1994       - Created by Jeremy Lea <reg@openpave.org>
		2002/01/23 - Modifications for use on a modern C++ compiler
		2007/09/11 - Modified to use realloc and placement new.

**************************************************************************/

#ifndef __SET_H
#define __SET_H

#include "event.h"
#include "mathplus.h"
#include "string.h"

#define DFLT_BLK	100

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
		return 1;
	}
	// The last element. Also nice for loops...
	inline const int end() const {
		return size;
	}
	// The length. Nice for lots of things...
	inline const int length() const {
		return size;
	}
	// Check if an index is within the set...
	inline bool inbounds(const int p) const {
		return (p >= 1 && p <= size ? true : false);
	}

protected:
	int size;                           // The size of the set...

	// Simple constructor...
	inline explicit set(const int s = 0,
			const int b = DFLT_BLK) {
		size = (s > 0 ? s : 0);
		block = (b > 1 ? b : 1);
	}
	inline explicit set(const set & s)
	 : size(s.size), block(s.block) {
	}
	inline ~set() {
	}

	// Calculate the buffer size.
	inline int bufsize(int s) {
		s = (s <= 0 ? 1 : s + 1);
		while (s > 10*block)
			block *= 10;
		while (100*s < block)
			block /= 10;
		return block*(s/block+(s%block?1:0));
	}

private:
	int block;                          // The minimum block size.
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
	  : set(s,b), value(0) {
		if (!allocate(size)) {
			size = 0;
			return;
		}
		for (int i = 0; v && i < size; i++)
			value[i+1] = v[i];
	}
	// Copy constructor.
	inline explicit fset(const fset<V> & v)
	  : set(v), value(0) {
		if (!allocate(size)) {
			size = 0;
			return;
		}
		for (int i = 0; i <= size; i++)
			value[i] = v.value[i];
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
		size = v.size;
		for (int i = 0; i <= size; i++)
			value[i] = v.value[i];
		return *this;
	}
	// Set the default element if constructor didn't do a good job.
	void setdefault(const V & d) {
		value[0] = d;
	}
	// Behave like an array. Zero indexed.
	V & operator[] (const int p) const {
		return (inbounds(p+1) ? value[p+1] : value[0]);
	}

protected:
	V * value;							// The buffer.
	struct _V {                         // Placement new wrapper
		V _v;
		explicit _V() : _v() {}
		explicit _V(const V & v) : _v(v) {}
		void *operator new(size_t, void * p) {
			return p;
		} 
	};

	// Since this is a fixed size set, hide the allocator...
	bool allocate(const int s) {
		V * temp = static_cast<V *>(realloc(value,bufsize(s)*sizeof(V)));
		if (temp == 0) {
			event_msg(EVENT_ERROR,"Out of memory in fset::allocate()!");
			return false;
		}
		if (value == 0)
			new(&temp[0]) _V();
		value = temp;
		return true;
	}
	void deallocate() {
		if (value) {
			for (int i = 0; i <= size; i++)
				value[i].~V();
			free(value);
			value = 0;
		}
		size = 0;
	}
	// Also hide the null constructor.
	inline explicit fset()
	  : set(), value(0) {
		// If it fails it fails...
		allocate(size);
	}
};

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
		return add(this->size,&v.value[1],v.size);
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
		return add(p,&v.value[1],v.size);
	}
	// Add an array at position p. (Actually do the work too).
	bool add(int p, const V * v, const int s = 1) {
		int i;
		if (++p <= 0 || s <= 0 || v == 0)
			return false;
		if (p > this->size) {
			if (!this->allocate(p+s-1))
				return false;
			for (i = this->size+1; i < p; i++)
				this->value[i] = this->value[0];
			this->size = p+s-1;
		} else {
			if (!allocate(this->size+s))
				return false;
			memmove(&this->value[p+s],&this->value[p],
					(this->size-p+1)*sizeof(V));
			this->size += s;
		}
		for (i = 0; i < s; i++)
			this->value[p+i] = v[i];
		return true;
	}
	// Remove the last element.
	inline bool remove() {
		return remove(this->size-1,1);
	}
	// Remove s elements, starting at position p.
	bool remove(int p, int s = 1) {
		if (!this->inbounds(++p) || s == 0)
			return false;
		s = (p+s-1 > this->size ? this->size-p+1 : s);
		for (int i = 0; i < s; i++)
			this->value[p+i].~V();
		if ((this->size-p+1) > s)
			memmove(&this->value[p],&this->value[p+s],
				(this->size-p+1-s)*sizeof(V));
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
		qsort(1,this->size);
	}

protected:
	// A highly optimised quick sort. Don't touch...
	void qsort(const int l, const int r) {
		if (r > l) {
			int i, j, k, p = (l+r)/2;
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
				isort(l, p-1);
				isort(p+1, r);
			} else {
				qsort(l, p-1);
				qsort(p+1, r);
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
		this->qsort(1,this->size);
		compact();
	}

protected:
	void compact() {
		int i, j, s;
		for (i = 2, j = 1; this->size > 1 && i <= this->size; i++) {
			if (this->value[j] != this->value[i])
				j++;
			if (j != i) {
				do
					this->value[i].~V();
				while (i++ < this->size
				 && this->value[j] == this->value[i]);
				s = i;
				while (s < this->size
				 && this->value[s] != this->value[s+1])
					s++;
				if (i <= this->size) {
					memmove(&this->value[j+1],&this->value[i],
						(s-i+1)*sizeof(V));
					j = s-(i-(j+1)); i = s;
				}
			}
		}
		this->size = j;
		allocate(this->size);
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
 *
 * XXX: These classes are almost identical to the fset based classes,
 *      but I'm too stupid to work out the correct heirarchy - JDL.
 */
template <class K, class V>
class kfset : public set {
public:
	// Simple constructor.
	inline explicit kfset(const int s, const int b = DFLT_BLK,
			const V * v = 0)
	  : set(s,b), value(0) {
		int i, p;
		if (!allocate(size)) {
			size = 0;
			return;
		}
		for (i = 0, size = 0; v && i < s; i++) {
			if ((p = haskey(v[i])) != 0)
				value[p] = v[i];
			else
				value[++size] = v[i];
		}
		allocate(size);
	}
	// Copy constuctor.
	inline explicit kfset(const kfset<K,V> & v)
	  : set(v), value(0) {
		if (!allocate(size)) {
			size = 0;
			return;
		}
		for (int i = 0; i <= size; i++)
			value[i] = v.value[i];
	}
	// Clean up.
	inline ~kfset() {
		deallocate();
	}

	// Do a key lookup, and return zero if the key is not found.
	inline int haskey(const K & k) const {
		for (int i = 0; i < size; i++) {
			if (static_cast<K &>(value[i+1]) == k)
				return i+1;
		}
		return 0;
	}
	// Assignment operator.
	inline kfset<K,V> & operator= (const kfset<K,V> & v) {
		deallocate();
		if (!allocate(v.size))
			return *this;
		size = v.size;
		for (int i = 0; i <= size; i++)
			value[i] = v.value[i];
		return *this;
	}
	// Set our default elelment.
	inline void setdefault(const V & d) {
		value[0] = d;
	}
	// Return data based on a key lookup. The default value is
	// returned if the key is not found.
	inline V & operator[] (const K & k) const {
		return value[haskey(k)];
	}
	// Allow integer keys.
	inline V & operator[] (const int p) const {
		return (inbounds(p+1) ? value[p+1] : value[0]);
	}

protected:
	V * value;                          // The data.
	struct _V {                         // Placement new wrapper
		V _v;
		explicit _V() : _v() {}
		explicit _V(const V & v) : _v(v) {}
		void *operator new(size_t, void * p) {
			return p;
		} 
	};

	// Hide the allocation function.
	bool allocate(const int s) {
		V * temp = static_cast<V *>(realloc(value,bufsize(s)*sizeof(V)));
		if (temp == 0) {
			event_msg(EVENT_ERROR,"Out of memory in kfset::allocate()!");
			return false;
		}
		if (value == 0)
			new(&temp[0]) _V();
		value = temp;
		return true;
	}
	void deallocate() {
		if (value) {
			for (int i = 0; i <= size; i++)
				value[i].~V();
			free(value);
			value = 0;
		}
		size = 0;
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
		return add(&v.value[1],v.size);
	}
	// Add an array of values... There is no point in
	// a position based addition. 
	bool add(const V * v, const int s = 1) {
		if (s <= 0 || v == 0)
			return false;
		if (!allocate(this->size+s))
			return false;
		for (int i = 0, p; i < s; i++) {
			if ((p = haskey(v[i])) != 0)
				this->value[p] = v[i];
			else
				this->value[++(this->size)] = v[i];
		}
		return allocate(this->size);
	}
	// Remove based on key.
	bool remove(const K & k) {
		int p = haskey(k);
		if (p == 0)
			return false;
		this->value[p].~V();
		if (p < this->size)
			memmove(&this->value[p],&this->value[p+1],
				(this->size-p)*sizeof(V));
		return allocate(--(this->size));
	}
	// Replace key/value with another.
	inline bool replace(const K & ko, const V & v) {
		int p = haskey(ko);
		if (p == 0)
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
 * class aoset - Keyed ordered set
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
	}
	inline explicit koset(const kfset<K,V> & v)
	  : ksset<K,V>(v) {
	}
	inline ~koset() {
	}

	inline void sort() {
		qsort(1,this->size);
	}
protected:
	void qsort(const int l, const int r) {
		if (r > l) {
			int i, j, k, p = (l+r)/2;
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
				isort(l, p-1);
				isort(p+1, r);
			} else {
				qsort(l, p-1);
				qsort(p+1, r);
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
	  : set(s,b), key(0), value(0) {
		int i, p;
		if (!allocate(size)) {
			size = 0;
			return;
		}
		for (i = 0, size = 0; k && i < s; i++) {
			if ((p = haskey(k[i])) != 0) {
				value[p] = (v ? v[i] : value[0]);
			} else {
				key[++size] = k[i];
				value[size] = (v ? v[i] : value[0]);
			}
		}
		allocate(size);
	}
	// Copy one...
	inline explicit afset(const afset<K,V> & v)
	  : set(v), key(0), value(0) {
		if (!allocate(size)) {
			size = 0;
			return;
		}
		for (int i = 0; i <= size; i++) {
			key[i] = v.key[i];
			value[i] = v.value[i];
		}
	}
	// Kill one...
	inline ~afset() {
		deallocate();
	}

	// Return an index for a key.
	inline int haskey(const K & k) const {
		for (int i = 0; i < size; i++) {
			if (key[i+1] == k)
				return i+1;
		}
		return 0;
	}
	// Assignement operator.
	inline afset<K,V> & operator= (const afset<K,V> & v) {
		deallocate();
		if (!allocate(v.size))
			return *this;
		size = v.size;
		for (int i = 0; i <= size; i++) {
			key[i] = v.key[i];
			value[i] = v.value[i];
		}
		return *this;
	}
	// Set our defaults.
	inline void setdefault(const K & k, const V & v) {
		key[0] = k;
		value[0] = v;
	}
	// Behave like an indexed array...
	inline V & operator[] (const K & k) const {
		return value[haskey(k)];
	}
	// Linear access.
	inline K & getkey(const int p) const {
		return (inbounds(p+1) ? key[p+1] : key[0]);
	}
	// More linear access.
	inline V & getvalue(const int p) const {
		return (inbounds(p+1) ? value[p+1] : value[0]);
	}
protected:
	K * key;                            // The keys.
	V * value;                          // Take a guess...
	struct _K {                         // Placement new wrapper
		K _k;
		explicit _K() : _k() {}
		explicit _K(const K & k) : _k(k) {}
		void *operator new(size_t, void * p) {
			return p;
		} 
	};
	struct _V {                         // Placement new wrapper
		V _v;
		explicit _V() : _v() {}
		explicit _V(const V & v) : _v(v) {}
		void *operator new(size_t, void * p) {
			return p;
		} 
	};

	// Make some space...
	bool allocate(const int s) {
		K * ktemp = static_cast<K *>(realloc(key,bufsize(s)*sizeof(K)));
		V * vtemp = static_cast<V *>(realloc(value,bufsize(s)*sizeof(V)));
		if (ktemp == 0 || vtemp == 0) {
			event_msg(EVENT_ERROR,"Out of memory in kfset::allocate()!");
			return false;
		}
		if (key == 0)
			new(&ktemp[0]) _K();
		if (value == 0)
			new(&vtemp[0]) _V();
		key = ktemp; value = vtemp;
		return true;

	}
	void deallocate() {
		if (key) {
			for (int i = 0; i <= size; i++)
				key[i].~K();
			free(key);
			key = 0;
		}
		if (value) {
			for (int i = 0; i <= size; i++)
				value[i].~V();
			free(value);
			value = 0;
		}
		size = 0;
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
		return add(&v.key[1],&v.value[1],v.size);
	}
	// Add a whole bunch...
	bool add(const K * k, const V * v = 0, const int s = 1) {
		if (s <= 0 || k == 0)
			return false;
		if (!allocate(this->size+s))
			return false;
		for (int i = 0, p; i < s; i++) {
			if ((p = haskey(k[i])) != 0) {
				this->value[p] = (v ? v[i] : this->value[0]);
			} else {
				this->key[++(this->size)] = k[i];
				this->value[this->size] = (v ? v[i] : this->value[0]);
			}
		}
		return allocate(this->size);
	}
	// Now start removing them...
	bool remove(const K & k) {
		int p = haskey(k);
		if (p == 0)
			return false;
		this->key[p].~K();
		this->value[p].~V();
		if (p < this->size) {
			memmove(&this->key[p],&this->key[p+1],
					(this->size-p)*sizeof(K));
			memmove(&this->value[p],&this->value[p+1],
					(this->size-p)*sizeof(V));
		}
		return allocate(--(this->size));
	}
	// Or replacing them...
	inline bool replace(const K & ko, const K & k, const V & v) {
		int p = haskey(ko);
		if (p == 0)
			return false;
		this->key[p] = k;
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
 * Ordered by key set. Does anyone want to order by value?
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
	}
	inline explicit aoset(const afset<K,V> & v)
	  : asset<K,V>(v) {
	}
	inline ~aoset() {
	}

	inline void sort() {
		qsort(1,this->size);
	}
protected:
	void qsort(const int l, const int r) {
		if (r > l) {
			int i, j, k, p = (l+r)/2;
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
				isort(l, p-1);
				isort(p+1, r);
			} else {
				qsort(l, p-1);
				qsort(p+1, r);
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
	}
	inline explicit avoset(const afset<K,V> & v)
	  : asset<K,V>(v) {
	}
	inline ~avoset() {
	}

	inline void sort() {
		qsort(1,this->size);
	}
protected:
	void qsort(const int l, const int r) {
		if (r > l) {
			int i, j, k, p = (l+r)/2;
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
				isort(l, p-1);
				isort(p+1, r);
			} else {
				qsort(l, p-1);
				qsort(p+1, r);
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
