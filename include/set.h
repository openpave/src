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
		be used as a default. However, the operator [] functions are
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
		The code is really niave, and will do a lot of moving of objects
		(probably using their assigment operators which might new/delete,
		etc.).  It should be redone to use realloc/free, placement new,
		and memcpy.  It should also learn to throw expections...

	History:
		1994       - Created by Jeremy Lea <reg@openpave.org>
		2002/01/23 - Modifications for use on a modern C++ compiler

**************************************************************************/

#ifndef __SET_H
#define __SET_H

#include "event.h"
#include "mathplus.h"

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
	int size;							// The size of the set...
	int block;							// The minimum block size.

	// Simple constructor...
	inline set(const int s = 0, const int b = DFLT_BLK) {
		size = (s > 0 ? s : 0);
		block = (b > 1 ? b : 1);
	}
	inline ~set() {
	}

	// Calculate the buffer size.
	inline int buffer(const int s) const {
		if (s > 0)
			return (block == 0 ? s+1 : block*(s/block+(s%block?1:0))+1);
		else
			return (block == 0 ? 1 : block);
	}
	// Check if we need to reallocate the buffer.
	inline bool reallocate(const int s) const {
		return (buffer(size) != buffer(s) ? true : false);
	}
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
	fset(const int s, const int b = DFLT_BLK, const V * v = 0)
			: set(s,b), value(0) {
		if (!allocate(size))
			return;
		for (int i = 0; v && i < size; i++)
			value[i+1] = v[i];
	}
	// Copy constructor.
	fset(const fset<V> & v) : set(v.size,v.block), value(0) {
		if (!allocate(size))
			return;
		for (int i = 0; i <= size; i++)
			value[i] = v.value[i];
	}
	// Wow, a destructor...
	~fset() {
		if (value != 0)
			delete [] value;
		value = 0;
	}

	// Assignment operator.
	fset<V> & operator = (const fset<V> & v) {
		if (reallocate(v.size)) {
			delete [] value;
			value = 0;
			if (!allocate(v.size)) {
				size = 0;
				return *this;
			}
		}
		size = v.size;
		block = v.block;
		for (int i = 0; i <= size; i++)
			value[i] = v.value[i];
		return *this;
	}
	// Set the default element if constructor didn't do a good job.
	void setdefault(const V & d) {
		value[0] = d;
	}
	// Return an element, or the default. One indexed.
	V & data(const int p) const {
		return (inbounds(p) ? value[p] : value[0]);
	}
	// Behave like an array. Zero indexed.
	V & operator [] (const int p) const {
		return value[p+1];
	}

protected:
	V * value;							// The buffer.

	// Since this is a fixed size set, hide the allocator...
	bool allocate(const int s) {
		V * tmp = value;
		if (s > 10*block)
			block *= 10;
		value = new V[buffer(s)];
		if (value == 0) {
			value = tmp;
			event_msg(EVENT_ERROR,"Out of memory in fset::allocate()!");
			return false;
		}
		return true;
	}
	// Also hide the null constructor.
	inline fset() : set(), value(0) {
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
	sset() : fset<V>() {
	}
	// Simple constructor.
	sset(const int s, const int b = DFLT_BLK, const V * v = 0)
		: fset<V>(s,b,v) {
	}
	// Copy constructor.
	sset(const fset<V> & v) : fset<V>(v) {
	}
	// We let someone else clean up...
	~sset() {
	}

	// Add one element, at the end.
	inline bool add(const V & v) {
		return add(this->size+1,&v,1);
	}
	// Add a whole set, at the end.
	inline bool add(const fset<V> & v) {
		return add(this->size+1,&v.value[1],v.size);
	}
	// Add an array, at the end.
	inline bool add(const V * v, const int s = 1) {
		return add(this->size+1,v,s);
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
	bool add(const int p, const V * v, const int s = 1) {
		int i, j;
		if (p < 1 || s < 1)
			return false;
		if (p > this->size) {
			if (this->reallocate(p+s-1)) {
				V * temp = this->value;
				if (!this->allocate(p+s-1))
					return false;
				for (i = 0; i <= this->size; i++)
					this->value[i] = temp[i];
				delete [] temp;
			}
			for (i = this->size+1; i < p; i++)
				this->value[i] = this->value[0];
			this->size = p+s-1;
			for (i = 0; i < s; i++)
				this->value[p+i] = v[i];
		} else {
			if (this->reallocate(this->size+s)) {
				V * temp = this->value;
				if (!this->allocate(this->size+s))
					return false;
				this->size += s;
				for (i = 0, j = i; i <= this->size; i++)
					this->value[i] = (i >= p && i < p+s ? v[i-p] : temp[j++]);
				delete [] temp;
			} else {
				this->size += s;
				for (i = this->size, j = i - 1; i >= p; i--)
					this->value[i] = (i < p+s ? v[i-p] : this->value[j--]);
			}
		}
		return true;
	}
	// Remove the last element.
	inline bool remove() {
		return remove(this->size,1);
	}
	// Remove s elements, starting at position p.
	bool remove(const int p, const int s = 1) {
		if (this->inbounds(p) && s > 0) {
			int as = (p+s-1 > this->size ? this->size-p+1 : s);
			if (reallocate(this->size-as)) {
				V * temp = this->value;
				if (!allocate(this->size-as))
					return false;
				this->size -= as;
				for (int i = 0; i <= this->size; i++)
					this->value[i] = (i < p ? temp[i] : temp[i+as]);
				delete [] temp;
			} else {
				this->size -= as;
				for (int i = p; i <= this->size; i++)
					this->value[i] = this->value[i+as];
			}
		} else {
			return false;
		}
		return true;
	}
	bool empty() {
		if (this->reallocate(0)) {
			delete [] this->value;
			this->value = 0;
			if (!this->allocate(0))
				return false;
		}
		this->size = 0;
		return true;
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
	oset() : sset<V>() {
	}
	// Simple constructor.
	oset(const int s, const int b = DFLT_BLK, const V * v = 0)
		: sset<V>(s,b,v) {
		if (v)
			sort();
	}
	// Copy constructor.
	oset(const fset<V> & v) : sset<V>(v) {
		sort();
	}
	// Destructor.
	~oset() {
	}

	// Guess what?
	void sort() {
		qsort(1,this->size);
	}

protected:
	// A highly optimised quick sort. Don't touch...
	void qsort(const int l, const int r) {
		if (r > l) {
			int i, j, k, p = (l+r)/2;
			if (this->value[p] <= this->value[l])
				swap(this->value[l],this->value[p]);
			if (this->value[r] <= this->value[l])
				swap(this->value[l],this->value[r]);
			if (this->value[r] <= this->value[p])
				swap(this->value[p],this->value[r]);
			for (i = l, j = r, k = 0; ;
								p = (p==i?j++:(p==j?i--:p)), k++) {
				while (++i < p && this->value[i] <= this->value[p]);
				while (p < --j && this->value[p] <= this->value[j]);
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
		for (int i = l+1; l+1 < r && i <= r; i++) {
			for (int j = i; j > l && this->value[j] <= this->value[j-1]; j--)
				swap(this->value[j],this->value[j-1]);
		}
	}
};

/*
 * class cset - An compact ordered sizeable set.
 *
 * Build on the ordered set to make a set which take unique
 * values only. This is also lazy.
 *
 * Custom classes should define operator !=.
 */
template <class V>
class cset : public oset<V> {
public:
	// Getting the hang of this yet?
	cset() : oset<V>() {
	}
	cset(const int s, const int b = DFLT_BLK, const V * v = 0)
		: oset<V>(s,b,v) {
		if (v)
			compact();
	}
	cset(const fset<V> & v) : oset<V>(v) {
		compact();
	}
	~cset() {
	}

	// Sort then compact the set.
	void sort() {
		this->qsort(1,this->size);
		compact();
	}

protected:
	void compact() {
		int i, j;
		if (this->size < 2)
			return;
		for (i = 2, j = 1; this->size > 1 && i <= this->size; i++) {
			if (this->value[j] != this->value[i] && ++j != i)
				this->value[j] = this->value[i];
		}
		this->size = j;
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
	kfset(const int s, const int b = DFLT_BLK, const V * v = 0)
		: set(s,b) {
		int i, p;
		if (!allocate(this->size))
			return;
		for (i = 0, this->size = 0; v && i < s; i++) {
			if (p = haskey(v[i])) {
				this->value[p] = v[i];
			} else {
				this->value[++(this->size)] = v[i];
			}
		}
	}
	// Copy constuctor.
	kfset(const kfset<K,V> & v) : set(v.size,v.block), value(0) {
		if (!allocate(this->size))
			return;
		for (int i = 0; i <= this->size; i++)
			this->value[i] = v.value[i];
	}
	// Clean up.
	~kfset() {
		if (this->value)
			delete [] this->value;
		this->value = 0;
	}

	// Do a key lookup, and return zero if the key is not found.
	int haskey(const K & k) const {
		for (int i = 1; this->size >= 1 && i <= this->size; i++) {
			if (static_cast<K &>(this->value[i]) == k)
				return i;
		}
		return 0;
	}
	// Assignment operator.
	kfset<K,V> & operator = (const kfset<K,V> & v) {
		if (reallocate(v.size)) {
			delete [] this->value;
			this->value = 0;
			if (!allocate(v.size)) {
				this->size = 0;
				return *this;
			}
		}
		this->size = v.size;
		block = v.block;
		for (int i = 0; i <= this->size; i++)
			this->value[i] = v.value[i];
		return *this;
	}
	// Set our default elelment.
	void setdefault(const V & v) {
		this->value[0] = v;
	}
	// Return data based on a key lookup. The default value is
	// returned if the key is not found.
	V & data(const K & k) const {
		return this->value[haskey(k)];
	}
	// Look and feel like an array.
	V & operator [] (const K & k) const {
		return data(k);
	}
	// Allow integer keys.
	V & operator [] (const int p) const {
		return this->value[p+1];
	}

protected:
	V * value;							// The data.

	// Hide the allocation function.
	bool allocate(const int s) {
		V * tmp = this->value;
		if (s > 10*block)
			block *= 10;
		this->value = new V[buffer(s)];
		if (this->value == 0) {
			this->value = tmp;
			event_msg(EVENT_ERROR,"Out of memory in kfset::allocate!");
			return false;
		}
		return true;
	}
	// And the null constructor.
	kfset() : set(), value(0) {
		allocate(this->size);
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
	ksset() : kfset<K,V>() {
	}
	ksset(const int s, const int b = DFLT_BLK, const V * v = 0)
		: kfset<K,V>(s,b,v) {
	}
	ksset(const kfset<K,V> & v) : kfset<K,V>(v) {
	}
	~ksset() {
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
		if (reallocate(this->size+s)) {
			V * tempval = this->value;
			if (!allocate(this->size+s))
				return false;
			for (int i = 0; i <= this->size; i++)
				this->value[i] = tempval[i];
			delete [] tempval;
		}
		for (int i = 0; i < s; i++) {
			if (int p = haskey(v[i])) {
				this->value[p] = v[i];
			} else {
				this->value[++(this->size)] = v[i];
			}
		}
		return true;
	}
	// Remove based on key.
	bool remove(const K & k) {
		if (int p = haskey(k)) {
			if (reallocate(this->size--)) {
				V * tempval = this->value;
				if (!allocate(this->size))
					return false;
				for (int i = 0, j = 0; j <= this->size; i++) {
					if (i != p)
						this->value[j++] = tempval[i];
				}
				delete [] tempval;
			} else {
				for (int i = p; i <= this->size; i++)
					this->value[i] = this->value[i+1];
			}
		} else {
			return false;
		}
		return true;
	}
	// Replace key/value with another.
	bool replace(const K & ko, const V & v) {
		if (int p = haskey(ko)) {
			this->value[p] = v;
			return true;
		}
		return false;
	}
	bool empty() {
		if (this->reallocate(0)) {
			delete [] this->value;
			this->value = 0;
			this->size = 0;
			if (!this->allocate(0))
				return false;
		}
		this->size = 0;
		return true;
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
	koset() : ksset<K,V>() {
	}
	koset(const int s, const int b = DFLT_BLK, const V * v = 0)
		: ksset<K,V>(s,b,v) {
	}
	koset(const kfset<K,V> & v) : ksset<K,V>(v) {
	}
	~koset() {
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
				while (++i < p && this->value[i] <= this->value[p]);
				while (p < --j && this->value[p] <= this->value[j]);
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
		for (int i = l+1; i <= r; i++) {
			for (int j = i; j > l && this->value[j-1] > this->value[j]; j--) {
				swap(this->value[j],this->value[j-1]);
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
	// Make empty one...
	afset(const int s, const int b = DFLT_BLK) : set(s,b) {
		allocate(size);
	}
	// Make one...
	afset(const int s, const int b, const K * k, const V * v = 0)
		: set(s,b) {
		int i, p;
		if (!allocate(size))
			return;
		for (i = 0, size = 0; k && i < s; i++) {
			if (p = haskey(k[i])) {
				if (v != 0)
					value[p] = v[i];
			} else {
				key[++size] = k[i];
				if (v != 0)
					value[size] = v[i];
			}
		}
	}
	// Copy one...
	afset(const afset<K,V> & v) : set(v.size,v.block) {
		if (!allocate(size))
			return;
		for (int i = 0; i <= size; i++) {
			key[i] = v.key[i];
			value[i] = v.value[i];
		}
	}
	// Kill one...
	~afset() {
		if (key)
			delete [] key;
		if (value)
			delete [] value;
		key = 0;
		value = 0;
	}

	// Return an index for a key.
	int haskey(const K & k) const {
		for (int i = 0; i < size; i++) {
			if (key[i+1] == k)
				return i+1;
		}
		return 0;
	}
	// Assignement operator.
	afset<K,V> & operator = (const afset<K,V> & v) {
		if (reallocate(v.size)) {
			delete [] key;
			delete [] value;
			key = value = 0;
			if (!allocate(v.size)) {
				size = 0;
				return *this;
			}
		}
		size = v.size;
		block = v.block;
		for (int i = 0; i <= size; i++) {
			key[i] = v.key[i];
			value[i] = v.value[i];
		}
		return *this;
	}
	// Set our defaults.
	void setdefault(const K & k, const V & v) {
		key[0] = k;
		value[0] = v;
	}
	// Return data...
	inline V & data(const K & k) const {
		return value[haskey(k)];
	}
	// Behave like an indexed array...
	inline V & operator [] (const K & k) const {
		return data(k);
	}
	// Linear access.
	inline K & getkey(const int p) const {
		return (inbounds(p) ? key[p] : key[0]);
	}
	// More linear access.
	inline V & getvalue(const int p) const {
		return (inbounds(p) ? value[p] : value[0]);
	}
protected:
	K * key;							// The keys.
	V * value;							// Take a guess...

	// Make some space...
	bool allocate(const int s) {
		K * ktmp = key;
		V * vtmp = value;
		if (s > 10*block)
			block *= 10;
		key = new K[buffer(s)];
		value = new V[buffer(s)];
		if (value == 0 || key == 0) {
			key = ktmp;
			value = vtmp;
			event_msg(EVENT_ERROR,"Out of memory in afset::allocate()!");
			return false;
		}
		return true;

	}
	// Don't allow yobos to make empty sets...
	afset() : set() {
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
	asset() : afset<K,V>() {
	}
	// Simple constructor.
	asset(const int s, const int b = DFLT_BLK, const K * k = 0,
		const V * v = 0) : afset<K,V>(s,b,k,v) {
	}
	// Copy constructor.
	asset(const afset<K,V> & v) : afset<K,V>(v) {
	}
	// Destructor.
	~asset() {
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
		if (s <= 0)
			return false;
		if (this->reallocate(this->size+s)) {
			K * tempkey = this->key;
			V * tempval = this->value;
			if (!allocate(this->size+s))
				return false;
			for (int i = 0; i <= this->size; i++) {
				this->key[i] = tempkey[i];
				this->value[i] = tempval[i];
			}
			delete [] tempkey;
			delete [] tempval;
		}
		for (int i = 0; i < s; i++) {
			if (int p = haskey(k[i])) {
				if (v != 0)
					this->value[p] = v[i];
			} else {
				this->size++;
				this->key[this->size] = k[i];
				if (v != 0)
					this->value[this->size] = v[i];
			}
		}
		return true;
	}
	// Now start removing them...
	bool remove(const K & k) {
		if (int p = haskey(k)) {
			if (this->reallocate(this->size--)) {
				K * tempkey = this->key;
				V * tempval = this->value;
				if (!allocate(this->size))
					return false;
				for (int i = 0,j = i; j <= this->size; i++) {
					if (i != p) {
						this->key[j] = tempkey[i];
						this->value[j++] = tempval[i];
					}
				}
				delete [] tempkey;
				delete [] tempval;
			} else {
				for (int i = p; i <= this->size; i++) {
					this->key[i] = this->key[i+1];
					this->value[i] = this->value[i+1];
				}
			}
		} else {
			return false;
		}
		return true;
	}
	// Or replacing them...
	bool replace(const K & ko, const K & k, const V & v) {
		if (int p = haskey(ko)) {
			this->key[p] = k;
			this->value[p] = v;
			return true;
		} else {
			return false;
		}
	}
	bool empty() {
		if (this->reallocate(0)) {
			delete [] this->key;
			delete [] this->value;
			this->key = this->value = 0;
			if (!this->allocate(0)) {
				this->size = 0;
				return false;
			}
		}
		this->size = 0;
		return true;
	}
};

/*
 * class aoset - Associative ordered set
 *
 * Ordered by vale set. Does anyone want to order by key?
 */
template <class K, class V>
class aoset : public asset<K,V> {
public:
	aoset() : asset<K,V>() {
	}
	aoset(const int s, const int b = DFLT_BLK, const K * k = 0,
		const V * v = 0) : asset<K,V>(s,b,k,v) {
	}
	aoset(const afset<K,V> & v) : asset<K,V>(v) {
	}
	~aoset() {
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
				while (++i < p && this->value[i] <= this->value[p]);
				while (p < --j && this->value[p] <= this->value[j]);
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
		for (int i = l+1; l+1 < r && i <= r; i++) {
			for (int j = i; j > l && this->value[j-1] > this->value[j]; j--) {
				swap(this->key[j],this->key[j-1]);
				swap(this->value[j],this->value[j-1]);
			}
		}
	}
};

#endif // SET_H
