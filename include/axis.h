/**************************************************************************

	AXIS.H - A type for nested multi-value sets

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
		This header implements classes for building multi-level 'axes',
		which resemble typical row headers one might see on a table in a
		document where the left-most column contains row headers in 'merged'
		cells.  The prior columns must also be an axis.  As elements are
		added and removed from the prior columns, the current axis listens
		for these and adjusts accordingly.
		
		+---+---+------
		| A | a | 1
		+---+---+------
		|   | a | 1
		| B +---+------
		|   | b | 1
		+---+---+------
		| C | b | 1
		+---+---+------
		|   |   | 1
		|   | a +------
		| D |   | 2
		|   +---+------
		|   | b | 1
		+---+---+------

		The class provides methods to iterate over the axis, and to get a
		tuple with the keys for a row.  The class automatically handles
		adds and removes to the prior axis, but this is quite expensive,
		so should not be used much in practice - it is better for the prior
		axis to be essentially constant.

	Design:
		These classes make use of many new C++11 features, including
		variadic templates and SFINAE.  This makes them possible, but also
		complex.
		
		The basic class is axis.  This contains a set which holds the keys
		for this axis, and optionally a reference to the prior axis.  The
		class listens for events on the prior axis, and sends events to any
		listening axes or tables, so they can auto adjust.  The axis does
		not auto expand for an added key in a prior axis, but does delete.

		The axis uses a tree internally to store the key and an array of
		index positions for each prior.  These are used to make sure we
		have a prior cell, and for sorting the current key.
		
	History:
		2016/03/05 - Created a basic implementation.

**************************************************************************/

#ifndef __AXIS_H
#define __AXIS_H

#include <functional>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include "listen.h"
#include "tree.h"
#include "hascompare.h"

enum class axis_message {
	add, remove, deleting
};

// Forward declare class axis as a variadic template
template<typename ...Ks> class axis {};

/*
 * Class axis - a multi-level set (for use as an "axis")
 *
 * This defines an "axis" which can have an number of prior axes.  This
 * level has key type K, and the parameter pack Ks has the prior levels.
 */
template<typename K, typename ...Ks>
class axis<K,Ks...>
  : public dispatcher<message<axis_message,unsigned>>, public listener
{
public:
	axis(axis<Ks...> & pr)
	  : prior(pr) {
		listen(pr,message<axis_message,unsigned>(
			[this](axis_message e, unsigned p){
				switch (e) {
				case axis_message::add:
					for (unsigned i = 0; i < internal.length(); i++) {
						axis_key & a = internal[i];
						if (a.pi >= p)
							a.pi++;
					}
					break;
				case axis_message::remove:
					// we remove in order
					for (unsigned i = 0; i < internal.length(); i++) {
						axis_key & a = internal.getatorder(i);
						if (a.pi == p) {
							dispatch(axis_message::remove,
								internal.getorderof(a));
							internal.remove(a);
							i--;
						} else if (a.pi > p)
							a.pi--;
					}
					break;
				case axis_message::deleting:
					throw std::runtime_error("Out of order destructors!");
				default:
					break;
				}
		}));
	}
	~axis() {
		dispatch(axis_message::deleting,UINT_MAX);
	}
	void add(const K & k, const Ks & ...ks) {
		axis_key a(prior,k,ks...);
		unsigned p = internal.getorderof(a);
		if (p != UINT_MAX)
			return; // XXX throw?
		internal.add(a);
		p = internal.getorderof(a);
		dispatch(axis_message::add,p);
	}
	void remove(const K & k, const Ks & ...ks) {
		axis_key a(prior,k,ks...);
		unsigned p = internal.getorderof(a);
		if (p == UINT_MAX)
			return; // XXX throw?
		internal.remove(a);
		dispatch(axis_message::remove,p);
	}
	unsigned length() const {
		return internal.length();
	}
	bool inbounds(const unsigned p) const {
		return internal.inbounds(p);
	}
	// Do a key lookup, and return UINT_MAX if the key is not found.
	bool haskey(const K & k, const Ks & ...ks) const {
		return internal.haskey(axis_key(prior,k,ks...));
	}
	std::tuple<const K &, const Ks &...>
	operator[] (const unsigned p) const {
		if (!inbounds(p))
			throw std::out_of_range("unordered index out of bounds!");
		const axis_key & a = internal.getatorder(p);
		return std::tuple_cat(std::tuple<const K &>(a.key),prior[a.pi]);
	}
	// Get the ordered position of an element in the sort.
	unsigned getorderof(const K & k, const Ks & ...ks) const {
		return internal.getorderof(axis_key(prior,k,ks...));
	}

private:
	struct axis_key {
		unsigned pi;
		K key;

		axis_key(const axis<Ks...> & pr, const K & k, const Ks &...ks)
		  : key(k) {
			pi = pr.getorderof(ks...);
		}
		template<typename T = K>
		typename std::enable_if<has_compare<T>::value,int>::type
		compare(const axis_key & k) const {
			if (pi < k.pi)
				return -1;
			if (pi > k.pi)
				return 1;
			return key.compare(k.key);
		}
		template<typename T = K>
		typename std::enable_if<!has_compare<T>::value,int>::type
		compare(const axis_key & k) const {
			if (pi < k.pi)
				return -1;
			if (pi > k.pi)
				return 1;
			if (key < k.key)
				return -1;
			if (key > k.key)
				return 1;
			return 0;
		}
	};
	ktree_avl<axis_key> internal;     // internal index
	const axis<Ks...> & prior;
};

/*
 * class axis<K> - specialized top level axis
 *
 * This specialization of axis does not have priors, so is actually the
 * first type that is created.  Since it does not have priors it is much
 * easier to implement - it is a simple wrapper around the set of keys.
 */
template<typename K>
class axis<K>
  : public dispatcher<message<axis_message,unsigned>>
{
public:
	axis() {
	}
	~axis() {
		dispatch(axis_message::deleting,UINT_MAX);
	}
	void add(const K & k) {
		axis_key a(k);
		unsigned p = internal.getorderof(a);
		if (p != UINT_MAX)
			return;
		internal.add(a);
		p = internal.getorderof(a);
		dispatch(axis_message::add,p);
	}
	void remove(const K & k) {
		axis_key a(k);
		unsigned p = internal.getorderof(a);
		if (p == UINT_MAX)
			return;
		internal.remove(a);
		dispatch(axis_message::remove,p);
	}
	unsigned length() const {
		return internal.length();
	}
	bool inbounds(const unsigned p) const {
		return internal.inbounds(p);
	}
	// Do a key lookup, and return UINT_MAX if the key is not found.
	bool haskey(const K & k) const {
		return internal.haskey(axis_key(k));
	}
	std::tuple<const K &> operator[] (const unsigned p) const {
		if (!inbounds(p))
			throw std::out_of_range("unordered index out of bounds!");
		const axis_key & k = internal.getatorder(p);
		return std::tuple<const K &>(k.key);
	}
	// Get the ordered position of an element in the sort.
	unsigned getorderof(const K & k) const {
		return internal.getorderof(axis_key(k));
	}

private:
	struct axis_key {
		K key;
		axis_key(const K & k) : key(k) {}
		template<typename T = K>
		typename std::enable_if<has_compare<T>::value,int>::type
		compare(const axis_key & k) const {
			return key.compare(k.key);
		}
		template<typename T = K>
		typename std::enable_if<!has_compare<T>::value,int>::type
		compare(const axis_key & k) const {
			if (key < k.key)
				return -1;
			if (key > k.key)
				return 1;
			return 0;
		}
	};
	ktree_avl<axis_key> internal;     // internal index
};

#endif // AXIS_H
