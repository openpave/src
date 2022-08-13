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

	Portions Copyright (C) 2016-2022 OpenPave.org.

	Contributor(s): Jeremy Lea <reg@openpave.org>.

	Purpose:
		This header implements classes for building multi-level 'axes',
		which resemble typical row headers one might see on a table in a
		document where the left-most column contains row headers in 'merged'
		cells.  The prior columns must also be an axis.  As elements are
		added and removed from the prior columns, the current axis listens
		for these and adjusts accordingly.  All axes are sorted and must
		be unique.

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

		The downside of variadic templates is that it means that the axes
		are in reverse order to their natural order in all the calls...

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

#if !defined(DFLT_BLK)
#define DFLT_BLK    64
#endif

#include <functional>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include "listen.h"
#include "tree.h"
#include "hascompare.h"

namespace OP {

enum class axis_message {
	add, remove, empty, deleting
};

// Forward declare class axis as a variadic template
template<typename...Ks> class axis {};

/*
 * Class axis - a multi-level set (for use as an "axis")
 *
 * This defines an "axis" which can have an number of prior axes.  This
 * level has key type K, and the parameter pack Ks has the prior levels.
 */
template<typename K, typename...Ks>
class axis<K,Ks...>
  : public dispatcher<message<axis_message,unsigned>>, public listener
{
public:
	typedef unsigned index_t;
	typedef std::tuple<const K, const Ks...> key_t;
	typedef std::tuple<const K &, const Ks &...> ref_t;

	axis(axis<Ks...> & pr)
	  : prior(pr) {
		listen(pr,message<axis_message,unsigned>(
			[this](axis_message e, unsigned p){
				switch (e) {
				case axis_message::add:
					for (unsigned i = 0; i < me.length(); i++) {
						axis_key & a = me[i];
						if (a.pi >= p)
							a.pi++;
					}
					break;
				case axis_message::remove:
					for (unsigned i = 0; i < me.length(); i++) {
						axis_key & a = me.getatorder(i);
						if (a.pi == p) {
							this->dispatch(axis_message::remove,i);
							me.remove(a);
							i--;
						} else if (a.pi > p)
							a.pi--;
					}
					break;
				case axis_message::empty:
					this->dispatch(axis_message::empty,p);
					me.empty();
					break;
				case axis_message::deleting:
					throw std::runtime_error("Out of order destructors!");
				}
		}));
	}
	~axis() {
		this->dispatch(axis_message::deleting,UINT_MAX);
	}
	void add(const K & k, const Ks &...ks) {
		if (!prior.haskey(ks...))
			throw std::runtime_error("attempting to insert key into axis without key in prior!");
		if (haskey(k,ks...))
			throw std::runtime_error("attempting to insert duplicate key into axis!");
		axis_key a(prior,k,ks...);
		me.add(a);
		this->dispatch(axis_message::add,me.getorderof(a));
	}
	void remove(const K & k, const Ks &...ks) {
		if (!haskey(k,ks...))
			throw std::runtime_error("removal key not found in axis!");
		axis_key a(prior,k,ks...);
		this->dispatch(axis_message::remove,me.getorderof(a));
		me.remove(a);
	}
	unsigned length() const noexcept {
		return me.length();
	}
	// Do a key lookup, and return UINT_MAX if the key is not found.
	bool haskey(const K & k, const Ks &...ks) const {
		if (!prior.haskey(ks...))
			return false;
		return me.haskey(axis_key(prior,k,ks...));
	}
	// Get the full key as a tuple
	template<unsigned N = sizeof...(Ks)>
	typename std::enable_if<N!=1,ref_t>::type
	operator [] (unsigned p) const {
		if (!me.inbounds(p))
			throw std::out_of_range("ordered index out of bounds!");
		const axis_key & a = me.getatorder(p);
		return std::tuple_cat(std::tuple<const K &>(a.key),prior[a.pi]);
	}
	// Special version for the second last since it does not get a tuple
	// from the first axis.
	template<unsigned N = sizeof...(Ks)>
	typename std::enable_if<N==1,ref_t>::type
	operator [] (unsigned p) const {
		if (!me.inbounds(p))
			throw std::out_of_range("ordered index out of bounds!");
		const axis_key & a = me.getatorder(p);
		return ref_t(a.key,prior[a.pi]);
	}
	// Get the ordered position of an element in the sort.
	unsigned getorderof(const K & k, const Ks &...ks) const {
		if (!haskey(k,ks...))
			throw std::runtime_error("key not found in axis!");
		return me.getorderof(axis_key(prior,k,ks...));
	}
	unsigned getorderof(const ref_t & t) const {
		return unpack_getorderof(t,typename idx<sizeof...(Ks)+1>::type());
	}
	// Remove all elements
	void empty() {
		this->dispatch(axis_message::empty,UINT_MAX);
		me.empty();
	}

private:
	// Tuple unpacker
	template<int...> struct seq {};
	template<int N, int...S> struct idx : idx<N-1,N-1,S...> {};
	template<int...S> struct idx<0,S...>{ typedef seq<S...> type; };
	template<int...S>
	unsigned unpack_getorderof(const ref_t & t, seq<S...>) const {
		return getorderof(std::get<S>(t)...);
	}
	struct axis_key {
		unsigned pi;
		K key;
		axis_key(const axis<Ks...> & pr, const K & k, const Ks &...ks)
		  : key(k) {
			pi = pr.getorderof(ks...);
			if (pi == UINT_MAX)
				throw std::runtime_error("key not found in prior axis!");
		}
		// use the compare function if it has one
		template<typename T = K>
		typename std::enable_if<has_compare<T>::value,int>::type
		compare(const axis_key & k) const {
			if (pi == k.pi)
				return key.compare(k.key);
			return (pi < k.pi ? -1 : 1);
		}
		// else resort to operators
		template<typename T = K>
		typename std::enable_if<!has_compare<T>::value,int>::type
		compare(const axis_key & k) const noexcept {
			if (pi == k.pi)
				return (key < k.key ? -1 : key == k.key ? 0 : 1);
			return (pi < k.pi ? -1 : 1);
		}
	};
	ktree_avl<axis_key> me;    // internal index
	const axis<Ks...> & prior; // the next axis up/to the left
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
	typedef unsigned index_t;
	typedef const K key_t;
	typedef const K & ref_t;

	axis() {
	}
	~axis() {
		this->dispatch(axis_message::deleting,UINT_MAX);
	}
	void add(const K & k) {
		if (haskey(k))
			throw std::runtime_error("attempting to insert duplicate key into axis!");
		axis_key a(k);
		me.add(a);
		this->dispatch(axis_message::add,me.getorderof(a));
	}
	void remove(const K & k) {
		if (!haskey(k))
			throw std::runtime_error("removal key not found in axis!");
		axis_key a(k);
		this->dispatch(axis_message::remove,me.getorderof(a));
		me.remove(a);
	}
	unsigned length() const noexcept {
		return me.length();
	}
	// Do a key lookup, and return UINT_MAX if the key is not found.
	bool haskey(const K & k) const {
		return me.haskey(axis_key(k));
	}
	// Get the key at some position in the sort order
	// This does not return a tuple to avoid making life complex.
	const K & operator [] (unsigned p) const {
		if (!me.inbounds(p))
			throw std::out_of_range("ordered index out of bounds!");
		return me.getatorder(p).key;
	}
	// Get the ordered position of an element in the sort.
	unsigned getorderof(const K & k) const {
		if (!haskey(k))
			throw std::runtime_error("key not found in axis!");
		return me.getorderof(axis_key(k));
	}
	// Remove all elements
	void empty() {
		this->dispatch(axis_message::empty,UINT_MAX);
		me.empty();
	}

private:
	struct axis_key {
		K key;
		axis_key(const K & k) noexcept : key(k) {}
		// use the compare function if it has one
		template<typename T = K>
		typename std::enable_if<has_compare<T>::value,int>::type
		compare(const axis_key & k) const noexcept {
			return key.compare(k.key);
		}
		// else resort to operators
		template<typename T = K>
		typename std::enable_if<!has_compare<T>::value,int>::type
		compare(const axis_key & k) const noexcept {
			return (key < k.key ? -1 : key == k.key ? 0 : 1);
		}
	};
	ktree_avl<axis_key> me;    // internal index
};

} // namespace OP

#endif // AXIS_H
