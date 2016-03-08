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
		cells.
		
		+----+----+------
		|    |    |
		+----+----+------
		|    |    |
		|    +----+------
		|    |    |
		+----+----+------
		|    |    |
		+----+----+------
		|    |    |
		|    |    +------
		|    |    |
		|    +----+------
		|    |    |
		+----+----+------

	Design:
		These classes make use of many new C++11 features, including
		variadic templates and SFINAE.  This makes them possible, but also
		complex.
		
		The basic class is axis.  This contains a set which holds the keys
		for this axis, and optionally a reference to the prior axis.
		
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

enum class axis_message {
	add, remove, deleting
};

/*
 * Class axis
 */
template<typename ...Ks> class axis {};

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
					break;
				case axis_message::remove:
					for (unsigned i = 0; i < internal.length(); i++) {
						axis_key & a = internal.getorder(i);
						if (a.priors[sizeof...(Ks)-1] == p) {
							dispatch(axis_message::remove,internal.getindex(a));
							internal.remove(a);
							i--;
						}
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
		unsigned p = internal.getindex(a);
		if (p != UINT_MAX)
			return;
		internal.add(a);
		p = internal.getindex(a);
		dispatch(axis_message::add,p);
	}
	void remove(const K & k, const Ks & ...ks) {
		axis_key a(prior,k,ks...);
		unsigned p = internal.getindex(a);
		if (p == UINT_MAX)
			return;
		internal.remove(a);
		dispatch(axis_message::remove,p);
	}
	inline unsigned length() const {
		return internal.length();
	}
	inline bool inbounds(const unsigned p) const {
		return internal.inbounds(p);
	}
	// Do a key lookup, and return UINT_MAX if the key is not found.
	inline bool haskey(const K & k, const Ks & ...ks) const {
		return internal.haskey(axis_key(prior,k,ks...));
	}
	inline K & operator[] (const unsigned p) const {
		if (!inbounds(p))
			throw std::out_of_range("unordered index out of bounds!");
		return internal[p];
	}
	// Do a key lookup, and return UINT_MAX if the key is not found.
	inline unsigned getposition(const K & k, const Ks & ...ks) const {
		return internal.getposition(axis_key(prior,k,ks...));
	}
	// Get the ordered position of an element in the sort.
	inline unsigned getindex(const K & k, const Ks & ...ks) const {
		return internal.getindex(axis_key(prior,k,ks...));
	}
	static const std::size_t levels = 1+sizeof...(Ks);
	void getlevels(unsigned (& p)[levels],
				   const K & k, const Ks & ...ks) const {
		axis_key a(prior,k,ks...);
		for (unsigned i = 0; i < sizeof...(Ks); i++)
			p[i] = a.priors[i];
		p[sizeof...(Ks)] = internal.getindex(a);
	}

private:
	struct axis_key {
		unsigned priors[sizeof...(Ks)];
		K key;

		axis_key(const axis<Ks...> & pr, const K & k, const Ks & ...ks)
		  : key(k) {
			pr.getlevels(priors,ks...);
		}
		int compare(const axis_key & k) const {
			for (unsigned i = 0; i < sizeof...(Ks); i++) {
				if (priors[i] < k.priors[i])
					return -1;
				if (priors[i] > k.priors[i])
					return 1;
			}
			return key.compare(k.key);
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
		unsigned p = internal.getindex(k);
		if (p != UINT_MAX)
			return;
		internal.add(k);
		p = internal.getindex(k);
		dispatch(axis_message::add,p);
	}
	void remove(const K & k) {
		unsigned p = internal.getindex(k);
		if (p == UINT_MAX)
			return;
		internal.remove(k);
		dispatch(axis_message::remove,p);
	}
	inline unsigned length() const {
		return internal.length();
	}
	inline bool inbounds(const unsigned p) const {
		return internal.inbounds(p);
	}
	// Do a key lookup, and return UINT_MAX if the key is not found.
	inline bool haskey(const K & k) const {
		return internal.haskey(k);
	}
	inline K & operator[] (const unsigned p) const {
		if (!inbounds(p))
			throw std::out_of_range("unordered index out of bounds!");
		return internal[p];
	}
	// Do a key lookup, and return UINT_MAX if the key is not found.
	inline unsigned getposition(const K & k) const {
		return internal.getposition(k);
	}
	// Get the ordered position of an element in the sort.
	inline unsigned getindex(const K & k) const {
		return internal.getindex(k);
	}
	// End of the line for recursive index look-up
	void getlevels(unsigned (& p)[1], const K & k) const {
		p[0] = internal.getindex(k);
	}

private:
	ktree_avl<K> internal;     // internal index
};

#endif // AXIS_H
