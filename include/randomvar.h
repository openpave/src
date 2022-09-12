/**************************************************************************

	RANDOMVAR.H - Random variables.

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
		This file contains classes that implement a random variable type,
		where realizations of a variable can be created that take on random
		values.  This is based on the old reliability.h file with more
		modern C++.

	Design:
		These classes a little like variants, in that they can implement
		a number of different random distributions.

		The central feature of the design is that random variables are
		created within a "house" (as the Monte Carlo casino), which
		generates all of the realizations of the random variables in one
		shot, so that they can be correlated.  The achieve this, the house
		stores a huge list of all of the realizations and their associated
		random variables, so that they can be linked together.

	Status:
		These need more types of distributions...

	History:
		2002/01/10 - Created by Jeremy Lea <reg@openpave.org>
		2016/10/19 - Created from from RELIABILITY.H

**************************************************************************/

#pragma once
#ifndef __RANDOMVAR_H
#define __RANDOMVAR_H

#include <cmath>
#include <stdexcept>
#include "mathplus.h"
#include "rng.h"
#include "statistics.h"
#include "set.h"

namespace OP {

#define RV_DIST_PARAM 3                // # of distribution parameters

// Possible types of distribution
enum class distribution {
	normal,
	lognormal
};

// forward declare
template<typename K>
struct realization;
struct random;
template<typename K>
struct randomvar;
template<typename K>
class house;

/*
 * struct realized - a store for a realization of a random variable.
 */
template<typename K>
class realized
{
	friend struct realization<K>;
	friend struct randomvar<K>;
	friend class house<K>;

	realized(const realized &) = delete;
	realized(realized &&) = delete;
	realized & operator = (const realized &) = delete;
	realized & operator = (realized &&) = delete;
	realized(house<K> &, const random &, float );
	~realized();
	// Get a new random number. XXX need covariance matrix.
	void rolldice() noexcept;

	unsigned trial;
	K sortkey;
	house<K> * dealer;
	const random * rv;
	float v;
	unsigned ref_cnt{0};
	bool rolled{false};
};

/*
 * struct realization - a handle to a realization of a random variable.
 */
template<typename K>
struct realization
{
	// Magic typedef for the variant class to know that we are a chameleon
	// that has a different realized type to the abstract variable.
	using cast_t = float;

	realization() noexcept
	  : me(nullptr) {
	}
	realization(const realization & r) noexcept
	  : me(r.me) {
		me->ref_cnt++;
	}
	realization(realization && r) noexcept
	  : me(r.me) {
		r.me = nullptr;
	}
	explicit realization(realized<K> * r) noexcept
	  : me(r) {
		me->ref_cnt++;
	}
	realization & operator = (const realization & r) noexcept {
		if (me && --me->ref_cnt == 0)
			delete me;
		me = r.me;
		me->ref_cnt++;
		return *this;
	}
	realization & operator = (realization && r) noexcept {
		std::swap(me, r.me);
		return *this;
	}
	~realization() {
		if (me && --me->ref_cnt == 0)
			delete me;
	}
	operator const float & () const noexcept {
		return me->v;
	}
	operator float () noexcept {
		return me->v;
	}
	void overwrite(float f) noexcept {
		me->v = f;
	}
	void set_house_sortkey(const unsigned t, const K & k) noexcept {
		me->trial = t;
		me->sortkey = k;
	}
	bool operator == (const realization & a) noexcept {
		return me->trial == a.me->trial && me->sortkey == a.me->sortkey;
	}
	bool operator != (const realization & a) noexcept {
		return me->trial != a.me->trial || me->sortkey != a.me->sortkey;
	}
	bool operator < (const realization & a) noexcept {
		return me->trial < a.me->trial
			|| (me->trial == a.me->trial && me->sortkey < a.me->sortkey);
	}
	bool operator > (const realization & a) noexcept {
		return me->trial > a.me->trial
			|| (me->trial == a.me->trial && me->sortkey > a.me->sortkey);
	}

private:
	friend class house<K>;

	realized<K> * me;
};

/*
 * class house - A centralized store for random variables.
 *
 * In any particular Monte Carlo simulation must have only one house.  All
 * created random variables and realizations store a reference to the house
 * and use that reference to ensure they are tracked.
 */
template<typename K>
class house
{
public:
	house() noexcept
	  : store(), dice() {
	}
	~house() = default;
	randomvar<K> make_rv(const random & r);
	randomvar<K> make_rv(distribution d, float m, float s);
	template<typename T>
	randomvar<K> make_rv(float m, float s);
	float rnd_stdnormal() noexcept {
		return static_cast<float>(dice.stdnormal());
	}
	void rolldice() noexcept;

private:
	friend class realized<K>;
	friend struct randomvar<K>;

	void add_realization(realized<K> * r) {
		store.add(realization<K>(r));
	}
	void rem_realization(realized<K> * r) {
		for (unsigned i = 0; i < store.length(); i++) {
			if (store[i].me == r)
				return store.remove(i);
		}
		throw std::runtime_error("Trying to remove unknown realization!");
	}

	oset<realization<K>> store;       // All the realized variables
	rng dice;                         // This is the source of randomness
};

/*
 * struct random - a generic random variable
 *
 * This wraps the specific distribution types below, and provides an
 * interface that can be used to interact with the variable, including
 * getting realizations.
 */
struct random
{
	random(random &&) = delete;
	random & operator = (const random &) = delete;
	random & operator = (random &&) = delete;
	virtual ~random() = default;
	// Distribution type.
	virtual distribution type() const = 0;
	// Mean value.
	virtual float mean() const noexcept = 0;
	// Standard Deviation.
	virtual float stddev() const noexcept = 0;
	// Get a new random number. XXX need covariance matrix.
	virtual float scale_stdnormal(float z) const noexcept = 0;
	// defined in statistics.cpp
	static random * make_rv(distribution d, float m, float s);
	static random * make_rv(const random & r);

protected:
	// Create a random variable.
	random() noexcept {
		for (size_t i = 0; i < RV_DIST_PARAM; i++)
			d[i] = NAN;
	}
	// copy a random variable.
	random(const random & r) noexcept {
		for (size_t i = 0; i < RV_DIST_PARAM; i++)
			d[i] = r.d[i];
	}
	// Get distribution parameter i
	float param(const size_t i) const {
		if (i > (RV_DIST_PARAM-1))
			throw std::runtime_error("Invalid distribution parameter");
		return d[i];
	}
	// Set distribution parameter i.
	float param(const size_t i, float dp) {
		if (i > (RV_DIST_PARAM-1))
			throw std::runtime_error("Invalid distribution parameter");
		return d[i] = dp;
	}

	float d[RV_DIST_PARAM];             // Four distribution parameters.

private:
	template<typename K>
	friend struct randomvar;

	volatile unsigned ref_cnt{0};       // for randomvar
};

/*
 * struct randomvar - A useable random variable.
 *
 * This is just a class holding two pointers: one to a random variable,
 * and another to a house.  This allows it to be used to realize instances
 * of the random variable.  The random variable is owned by this class.
 */
template<typename K>
struct randomvar
{
	// Magic typedef for the variant class to know that we are a chameleon
	// that has a different realized type to the abstract variable.
	using real_t = realization<K>;

	// Create a random variable.
	randomvar() = delete;
	randomvar(const randomvar & r) noexcept
		: dealer(r.dealer), rv(r.rv) {
		rv->ref_cnt++;
	}
	randomvar(randomvar && r) noexcept
		: dealer(r.dealer), rv(r.rv) {
		r.rv = nullptr;
	}
	randomvar(house<K> * h, random * r)
		: dealer(h), rv(r) {
		if (rv == nullptr)
			throw std::runtime_error("Cannot set a null random variable!");
		rv->ref_cnt++;
	}
	randomvar & operator = (const randomvar & r) noexcept {
		if (rv && --rv->ref_cnt == 0)
			delete rv;
		rv = r.rv;
		rv->ref_cnt++;
		return *this;
	}
	randomvar & operator = (randomvar && r) noexcept {
		std::swap(dealer, r.dealer);
		std::swap(rv, r.rv);
		return *this;
	}
	~randomvar() {
		if (rv && --rv->ref_cnt == 0)
			delete rv;
	}
	realization<K> realize() const {
		realized<K> * r = new realized<K>(*dealer, *rv, rv->mean());
		return realization<K>(r);
	}
	house<K> * get_random() const noexcept {
		return dealer;
	}
	const random & get_rv() const noexcept {
		return *rv;
	}

private:
	house<K> * dealer;
	random * rv;
};

template<distribution D>
struct rv_base
  : public random
{
	static constexpr const distribution distribution_t = D;

	distribution type() const noexcept final {
		return D;
	}

protected:
	rv_base() = default;
	rv_base(const random & r) noexcept
	  : random(r) {
	}
};

struct rv_normal
  : public rv_base<distribution::normal>
{
	rv_normal(float m, float s)
	  : rv_base() {
		param(0,m);
		param(1,s);
	}
	float mean() const noexcept final {
		return d[0];
	}
	float stddev() const noexcept final {
		return d[1];
	}
	float scale_stdnormal(float z) const noexcept final {
		return static_cast<float>(d[0]+d[1]*z);
	}

protected:
	friend struct random;

	rv_normal(const random & r) noexcept
	  : rv_base(r) {
	}
};

struct rv_lognormal
  : public rv_base<distribution::lognormal>
{
	rv_lognormal(float m, float s)
	  : rv_base() {
		param(1, std::isinf(m) ? NAN :
			std::isnan(m) || std::isnan(s) || m <= 0 || s < 0 ? NAN :
			sqrt(log(1 + s*s / m / m)));
		param(0, std::isinf(m) ? INFINITY :
			std::isnan(m) || std::isnan(s) || m <= 0 || s < 0 ? NAN :
			log(m) - d[1] * d[1] / 2);
	}
	float mean() const noexcept final {
		return std::isinf(d[0]) ? INFINITY :
		       std::isnan(d[0]) ? NAN : exp(d[0] + d[1]*d[1]/2);
	}
	float stddev() const noexcept final {
		return std::isinf(d[0]) || std::isnan(d[0]) ? NAN :
			std::isnan(d[1]) ? NAN : mean()*sqrt(exp(d[1]*d[1])-1);
	}
	float scale_stdnormal(float z) const noexcept final {
		return std::isinf(d[0]) ? INFINITY :
			std::isnan(d[0]) || std::isnan(d[1]) ? NAN :
			static_cast<float>(exp(d[0] + d[1] * z));
	}

protected:
	friend struct random;

	rv_lognormal(const random & r) noexcept
	  : rv_base(r) {
	}
};

// Out-of-line constructor to add ourselves to the house.
template<typename K>
inline
realized<K>::realized(house<K> & h, const random & r, float d)
  : trial(UINT_MAX), sortkey(), dealer(&h), rv(&r), v(d)
{
	if (dealer != nullptr)
		dealer->add_realization(this);
}

// Out-of-line constructor to add ourselves to the house.
template<typename K>
inline
realized<K>::~realized()
{
	try {
		if (dealer != nullptr)
			dealer->rem_realization(this);
	} catch (...) {}
}

template<typename K>
inline void
realized<K>::rolldice() noexcept
{
	if (!rolled)
		v = rv->scale_stdnormal(dealer->rnd_stdnormal());
	rolled = true;
}

template<typename K>
inline void
house<K>::rolldice() noexcept
{
	store.sort();
	for (unsigned i = 0; i < store.length(); i++) {
		store[i].me->rolldice();
	}
}

// Let the house make random variables too..
template<typename K>
inline randomvar<K>
house<K>::make_rv(const random & r)
{
	return randomvar<K>(this,random::make_rv(r));
}

template<typename K>
inline randomvar<K>
house<K>::make_rv(distribution d, float m, float s)
{
	return randomvar<K>(this,random::make_rv(d,m,s));
}

template<typename K>
template<typename T>
inline randomvar<K>
house<K>::make_rv(float m, float s)
{
	return randomvar<K>(this,random::make_rv(T::distribution_t,m,s));
}

#undef RV_DIST_PARAM

} // namespace OP

#endif // RANDOMVAR_H
