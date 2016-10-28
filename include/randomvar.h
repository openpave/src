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

	Portions Copyright (C) 2006-2008 OpenPave.org.

	Contributor(s): Jeremy Lea <reg@openpave.org>.

	Purpose:
		This file contains classes that implement a random variable type,
		where realizations of a variable can be created that take on random
		values.  This is based on the old reliability.h file with more modern
		C++.

	Design:
		These classes a little like variants, in that they can implement
		a number of different random distributions.

		The central feature of the design is that random variables are
		created within a "house" (as the Monte Carlo casino), which generates
		all of the realizations of the random variables in one shot, so that
		they can be correlated.  The achieve this, the house stores a huge
		list of all of the realizations and their associated random variables,
		so that they can be linked together.
	
	Status:
		These need more types of distributions...
  
	History:
		2002/01/10 - Created by Jeremy Lea <reg@openpave.org>
		2016/10/19 - Created from from RELIABILITY.H

**************************************************************************/

#ifndef __RANDOMVAR_H
#define __RANDOMVAR_H

#include <cmath>
#include <stdexcept>
#include <typeindex>
#include <tuple>
#include "mathplus.h"
#include "rng.h"
#include "statistics.h"
#include "tree.h"

namespace OP {

#define RV_DIST_PARAM		3			// # of distribution parameters

// Possible types of distribution
enum class distribution {
	normal,
	lognormal
};

// forward declare
struct random;
struct house;

/*
 * struct realization - a store for a realization of a random variable.
 */
struct realization
{
	realization(house & h, const random & r, double d);
	~realization();
	operator const double & () const {
		return v;
	}
	operator double () {
		return v;
	}
	house * dealer;
	const random * rv;
	double v;
};

/*
 * struct house - A centralized store for random variables.
 *
 * In any particular Monte Carlo simulation must have only one house.  All
 * created random variables and realizations store a reference to the house
 * and use that reference to ensure they are tracked.
 */
struct house
{
	house() :
		store(), vars(), dice() {
	}
	~house();
	random * make_rv(const random & r);
	random * make_rv(distribution d, double m, double s);
	template<typename T>
	random * make_rv(double m, double s);

private:
	friend struct realization;
	void add_realization(const realization * r) {
		store.add(info(r));
	}
	void rem_realization(const realization * r) {
		store.remove(key(r));
	}
	struct key {
		key(const realization * r) :
			x(r) {
		}
		int compare (const key & k) const {
			return (x==k.x ? 0 : SGN(x-k.x));
		}
		const realization * x;
	};
	struct info : public key {
		info(const realization * r) :
			key(r), X(r->rv) {
		}
		info & operator = (const info & i) {
			x = i.x;
			X = i.X;
			return *this;
		}
		const random * X;
	};
	ktree_avl<key,info> store;

	friend struct random;
	void add_variable(random * r) {
		vars.add(r);
	}
	void rem_variable(random * r) {
		vars.remove(r);
	}
	ktree_avl<random *> vars;

	rng dice;
};

/*
 * struct random - a generic random variable
 *
 * This wraps the specific distribution types below, and provides an interface
 * that can be used to interact with the variable, including getting
 * realizations.
 */
struct random
{
	// Magic typedef for the variant class to know that we are a cameleon that
	// has a different realized type to the abstract variable.
	using real = realization;

	virtual ~random() {
		if (dealer != nullptr)
			dealer->rem_variable(this);
	}
	// Distribution type.
	virtual distribution type() const = 0;
	// Mean value.
	virtual double mean() const = 0;
	// Standard Deviation.
	virtual double stddev() const = 0;
	// Get the expected value for use in variants
	realization expected() const {
		return realization(*dealer,*this,mean());
	}
	static random * make_rv(distribution d, double m, double s,
		house * h = nullptr);
	static random * make_rv(const random & r, house * h = nullptr);

protected:
	house * dealer;                     // The house - can be null.
	double d[RV_DIST_PARAM];			// Four distribution parmeters.

	// Create a random variable.
	random(house * h = nullptr) :
		dealer(h) {
		for (size_t i = 0; i < RV_DIST_PARAM; i++)
			d[i] = NAN;
		if (dealer != nullptr)
			dealer->add_variable(this);
	}
	// copy a random variable.
	random(const random & r, house * h = nullptr) :
		dealer(h == nullptr ? r.dealer : h) {
		for (size_t i = 0; i < RV_DIST_PARAM; i++)
			d[i] = r.d[i];
		if (dealer != nullptr)
			dealer->add_variable(this);
	}
	// Get distribution parmater i
	double param(size_t i) {
		if (i > (RV_DIST_PARAM-1))
			throw std::runtime_error("Invalid distribution parameter");
		return d[i];
	}
	// Set distribution parmater i.
	double param(size_t i, double dp) {
		return d[i] = dp;
	}
};

template<distribution D>
struct random_var :
	public random
{
	static constexpr const distribution distribution_t = D;

	virtual distribution type() const override final {
		return D;
	}
protected:
	using random::random;
};

struct rv_normal :
	public random_var<distribution::normal>
{
	virtual double mean() const override final {
		return d[0];
	}
	virtual double stddev() const override final {
		return d[1];
	}

	rv_normal(double m, double s, house * h = nullptr) :
		random_var(h) {
		param(0,m);
		param(1,s);
	}
protected:
	friend struct random;
	rv_normal(const random & r, house * h = nullptr) :
		random_var(r,h) {
	}
};

struct rv_lognormal :
	public random_var<distribution::lognormal>
{
	virtual double mean() const override final {
		return std::isinf(d[0]) ? INFINITY :
			   std::isnan(d[0]) ? NAN : exp(d[0] + d[1]*d[1]/2);
	}
	virtual double stddev() const override final {
		return std::isinf(d[0]) ? INFINITY :
			   std::isnan(d[0]) ? NAN : mean()*sqrt(exp(d[1]*d[1])-1);
	}

	rv_lognormal(double m, double s, house * h = nullptr) :
		random_var(h) {
		param(1,std::isinf(m) ? NAN :
			    m <= 0 ? NAN : sqrt(log(1+s*s/m/m)));
		param(0,std::isinf(m) ? INFINITY :
			    m <= 0 ? NAN : log(m) - d[1]*d[1]/2);
	}
protected:
	friend struct random;
	rv_lognormal(const random & r, house * h = nullptr) :
		random_var(r,h) {
	}
};

// Out-of-line constructor to add ourselves to the house.
inline
realization::realization(house & h, const random & r, double d) :
	dealer(&h), rv(&r), v(d) {
	if (dealer != nullptr)
		dealer->add_realization(this);
}

// Out-of-line constructor to add ourselves to the house.
inline
realization::~realization() {
	if (dealer != nullptr)
		dealer->rem_realization(this);
}

// Out-line-destructor to call ~random().
inline
house::~house() {
	for (unsigned i = 0; i < vars.length(); i++)
		delete vars[i];
}

// Let the house make random variables too..
inline random *
house::make_rv(const random & r)
{
	return random::make_rv(r,this);
}

inline random *
house::make_rv(distribution d, double m, double s)
{
	return random::make_rv(d,m,s,this);
}

template<typename T>
inline random *
house::make_rv(double m, double s)
{
	return random::make_rv(T::distribution_t,m,s,this);
}

#undef RV_DIST_PARAM

} // namespace OP

#endif // RANDOMVAR_H