/**************************************************************************

	NUMERICS.H - Simple numerical methods

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

	Portions Copyright (C) 2017 OpenPave.org.

	Contributor(s): Jeremy Lea <reg@openpave.org>.

	Purpose:
		This header contains prototypes for several basic numerical
		methods, implemented in C style functions.  Most are based on
		reference implementations from online sources.

	Design:
		These provide the basic building blocks for doing various
		non-linear numerical opertions, such as minimization, root
		finding, etc.

	History:
		2017/05/09 - Created by Jeremy Lea <reg@openpave.org>

**************************************************************************/

#ifndef __NUMERICS_H
#define __NUMERICS_H

#include <functional>
#include <stdexcept>
#include "mathplus.h"

namespace OP {

/*
 * Matlab's eps function, which is used for testing floating point error.
 */
template<typename T>
inline T
eps(T x = 1)
{
	x = std::abs(x);
	return std::nextafter(x,x+1) - x;
}

double fmin_brent(const std::function<double(double)> & f,
		double a = NAN, double b = NAN);
// Core root finding algorithm.
double fzero_core(const std::function<double(double)> & f, double b,
		double y, double ll, double lr);

/*
 * Simple wrapper for fzero_core() to set up initial values correctly.
 */
inline double
fzero(const std::function<double(double)> & f, double b = NAN,
		double ll = -INFINITY, double lr = INFINITY)
{
	double y;

	// Allow a more natural form of putting limit,value,limit.
	if (std::isfinite(b) && std::isfinite(ll) && ll > b) {
		if (std::isfinite(lr))
			std::swap(ll,b);
		else
			lr = ll, ll = b;
	}
	// If all limits are supplied then avoid the exact upper and lower
	// limits.  If one of them is the actual root, then don't call this
	// function!
	if (std::isfinite(b) && std::isfinite(ll) && std::isfinite(lr)
			&& ll < b && lr > b) {
		ll = std::nextafter(ll,b);
		lr = std::nextafter(lr,b); // actually next before...
	}
	// Reject some silly values.
	if (std::isnan(ll) || std::isnan(ll))
		throw std::runtime_error("Cannot use NaN's as a bracket in fzero()!");
	if (lr <= ll)
		throw std::runtime_error("Invalid bracket in fzero()!");
	// If no start value was supplied, then just use zero.
	if (std::isnan(b))
		b = (std::isfinite(ll) ? ll : 0.0);
	// Now guess another start value, if no bracket was supplied.
	if (std::isinf(lr)) {
		y = b+(std::abs(b) < 1 ? 0.05 : b/20);
	} else
		y = (b == lr ? ll : lr);
	return fzero_core(f,b,y,ll,lr);
}

} // namespace OP

#endif // NUMERICS_H
