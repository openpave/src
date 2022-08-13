/*************************************************************************

	NUMERICS.CPP - Implementation for various numerical methods.

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

	Portions Copyright (C) 2017-2022 OpenPave.org.

	Contributor(s): Jeremy Lea <reg@openpave.org>.

	See NUMERICS.H.

	History:
		2017/05/09 - Created by Jeremy Lea <reg@openpave.org>

*************************************************************************/

#include "../include/numerics.h"

namespace OP {

/*
 * Brent's method for 1D function minimization, with initial bracketing.
 * This is based around the code from "Numerical Recipes in C",
 * although the algorithms are in much more general use.
 *
 * We start with one value, and construct an interval around it.  Then
 * we use that interval to determine if we have bracketed the minimum.
 * Once we have bracketed the minimum we close in on the final minimum.
 */
double
fmin_brent(const std::function<double(double)> & f, double a, double b)
{
	static const double GC = (3.0-std::sqrt(5.0))/2.0;
	static const double GR = 1/(1-GC);

	double A, B, C, X, Y, Z, c = NAN, x, y, z;
	double r, q;
	unsigned i = 0;

	if (!std::isnan(a) && !std::isnan(b)) {
		// User supplied bracket.  Check it a little.
		if (a == b)
			// We can't start with an empty bracket.
			a = NAN;
		else {
			// Try the golden section mid-point to start.
			c = b, b = b+GC*(a-b);
			A = f(a), B = f(b), C = f(c);
			if (A <= B && B <= C) {
				// The result is not bracketed and going the wrong way...
				b = a, a = c, c = NAN;
				// Don't care about the extra evaluations here.
			}
		}
	}
	if (std::isnan(c)) {
		if (std::isnan(b))
			std::swap(a,b);
		if (std::isnan(b))
			b = 0.0;
		if (std::isnan(a))
			a = (b == 0.0 ? -0.05 : b-b/20);
		// To start with we only have two points.  We need three,
		// so dream up one more...  We think that the direction from a
		// to b is a descent direction, so test that.
		if ((B = f(b)) > (A = f(a)))
			c = a, C = A, a = b+(b-a)*GR, A = f(a);
		else
			c = b+(b-a)*GR, C = f(c);
	}
	// Now search until we bracket a minimum...
	while (A >= B && B >= C) {
		if (++i > 10)
			throw std::runtime_error("Failed to bracket minimum!");
		r = (b-a)*(B-C), q = (b-c)*(B-A);
		x = (r != q ? b-((b-a)*r-(b-c)*q)/std::abs(2*(q-r)) : a);
		y = c+(c-b)*GR;
		// x is the turning point of a parabola through A, B, C.
		if (x > std::min(b,c) && x < std::max(b,c)) {
			X = f(x);
			if (X < C) {
				// A > B > C && B > X < C so we've turned between
				// b and c.  Swap c & x, then let them be shifted so we
				// terminate normally.
				std::swap(x,c), std::swap(C,X);
			} else if (X > B) {
				// A > B > C && A > B < X so we've turned between
				// A and X.  Replace c with x and skip the final shift.
				c = x, C = X;
				continue;
			} else {
				// Inconceivable! The parabola turned, but the function
				// did not. Replace a with b and b with x and skip the
				// shift.
				a = b, b = x, A = B, B = X;
				continue;
			}
		} else if (x > std::min(c,y) && x < std::max(c,y)) {
			X = f(x);
			if (X < C) {
				// The function keeps going down, so take another step.
				b = c, c = x, x = y;
				B = C, C = X, X = f(x);
			}
		} else {
			// The parabola didn't work out, so take a golden section step.
			x = y, X = f(x);
		}
		a = b, b = c, c = x, A = B, B = C, C = X;
	}
	// Now that we know we have a minimum, hunt it down...
	// Since we've already evaluated the function at the end points,
	// use this to initialise Y and Z. A and C are still the outer
	// bracket, and B the current point.  Y and Z are A and C but ordered
	// so Y < Z.  X is the new test point.
	if (a > c)
		std::swap(a,c), std::swap(A,C);
	if (C < A)
		z = a, y = c, Z = A, Y = C;
	else
		z = c, y = a, Z = C, Y = A;
	// Continue tightening the bracket until we're as close as possible...
	while (c-a > 4*eps(b) && (Y > B || Z > B)) {
		if (++i > 100)
			throw std::runtime_error("failed to converge!");
		// Make a parabolic approximation of the minimum.
		r = (b-z)*(B-Y), q = (b-y)*(B-Z);
		x = (r != q ? b-((b-z)*r-(b-y)*q)/(2*(r-q)) : b);
		// if we failed or are outside the centre of the two intervals
		// then just take a conservative golden section step.
		if (r == q || 2*x < (a+b) || 2*x > (b+c))
			x = b+GC*(c-b < b-a ? a-b : c-b);
		X = f(x);
		if (X <= B) {
			// We found a new best point.  Update everything.
			(x < b ? c : a) = b;
			z = y, y = b, b = x, Z = Y, Y = B, B = X;
		} else {
			// Our new point wasn't better, but update what we can.
			(x < b ? a : c) = x;
			if (X <= Y || x == y)
				z = y, y = x, Z = Y, Y = X;
			else if (X <= Z || x == z || y == z)
				z = x, Z = X;
		}
	}
	return b;
}

/*
 * Brent's method for 1D root finding, with initial bracketing, with some
 * improvements. This is based around the Wikipedia code, along with some
 * numerical improvements from Matlab's fzero implementation.  In addition
 * this incorporates my own improvements to use Ridder's, Muller's or the
 * Secant method if they seem appropriate.
 *
 * In addition, this takes limits that will not be crossed in testing the
 * function.  We can start with one value, and construct an interval
 * around it.  This is done by using a doubling of the interval or the
 * secant method.  Without any starting point, zero is assumed.  If a
 * bracket is supplied, it will be used as limits.
 *
 * There are a lot of moving parts involved here.  (x,X) is the latest
 * point, (x,X) is the new test point.  (y,Y) and (z,Z) are the last two
 * function evaluations, sorted by abs(f).  These are used for the
 * approximations.  (x,X) and (c,C) are the bracket.  C always has the
 * opposite sign to X.
 */
double
fzero_core(const std::function<double(double)> & f, double x, double y,
		double ll, double lr)
{
	double c, C, z, X, Y, Z, e, d, m, tol, w = 1.0;
	bool fp = false;
	unsigned i = 0;

	// Initial evaluation
	X = f(x), Y = f(y), z = y, Z = Y;
	if (!std::isfinite(X) || !std::isfinite(Y))
		throw std::runtime_error("Function not sufficiently bounded in fzero()!");
	// Now start looking for a valid bracket.
	while (true) {
		if (std::abs(Y) < std::abs(X))
			std::swap(y,x), std::swap(Y,X);
		if (Y > 0 != X > 0)
			break;
		if (++i > 10)
			throw std::runtime_error("Failed to bracket minimum!");
		e = (x-y);
		if (x+2*e < ll)
			e = (ll-x)/2;
		if (x+2*e > lr)
			e = (lr-x)/2;
		// Linear interpolation with a little extra step.
		d = 1.05*(Y != X ? -X*e/(X-Y) : e);
		if (std::isfinite(d) && std::abs(d) > eps(e)
				&& std::abs(d) < 8*std::abs(e) && x+d > ll && x+d < lr)
			z = y, y = x, x += d;
		else
			z = y, y = x, x += 2*e;
		if (!std::isfinite(x) || x < ll || x > lr)
			throw std::runtime_error("Unable to automatically bracket function in fzero()!");
		Z = Y, Y = X, X = f(x);
	}
	// Now we can actually start looking for a root!
	// (c,C) is always on the other side from (x,X).
	c = y, C = Y, e = d = x-c;
	while (X != 0 && d != 0) {
		// Anderson & Bjork's tweak for false position method.
		if (fp) {
			w = 1-X/Y, fp = false;
			if (w <= 0.1) // zero seems to risky...
				w = 0.5;
		} else
			w = 1.0;
		if (++i > 100)
			throw std::runtime_error("failed to converge on root!");
		// m is the working bracket width.
		m = c-x;
		// tol is the tolerance.  Since most people don't know the right
		// answer this is not user settable.  We use Matlab's eps()
		// function to determine the minimum possible error.
		tol = 2*eps(std::max(std::abs(x),1e-6));
		// Actual loop termination condition: is the bracket small enough.
		if (std::abs(m) <= 2*tol)
			break;
		// Choose a fitter!  These leave p and q as the -numerator and
		// denominator of the step, to be tested later.  If left as zero
		// then other methods are attempted.  We don't do Brent's test for
		// small steps here, or increasing function values.  If the
		// function value increases, one of the fits might still be very
		// good.
		double p = 0.0, q, r, s;
		if (y != z && y != x && x != z) { // parabolic fit possible.
			// Ridder's method.  This fits a parabola in exponential
			// space.  This has quadratic convergence, but theoretically
			// requires two function evaluations, with one at the
			// mid-point.  Since we are often at the mid-point already,
			// test for that condition.
			if (x == y+0.5*(z-y) && X*X > Z*Y) {
				q = std::sqrt(X*X-Z*Y);
				p = (x-y)*std::copysign(X,Z-Y);
			}
			// Check if x and y were swapped.
			if (y == x+0.5*(z-x) && Y*Y > Z*X) {
				q = std::sqrt(Y*Y-Z*X);
				p = (y-x)*(std::copysign(Y,Z-X)-q);
			}
			// Muller's method (quadratic interpolation), if the function
			// is more 'vertical' at this point.
			if (p == 0.0 && std::abs(X-C) > std::abs(m)) {
				p = (Z-Y)/(z-y), q = (Y-X)/(y-x);
				r = (p-q)/(z-x), s = 0.5*(q+(x-y)*r);
				if (s*s-X*r >= 0) {
					p = X;
					q = s+std::copysign(std::sqrt(s*s-X*r),s);
				} else
					p = 0.0;
			}
			// The classic inverse parabolic interpolation from Brent's
			// method.
			if (p == 0.0) {
				s = X/Y, q = Y/Z, r = X/Z;
				p = s*((z-x)*q*(q-r)-(x-y)*(r-1));
				q = (q-1)*(r-1)*(s-1);
			}
		} else {
			// Try the secant method first.
			s = X/Y, p = (y-x)*s, q = 1-s;
			// Resort to false position.
			if (p == 0 || 2*std::abs(p) > std::abs(e*q))
				s = X/(w*C), p = m*s, q = 1-s, fp = true;
		}
		(p > 0 ? q : p) = -(p > 0 ? q : p);
		if (p != 0.0 && q != 0.0 && 2*p < (1.5*m*q-std::abs(tol*q))
				&& 2*p < std::abs(e*q))
			e = d, d = p/q;
		else
			e = d = 0.5*m, fp = false;
		z = y, Z = Y, y = x, Y = X, x += d, X = f(x);
		if (std::abs(d) < tol && std::abs(X) > eps(C))
			// We took a very small step.  If we didn't get a really good
			// result just bisect, and eat the cost of an extra
			// evaluation.  This replaces Brent's small step test.
			e = d = 0.5*m, x = y+d, X = f(x), fp = false;
		// Make sure C is on the other side from X.
		if (X > 0 == C > 0) {
			// check if Y is better than X.
			if (std::abs(Y) < std::abs(X))
				c = x, std::swap(x,y), C = X, std::swap(X,Y);
			else
				c = y, C = Y;
			e = d = c-x, fp = false;
		}
		if (X > 0 == C > 0)
			throw std::runtime_error("oops.");
	}
	return x;
}

} // namespace OP
