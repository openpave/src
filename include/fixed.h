/**************************************************************************

	FIXED.H - A fixed point type

	$OpenPave$

	Purpose:
		This header implements fixed point type.  This is suitable for
		number which need to be 'almost' floating point, but with no round
		off error.

		This does not implement a 'fixed precision floating point' type,
		which would only store a limited number of places after the decimal
		point.  That would require working in a decimal space not in a
		binary (power of two) space.  Just use integers for that...

	Design:
		Based on code from various places on the Internet...  

	History:
		2006/07/20 - Created by Jeremy Lea <jdlea@ucdavis.edu>

**************************************************************************/

#ifndef __FIXED_H
#define __FIXED_H

#include "config.h"

/*
 * class fixed - A fixed point type.
 *
 * The precision template arguement is the number of bits to the left of
 * the decimal point.
 */
template <unsigned int P>
class fixed {
  public:
	// Some constructors for common cases...
	inline fixed() {
	};
	inline fixed(const fixed & f) {
		value = f.value;
	};
	inline fixed(const int i) {
		value = i << P;
	};
	inline fixed(const double d) {
		value = (d >= 0 ? int(d*(1 << P) + 0.5)
			 : -int(0.5 - d*(1 << P)) );
	};
	inline ~fixed() {
	};
	
	// Some assignment operators...
	inline fixed & operator = (const fixed & f) {
		value = f.value;
		return *this;
	};
	inline fixed & operator += (const fixed & f) {
		value += f.value;
		return *this;
	};
	inline fixed & operator -= (const fixed & f) {
		value -= f.value;
		return *this;
	};
	inline fixed & operator *= (const fixed & f) {
		value *= f.value;
		return *this;
	};
	inline fixed & operator /= (const fixed & f) {
		value /= f.value;
		return *this;
	};
	
	// Some conversion operators...
	inline fixed<P> & operator = (const int i) {
		value = i << P;
		return *this;
	};
	inline fixed<P> & operator = (const double d) {
		value = (d >= 0 ? int(d*(1 << P) + 0.5)
			 : -int(0.5 - d*(1 << P)) );
		return *this;
	};
	inline operator int () const {
		return value >> P;
	};
	inline operator double () const {
		return (double(value))/(1 << P);
	};
  
	// Comparison operators...
	inline bool operator == (const fixed & f) const {
		return (value == f.value);
	};
	inline bool operator != (const fixed & f) const {
		return (value != f.value);
	};
	inline bool operator  < (const fixed & f) const {
		return (value < f.value);
	};
	inline bool operator <= (const fixed & f) const {
		return (value <= f.value);
	};
	inline bool operator  > (const fixed & f) const {
		return (value > f.value);
	};
	inline bool operator >= (const fixed & f) const {
		return (value >= f.value);
	};
	
	// Some math operators...
	inline fixed operator - () const {
		fixed t;
		t.value = -value;
		return t;
	};
	inline fixed operator + (const fixed<P> & f) const {
		fixed t;
		t.value = value + f.value;
		return t;
	};
	inline fixed operator - (const fixed<P> & f) const {
		fixed t;
		t.value = value - f.value;
		return t;
	};
	inline fixed operator * (const fixed<P> & f) const {
		fixed t;
		t.value = value * f.value;
		return t;
	};
	inline fixed operator / (const fixed<P> & f) const {
		fixed t;
		t.value = value / f.value;
		return t;
	};

	// The prefix and postfix operators, for convienience.
	inline fixed & operator ++ () {
		value += (1 << P);
		return *this;
	};
	inline fixed & operator -- () {
		value -= (1 << P);
		return *this;
	};
	inline fixed operator ++ (int) {
		fixed t(*this);
		value += (1 << P);
		return t;
	};
	inline fixed operator -- (int) {
		fixed t(*this);
		value -= (1 << P);
		return t;
	};

  private:
	int value;
};


#endif // FIXED_H

