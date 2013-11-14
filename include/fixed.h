/**************************************************************************

	FIXED.H - A fixed point type

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

	The Initial Developer of the Original Software is Jeremy Lea.

	Portions Copyright (C) 2006-2008 OpenPave.org.

	Contributor(s): Jeremy Lea <reg@openpave.org>.

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
		2006/07/20 - Created by Jeremy Lea <reg@openpave.org>

**************************************************************************/

#ifndef __FIXED_H
#define __FIXED_H

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
	}
	inline fixed(const fixed & f) {
		value = f.value;
	}
	inline fixed(const int i) {
		value = i << P;
	}
	inline fixed(const unsigned i) {
		value = i << P;
	}
	inline fixed(const float f) {
		value = (f >= 0 ? int(f*(1 << P) + 0.5)
			 : -int(0.5 - f*(1 << P)) );
	}
	inline fixed(const double d) {
		value = (d >= 0 ? int(d*(1 << P) + 0.5)
			 : -int(0.5 - d*(1 << P)) );
	}
	inline ~fixed() {
	}
	
	// Some assignment operators...
	inline fixed & operator = (const fixed & f) {
		value = f.value;
		return *this;
	}
	inline fixed & operator += (const fixed & f) {
		value += f.value;
		return *this;
	}
	inline fixed & operator -= (const fixed & f) {
		value -= f.value;
		return *this;
	}
	inline fixed & operator *= (const fixed & f) {
		*this = double(*this) * double(f);
		return *this;
	}
	inline fixed & operator /= (const fixed & f) {
		*this = double(*this) / double(f);
		return *this;
	}
	
	// Some conversion operators...
	inline fixed<P> & operator = (const int i) {
		value = i << P;
		return *this;
	}
	inline fixed<P> & operator = (const double d) {
		value = (d >= 0 ? int(d*(1 << P) + 0.5)
			 : -int(0.5 - d*(1 << P)) );
		return *this;
	}
	inline operator int () const {
		return value >> P;
	}
	inline operator float () const {
		return (float(value))/(1 << P);
	}
	inline operator double () const {
		return (double(value))/(1 << P);
	}
  
	// Comparison operators...
	inline bool operator == (const fixed & f) const {
		return (value == f.value);
	}
	inline bool operator != (const fixed & f) const {
		return (value != f.value);
	}
	inline bool operator  < (const fixed & f) const {
		return (value < f.value);
	}
	inline bool operator <= (const fixed & f) const {
		return (value <= f.value);
	}
	inline bool operator  > (const fixed & f) const {
		return (value > f.value);
	}
	inline bool operator >= (const fixed & f) const {
		return (value >= f.value);
	}
	
	// Some math operators...
	inline fixed operator - () const {
		fixed t;
		t.value = -value;
		return t;
	}
	inline fixed operator + (const fixed<P> & f) const {
		fixed t;
		t.value = value + f.value;
		return t;
	}
	inline fixed operator - (const fixed<P> & f) const {
		fixed t;
		t.value = value - f.value;
		return t;
	}
	inline fixed operator * (const fixed<P> & f) const {
		fixed t = double(*this) * double(f);
		return t;
	}
	inline fixed operator / (const fixed<P> & f) const {
		fixed t = double(*this) / double(f);
		return t;
	}

	// The prefix and postfix operators, for convienience.
	inline fixed & operator ++ () {
		value += (1 << P);
		return *this;
	}
	inline fixed & operator -- () {
		value -= (1 << P);
		return *this;
	}
	inline fixed operator ++ (int) {
		fixed t(*this);
		value += (1 << P);
		return t;
	}
	inline fixed operator -- (int) {
		fixed t(*this);
		value -= (1 << P);
		return t;
	}

  private:
	template <unsigned int P1>
	friend fixed<P1> fabs(const fixed<P1> & f);
	int value;
};

template<unsigned int P>
inline fixed<P> fabs(const fixed<P> & f) {
	fixed<P> t;
	t.value = abs(f.value);
	return t;
}

#endif // FIXED_H
