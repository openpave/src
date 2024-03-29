/**************************************************************************

	FIXED.H - A fixed point type

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

#pragma once
#ifndef __FIXED_H
#define __FIXED_H

namespace OP {

/*
 * class fixed - A fixed point type.
 *
 * The precision template argument is the number of bits to the left of
 * the decimal point.
 */
template<unsigned int P>
class fixed {
public:
	// Some constructors for common cases...
	fixed() noexcept {
	}
	fixed(const fixed & f) noexcept {
		value = f.value;
	}
	fixed(int i) noexcept {
		value = i << P;
	}
	fixed(unsigned i) noexcept {
		value = i << P;
	}
	fixed(float f) noexcept {
		value = (f >= 0 ? int(f*(1 << P) + 0.5)
			 : -int(0.5 - f*(1 << P)) );
	}
	fixed(double d) noexcept {
		value = (d >= 0 ? int(d*(1 << P) + 0.5)
			 : -int(0.5 - d*(1 << P)) );
	}
	~fixed() {
	}

	// Some assignment operators...
	fixed & operator = (const fixed & f) noexcept {
		value = f.value;
		return *this;
	}
	fixed & operator += (const fixed & f) noexcept {
		value += f.value;
		return *this;
	}
	fixed & operator -= (const fixed & f) noexcept {
		value -= f.value;
		return *this;
	}
	fixed & operator *= (const fixed & f) noexcept {
		*this = double(*this) * double(f);
		return *this;
	}
	fixed & operator /= (const fixed & f) noexcept {
		*this = double(*this) / double(f);
		return *this;
	}

	// Some conversion operators...
	fixed<P> & operator = (int i) noexcept {
		value = i << P;
		return *this;
	}
	fixed<P> & operator = (double d) noexcept {
		value = (d >= 0 ? int(d*(1LL << P) + 0.5)
			 : -int(0.5 - d*(1LL << P)) );
		return *this;
	}
	operator int () const noexcept {
		return value >> P;
	}
	operator float () const noexcept {
		return (float(value))/(1LL << P);
	}
	operator double () const noexcept {
		return (double(value))/(1LL << P);
	}

	// Comparison operators...
	bool operator == (const fixed & f) const noexcept {
		return (value == f.value);
	}
	bool operator != (const fixed & f) const noexcept {
		return (value != f.value);
	}
	bool operator < (const fixed & f) const noexcept {
		return (value < f.value);
	}
	bool operator <= (const fixed & f) const noexcept {
		return (value <= f.value);
	}
	bool operator > (const fixed & f) const noexcept {
		return (value > f.value);
	}
	bool operator >= (const fixed & f) const noexcept {
		return (value >= f.value);
	}

	// Some math operators...
	fixed operator - () const noexcept {
		fixed t;
		t.value = -value;
		return t;
	}
	fixed operator + (const fixed<P> & f) const noexcept {
		fixed t;
		t.value = value + f.value;
		return t;
	}
	fixed operator - (const fixed<P> & f) const noexcept {
		fixed t;
		t.value = value - f.value;
		return t;
	}
	fixed operator * (const fixed<P> & f) const noexcept {
		fixed t = double(*this) * double(f);
		return t;
	}
	fixed operator / (const fixed<P> & f) const noexcept {
		fixed t = double(*this) / double(f);
		return t;
	}

	// The prefix and postfix operators, for convenience.
	fixed & operator ++ () noexcept {
		value += (1 << P);
		return *this;
	}
	fixed & operator -- () noexcept {
		value -= (1 << P);
		return *this;
	}
	fixed operator ++ (int) noexcept {
		fixed t(*this);
		value += (1 << P);
		return t;
	}
	fixed operator -- (int) noexcept {
		fixed t(*this);
		value -= (1 << P);
		return t;
	}

private:
	template<unsigned int P1>
	friend fixed<P1> fabs(const fixed<P1> & f) noexcept;

	int value;
};

template<unsigned int P>
inline fixed<P>
fabs(const fixed<P> & f) noexcept
{
	fixed<P> t;
	t.value = abs(f.value);
	return t;
}

} // namespace OP

#endif // FIXED_H
