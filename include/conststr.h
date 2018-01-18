/**************************************************************************

	CONSTSTR.H - Compile-time constant strings

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
		This header implements a basic string class for compile time
		string constants, including compile time hashing.

	History:
		2017/09/06 - Created a basic implementation.

		Based on CTTI <https://github.com/Manu343726/ctti> and other
		on-line sources.

**************************************************************************/

/* CTTI License

The MIT License

Copyright (c) 2015 Manuel Sánchez

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/

#ifndef __CONSTSTR_H
#define __CONSTSTR_H

#include <algorithm>
#include <cstring>
#include <string>

namespace OP {

/*
 * class conststr - A compile time constant string
 *
 * This handles compile time constant strings.  These
 * can only be static const values in the code (e.g.
 * explicit strings).  The strings need not be null
 * terminated, but you cannot mix them up, and you
 * need to be sure they are null terminated before
 * getting a char * for printing.
 */
class conststr
{
public:
	typedef std::size_t size_t;
	typedef std::size_t hash_t;

	// This is the default constructor within most code.
	// Template deduction is used to get the length.
	template<size_t N>
	constexpr conststr(const char (&s)[N]) :
		str{s}, len{N-1} {}
	// Explicit construction with known length.
	constexpr conststr(const char * const s, size_t l) :
		str{s}, len{l} {}
	// Explicit construction from a pointer.  Must be null
	// terminated.
	explicit constexpr conststr(const char * const s) :
		conststr{s,calclen(s)} {}
	// Explicit construction from a std::string.  Mostly
	// should only be used for run time tests.
	explicit conststr(const std::string & s) :
		str{s.c_str()}, len{s.length()} {}
	// Explicit constructor from begin and end pointers.
	constexpr conststr(const char * const b, const char * const e) :
		conststr{b,static_cast<size_t>(e-b)} {}
	// Get the length of the string  (without null).
	constexpr size_t length() const {
		return len;
	}
	// Get an explicit hash from the string.
	constexpr hash_t hash() const {
		return fnv1a_hash(len,str);
	}
	// Return a std::string for printing, etc.
	operator std::string () const {
		return {str,str+len};
	}
	// only call this if you know the strings are null terminated.
	constexpr operator const char * () const {
		return str;
	}
	// Begin and end for iterators.
	constexpr const char * begin() const {
		return str;
	}
	constexpr const char * end() const {
		return str+len;
	}
	// Get chars from the string.
	constexpr char operator [] (size_t i) const {
		return str[i];
	}
	// Get a chunk from a const string.
	constexpr conststr operator ()(size_t b, size_t e) const {
		return {str+b,str+e};
	}
	// Compares need to be done carefully for non-null terminated
	// strings.
	int compare(const conststr & k) const {
		int c = strncmp(str,k.str,std::min(len,k.len));
		return len == k.len ? c : (c == 0 ? (len < k.len ? -1 : 1) : c);
	}
	bool operator == (const conststr & k) const {
		return compare(k) == 0;
	}
	bool operator != (const conststr & k) const {
		return compare(k) != 0;
	}
	bool operator < (const conststr & k) const {
		return compare(k) < 0;
	}
	bool operator <= (const conststr & k) const {
		return compare(k) <= 0;
	}
	bool operator > (const conststr & k) const {
		return compare(k) > 0;
	}
	bool operator >= (const conststr & k) const {
		return compare(k) >= 0;
	}

private:
	// The actual string pointer.
	const char * const str;
	// The length without a possible terminating null.
	size_t len;

	// The FNV1A hash function.
	static constexpr const hash_t basis =
		sizeof(hash_t) >= 8 ? 0xcbf29ce484222325 : 0x811c9dc5;
	static constexpr const hash_t prime =
		sizeof(hash_t) >= 8 ? 0x00000100000001b3 : 0x01000193;
	static constexpr hash_t
	fnv1a_hash(size_t n, const char * const str, hash_t hash = basis) {
		return n > 0 ? fnv1a_hash(n-1,str+1,
			(hash^static_cast<hash_t>(*str))*prime) : hash;
	}
	// Static strlen for explicit constructor.
	static constexpr size_t
	calclen(const char * const s) {
		return *s ? 1+calclen(s+1) : 0;
	}
};

} // namespace OP

#endif // CONSTSTR_H
