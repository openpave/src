/**************************************************************************

	CTTI.H - Compile-time type information

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
		This header implements a basic compile time typing class, to
		avoid the poor design of the C++ RTTI classes.

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

#ifndef __CTTI_H
#define __CTTI_H

#include <functional>
#include <type_traits>
#include "conststr.h"

namespace OP {

/*
 * class type_info - Compile time type information
 *
 * This follows std::type_info as closely as possible, given the limits
 * of having to be compile time constant.
 */
class type_info
{
public:
	typedef OP::conststr::hash_t hash_t;

	// Do not allow default construction.
	type_info() = delete;
	// Only construct directly from a name.
	constexpr type_info(const OP::conststr & name) :
		magic{name} {}
	// Destructor cannot be virtual
	~type_info() = default;
	// Allow copy assignment unlike std::type_info.
	constexpr type_info(const OP::type_info &) = default;
	constexpr type_info(OP::type_info &&) = default;
	OP::type_info & operator = (const OP::type_info &) = default;
	OP::type_info & operator = (OP::type_info &&) = default;
	// Return our hash
	constexpr hash_t hash_code() const {
		return magic.hash_code();
	}
	// Note that the name is probably not null terminated.
	constexpr const OP::conststr & name() const {
		return magic;
	}
	bool before(const type_info & r) const {
		return magic < r.name();
	}
	constexpr bool operator == (const type_info & r) const {
		return hash_code() == r.hash_code();
	}
	constexpr bool operator != (const type_info & r) const {
		return hash_code() != r.hash_code();
	}

private:
	OP::conststr magic;
};

/*
 * class type_index - Compile time type index (hash)
 *
 * This follows std::type_index as closely as possible, given the limits
 * of having to be compile time constant.
 */
class type_index
{
public:
	typedef OP::conststr::hash_t hash_t;

	type_index() = delete;
	// Default constructor
	constexpr type_index(const OP::type_info & info) :
		type{&info}, hash{info.hash_code()} {}
	// Default copy and move constructors and assignments
	constexpr type_index(const OP::type_index &) = default;
	constexpr type_index(OP::type_index &&) = default;
	OP::type_index & operator = (const OP::type_index &) = default;
	OP::type_index & operator = (OP::type_index &&) = default;
	// Get the hash directly.
	constexpr hash_t hash_code() const {
		return hash;
	}
	// Note that the name is probably not null terminated.
	constexpr const OP::conststr & name() const {
		return type->name();
	}
	constexpr bool operator == (const type_index & r) const {
		return hash == r.hash_code();
	}
	constexpr bool operator != (const type_index & r) const {
		return hash != r.hash_code();
	}
	constexpr bool operator < (const type_index & r) const {
		return hash < r.hash_code();
	}
	constexpr bool operator <= (const type_index & r) const {
		return hash <= r.hash_code();
	}
	constexpr bool operator > (const type_index & r) const {
		return hash > r.hash_code();
	}
	constexpr bool operator >= (const type_index & r) const {
		return hash >= r.hash_code();
	}

private:
	const OP::type_info * type;
	hash_t hash;
};

// The actual magic happens here.
#if defined(__clang__)
	#define CTTI __PRETTY_FUNCTION__
	#define CTTI_PREFIX "OP::type_info OP::type_deducer() [T = "
	#define CTTI_SUFFIX "]"
#elif defined(__GNUC__) && !defined(__clang__)
	#define CTTI __PRETTY_FUNCTION__
	#define CTTI_PREFIX "constexpr OP::type_info OP::type_deducer() [with T = "
	#define CTTI_SUFFIX "]"
#elif defined(_MSC_VER)
	#define CTTI __FUNCSIG__
	#define CTTI_PREFIX "class OP::type_info __cdecl OP::type_deducer<"
	#define CTTI_SUFFIX ">(void)"
#else
	#error "No support for this compiler."
#endif
#define CTTI_BEGIN (sizeof(CTTI_PREFIX)-1)
#define CTTI_END ((sizeof(CTTI)-1)-(sizeof(CTTI_SUFFIX)-1))

// Internal only function, to infer the type
template<typename T>
constexpr OP::type_info type_deducer()
{
	// This constructs a conststr using the compiler macro, then chops a bit
	// out of it - the deduced template parameter(not null terminated) and then
	// constructs a type_info object from that.
	return { OP::conststr{CTTI}(CTTI_BEGIN,CTTI_END) };
};

#undef CTTI_END
#undef CTTI_BEGIN
#undef CTTI_SUFFIX
#undef CTTI_PREFIX
#undef CTTI

// Compile-time equivalent of the typeid() keyword.
template<typename T>
constexpr type_info
type_id(T&&)
{
	return type_deducer<typename std::decay<T>::type>();
}
template<typename T>
constexpr type_info
type_id()
{
	return type_deducer<T>();
}

} // namespace OP

// Inject the hashing into the standard name space, like type_info and
// type_index.
namespace std {

template<>
struct hash<OP::type_info>
{
	constexpr std::size_t operator()(const OP::type_info & id) const {
		return id.hash_code();
	}
};
template<>
struct hash<OP::type_index>
{
	constexpr std::size_t operator()(const OP::type_index & id) const {
		return id.hash_code();
	}
};

}

#endif // CTTI_H
