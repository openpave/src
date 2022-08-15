/**************************************************************************

	HASCOMPARE.H - Simple SFINAE class to check if we have compare()

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

	Portions Copyright (C) 2016-2022 OpenPave.org.

	Contributor(s): Jeremy Lea <reg@openpave.org>.

	Purpose:
		This header implements templated C++ classes for SFINAE
		(Substitution failure is not an error) checking if a class has
		a compare() method that returns -1/0/1.  It is used to specialise
		various classes that need to to order comparisons.

	History:
		2016/03/08 - Created a basic implementation.

**************************************************************************/

#pragma once
#ifndef __HASCOMPARE_H
#define __HASCOMPARE_H

#include <type_traits>
#include <tuple>

namespace OP {

// Macroized class that creates a template meta-program that tests if a
// class has a member function with a certain signature (name, return type
// and arguements).
#define HAS_MEMBER_FUNCTION(NAME, RET, ...)                               \
template<typename T>                                                      \
class has_##NAME                                                          \
{                                                                         \
	template<typename C>                                                  \
	static typename std::is_same<decltype(std::declval<C>().NAME(         \
			__VA_ARGS__)),RET>::type                                      \
	check(C *) {                                                          \
		return {};                                                        \
	}                                                                     \
	template<typename>                                                    \
	static std::false_type                                                \
	check(...) {                                                          \
		return {};                                                        \
	}                                                                     \
	using type = decltype(check<T>(nullptr));                             \
public:                                                                   \
	static constexpr bool value = type::value;                            \
};

// Macroized class that creates a template meta-program that tests if a
// class has a nested class with a certain name.
#define HAS_NESTED_CLASS(NAME)                                            \
template<typename T>                                                      \
class has_##NAME                                                          \
{                                                                         \
	template<typename C>                                                  \
	static typename std::is_class<typename C::NAME>::type                 \
	check(C *) {                                                          \
		return {};                                                        \
	}                                                                     \
	template<typename>                                                    \
	static std::false_type                                                \
	check(...) {                                                          \
		return {};                                                        \
	}                                                                     \
	using type = decltype(check<T>(nullptr));                             \
public:                                                                   \
	static constexpr bool value = type::value;                            \
};

// Create a has_compare meta-program to check if a class has a compare()
// method that returns -1/0/1.  It is used to specialise various classes
// that need to do order comparisons.
HAS_MEMBER_FUNCTION(compare,int,std::declval<C>())

/*
 * is_callable - check if a type can be called
 *
 * Template meta-program using SFINAE to test if a type can be called as a
 * function.  This implies no parameters or return value.  Replace once
 * C++17 version is available.
 */
template<typename F>
class is_callable
{
	// base of F
	using T = typename std::remove_cv<
			typename std::remove_reference<F>::type>::type;

	// check1 tests for a callable type
	template<typename C>
	static auto check1(C * c) -> decltype((*c)(),std::true_type()) {
		return {};
	}
	template<typename>
	static std::false_type check1(...) {
		return {};
	}
	// check2 tests for operator () (a functor)
	template<typename C>
	static auto check2(C *)	-> decltype(&C::operator(),std::true_type()) {
		return {};
	}
	template<typename>
	static std::false_type check2(...) {
		return {};
	}
public:
	static constexpr bool value = decltype(check1<T>(nullptr))::value
	                           || decltype(check2<T>(nullptr))::value;
};

/*
 * class function_traits - get various details about a function
 *
 * This looks at callable objects and extracts the return type and
 * parameters.  We really only need arg<0> from this...
 */
// Forward declare
template<typename>
struct function_traits;
// Specialize for functors and lambdas.
template<typename F>
struct function_traits
  : public function_traits<decltype(&F::operator())>
{
};
// Real meta-program, using template deduction to split type
template<typename C, typename R, typename... Ts>
struct function_traits<R(C::*)(Ts...) const>
{
	// Get the number of function arguments
	static constexpr std::size_t arg_count = sizeof...(Ts);
	// Get the return type of the function
	using result_type = R;
	// Get the type of the Nth argument
	template<std::size_t N>
	struct arg {
		static_assert(N < arg_count,"error: index exceeds argument count");
		using type = typename std::tuple_element<N,std::tuple<Ts...>>::type;
	};
};
// Specialize for references
template<typename F>
struct function_traits<F &>
  : public function_traits<F>
{
};
// Specialize for const references
template<typename F>
struct function_traits<const F &>
  : public function_traits<F>
{
};
// Specialize for r-value references
template<typename F>
struct function_traits<F &&>
  : public function_traits<F>
{
};

} // namespace OP

#endif // HASCOMPARE_H
