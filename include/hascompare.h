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

namespace OP {

#define HAS_MEMBER_FUNCTION(NAME,RET, ...)                                \
template<typename T>                                                      \
class has_##NAME                                                          \
{                                                                         \
	template<typename C>                                                  \
	static typename std::is_same<                                         \
		decltype(std::declval<C>().NAME(__VA_ARGS__)),RET>::type          \
		check(C *) { return std::true_type(); }                           \
	template<typename>                                                    \
	static std::false_type check(...) { return std::false_type(); }       \
	typedef decltype(check<T>(0)) type;                                   \
public:                                                                   \
	static constexpr bool value = type::value;                            \
	~has_##NAME() {}  /* Keep GCC happy */                                \
};

#define HAS_NESTED_CLASS(NAME)                                            \
template<typename T>                                                      \
class has_##NAME                                                          \
{                                                                         \
	template<typename C>                                                  \
	static typename std::is_class<typename C::NAME>::type                 \
		check(C *) { return std::true_type(); }                           \
	template<typename>                                                    \
	static std::false_type check(...) { return std::false_type(); }       \
	typedef decltype(check<T>(0)) type;                                   \
public:                                                                   \
	static constexpr bool value = type::value;                            \
	~has_##NAME() {}  /* Keep GCC happy */                                \
};

HAS_MEMBER_FUNCTION(compare,int,std::declval<C>())

// Replace once C++17 version is available.
template<typename F>
class is_callable
{
	using T = typename std::remove_cv<
		typename std::remove_reference<F>::type>::type;

	template<typename C>
	static auto check1(C * c) -> decltype((*c)(),void(),std::true_type())
	{ return std::true_type(); }
	template<typename>
	static auto check1(...) -> decltype(std::false_type())
	{ return std::false_type(); }
	template<typename C>
	static auto check2(C *) -> decltype(&C::operator(),void(),std::true_type())
	{ return std::true_type(); }
	template<typename>
	static auto check2(...) -> decltype(std::false_type())
	{ return std::false_type(); }
public:
	static constexpr bool value = decltype(check1<T>(nullptr))::value
	                           || decltype(check2<T>(nullptr))::value;
};

// We really only need arg<0> from this...
template<typename>
struct function_traits;
template<typename F>
struct function_traits :
	public function_traits<decltype(&F::operator())>
{};
template<typename C, typename R, typename...Ts>
struct function_traits<R(C::*)(Ts...) const>
{
	static constexpr std::size_t arg_count = sizeof...(Ts);
	using result_type = R;
	template<std::size_t N>
	struct arg {
		static_assert(N < arg_count,"error: index exceeds arguement count");
		using type = typename std::tuple_element<N,std::tuple<Ts...>>::type;
	};
};
/*
template<typename R, typename...Ts>
struct function_traits<R(*)(Ts...)> :
	public function_traits<R(Ts...)>
{};
// member function pointer
template<typename C, typename R, typename...Ts>
struct function_traits<R(C::*)(Ts...)> :
	public function_traits<R(C&,Ts...)>
{};
// const member function pointer
template<typename C, typename R, typename...Ts>
struct function_traits<R(C::*)(Ts...) const> :
	public function_traits<R(C&,Ts...)>
{};
template<typename F>
struct function_traits<std::function<F>> :
	public function_traits<F>
{};
// member object pointer
template<typename C, typename R>
struct function_traits<R(C::*)> :
	public function_traits<R(C&)>
{};
*/
template<typename F> struct function_traits<F&> :
	public function_traits<F>
{};
template<typename F>
struct function_traits<const F&> :
	public function_traits<F>
{};
template<typename F>
struct function_traits<F&&> :
	public function_traits<F>
{};

} // namespace OP

#endif // HASCOMPARE_H
