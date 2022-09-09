/**************************************************************************

	VARIANT.H - Variant types

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
		This header implements three variant data types that can store one
		of the types in the template parameters.  These classes are a
		little different to the classic variant class, in that we always
		know in code what we are putting in, and what type we are expecting,
		and the type never changes.  This is not quite as flexible as
		some designs, but is still useful.

		The first type "variant" is a basic variant, that stores just a
		value, without much magic.  You put values in and can get them out.

		The second type, which is more complex, is a variant functor or
		"vunctor" that store the type and aa callback function, so the
		value could be the result of a functor of some form.

		The third type is a "validator" that takes a value from a vunctor
		and tests if against a function to determine if it is valid.  This
		allows one to set ranges on values and automatically test before
		use, so you can throw immediately if code attempts to use an
		invalid value.

	Design:
		This uses a rather complex nested set of variadic union templates,
		rather than the more traditional union approach, to make the
		class much more friendly to use.

		There is significant SFINAE and meta-template magic involved, which
		makes this class tricky, and a little fragile.  In addition, the
		class distiguishes between values and pointers - if you store a
		pointer in the class it will magically become a unique_ptr like
		implementation that manages the lifetime of the memory.

		To make things more comlex, if a contained type exposes a real_t
		type, and has a realize() function, the class allows one to create
		realized versions of the contained value.  And if the contained
		type has a cast_t type, then the variant will allow direct casting
		to that type as a valid form of extraction.

		The drawbacks to the design are that one must explicitly cast to
		extract a value in most cases, and that one cannot return a
		reference because the vunctor requires the extraction to run the
		callback.

		There is a lot of replication between the classes, but there does
		not seem to a be an easy way around that.  They are just too
		complex to try to merge.

	History:
		2016/07/20 - Created by Jeremy Lea <reg@openpave.org>

**************************************************************************/

#pragma once
#ifndef __VARIANT_H
#define __VARIANT_H

#include <functional>
#include <type_traits>
#include <utility>
#include "ctti.h"
#include "meta.h"

namespace OP {

// Check if the types being stored in the variant have a magic nested typedef
// for a realization type, in which case define a new variant type for the
// actual realizations of the variant.
//
// First fake std::void_t (remove when we move to C++17)
template<typename...Ts> struct make_void { typedef void type;};
template<typename...Ts> using void_t = typename make_void<Ts...>::type;
// Then set up the wrapper type to handle pointers...
template<typename T>
struct wrap_type
{
	wrap_type() noexcept : t() {
	}
	explicit wrap_type(T n) noexcept : t(n) {
	}
	wrap_type(const wrap_type &) = delete;
	wrap_type(wrap_type &&) = delete;
	wrap_type & operator = (const wrap_type &) = delete;
	wrap_type & operator = (wrap_type &&) = delete;
	template<typename U = T>
	typename std::enable_if<std::is_pointer<U>::value>::type
	set(T n) noexcept {
		delete t;
		t = n;
	}
	template<typename U = T>
	typename std::enable_if<!std::is_pointer<U>::value>::type
	set(T n) noexcept {
		t = n;
	}
	template<typename U = T>
	typename std::enable_if<std::is_pointer<U>::value>::type
	move(T && n) noexcept {
		delete t;
		t = n;
		n = nullptr;
	}
	template<typename U = T>
	typename std::enable_if<!std::is_pointer<U>::value>::type
	move(T && n) noexcept {
		t = std::move(n);
	}
	template<typename U = T>
	typename std::enable_if<std::is_pointer<U>::value>::type
	clear() noexcept {
		delete t;
		t = nullptr;
	}
	template<typename U = T>
	typename std::enable_if<!std::is_pointer<U>::value>::type
	clear() noexcept {
		t.~T();
	}
	operator const T & () const noexcept {
		return t;
	}
	operator T & () noexcept {
		return t;
	}
	T t;
};
// Now define a class to get the realized type, handling pointers.
// This is the default case with no special types.
template<typename T, typename = void_t<>>
struct real_type
{
	using real_t = typename std::remove_pointer<T>::type;
	template<typename U = T>
	static typename std::enable_if<std::is_pointer<U>::value,real_t>::type
	realize(const U & t) noexcept {
		return *t;
	}
	template<typename U = T>
	static typename std::enable_if<!std::is_pointer<U>::value,real_t>::type
	realize(const U & t) noexcept {
		return t;
	}
};
// This is the case for types with a "real" member type.
template<typename T>
struct real_type<T, void_t<typename std::remove_pointer<T>::type::real_t>>
{
	using real_t = typename std::remove_pointer<T>::type::real_t;
	template<typename U = T>
	static typename std::enable_if<std::is_pointer<U>::value,real_t>::type
	realize(const U & t) {
		return t->realize();
	}
	template<typename U = T>
	static typename std::enable_if<!std::is_pointer<U>::value,real_t>::type
	realize(const U & t) {
		return t.realize();
	}
};

// Template meta-program to compare two types to see if they are the same
// at a fundamental level (can be directly cast, etc.).
template<typename T>
struct clean_type
{
	using type = typename std::remove_cv<
		typename std::remove_reference<T>::type>::type;
};
template<typename T, typename = void_t<>>
struct cast_type
{
	using type = typename clean_type<T>::type;
};
template<typename T>
struct cast_type<T, void_t<typename clean_type<T>::type::cast_t>>
{
	using cast_t = typename clean_type<T>::type::cast_t;
	using type = typename clean_type<cast_t>::type;
};
template<typename U, typename T>
struct compare_v
{
	using setable_t = typename std::conditional<
	std::is_same<U,T>::value
	|| (std::is_base_of<typename clean_type<T>::type,
			            typename clean_type<U>::type>::value)
	,std::true_type,std::false_type>::type;
	static constexpr const bool setable = setable_t::value;
	using castable_t = typename std::conditional<
	std::is_same<U,T>::value
	|| std::is_same<typename clean_type<U>::type,
		            typename cast_type<T>::type>::value
	|| (std::is_base_of<typename clean_type<U>::type,
			            typename clean_type<T>::type>::value)
	,std::true_type,std::false_type>::type;
	static constexpr const bool castable = castable_t::value;
};

/*
 * class vunctor - A variant functor type.
 */
// Tail for the recursion.
template<typename...Ts>
class vunctor {
	struct store {
		store() {}
		store(const store &) = delete;
		store(store &&) = delete;
		store & operator = (const store &) = delete;
		store & operator = (store &&) = delete;
		template<typename U>
		U get(const OP::type_index &) const {
			throw std::runtime_error("Attempting to get invalid type from variant!");
		}
		void set(const OP::type_index &, const store &) {
			throw std::runtime_error("Attempting to store invalid type in variant!");
		}
		void set(const OP::type_index &, store &&) {
			throw std::runtime_error("Attempting to store invalid type in variant!");
		}
		void chain(const OP::type_index &, const store &) {
			throw std::runtime_error("Attempting to chain invalid type in variant!");
		}
		void clear(const OP::type_index &) noexcept {}
		// Set value from a different type of store
		template<typename...Vs>
		OP::type_index set_r(const OP::type_index &, const typename vunctor<Vs...>::store &) const {
			throw std::runtime_error("Attempting to realize invalid value type in variant!");
		}
	};
	template<typename...S>
	friend class vunctor;
	template<typename...S>
	friend class variant;
};

/*
 * class validator - A variant functor to validate a variant value.
 */
// Tail for the recursion.
template<typename...Ts>
class validator {
	struct store {
		store() = delete;
		store(const store &) = delete;
		store(store &&) = delete;
		store & operator = (const store &) = delete;
		store & operator = (store &&) = delete;
		void set(const OP::type_index &, const store &) {
			throw std::runtime_error("Attempting to store invalid type in validator!");
		}
		void set(const OP::type_index &, store &&) {
			throw std::runtime_error("Attempting to store invalid type in validator!");
		}
		void clear(const OP::type_index &) noexcept {}
		template<typename U>
		bool check_v(const OP::type_index &, const U &) const {
			throw std::runtime_error("Trying to validate incorrect type from variant!");
		}
		template<typename U>
		bool check(const OP::type_index &, const U &) const {
			throw std::runtime_error("Trying to validate incorrect type from variant!");
		}
	};
	template<typename...S>
	friend class validator;
};

/*
 * class variant - A simple variant type, without chaining.
 */
// Tail for the recursion.
template<typename...Ts>
class variant {
	struct store {
		store() {}
		store(const store &) = delete;
		store(store &&) = delete;
		store & operator = (const store &) = delete;
		store & operator = (store &&) = delete;
		void set(const OP::type_index &, const store &) {
			throw std::runtime_error("Attempting to store invalid type in variant!");
		}
		void set(const OP::type_index &, store &&) {
			throw std::runtime_error("Attempting to store invalid type in variant!");
		}
		void clear(const OP::type_index &) noexcept {}
		// Set value from a different type of store
		template<typename...Vs>
		OP::type_index set_r(const OP::type_index &, const typename vunctor<Vs...>::store &) const {
			throw std::runtime_error("Attempting to realize invalid value type in variant!");
		}
	};
	template<typename...S>
	friend class variant;
	template<typename...S>
	friend class vunctor;
};

// Actual variant functor for type T
// This is a little complex.  This class is only instatiated for the full
// template list <T,Ts...>, but creates a recursive definition only on the
// contained store.  So it is not actually recursively constructed.
template<typename T, typename...Ts>
class vunctor<T,Ts...> {
	using base_t = vunctor<Ts...>;

	OP::type_index k;
	union store {
		// A callback that can be used to get a value.
		typedef std::function<T()> callback;

		store()
		  : f(), t() {
		}
		store(const store &) = delete;
		store(store &&) = delete;
		store & operator = (const store &) = delete;
		store & operator = (store &&) = delete;
		store(const OP::type_index & d, const store & v)
		  : f(), t() {
			set(d,v);
		}
		store(const OP::type_index & d, store && v)
		  : f(), t() {
			set(d,std::move(v));
		}
		// U templates are conditioned on the union type for this slot.
		template<typename U>
		store(U u, OP::type_index * d, typename std::enable_if<
				compare_v<U,T>::setable>::type * = nullptr) noexcept
		  : f(), t(u) {
			*d = OP::type_index(OP::type_id<T>());
		}
		template<typename U>
		store(U u, OP::type_index * d, typename std::enable_if<
				!compare_v<U,T>::setable>::type * = nullptr) noexcept
		  : b(u,d) {
		}
		~store() {
			// clear() must be called outside...
		}
		// Get the function or value
		template<typename U>
		typename std::enable_if<compare_v<U,T>::castable,U>::type
		get(const OP::type_index & d) const {
			if (d != OP::type_index(OP::type_id<T>()))
				throw std::runtime_error("Attempting to get wrong type from variant!");
			return f ? f() : T(t);
		}
		template<typename U>
		typename std::enable_if<!compare_v<U,T>::castable,U>::type
		get(const OP::type_index & d) const {
			return b.template get<U>(d);
		}
		// Set from a value
		template<typename U>
		typename std::enable_if<compare_v<U,T>::setable,void>::type
		set_t(const OP::type_index & d, const U & u) {
			if (d != OP::type_index(OP::type_id<T>()))
				throw std::runtime_error("Attempting to set wrong type to variant!");
			t.set(u);
			f = nullptr;
		}
		template<typename U>
		typename std::enable_if<!compare_v<U,T>::setable,void>::type
		set_t(const OP::type_index & d, const U & u) {
			b.template set_t<U>(d,u);
		}
		// Set from a functional
		template<typename U>
		typename std::enable_if<compare_v<U,T>::setable,void>::type
		set_f(const OP::type_index & d, std::function<U()> && u) {
			if (d != OP::type_index(OP::type_id<T>()))
				throw std::runtime_error("Attempting to set wrong type of function to variant!");
			f = std::move(u);
		}
		template<typename U>
		typename std::enable_if<!compare_v<U,T>::setable,void>::type
		set_f(const OP::type_index & d, std::function<U()> && u) {
			b.template set_f<U>(d,std::move(u));
		}
		// Set directly based on another store
		void set(const OP::type_index & d, const store & v) {
			if (d == OP::type_index(OP::type_id<T>())) {
				t.set(v.get<T>(d));
				f = nullptr;
			} else
				b.set(d,v.b);
		}
		void set(const OP::type_index & d, store && v) {
			if (d == OP::type_index(OP::type_id<T>())) {
				t.move(std::move(T(v.t)));
				f = std::move(v.f);
			} else
				b.set(d,std::move(v.b));
		}
		// Chain our value to that of another variant
		void chain(const OP::type_index & d, const store & v) {
			if (d == OP::type_index(OP::type_id<T>())) {
				f = [&,d]() -> T { return v.get<T>(d); };
			} else
				b.chain(d,v.b);
		}
		// Clear if we are the tagged type, else pass
		void clear(const OP::type_index & d) noexcept {
			if (d == OP::type_index(OP::type_id<T>())) {
				t.clear();
				f.~callback();
			} else
				b.clear(d);
		}
		// Set value from a different type of store
		template<typename V, typename...Vs>
		OP::type_index set_r(const OP::type_index & d, const typename vunctor<V,Vs...>::store & v) {
			if (d == OP::type_index(OP::type_id<V>())) {
				OP::type_index i = OP::type_index(OP::type_id<typename real_type<V>::real_t>());
				this->template set_t<typename real_type<V>::real_t>(i,
					real_type<V>::realize(v.template get<V>(d)));
				return i;
			} else
				return b.template set_r<Vs...>(d,v.b);
		}
		// The rest of the union is in b
		typename base_t::store b;
		// Anonymous struct for value and functional
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable: 4201)
#endif
		struct {
			callback f;
			wrap_type<T> t;
		};
#if defined(_MSC_VER)
#pragma warning(pop)
#endif
	} s;
	// Allow all other variants to mess with us
	template<typename...S>
	friend class vunctor;
	// And variants for fixate.
	template<typename...S>
	friend class variant;
	// And for the validator to check types.
	template<typename...S>
	friend class validator;

	// Allow only private access to create an empty variant
	vunctor(const OP::type_index & d)
	  : k(d), s() {
	}
	// Convoluted templates to handle assignments from lambda functions.
	template<typename V>
	void copy_f(std::function<V()> && v) {
		s.template set_f<V>(k,std::move(v));
	}
	template<typename V>
	typename std::enable_if<!is_callable<V>::value,void>::type
	copy_v(V && v) {
		s.template set_t<V>(k,std::forward<V>(v));
	}
	template<typename V>
	typename std::enable_if<is_callable<V>::value,void>::type
	copy_v(V && v) {
		copy_f<decltype(v())>(std::move(v));
	}
	template<typename V>
	typename std::enable_if<!is_callable<V>::value,void>::type
	copy_c(const V & v) {
		s.template set_t<V>(k,v);
	}
	template<typename V>
	typename std::enable_if<is_callable<V>::value,void>::type
	copy_c(const V & v) {
		copy_f<decltype(v())>(v);
	}
	// Set value from a different type of store
	template<typename...Vs>
	void set_r(const OP::type_index & d, const vunctor<Vs...> & v) {
		 k = std::move(s.template set_r<Vs...>(d,v.s));
	}

public:
	// This is the type of variant we will be when realized.
	using real_t = vunctor<typename real_type<T>::real_t,typename real_type<Ts>::real_t...>;
	// This is the type of our validator.
	using validator_t = validator<typename cast_type<T>::type,typename cast_type<Ts>::type...>;
	// This is the type of our simple variant.
	using variant_t = variant<typename cast_type<T>::type,typename cast_type<Ts>::type...>;

	vunctor() = delete;
	vunctor(const vunctor & v)
	  : k(v.k), s(k,v.s) {
	}
	vunctor(vunctor && v)
	  : k(std::move(v.k)), s(k,std::move(v.s)) {
	}
	// Create a variant, fixing the type.
	// fake out the constructor with V.  It will actually be replaced.
	template<typename V>
	vunctor(V v) noexcept
	  : k(OP::type_index(OP::type_id<V>())), s(v,&k) {
	}
	// Destruct depending on type
	~vunctor() {
		s.clear(k);
	}
	vunctor & operator = (const vunctor & v) {
		if (k != v.k)
			throw std::runtime_error("Trying to set incorrect type into variant!");
		s.set(k,v.s);
		return *this;
	}
	vunctor & operator = (vunctor && v) {
		if (k != v.k)
			throw std::runtime_error("Trying to set incorrect type into variant!");
		s.set(k,std::move(v.s));
		return *this;
	}
	// Assign from value (of the correct type).  This cannot change the type.
	template<typename V>
	vunctor & operator = (const V & v) {
		copy_c(v);
		return *this;
	}
	template<typename V>
	vunctor & operator = (V && v) {
		copy_v(std::forward<V>(v));
		return *this;
	}
	vunctor & chain(const vunctor & v) {
		if (k != v.k)
			throw std::runtime_error("Trying to chain incorrect type with variant!");
		s.chain(k,v.s);
		return *this;
	}
	// Return the contained value.
	template<typename V>
	V get() const {
		return s.template get<V>(k);
	}
	// Also allow casting in some applications.
	template<typename V>
	operator V () const {
		return get<V>();
	}
	// Get a realized value for this variant.
	real_t realize() const {
		real_t rv(k); // k is fake here.  We needed k's real type, which set_r sets.
		rv.template set_r<T,Ts...>(k,*this);
		return rv;
	}
	// Get a fixated (non-vunctor) value for this variant.
	variant_t fixate() const;
	// Get the actual type index for this variant.
	const OP::type_index & get_type() const noexcept {
		return k;
	}
};

// Actual validator for type T
template<typename T, typename...Ts>
class validator<T,Ts...> {
	typedef validator<Ts...> base_t;

	OP::type_index k;
	union store {
		typedef std::function<bool(const typename cast_type<T>::type &)> callback;

		store() = delete;
		store(const store &) = delete;
		store(store &&) = delete;
		store & operator = (const store &) = delete;
		store & operator = (store &&) = delete;
		store(const OP::type_index & d, const store & v) {
			set(d,v);
		}
		store(const OP::type_index & d, store && v) {
			set(d,std::move(v));
		}
		// U templates are conditioned on the union type for this slot.
		template<typename U>
		store(std::function<bool(const U &)> && u,
			  OP::type_index * k,
			  typename std::enable_if<compare_v<U,T>::castable>::type * = nullptr)
		  : f(u) {
			*k = OP::type_index(OP::type_id<T>());
		}
		template<typename U>
		store(std::function<bool(const U &)> && u,
			  OP::type_index * k,
			  typename std::enable_if<!compare_v<U,T>::castable>::type * = nullptr)
		  : b(std::move(u),k) {
		}
		~store() {
			// clear() must be called outside...
		}
		// Get the function or value
		template<typename U>
		typename std::enable_if<compare_v<T,U>::castable,bool>::type
		check_v(const OP::type_index & d, const U & u) const {
			if (d != OP::type_index(OP::type_id<T>()))
				throw std::runtime_error("Attempting to check wrong type of variant!");
			return f(u);
		}
		template<typename U>
		typename std::enable_if<!compare_v<T,U>::castable,bool>::type
		check_v(const OP::type_index & d, const U & u) const {
			return b.template check_v<U>(d,u);
		}
		template<typename U>
		bool check(const OP::type_index & d, const U & u) const {
			if (d == OP::type_index(OP::type_id<T>()))
				return f(u);
			else
				return b.check(d,u);
		}
		// Set from a functional
		template<typename U>
		typename std::enable_if<compare_v<U,T>::castable,void>::type
		set_f(const OP::type_index & d, std::function<bool(const U &)> && u) {
			if (d != OP::type_index(OP::type_id<T>()))
				throw std::runtime_error("Attempting to assign wrong function type to validator!");
			f = std::move(u);
		}
		template<typename U>
		typename std::enable_if<!compare_v<U,T>::castable,void>::type
		set_f(const OP::type_index & d, std::function<bool(const U &)> && u) {
			b.template set_f<U>(d,std::move(u));
		}
		void set(const OP::type_index & d, const store & v) {
			if (d == OP::type_index(OP::type_id<T>())) {
				f = v.f;
			} else
				b.set(d,v.b);
		}
		void set(const OP::type_index & d, store && v) {
			if (d == OP::type_index(OP::type_id<T>())) {
				f = std::move(v.f);
			} else
				b.set(d,std::move(v.b));
		}
		// Clear if we are the tagged type, else pass
		void clear(const OP::type_index & d) noexcept {
			if (d == OP::type_index(OP::type_id<T>()))
				f.~callback();
			else
				b.clear(d);
		}
		// The rest of the union is in b
		typename base_t::store b;
		// The validation callback functional
		callback f;
	} s;
	// Allow all other variants to mess with us
	template<typename...S>
	friend class validator;

	// Convoluted templates to handle assignments from lambda functions.
	template<typename V>
	void copy_f(std::function<bool(const V &)> && v) {
		s.set_f(std::move(v));
	}

public:
	// Create a variant, fixing the type.
	validator() = delete;
	validator(const validator & v)
	  : k(v.k), s(k,v.s) {
	}
	validator(validator && v)
	  : k(std::move(v.k)), s(k,std::move(v.s)) {
	}
	template<typename F,
		typename = typename std::enable_if<is_callable<F>::value>::type,
		typename V = typename function_traits<F>::template arg<0>::type>
	validator(F && v)
	  : k(OP::type_index(OP::type_id<V>())),
		s(std::function<bool(const V &)>(std::forward<F>(v)),&k) {
	}
	// Destruct depending on type
	~validator() {
		s.clear(k);
	}
	validator & operator = (const validator & v) {
		if (k != v.k)
			throw std::runtime_error("Trying to set incorrect type into validator!");
		s.set(k,v.s);
		return *this;
	}
	validator & operator = (validator && v) {
		if (k != v.k)
			throw std::runtime_error("Trying to set incorrect type into validator!");
		s.set(k,std::move(v.s));
		return *this;
	}
	template<typename V>
	validator & operator = (const V & v) {
		copy_f<V>(v);
		return *this;
	}
	template<typename V>
	validator & operator = (V && v) {
		copy_f<V>(std::forward<V>(v));
		return *this;
	}
	// Validate a value of some type.
	template<typename V>
	bool validate(const V & v) const {
		return s.template check_v<V>(k,v);
	}
	bool validate(const vunctor<T,Ts...> & t) const {
		if (k != t.k)
			throw std::runtime_error("Trying to validate incorrect type!");
		return s.check(k,t);
	}
};

// Actual variant for type T
// This is like a vunctor but does not support realize, chain or validate.
template<typename T, typename...Ts>
class variant<T,Ts...> {
	using base_t = variant<Ts...>;

	OP::type_index k;
	union store {
		// A callback that can be used to get a value.
		typedef std::function<T()> callback;

		store() {
		}
		store(const store &) = delete;
		store(store &&) = delete;
		store & operator = (const store &) = delete;
		store & operator = (store &&) = delete;
		store(const OP::type_index & d, const store & v) {
			set(d,v);
		}
		store(const OP::type_index & d, store && v) {
			set(d,std::move(v));
		}
		template<typename U>
		store(U u, OP::type_index * k, typename std::enable_if<
				compare_v<U,T>::setable>::type * = 0)
		  : t(u) {
			*k = OP::type_index(OP::type_id<T>());
		}
		template<typename U>
		store(U u, OP::type_index * k, typename std::enable_if<
				!compare_v<U,T>::setable>::type * = 0)
		  : b(u,k) {
		}
		~store() {
			// clear() must be called outside...
		}
		// Get the function or value
		template<typename U>
		typename std::enable_if<compare_v<U,T>::castable,U>::type
		get(const OP::type_index & d) const {
			if (d != OP::type_index(OP::type_id<T>()))
				throw std::runtime_error("Attempting to get wrong type from variant!");
			return T(t);
		}
		template<typename U>
		typename std::enable_if<!compare_v<U,T>::castable,U>::type
		get(const OP::type_index & d) const {
			return b.template get<U>(d);
		}
		// Set from a value
		template<typename U>
		typename std::enable_if<compare_v<U,T>::setable,void>::type
		set_t(const OP::type_index & d, const U & u) {
			if (d != OP::type_index(OP::type_id<T>()))
				throw std::runtime_error("Attempting to set wrong type to variant!");
			t.set(u);
		}
		template<typename U>
		typename std::enable_if<!compare_v<U,T>::setable,void>::type
		set_t(const OP::type_index & d, const U & u) {
			b.template set_t<U>(d,u);
		}
		// Set directly based on another store
		void set(const OP::type_index & d, const store & v) {
			if (d == OP::type_index(OP::type_id<T>()))
				t.set(v.get<T>(d));
			else
				b.set(d,v.b);
		}
		void set(const OP::type_index & d, store && v) {
			if (d == OP::type_index(OP::type_id<T>()))
				t.move(std::move(T(v.t)));
			else
				b.set(d,std::move(v.b));
		}
		// Clear if we are the tagged type, else pass
		void clear(const OP::type_index & d) noexcept {
			if (d == OP::type_index(OP::type_id<T>()))
				t.clear();
			else
				b.clear(d);
		}
		// Set value from a different type of store
		template<typename V, typename...Vs>
		OP::type_index set_r(const OP::type_index & d, const typename vunctor<V,Vs...>::store & v) {
			if (d == OP::type_index(OP::type_id<V>())) {
				OP::type_index i = OP::type_index(OP::type_id<typename cast_type<V>::type>());
				this->template set_t<typename cast_type<V>::type>(i,
					v.template get<V>(d));
				return i;
			} else
				return b.template set_r<Vs...>(d,v.b);
		}
		// The rest of the union is in b
		typename base_t::store b;
		wrap_type<T> t;
	} s;
	// Allow all other variants to mess with us
	template<typename...S>
	friend class variant;
	// And also allow vunctors.
	template<typename...S>
	friend class vunctor;

	// Allow only private access to create an empty variant
	variant(const OP::type_index & d)
	  : k(d), s() {
	}
	// Set value from a different type of store
	template<typename...Vs>
	void set_r(const OP::type_index & d, const vunctor<Vs...> & v) {
		k = std::move(s.template set_r<Vs...>(d,v.s));
	}

public:
	variant() = delete;
	variant(const vunctor<T,Ts...> & v)
	  : k(v.k), s(k,v.s) {
	}
	variant(const variant & v)
	  : k(v.k), s(k,v.s) {
	}
	variant(variant && v)
	  : k(std::move(v.k)), s(k,std::move(v.s)) {
	}
	template<typename V>
	variant(V v)
	  : k(OP::type_index(OP::type_id<V>())), s(v,&k) {
	}
	~variant() {
		s.clear(k);
	}
	variant & operator = (const vunctor<T,Ts...> & v) {
		if (k != v.k)
			throw std::runtime_error("Trying to set incorrect type into variant!");
		s.set(k,v.s);
		return *this;
	}
	variant & operator = (const variant & v) {
		if (k != v.k)
			throw std::runtime_error("Trying to set incorrect type into variant!");
		s.set(v.k,v.s);
		return *this;
	}
	variant & operator = (variant && v) {
		if (k != v.k)
			throw std::runtime_error("Trying to set incorrect type into variant!");
		s.set(std::move(v.k),std::move(v.s));
		return *this;
	}
	template<typename V>
	variant & operator = (const V & v) {
		s.template set_t<V>(k,v);
		return *this;
	}
	template<typename V>
	variant & operator = (V && v) {
		s.template set_t<V>(k,std::forward<V>(v));
		return *this;
	}
	template<typename V>
	V get() const {
		return s.template get<V>(k);
	}
	template<typename V>
	operator V () const {
		return get<V>();
	}
	const OP::type_index & get_type() const noexcept {
		return k;
	}
};

template<typename T, typename...Ts>
inline typename vunctor<T,Ts...>::variant_t
vunctor<T,Ts...>::fixate() const {
	variant_t rv(k); // k is fake here.
	rv.template set_r<T,Ts...>(k,*this);
	return rv;
}

} // namespace OP

#endif // VARIANT_H
