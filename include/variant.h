/**************************************************************************

	VARIANT.H - A variant type

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
		This header implements a variant data type that can store one of
		the types in the template parameters.  This class is a little
		different to the classic variant class, in that we always know
		in code what we are putting in, and what type we are expecting,
		and the type never changes.  This is not quite as flexible as
		some designs, but is still useful.

		In addition, the type stores a callback function, so the value
		could be the result of a functor of some form.

	Design:
		This uses a rather complex nested set of variadic templates,
		rather than the more traditional union approach, to make the
		class much more friendly to use.

	History:
		2016/07/20 - Created by Jeremy Lea <reg@openpave.org>

**************************************************************************/

#ifndef __VARIANT_H
#define __VARIANT_H

#include <functional>
#include <typeindex>
#include <type_traits>
#include <utility>
#include "hascompare.h"

namespace OP {

namespace {
// Check if the types being stored in the variant have a magic nested typedef
// for a realization type, in which case define a new variant type for the
// actual realizations of the variant.
//
// First fake std::void_t
template<typename...Ts> struct make_void { typedef void type;};
template<typename...Ts> using void_t = typename make_void<Ts...>::type;
// Then set up the cleaner...
template<typename T>
struct clear_type
{
	template<typename U = T>
	static typename std::enable_if<std::is_pointer<U>::value,void>::type
	clear(U & t) {
		delete t;
		t = nullptr;
	}
	template<typename U = T>
	static typename std::enable_if<!std::is_pointer<U>::value,void>::type
	clear(U & t) {
		t.~U();
	}
};
// This is the default case with no special types.
template<typename T, typename = void_t<>>
struct real_type
{
	using real_t = typename std::remove_pointer<T>::type;
	template<typename U = T>
	static typename std::enable_if<std::is_pointer<U>::value,real_t>::type
	value(const U & t) {
		return *t;
	}
	template<typename U = T>
	static typename std::enable_if<!std::is_pointer<U>::value,real_t>::type
	value(const U & t) {
		return t;
	}
};
// This is the case for types with a "real" member type.
template<typename T>
struct real_type<T, void_t<typename std::remove_pointer<T>::type::real>>
{
	using real_t = typename std::remove_pointer<T>::type::real;
	template<typename U = T>
	static typename std::enable_if<std::is_pointer<U>::value,real_t>::type
	value(const U & t) {
		return t->realize();
	}
	template<typename U = T>
	static typename std::enable_if<!std::is_pointer<U>::value,real_t>::type
	value(const U & t) {
		return t.realize();
	}
};

// Template meta-program to compare two types to see if they are the same
// at a fundamental level (can be directly cast, etc.).
template<typename T, typename = void_t<>>
struct unwrap_type
{
	using type = typename std::remove_reference<T>::type;
};
template<typename T>
struct unwrap_type<T, void_t<typename T::type>>
{
	using type = typename std::remove_reference<typename T::type>::type;
};
template<typename U, typename T>
struct compare_v
{
	using type = typename std::conditional<
	std::is_same<U,T>::value
	|| std::is_same<typename std::remove_cv<U>::type,
		            typename unwrap_type<T>::type>::value
	|| (std::is_base_of<typename unwrap_type<T>::type,
			            typename std::remove_cv<U>::type>::value)
	,std::true_type,std::false_type>::type;
	static constexpr bool value = type::value;
};

}

/*
 * class variant - A variant type.
 */
// Tail for the recursion.
template<typename...Ts>
class variant {
	struct store {
		store() {}
		store(const store &) = delete;
		store(store &&) = delete;
		store & operator= (const store &) = delete;
		store & operator= (store &&) = delete;
		void set(std::type_index, const store &) {
			throw std::runtime_error("Attempting to store invalid type in variant!");
		};
		void set(std::type_index, store &&) {
			throw std::runtime_error("Attempting to store invalid type in variant!");
		};
		void clear(std::type_index) {};
		// Set value from a different type of store
		template<typename...Vs>
		std::type_index set_e(const std::type_index, const typename variant<Vs...>::store &) const {
			throw std::runtime_error("Attempting to store invalid value type in variant!");
		}
	};
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
		store & operator= (const store &) = delete;
		store & operator= (store &&) = delete;
		void set(std::type_index, const store &) {
			throw std::runtime_error("Attempting to store invalid type in validator!");
		}
		void set(std::type_index, store &&) {
			throw std::runtime_error("Attempting to store invalid type in validator!");
		}
		void clear(std::type_index) {}
		template<typename U>
		bool check(std::type_index, const U &) const {
			throw std::runtime_error("Trying to validate incorrect type from variant!");
		}
	};
	template<typename...S>
	friend class validator;
};

// Actual variant for type T
template<typename T, typename...Ts>
class variant<T,Ts...> {
	using base_t = variant<Ts...>;

	std::type_index k;
	union store {
		typedef std::function<T()> callback;

		store() :
			f() {
		};
		store(const store &) = delete;
		store(store &&) = delete;
		store & operator= (const store &) = delete;
		store & operator= (store &&) = delete;
		store(std::type_index d, const store & v) :
			f() {
			set(d,v);
		}
		store(std::type_index d, store && v)  :
			f() {
			set(d,std::move(v));
		}
		// U templates are conditioned on the union type for this slot.
		template<typename U>
		store(U u, typename std::enable_if<compare_v<U,T>::value>::type * = 0) :
			f(), t(u) {
		}
		template<typename U>
		store(U u, typename std::enable_if<!compare_v<U,T>::value>::type * = 0) :
			b(u) {
		}
		~store() {
			// clear() must be called outside...
		}
		// Get the function or value
		template<typename U>
		typename std::enable_if<compare_v<U,T>::value,typename real_type<U>::real_t>::type
		get() const {
			return real_type<U>::value(f ? t = f() : t);
		}
		template<typename U>
		typename std::enable_if<!compare_v<U,T>::value,typename real_type<U>::real_t>::type
		get() const {
			return b.template get<U>();
		}
		// Set from a value
		template<typename U>
		typename std::enable_if<compare_v<U,T>::value,void>::type
		set_t(const U & u) {
			t = u;
			f = nullptr;
		}
		template<typename U>
		typename std::enable_if<!compare_v<U,T>::value,void>::type
		set_t(const U & u) {
			b.template set_t<U>(u);
		}
		// Set from a functional
		template<typename U>
		typename std::enable_if<compare_v<U,T>::value,void>::type
		set_f(std::function<U()> && u) {
			f = std::move(u);
		}
		template<typename U>
		typename std::enable_if<!compare_v<U,T>::value,void>::type
		set_f(std::function<U()> && u) {
			b.template set_f<U>(std::move(u));
		}
		// Set directly based on another store
		void set(std::type_index d, const store & v) {
			if (d == std::type_index(typeid(T))) {
				t = v.get<T>();
				f = nullptr;
			} else
				b.set(d,v.b);
		}
		void set(std::type_index d, store && v) {
			if (d == std::type_index(typeid(T))) {
				t = std::move(v.t);
				f = std::move(v.f);
			} else
				b.set(d,std::move(v.b));
		}
		// Chain our value to that of another variant
		void chain(std::type_index d, const store & v) {
			if (d == std::type_index(typeid(T))) {
				f = [&]() { return v.get<T>(); };
			} else
				b.set(d,v.b);
		}
		// Clear if we are the tagged type, else pass
		void clear(std::type_index d) {
			if (d == std::type_index(typeid(T))) {
				clear_type<T>::clear(t);
				f.~callback();
			} else
				b.clear(d);
		}
		// Set value from a different type of store
		template<typename V, typename...Vs>
		std::type_index set_e(const std::type_index d, const typename variant<V,Vs...>::store & v) {
			if (d == std::type_index(typeid(V))) {
				this->template set_t<typename real_type<V>::real_t>(
					real_type<V>::value(v.t));
				return std::type_index(typeid(typename real_type<V>::real_t));
			} else
				return b.template set_e<Vs...>(d,v.b);
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
			mutable T t;
		};
#if defined(_MSC_VER)
#pragma warning(pop)
#endif
	} s;
	// Allow all other variants to mess with us
	template<typename...S>
	friend class variant;

	// Allow only private access to create an empty variant
	variant(const std::type_index d) :
		k(d), s() {
	}
	// Convoluted templates to handle assignments from lambda functions.
	template<typename V>
    void copy_f(std::function<V()> && v) {
		if (k != std::type_index(typeid(V)))
			throw std::runtime_error("Trying to set incorrect type into variant!");
		s.template set_f<V>(std::move(v));
    }
	template<typename V>
	typename std::enable_if<!is_callable<V>::value,void>::type
	copy_v(V && v) {
		if (k != std::type_index(typeid(V)))
			throw std::runtime_error("Trying to set incorrect type into variant!");
		s.template set_t<V>(std::forward<V>(v));
    }
	template<typename V>
	typename std::enable_if<is_callable<V>::value,void>::type
	copy_v(V && v) {
		copy_f<decltype(v())>(std::move(v));
    }
	template<typename V>
	typename std::enable_if<!is_callable<V>::value,void>::type
	copy_c(const V & v) {
		if (k != std::type_index(typeid(V)))
			throw std::runtime_error("Trying to set incorrect type into variant!");
		s.template set_t<V>(v);
    }
	template<typename V>
	typename std::enable_if<is_callable<V>::value,void>::type
	copy_c(const V & v) {
		copy_f<decltype(v())>(v);
    }
	// Set value from a different type of store
	template<typename...Vs>
	void set_e(const std::type_index d, const variant<Vs...> & v) {
		 k = s.template set_e<Vs...>(d,v.s);
	}

public:
	// This is the type of variant we will be when realized.
	using real_t = variant<typename real_type<T>::real_t,typename real_type<Ts>::real_t...>;
	// This is the type of out validator.
	using validator_t = validator<T,Ts...>;

	variant() = delete;
	variant(const variant & v) :
		k(v.k), s(k,v.s) {
	}
	variant(variant && v) :
		k(std::move(v.k)), s(k,std::move(v.s)) {
	}
	// Create a variant, fixing the type.
	template<typename V>
	variant(V v) :
		k(std::type_index(typeid(V))), s(v) {
	}
	// Destruct depending on type
	~variant() {
		s.clear(k);
	}
	variant & operator = (const variant & v) {
		if (k != v.k)
			throw std::runtime_error("Trying to set incorrect type into variant!");
		s.set(k,v.s);
		return *this;
	}
	variant & operator = (variant && v) {
		if (k != v.k)
			throw std::runtime_error("Trying to set incorrect type into variant!");
		s.set(k,std::move(v.s));
		return *this;
	}
	// Assign from value (of the correct type).  This cannot change the type.
	template<typename V>
	variant & operator = (const V & v) {
		copy_c(v);
		return *this;
	}
	template<typename V>
	variant & operator = (V && v) {
		copy_v(std::forward<V>(v));
		return *this;
	}
	variant & chain(const variant & v) {
		if (k != v.k)
			throw std::runtime_error("Trying to chain incorrect type with variant!");
		s.chain(k,v.s);
		return *this;
	}
	// Return the contained value.
	template<typename V>
	operator V () const {
		if (k != std::type_index(typeid(V)))
			throw std::runtime_error("Trying to get incorrect type from variant!");
		return s.template get<V>();
	}
	// Get a realized value for this variant.
	real_t realize() const {
		real_t rv(k);
		rv.template set_e<T,Ts...>(k,*this);
		return rv;
	}
	std::type_index get_type() const {
		return k;
	}
};

// Actual validator for type T
template<typename T, typename...Ts>
class validator<T,Ts...> {
	typedef validator<Ts...> base_t;

	std::type_index k;
	union store {
		typedef std::function<bool(const T &)> callback;

		store() = delete;
		store(const store &) = delete;
		store(store &&) = delete;
		store & operator= (const store &) = delete;
		store & operator= (store &&) = delete;
		store(std::type_index d, const store & v) {
			set(d,v);
		}
		store(std::type_index d, store && v) {
			set(d,std::move(v));
		}
		// U templates are conditioned on the union type for this slot.
		template<typename U>
		store(std::function<bool(const U &)> && u,
			  typename std::enable_if<compare_v<U,T>::value>::type * = 0) :
			f(u) {
		}
		template<typename U>
		store(std::function<bool(const U &)> && u,
			  typename std::enable_if<!compare_v<U,T>::value>::type * = 0) :
			b(std::move(u)) {
		}
		~store() {
			// clear() must be called outside...
		}
		// Get the function or value
		template<typename U>
		typename std::enable_if<compare_v<U,T>::value,bool>::type
		check(const U & u) const {
			return f(u);
		}
		template<typename U>
		typename std::enable_if<!compare_v<U,T>::value,bool>::type
		check(const U & u) const {
			return b.template check<U>(u);
		}
		template<typename U>
		bool check(std::type_index d, const U & u) const {
			if (d == std::type_index(typeid(T)))
				return f(static_cast<T>(u));
			else
				return b.check(d,u);
		}
		// Set from a functional
		template<typename U>
		typename std::enable_if<compare_v<U,T>::value,void>::type
		set_f(std::function<bool(const U &)> && u) {
			f = std::move(u);
		}
		template<typename U>
		typename std::enable_if<!compare_v<U,T>::value,void>::type
		set_f(std::function<bool(const U &)> && u) {
			b.template set_f<U>(std::move(u));
		}
		void set(std::type_index d, const store & v) {
			if (d == std::type_index(typeid(T))) {
				f = v.f;
			} else
				b.set(d,v.b);
		}
		void set(std::type_index d, store && v) {
			if (d == std::type_index(typeid(T))) {
				f = std::move(v.f);
			} else
				b.set(d,std::move(v.b));
		}
		// Clear if we are the tagged type, else pass
		void clear(std::type_index d) {
			if (d == std::type_index(typeid(T)))
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
		if (k != std::type_index(typeid(V)))
			throw std::runtime_error("Trying to set incorrect type into validator!");
		s.set_f(std::move(v));
    }

public:
	// Create a variant, fixing the type.
	validator() = delete;
	validator(const validator & v) :
		k(v.k), s(k,v.s) {
	}
	validator(validator && v) :
		k(std::move(v.k)), s(k,std::move(v.s)) {
	}
	template<typename F,
		typename = typename std::enable_if<is_callable<F>::value>::type,
		typename V = typename function_traits<F>::template arg<0>::type>
	validator(F && v) :
		k(std::type_index(typeid(V))),
		s(std::function<bool(const V &)>(std::forward<F>(v))) {
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
	// Return the contained value.
	template<typename V>
	bool validate(const V & v) const {
		if (k != std::type_index(typeid(V)))
			throw std::runtime_error("Trying to validate incorrect type!");
		return s.template check<V>(v);
	}
	bool validate(const variant<T,Ts...> & t) const {
		return s.check(k,t);
	}
};

} // namespace OP

#endif // VARIANT_H
