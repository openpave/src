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

/*
 * class variant - A variant type.
 */
// Tail for the recursion.
template<typename ...Ts>
class variant {
	struct store {
		store() = delete;
		store(const store &) = delete;
		store(store &&) = delete;
		store & operator= (const store &) = delete;
		store & operator= (store &&) = delete;
		void set(std::type_index, const store &) {};
		void set(std::type_index, store &&) {};
		void clear(std::type_index) {};
	};
	template<typename ...S>
	friend class variant;
};

// Actual variant for type T
template<typename T, typename ...Ts>
class variant<T,Ts...> {
	typedef variant<Ts...> base_t;

	std::type_index k;
	union store {
		typedef std::function<T()> callback;

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
		store(U u, typename std::enable_if<std::is_same<U,T>::value>::type * = 0) :
			t(u), f() {
		}
		template<typename U>
		store(U u, typename std::enable_if<!std::is_same<U,T>::value>::type * = 0) :
			b(u) {
		}
		~store() {
			// clear() must be called outside...
		}
		// Get the function or value
		template<typename U>
		typename std::enable_if<std::is_same<U,T>::value,U>::type
		get() const {
			return (f ? f() : t);
		}
		template<typename U>
		typename std::enable_if<!std::is_same<U,T>::value,U>::type
		get() const {
			return b.template get<U>();
		}
		// Set from a value
		template<typename U>
		typename std::enable_if<std::is_same<U,T>::value,void>::type
		set_t(const U & u) {
			t = u;
			f = nullptr;
		}
		template<typename U>
		typename std::enable_if<!std::is_same<U,T>::value,void>::type
		set_t(const U & u) {
			b.template set_t<U>(u);
		}
		// Set from a functional
		template<typename U>
		typename std::enable_if<std::is_same<U,T>::value,void>::type
		set_f(std::function<U()> && u) {
			f = std::move(u);
		}
		template<typename U>
		typename std::enable_if<!std::is_same<U,T>::value,void>::type
		set_f(std::function<U()> && u) {
			b.template set_f<U>(std::move(u));
		}
		// Set directly based on another store
		void set(std::type_index d, const store & v) {
			if (d == std::type_index(typeid(T))) {
				t = v.t;
				f = v.f;
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
				f = [&]() { return v.f(); };
			} else
				b.set(d,v.b);
		}
		// Clear if we are the tagged type, else pass
		void clear(std::type_index d) {
			if (d == std::type_index(typeid(T))) {
				t.~T();
				f.~callback();
			} else
				b.clear(d);
		}
		// The rest of the union is in b
		typename base_t::store b;
		// Anonymous struct for value and functional
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable: 4201)
#endif
		struct {
			T t;
			callback f;
		};
#if defined(_MSC_VER)
#pragma warning(pop)
#endif
	} s;
	// Allow all other variants to mess with us
	template<typename ...S>
	friend class variant;

	// Convoluted templates to handle assignments from lambda functions.
	template<typename V>
    void copy_f(std::function<V()> && v) {
		if (k != std::type_index(typeid(V)))
			throw std::runtime_error("Trying to set incorrect type into variant!");
		s.set_f(std::move(v));
    }
	template<typename V>
	typename std::enable_if<!is_callable<V>::value,void>::type
	copy_v(V && v) {
		if (k != std::type_index(typeid(V)))
			throw std::runtime_error("Trying to set incorrect type into variant!");
		s.set_t(std::move(v));
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
		s.set_t(v);
    }
	template<typename V>
	typename std::enable_if<is_callable<V>::value,void>::type
	copy_c(const V & v) {
		copy_f<decltype(v())>(v);
    }

public:
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
		copy_v(std::move(v));
		return *this;
	}
	// Return the contained value.
	template<typename V>
	operator V () const {
		if (k != std::type_index(typeid(V)))
			throw std::runtime_error("Trying to get incorrect type from variant!");
		return s.template get<V>();
	}
	variant & chain(const variant & v) {
		if (k != v.k)
			throw std::runtime_error("Trying to chain incorrect type with variant!");
		s.chain(k,v.s);
		return *this;
	}
	std::type_index get_type() const {
		return k;
	}
};

/*
 * class validator - A variant functor to validate a variant value.
 */
// Tail for the recursion.
template<typename ...Ts>
class validator {
	struct store {
		store() = delete;
		store(const store &) = delete;
		store(store &&) = delete;
		store & operator= (const store &) = delete;
		store & operator= (store &&) = delete;
		void set(std::type_index, const store &) {}
		void set(std::type_index, store &&) {}
		void clear(std::type_index) {}
		template<typename U>
		bool check(std::type_index, const U &) const {
			throw std::runtime_error("Trying to validate incorrect type from variant!");
		}
	};
	template<typename ...S>
	friend class validator;
};

// Actual validator for type T
template<typename T, typename ...Ts>
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
			  typename std::enable_if<std::is_same<U,T>::value>::type * = 0) :
			f(u) {
		}
		template<typename U>
		store(std::function<bool(const U &)> && u,
			  typename std::enable_if<!std::is_same<U,T>::value>::type * = 0) :
			b(std::move(u)) {
		}
		~store() {
			// clear() must be called outside...
		}
		// Get the function or value
		template<typename U>
		typename std::enable_if<std::is_same<U,T>::value,bool>::type
		check(const U & u) const {
			return f(u);
		}
		template<typename U>
		typename std::enable_if<!std::is_same<U,T>::value,bool>::type
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
		typename std::enable_if<std::is_same<U,T>::value,void>::type
		set_f(std::function<bool(const U &)> && u) {
			f = std::move(u);
		}
		template<typename U>
		typename std::enable_if<!std::is_same<U,T>::value,void>::type
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
	template<typename ...S>
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
		s(std::function<bool(const V &)>(std::move(v))) {
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
		copy_f<V>(std::move(v));
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
