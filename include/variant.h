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

namespace OP {

namespace {

// Replace once C++17 version is available.
template<typename T>
class is_callable
{
	template<typename C>
	static auto check(C * c) -> decltype((*c)(),void(),std::true_type());
	template<typename>
	static auto check(...) -> decltype(std::false_type());
public:
	static constexpr bool value = decltype(check<T>(nullptr))::value;
};

}
            
/*
 * class variant - A variant type.
 *
 */
// Tail for the recursion.
template<typename ...Ts>
class variant {
	struct store {
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
		 
		// U templates are conditioned on the union type for this slot.
		template<typename U>
		store(U u, typename std::enable_if<std::is_same<U,T>::value>::type * = 0)
		  : t(u), f() {
		}
		template<typename U>
		store(U u, typename std::enable_if<!std::is_same<U,T>::value>::type * = 0)
		  : b(u) {
		}
		~store() {
			// clear() must be called outside...
		}
		// Get the function or value
		template<typename U>
		typename std::enable_if<std::is_same<U,T>::value,U>::type
		get() {
			return (f ? f() : t);
		}
		template<typename U>
		typename std::enable_if<!std::is_same<U,T>::value,U>::type
		get() {
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
		// Clear if we are the tagged type, else pass
		void clear(std::type_index d) {
			if (d == std::type_index(typeid(T))) {
				t.~T();
				if (f)
					f.~callback();
			} else
				b.clear(d);
		}
		// The rest of the union is in b
		typename base_t::store b;
		// Anonymous struct for value and functional
		struct {
			T t;
			callback f;
		};
	} s;
	// Allow all other variants to mess with us
	template<typename ...S>
	friend class variant;
	
	// Convoluted templates to handle assignments from lambda functions.
	template<typename V>
    void copy_f(std::function<V()> && v) {
		if (k != std::type_index(typeid(V)))
			throw std::runtime_error("Trying to set incorrect type into variant!");
		s.template set_f(std::move(v));
    }
	template<typename V>
	typename std::enable_if<!is_callable<V>::value,void>::type
	copy_v(V && v) {
		if (k != std::type_index(typeid(V)))
			throw std::runtime_error("Trying to set incorrect type into variant!");
		s.template set_t(std::move(v));
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
		s.template set_t(v);
    }
	template<typename V>
	typename std::enable_if<is_callable<V>::value,void>::type
	copy_c(const V & v) {
		copy_f<decltype(v())>(v);
    }

public:
	// Create a variant, fixing the type.
	template<typename V>
	variant(V v)
	  : k(std::type_index(typeid(V))), s(v) {
	}
	// Destruct depending on type
	~variant() {
		s.clear(k);
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
	operator V () {
		if (k != std::type_index(typeid(V)))
			throw std::runtime_error("Trying to get incorrect type from variant!");
		return s.template get<V>();
	}
	// Allow this to act as a function
	variant & operator () () {
		return *this;
	}
};

} // namespace OP

#endif // VARIANT_H
