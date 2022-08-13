/**************************************************************************

	LIST.H - Simple list template classes

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
		This header implements a number of templated C++ classes for
		various types of lists.

	Design:
		There are three types of list elements, and obviously the lists
		require that you use the right types.  There are singly linked,
		doubly linked and owned elements.  Owned elements are the only
		strange ones, since they require that you define the type of
		class that can own this list type.  This is very useful for a
		number of engineering problems, where you want the list class
		to be able to do things with the logical unit above, plus you can
		have lists which manage themselves entirely in their constructors
		and destructors.

	Status:
		At the moment this only does owned doubly linked lists, since
		that's all I need.

	History:
		2002/01/23 - Created by Jeremy Lea <reg@openpave.org>

**************************************************************************/

#ifndef __LIST_H
#define __LIST_H

#include <stdexcept>

namespace OP {

// Forward declare some classes...
template<typename T>
class list_single;
template<typename T>
class list_double;
template<typename O, typename T>
class list_owned;

/*
 * class listelement_s - Templated linked list element
 *
 * Use this as a base class for the data you want to store in the list.
 */
template<typename T>
class listelement_s {
protected:
	T * next;                   // The next element.

	listelement_s(T * n = nullptr)
	  : next(n) {
	}
	listelement_s(const listelement_s &) = delete;
	listelement_s & operator = (const listelement_s &) = delete;
	listelement_s(listelement_s &&) = delete;
	listelement_s & operator = (listelement_s &&) = delete;
	~listelement_s() {
	}

	friend class list_single<T>;
};

/*
 * class list_single - Templated list management class
 *
 * This class manages a singally linked list of type listelement_s<T>.
 */
template<typename T>
class list_single {
protected:
	T * next;                   // The head of the list.

	// All lists start empty...
	list_single()
	  : next(nullptr) {
	}
	list_single(const list_single &) = delete;
	list_single & operator = (const list_single &) = delete;
	// We just need to move the head to steal the chain.
	list_single(list_single && e)
	  : next(e.next) {
		e.next = nullptr;
	}
	list_single & operator = (list_single && e) {
		next = e.next;
		e.next = nullptr;
		return *this;
	}
	~list_single() {
		empty();
	}

	// Insert an element at the end.
	T * insert(T * e) {
		if (e == nullptr)
			return e;          // Pass-through allocation failure
		T ** t = &next, ** l = &next;
		while (*t != nullptr) {
			if (e->next != nullptr && e->next == *t) {
				(*l)->next = e;
				return e;
			} else
				l = t, t = &((*t)->next);
		};
		return *t = e;
	}
	// Remove an element from the list (but don't delete it)
	T * remove(T * e = nullptr) {
		T ** t = &next, ** l = &next;
		while (*t != nullptr) {
			if (e != nullptr && *t == e) {
				*l = e->next;
				e->next = nullptr;
				return e;
			} else
				l = t, t = &((*t)->next);
		}
		return *l;
	}
	// Push a new element as the head
	T * push(T * e) {
		if (e == nullptr)
			return e;          // Pass-through allocation failure
		T ** t = &(e->next);
		while (*t != nullptr)
			t = &((*t)->next);
		*t = next, next = e;
		return e;
	}
	// Remove the head element (but don't delete it)
	T * pop() {
		T * e = next;
		next = e->next;
		e->next = nullptr;
		return e;
	}
	// Check in the list is empty.
	bool isempty() const {
		return (next == nullptr ? true : false);
	}
	// Figure out the length of the list.
	unsigned length() const {
		unsigned s = 0;
		T * t = next;
		while (t != nullptr)
			s++, t = t->next;
		return s;
	}
	// Sometimes you just need to start a new list.
	void empty() {
		while (!isempty())
			delete pop();
	}
};

/*
 * class listelement_d - Templated doubly linked list element.
 *
 * Use this class as a base class for the data you want to store in the
 * list.
 */
template<typename T>
class listelement_d {
protected:
	T * next;                   // The next element.
	T * prev;                   // The previous element.

	// Create a new list element, with an optional previous element.
	// We force the consumer to use nullptr so that they have to think.
	listelement_d(T * p, T * n = nullptr)
	  : next(n), prev(p) {
		if (next != nullptr) {
			if (next->prev != nullptr) {
				if (prev != nullptr && prev != next->prev)
					throw std::runtime_error("Linked list pointer mismatch!");
				prev = next->prev;
				prev->next = static_cast<T *>(this);
			}
			next->prev = static_cast<T *>(this);
		} else if (prev != nullptr) {
			if (prev->next != nullptr) {
				next = prev->next;
				next->prev = static_cast<T *>(this);
			}
			prev->next = static_cast<T *>(this);
		}
	}
	// Copying or moving list elements is inherently unsafe.
	listelement_d(const listelement_d &) = delete;
	listelement_d & operator = (const listelement_d &) = delete;
	listelement_d(listelement_d && e) = delete;
	listelement_d & operator = (listelement_d && e) = delete;
	// Unlink ourselves from the list before we die...
	~listelement_d() {
		if (prev != nullptr)
			prev->next = next;
		if (next != nullptr)
			next->prev = prev;
	}

	friend class list_double<T>;
};

/*
 * class list_double - Templated list management class.
 *
 * This class manages a doubly linked list of type listelement_d<T>.
 * It has both a head and a tail pointer.
 */
template<typename T>
class list_double {
protected:
	T * first;                  // The head of the list.
	T * last;                   // The tail of the list.

	// All lists start empty...
	list_double()
	  : first(nullptr), last(nullptr) {
	}
	list_double(const list_double &) = delete;
	list_double & operator = (const list_double &) = delete;
	list_double(list_double && e)
	  : first(e.first), last(e.last) {
		e.first = e.last = nullptr;
	}
	list_double & operator = (list_double && e) {
		first = e.first, last = e.last;
		e.first = e.last = nullptr;
		return *this;
	}
	~list_double() {
		empty();
	}

	// Insert an element.  If the element is not already part of a list
	// it is inserted at the end.
	T * insert(T * e) {
		if (e == nullptr)
			return e;
		if (e->prev == nullptr && e->next == nullptr) {
			e->prev = this->last;
			if (e->prev != nullptr)
				e->prev->next = e;
		}
		if (e->prev == nullptr) {
			if (this->first != e->next)
				throw std::runtime_error("Linked list pointer mismatch!");
			this->first = e;
			if (e->next != nullptr)
				e->next->prev = e;
		}
		if (e->next == nullptr) {
			if (this->last != e->prev)
				throw std::runtime_error("Linked list pointer mismatch!");
			this->last = e;
			if (e->prev != nullptr)
				e->prev->next = e;
		}
		return e;
	}
	// Remove an element from the list (but don't delete it).
	T * remove(T * e = nullptr) {
		if (e == nullptr)
			e = last;
		if (e == nullptr)
			return e;
		if (e->prev == nullptr)
			this->first = e->next;
		else
			e->prev->next = this->next;
		if (e->next == nullptr)
			this->last = e->prev;
		else
			e->next->prev = this->prev;
		e->next = nullptr;
		e->prev = nullptr;
		return e;
	}
	// Push an element onto the head of the list (or just insert)
	T * push(T * e) {
		if (e == nullptr)
			return e;
		if (e->prev == nullptr && e->next == nullptr) {
			e->next = this->first;
			if (e->next != nullptr)
				e->next->prev = e;
		}
		return insert(e);
	}
	// Pop the first element off the list
	T * pop() {
		return remove(first);
	}
	// Check in the list is empty.
	bool isempty() const {
		return (first == nullptr ? true : false);
	}
	// Figure out the length of the list.
	unsigned length() const {
		unsigned s = 0;
		T * t = first;
		while (t != nullptr)
			s++, t = t->next;
		return s;
	}
	// Sometimes you just need to start a new list.
	void empty() {
		while (first != nullptr)
			delete first;
	}
};

/*
 * class listelement_o - Templated owned list element.
 *
 * List is like listelemt_d, but requires an owner of type O.
 * The owner class must be derived from class list_double.
 */
template<typename O, typename T>
class listelement_o : public listelement_d<T> {
protected:
	O * owner;                  // Our owner.

	// Create an element.
	listelement_o(O * o, T * p = nullptr, T * n = nullptr)
	  : listelement_d<T>(p,n), owner(o) {
		if (owner == nullptr)
			throw std::invalid_argument("Owner cannot be null!");
		if (this->prev == nullptr && this->next == nullptr) {
			this->prev = owner->last;
			if (this->prev != nullptr)
				this->prev->next = static_cast<T *>(this);
		}
		if (this->prev == nullptr) {
			if (owner->first != this->next)
				throw std::runtime_error("Linked list pointer mismatch!");
			owner->first = static_cast<T *>(this);
		}
		if (this->next == nullptr) {
			if (owner->last != this->prev)
				throw std::runtime_error("Linked list pointer mismatch!");
			owner->last = static_cast<T *>(this);
		}
	}
	// Copying or moving list elements is inherently unsafe.
	listelement_o(const listelement_o &) = delete;
	listelement_o & operator = (const listelement_o &) = delete;
	listelement_o(listelement_o &&) = delete;
	listelement_o & operator = (listelement_o &&) = delete;
	// Also manage our owner's pointers.
	~listelement_o() {
		if (owner == nullptr)
			return;
		if (this->prev == nullptr)
			owner->first = this->next;
		if (this->next == nullptr)
			owner->last = this->prev;
	}
	// Figure out the position of this element in the list.
	unsigned position() const {
		unsigned s = 0;
		T * t = owner->first;
		while (static_cast<listelement_o *>(t) != this)
			s++, t = t->next;
		return s;
	}

	friend class list_owned<O,T>;
};

/*
 * class list_owned - Templated list owner class.
 *
 * This is the base class for the list owner.
 */
template<typename O, typename T>
class list_owned : public list_double<T> {
protected:
	using list_double<T>::first;

	list_owned()
	  : list_double<T>() {
	}
	list_owned(const list_owned &) = delete;
	list_owned & operator = (const list_owned &) = delete;
	list_owned(list_owned && l)
	  : list_double<T>(std::move(l)) {
		listelement_o<O,T> * t = first;
		while (t != nullptr) {
			t->owner = static_cast<O *>(this);
			t = t->next;
		}
	}
	list_owned & operator = (list_owned && l) {
		list_double<T>::operator=(std::move(l));
		listelement_o<O,T> * t = first;
		while (t != nullptr) {
			t->owner = static_cast<O *>(this);
			t = t->next;
		}
	}
	~list_owned() {
	}

	friend class listelement_o<O,T>;
};

} // namespace OP

#endif // LIST_H
