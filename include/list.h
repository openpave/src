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

	Portions Copyright (C) 2006-2008 OpenPave.org.

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

// Forward declare some classes...
template <class T>
class list_double;
template <class O, class T>
class list_owned;

/*
 * class listelement_d - Templated doubly linked list element.
 *
 * Use this class as a base class for the data you want to store in the
 * list.
 */
template <class T>
class listelement_d {
protected:
	T * next;                   // The next element.
	T * prev;                   // The previous element.

	// Create a new list element, with an optional previous element.
	// We force the consumer to use 0 so that they have to think.
	listelement_d(T * p)
	  : next(0), prev(p) {
		if (prev != 0) {
			if (prev->next != 0) {
				next = prev->next;
				next->prev = static_cast<T *>(this);
			}
			prev->next = static_cast<T *>(this);
		}
	}
	// Unlink ourselves from the list before we die...
	~listelement_d() {
		if (prev != 0)
			prev->next = next;
		if (next != 0)
			next->prev = prev;
	}
	// These are so the other classes can access our points.
#if defined(_MSC_VER) && (_MSC_VER <= 1200)
	friend class listelement_d<T>;
#endif
	friend class list_double<T>;
};

/*
 * class listelement_o - Templated owned list element.
 *
 * List is like listelemt_d, but requires an owner of type O.
 * The owner class must be derived from class list_double.
 */
template <class O, class T>
class listelement_o : public listelement_d<T> {
protected:
	O * owner;                  // Our owner.

	// Create an element.
	listelement_o(O * o, T * p = 0)
	  : listelement_d<T>(p), owner(o) {
		if (owner != 0) {
			if (this->prev == 0) {
				this->prev = owner->last;
				if (this->prev != 0)
					this->prev->next = static_cast<T *>(this);
				else
					owner->first = static_cast<T *>(this);
			}
			if (this->next == 0)
				owner->last = static_cast<T *>(this);
		}
	}
	// Also manage our owner's pointers.
	~listelement_o() {
		if (owner != 0) {
			if (this->prev == 0)
				owner->first = this->next;
			if (this->next == 0)
				owner->last = this->prev;
		}
	}

	//friend class O;
	friend class list_owned<O,T>;
};

/*
 * class list_double - Templated list management class.
 *
 * This class manages a doubly linked list of type listelement_d<T>.
 * It has both a head and a tail pointer.
 */
template <class T>
class list_double {
protected:
	T * first;                  // The head of the list.
	T * last;                   // The tail of the list.

	// All lists start empty...
	list_double()
	  : first(0), last(0) {
	}
	~list_double() {
		empty();
	}

	// Check in the list is empty.
	bool inline isempty() const {
		return (first == 0 ? true : false);
	}
	// Figure out the length of the list.
	unsigned length() const {
		unsigned s = 0;
		T * t = first;
		while (t != 0) {
			t = t->next;
			s++;
		}
		return s;
	}
	// Sometimes you just need to start a new list.
	void empty() {
		if (first != 0) {
			while (first->next != 0)
				delete first->next;
			delete first;
		}
		first = 0;
		last = 0;
	}
};

/*
 * class list_owned - Templated list owner class.
 *
 * This is the base class for the list owner.
 */
template <class O, class T>
class list_owned : public list_double<T> {
protected:
	list_owned()
	  : list_double<T>() {
	}
	~list_owned() {
	}

	friend class listelement_o<O,T>;
	//friend class O;
};

#endif // LIST_H
