/**************************************************************************

	LIST.H - Simple list template classes

	$OpenPave$

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

	The Original Code is OpenPave.org Core Libraries.

	The Initial Developer of the Original Code is OpenPave.org.

	Portions Copyright (C) 2006 OpenPave.org.

	Contributor(s): Jeremy Lea <reg@openpave.org>.

	Purpose:
		This header implements a number of templated C++ classes for
		various types of lists.

	Design:
		There are three types of list elements, and obviously the lists
		require that you use the right types.  There are singly linked,
		doubally linked and owned elements.  Owned elements are the only
		strange ones, since they require that you define the type of
		class that can own this list type.  This is very useful for a
		number of engineering problems, where you want the list class
		to be able to do things with the logical unit above, plus you can
		have lists which manage themselves entirely in their constructors
		and destructors.

	Status:
		At the moment this only does owned doubally linked lists, since
		that's all I need.

	History:
		2002/01/23 - Created by Jeremy Lea <reg@openpave.org>

**************************************************************************/

#ifndef __LIST_H
#define __LIST_H

#include "config.h"

// Forward declare some classes...
template <class T>
class list_double;
template <class O, class T>
class list_owned;

/*
 * class listelement_d - Templated doublally linked list element.
 *
 * Use this class as a base class for the data you want to store in the
 * list.
 */
template <class T>
class listelement_d {
protected:
	T *next;							// The next element.
	T *prev;							// The previous element.

	// Create a new list element, with an optional previous element.
	// We force the consumer to use NULL so that they have to think.
	listelement_d(T *p) {
		prev = p;
		next = NULL;
		if (prev != NULL) {
			if (prev->next != NULL) {
				next = prev->next;
				next->prev = (T *)this;
			}
			prev->next = (T *)this;
		}
	};
	// Unlink ourselves from the list before we die...
	virtual ~listelement_d() {
		if (prev != NULL)
			prev->next = next;
		if (next != NULL)
			next->prev = prev;
	};
	// These are so the other classes can access our points.
	//friend class listelement_d<T>;
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
	O *owner;							// Our owner.

	// Create an element.
	listelement_o(O *o, T *p) : listelement_d<T>(p) {
		owner = o;
		if (owner != NULL) {
			if (this->prev == NULL) {
				this->prev = owner->last;
				if (this->prev != NULL)
					this->prev->next = (T *)this;
				else
					owner->first = (T *)this;
			}
			if (this->next == NULL)
				owner->last = (T *)this;
		}
	};
	// Also manage our owner's pointers.
	virtual ~listelement_o() {
		if (owner != NULL) {
			if (this->prev == NULL)
				owner->first = this->next;
			if (this->next == NULL)
				owner->last = this->prev;
		}
	};
	// Be a bit social.
	//friend class O;
	friend class list_owned<O,T>;
};

/*
 * class list_double - Templated list management class.
 *
 * This class manages a doubally linked list of type listelement_d<T>.
 * It has both a head and a tail pointer.
 */
template <class T>
class list_double {
protected:
	T *first;							// The head of the list.
	T *last;							// The tail of the list.

	// All lists start empty...
	list_double() {
		first = NULL;
		last = NULL;
	};
	virtual ~list_double() {
		empty();
	};

	// Check in the list is empty.
	bool inline isempty() const {
		return (first == NULL ? true : false);
	};
	// Figure out the length of the list.
	int length() const {
		int s = 0;
		T *t = first;
		while (t != NULL) {
			t = t->next;
			s++;
		}
		return s;
	};
	// Sometimes you just need to start a new list.
	void empty() {
		if (first != NULL) {
			while (first->next != NULL)
				delete first->next;
			delete first;
		}
		first = NULL;
		last = NULL;
	};
};

/*
 * class list_owned - Templated list owner class.
 *
 * This is the base class for the list owner.
 */
template <class O, class T>
class list_owned : public list_double<T> {
protected:
	list_owned() : list_double<T>() {
	};
	virtual ~list_owned() {
	};

	friend class listelement_o<O,T>;
	//friend class O;
};

#endif // LIST_H
