/**************************************************************************

	LISTEN.H - Simple event dispatch/listener template classes

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
		This header implements a number of templated C++ classes for
		a simple event dispatcher/listener model.

	Design:
		There are two sets of two classes: a listener item and dispatcher
		that stores these, and a listener object with a list of dispatchers.
		The dispatcher is created first, and listeners register with it.

	Status:
		This is not intended for an async event handling model with objects
		that might go away.

	History:
		2016/02/22 - Created by Jeremy Lea <reg@openpave.org>

**************************************************************************/

#pragma once
#ifndef __LISTEN_H
#define __LISTEN_H

#include <cassert>
#include <functional>
#include <memory>

namespace OP {

template<typename M> class dispatcher;

/*
 * class message
 */
template<typename...Ts>
struct message {
	message(std::function<void(Ts...)> && f)
	  : handler(std::move(f)) {
	}
	message(message && m)
	  : handler(std::move(m.handler)) {
	}
	message & operator = (message && m) {
		handler = std::move(m.handler);
		return *this;
	}
	message(const message &) = delete;
	message & operator = (const message &) = delete;
	~message() {
	}
private:
	friend class dispatcher<message<Ts...>>;

	void dispatch(Ts...args) {
		handler(args...);
	}
	std::function<void(Ts...)> handler;
};

/*
 * class dispatcher
 */
template<typename...Ts>
class dispatcher<message<Ts...>> {
public:
	std::function<void(void)> attach(message<Ts...> && m) {
		// Attach to end so they are dispatched in the order attached.
		sink ** s = &head, * t;
		while (*s != nullptr)
			s = &((*s)->next);
		t = *s = new sink(std::move(m));
		return std::function<void(void)>([this,t](){
			sink ** c = &head;
			while (*c != nullptr) {
				if (*c == t) {
					*c = t->next;
					delete t;
					return;
				}
				c = &(*c)->next;
			}
		});
	}

protected:
	dispatcher()
	  : head(nullptr) {
	}
	~dispatcher() {
		assert(head == nullptr);
	}
	void dispatch(const Ts &...args) {
		sink * s = head;
		while (s != nullptr) {
			s->dispatch(args...);
			s = s->next;
		}
	}

private:
	struct sink : message<Ts...> {
		sink(message<Ts...> && m)
		  : message<Ts...>(std::move(m)), next(nullptr) {
		}
		~sink() {
		}
		sink * next;           // The next listener
	} * head;                  // The list head
};

/*
 * class listener
 */
class listener {
protected:
	listener()
	  : head(nullptr) {
	}
	~listener() {
		source * s;
		while ((s = head) != nullptr) {
			head = s->next;
			delete s;
		}
	}
	template<typename...Ts>
	void listen(dispatcher<message<Ts...>> & d, message<Ts...> && m) {
		source ** s = &head;
		while (*s != nullptr)
			s = &((*s)->next);
		*s = new source(d.attach(std::move(m)));
	}

private:
	// Simple linked list of callbacks to remove ourself from dispatcher
	struct source {
		source(std::function<void(void)> && r)
		  : remover(r), next(nullptr) {
		}
		~source() {
			// Ignore any failures here.
			try {
				remover();
			} catch (...) {}
		}
		std::function<void(void)> remover;
		source * next;         // The next source
	} * head;                  // The head of the list
};

} // namespace OP

#endif // LISTEN_H
