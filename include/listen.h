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

	Portions Copyright (C) 2016 OpenPave.org.

	Contributor(s): Jeremy Lea <reg@openpave.org>.

	Purpose:
		This header implements a number of templated C++ classes for
		a simple event dispatcher/listener model.

	Design:
		There are two sets of two classes: a listener item and dispatcher
		that stores these, and a listener object with a list of dispatchers. 
		The dispatcher is created first, and listener's register with it.

	Status:
		This is not intended for an async event handling model with objects
		that might go away.

	History:
		2016/02/22 - Created by Jeremy Lea <reg@openpave.org>

**************************************************************************/

#ifndef __LISTEN_H
#define __LISTEN_H

#include <memory>
#include <functional>

template<typename M> class dispatcher;

/*
 * class message
 */
template<typename T, typename ...Ts>
struct message {
	message(std::function<void(const T &, Ts...)> && f)
	  : handler(std::move(f)) {
	}
	~message() {
	}
private:
	friend class dispatcher<message<T,Ts...>>;

	void dispatch(const T & t, Ts...args) {
		handler(t,args...);
	}
	std::function<void(const T &, Ts...)> handler;
};

/*
 * class dispatcher
 */
template<typename T, typename ...Ts>
class dispatcher<message<T,Ts...>> {
public:
	std::function<void(void)> listen(message<T,Ts...> && m) {
		sink ** s = &head;
		while (*s != nullptr)
			s = &((*s)->next);
		*s = new sink(std::move(m));
		return std::function<void(void)>([this,s](){
			sink ** c = &head, ** l = &head;
			while (*c != nullptr) {
				if (*c == *s) {
					*l = (*c)->next;
					delete *c;
					c = l;
				} else
					l = c, c = &((*c)->next);
			}
		});
	}

protected:
	dispatcher()
	  : head(nullptr) {
	}
	~dispatcher() {
		sink * s;
		while ((s = head) != nullptr) {
			head = s->next;
			delete s;
		}
	}
	void dispatch(const Ts &...args) {
		sink * s = head;
		while (s != nullptr) {
			s->dispatch(static_cast<const T &>(*this),args...);
			s = s->next;
		}
	}

private:
	struct sink : message<T,Ts...> {
		sink(message<T,Ts...> && m)
		  : message(std::move(m)), next(nullptr) {
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
			s->remover();
			delete s;
		}
	}
	template<typename T, typename ...Ts>
	void listen(dispatcher<message<T,Ts...>> & d, message<T,Ts...> && m) {
		source ** s = &head;
		while (*s != nullptr)
			s = &((*s)->next);
		*s = new source(d.listen(std::move(m)));
	}

private:
	// Simple linked list of callbacks to remove ourself from dispatcher
	struct source {
		source(std::function<void(void)> & r)
		  : remover(r), next(nullptr) {
		}
		~source() {
		}
		std::function<void(void)> remover;
		source * next;         // The next source
	} * head;                  // The head of the list
};

#endif // LISTEN_H
