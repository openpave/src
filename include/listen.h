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
 * class channel - Base class for channels
 *
 * The template parameter pack Ts defines the message.
 */
template<typename... Ts>
struct channel
{
	// Construct a channel with a callback
	explicit channel(std::function<void(Ts...)> && f) noexcept
	  : handler(std::move(f)) {
	}
	channel(channel && c) noexcept
	  : handler(std::move(c.handler)) {
	}
	channel & operator = (channel && c) noexcept {
		handler = std::move(c.handler);
		return *this;
	}
	channel(const channel &) = delete;
	channel & operator = (const channel &) = delete;
	~channel() = default;

private:
	friend class dispatcher<channel<Ts...>>;

	// Call the message handler provided by the listener.
	void dispatch(Ts... args) {
		handler(args...);
	}
	std::function<void(Ts...)> handler;
};

/*
 * class dispatcher - Base class for classes that send messages
 *
 * Maintains a list of channels to various listeners.  These are registered
 * via attach().  Each dispatcher only sends one type of message, but you
 * can inherit multiple times to send different types of message.
 */
template<typename... Ts>
class dispatcher<channel<Ts...>>
{
public:
	// Disallow copying
	dispatcher(const dispatcher &) = delete;
	dispatcher & operator = (const dispatcher &) = delete;
	// Attach a message handler, return a remover
	std::function<void(void)> attach(channel<Ts...> && c) {
		// Attach to end so they are dispatched in the order attached.
		sink ** s = &head, * t;
		while (*s != nullptr)
			s = &((*s)->next);
		t = *s = new sink(std::move(c));
		return std::function<void(void)>([this,t] () noexcept {
			sink ** c = &head;
			while (*c != nullptr) {
				if (*c == t) {
					*c = t->next;
					delete t;
					return;
				}
				c = &((*c)->next);
			}
		});
	}

protected:
	dispatcher() noexcept = default;
	dispatcher(dispatcher &&) noexcept = default;
	dispatcher & operator = (dispatcher &&) noexcept = default;
	~dispatcher() {
		// Check that we are deleted last
		assert(head == nullptr);
	}
	// Only derived classes can dispatch messages
	void dispatch(const Ts &... args) {
		sink * s = head;
		while (s != nullptr) {
			s->dispatch(args...);
			s = s->next;
		}
	}

private:
	// Private internal class of sinks that are listening for dispatched
	// messages.  Simple singularly linked list of message classes.
	struct sink
	  : channel<Ts...> {
		explicit sink(channel<Ts...> && c) noexcept
		  : channel<Ts...>(std::move(c)) {
		}
		sink(const sink &) = delete;
		sink(sink &&) = delete;
		sink & operator = (const sink &) = delete;
		sink & operator = (sink &&) = delete;
		~sink() = default;
		sink * next{nullptr};  // The next listener
	} * head{nullptr};         // The list head
};

/*
 * class listener - Base class for listening classes
 *
 * Maintains a list of callbacks to remove ourselves from any dispatcher
 * we have attached to.  A listener can listen to multiple different types
 * of message from different dispatcher sources.
 */
class listener
{
public:
	// Disallow copying (else who gets the messages?)
	listener(const listener &) = delete;
	listener & operator = (const listener &) = delete;

protected:
	listener() noexcept = default;
	listener(listener &&) noexcept = default;
	listener & operator = (listener &&) noexcept = default;
	~listener() {
		// Remove ourselves from any sources
		source * s;
		while ((s = head) != nullptr) {
			head = s->next;
			delete s;
		}
	}
	// Listen to a dispatcher for a type of message.
	template<typename... Ts>
	void listen(dispatcher<channel<Ts...>> & d, channel<Ts...> && c) {
		source ** s = &head;
		while (*s != nullptr)
			s = &((*s)->next);
		*s = new source(d.attach(std::move(c)));
	}

private:
	// Simple linked list of callbacks to remove ourself from dispatcher
	struct source {
		source(std::function<void(void)> && r) noexcept
		  : remover(std::move(r)) {
		}
		source(const source &) = delete;
		source(source &&) = delete;
		source & operator = (const source &) = delete;
		source & operator = (source &&) = delete;
		~source() {
			// Ignore any failures here.
			try {
				remover();
			} catch (...) {}
		}
		std::function<void(void)> remover;
		source * next{nullptr}; // The next source
	} * head{nullptr};          // The head of the list
};

} // namespace OP

#endif // LISTEN_H
