/**************************************************************************

	CACHE.H - Fast caches to avoid the overhead of new/delete

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

	Portions Copyright (C) 2017-2022 OpenPave.org.

	Contributor(s): Jeremy Lea <reg@openpave.org>.

	Purpose:
		This header implements classes for fast caching of memory, to
		avoid firstly the overhead of new/delete on repeated calls, but
		also so that values from prior computation can be reused.  The
		implementation is rather fragile and relies of the user to make
		sure that the data is valid.

	History:
		2017/09/06 - Created a basic implementation.

**************************************************************************/

#pragma once
#ifndef __CACHE_H
#define __CACHE_H

#include <algorithm>
#include <cassert>
#include <cstddef>

namespace OP {

/*
 * class slab_cache - A hackish slab allocator.
 *
 * The result cache is a very hackish slab allocator which stores
 * intermediate calculation values.  It is implemented via a single pointer.
 * That is the head of a linked list of malloc()'ed pages where memory is
 * handed out from.  The pointer to the next page is stored in the first
 * slot in the slab, and the second slot stores the fill level of the slab.
 * If this points to the same address as the cache pointer then we are in a
 * reset state and will not allocate memory but will give back pointers to
 * already allocated space.  The LEsystem classes know if they can just
 * reuse these arrays without recalculating the contents.
 *
 * The pointers are always aligned to sixteen byte boundaries.  On some
 * platforms this could be assisted by allocating aligned space.
 */
class slab_cache
{
public:
	void* allocate(size_t len) {
		const size_t s = align(std::max(len,page_size))
				+ align(sizeof(slab_cache));
		slab_cache* c = this;
		void* p = this, * q = this->mark;
		ptrdiff_t l;

		assert(this->mark != nullptr);
		while (p != nullptr) {
			// now get the new mark
			c = static_cast<slab_cache*>(p);
			if (c->mark == c)
				// The cache is empty or has been reset
				c->mark = align(static_cast<char*>(static_cast<void*>(c))
					+ sizeof(slab_cache));
			q = c->mark;
			// get the amount of space left in this slab.
			l = static_cast<char*>(p) + s - static_cast<char*>(q);
			// check if there is enough space for this request
			if (l-static_cast<ptrdiff_t>(align(len)) >= 0)
				break;
			// get the next slab
			p = c->next;
		}
		// Allocate a new slab
		if (p == nullptr) {
			// First mark the current slab as full.
			c->mark = align(static_cast<char*>(static_cast<void*>(c)) + s);
			p = malloc(s);
			if (p == nullptr)
				throw std::bad_alloc();
			// Connect to current slab
			c->next = static_cast<slab_cache*>(p);
			// Create the new entry at the start of the slab
			c = new(p) slab_cache();
			// set the new mark
			q = c->mark = align(static_cast<char*>(static_cast<void*>(c))
				+ sizeof(slab_cache));
		}
		// Finally increment the mark to cover the space
		c->mark = align(static_cast<char*>(q) + len);
		// Return the original mark to the caller
		return q;
	}
	void reset() noexcept {
		slab_cache* c = this;

		while (c != nullptr) {
			c->mark = c;
			c = c->next;
		}
	}
	static slab_cache* create() {
		void* p;
		slab_cache* c;

		p = malloc(page_size);
		if (p == nullptr)
			throw std::bad_alloc();
		c = new(p) slab_cache();
		c->mark = c;
		return c;
	}
	static void destroy(slab_cache* cache) noexcept {
		void* p;
		slab_cache* c = cache;

		if (c == nullptr)
			return;
		p = c->next;
		while (p != nullptr) {
			c = static_cast<slab_cache*>(p);
			p = c->next;
			c->~slab_cache();
			free(c);
		}
		cache->~slab_cache();
		free(cache);
	}

private:
	// always under allocate a little to help malloc()
	// 6*sizeof(void *) ~= 32+sizeof(this);
	static constexpr const size_t page_size{0x400000 - 6 * sizeof(void*)};

	slab_cache* next;                  // the next slab
	void* mark;                        // the fill level for this slab

	slab_cache() noexcept
		: next(nullptr), mark(nullptr) {
	}
	slab_cache(const slab_cache&) = delete;
	slab_cache(slab_cache&&) = delete;
	slab_cache& operator = (const slab_cache& ) = delete;
	slab_cache& operator = (slab_cache&& ) = delete;
	~slab_cache() {
	}

	static constexpr uintptr_t align(uintptr_t s) noexcept {
		return s + (16 - s % 16) % 16;
	}
	static void* align(const void* p) noexcept {
		const uintptr_t s = reinterpret_cast<uintptr_t>(p);
		return reinterpret_cast<void*>(align(s));
	}
	void* operator new (size_t, void* p) noexcept {
		return p;
	}
	void operator delete (void*, void*) noexcept {
	}
};

}

#endif // CACHE_H
