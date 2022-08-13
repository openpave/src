/**************************************************************************

	TASK.H - A basic producer/consumer task queue

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
		This header implements a basic producer/consumer model task queue,
		with appropriate locking and other features.

	Design:
		The design follows many similar designs from on-line tutorials.
		The class has four main components: a pool of worker threads that
		run and work on tasks as they are assigned, a queue of tasks that
		can be added to, a queue of results that are reported back to the
		caller, and the associated locking primitives to make sure we
		don't get in a tangle.

		Tasks cannot take any arguments, but can return values or throw. If
		a task throws, then thread processing is aborted.  Tasks will
		be processed in the order they are enqueued and the callbacks are
		likely to also be ordered.

		The class makes use of the standard C++ threading tools as much as
		possible.  A helper class takes care of the special processing
		needed for void return values.

	History:
		2017/12/20 - Created a basic implementation.

**************************************************************************/

#ifndef __TASK_H
#define __TASK_H

#include <algorithm>
#include <future>
#include <queue>
#include <thread>

namespace OP {

// Wrap the types so we can specialize for void.
template<typename T>
struct task_helper {
	using function_t = std::packaged_task<T()>; // The actual task
	using callback_t = std::function<void(T)>;  // The callback
	using task_t = std::tuple<function_t,callback_t>;
	using result_t = std::tuple<std::future<T>,callback_t>;
	// Do the callback with the result
	static void callback(result_t & r) {
		std::get<1>(r)(std::get<0>(r).get());
	}
	// Fallback for callback
	static void fallback(T) {
	}
	// Direct call without threading
	static void direct(task_t & t) {
		std::get<0>(t)();
		std::get<1>(t)(std::get<0>(t).get_future().get());
	}
};
// and then for void...
template<>
struct task_helper<void> {
	using function_t = std::packaged_task<void()>; // The actual task
	using callback_t = std::function<void()>;      // The callback
	using task_t = std::tuple<function_t,callback_t>;
	using result_t = std::tuple<std::future<void>,callback_t>;
	static void callback(result_t & r) {
		std::get<0>(r).get();
		std::get<1>(r)();
	}
	static void fallback() {
	}
	static void direct(task_t & t) {
		try {
			std::get<0>(t)();
			std::get<0>(t).get_future().get();
			std::get<1>(t)();
		} catch (...) {
			throw;
		}
	}
};

/*
 * class task_queue - a basic task queue
 */
template<typename T = void>
class task_queue
{
public:
	using callback_t = typename task_helper<T>::callback_t;

	// basic constructor.
	task_queue(std::size_t s = 0) {
		s = s == 0 ? std::thread::hardware_concurrency()+2 : s;
		s = std::max(std::size_t(1),s);
		for (std::size_t i = 0; i < s; i++)
			// implicit call to std::thread constructor.
			workers.emplace_back(std::bind(&task_queue::worker, this));
	}
	// basic destructor.
	~task_queue() {
		bool done = false;

		if (exiting.exchange(true))
			return; // already exiting?
		while (!done) {
			std::unique_lock<std::mutex> lock{mtx};

			if (!havework()) {
				done = true;
				if (empty.exchange(true))
					return; // already empty?
			}
			cv.notify_all(); // wake everyone up to exit.
			lock.unlock(); // unlock so they can work.
			std::this_thread::yield();
		}
		for (auto & thread : workers)
			thread.join();
		drain();
	}
	void drain() {
		bool done = false;

		while (!done) {
			std::unique_lock<std::mutex> lock{mtx};

			if (!havework() && !haveresults())
				done = true;
			if (!empty && havework())
				cv.notify_one();
			lock.unlock(); // unlock so they can work.
			std::this_thread::yield();
			result(true); // wait and return;
		}
	}
	// enqueue some work into the task queue and try to wake up a worker
	template <typename F, typename C = callback_t>
	void enqueue(F && f, C && cb = callback_t(task_helper<T>::fallback)) {
		std::unique_lock<std::mutex> lock{mtx};

		// Don't add more work if we are exiting
		if (exiting || empty || abort)
			return; // throw?
		if (single_thread) {
			task_t t{std::forward<F>(f),std::forward<C>(cb)};
			task_helper<T>::direct(t);
			return;
		}
		tasks.emplace(std::forward<F>(f),std::forward<C>(cb));
		cv.notify_one();
		bool rv = haveresults();
		lock.unlock(); // unlock so they can work.
		if (rv)
			result(false); // try to return a result.
	}
	void force_single_threaded() {
		single_thread = true;
		drain();
	}

private:
	using task_t = typename task_helper<T>::task_t;
	using result_t = typename task_helper<T>::result_t;

	// prevent move and copy
	task_queue(const task_queue &) = delete;
	task_queue(task_queue &&) = delete;
	task_queue & operator = (const task_queue &) = delete;
	task_queue & operator = (task_queue &&) = delete;

	// This runs in the threads but in the same context as the main thread.
	void worker() {
		task_t task;

		while (true) try {
			std::unique_lock<std::mutex> lock{mtx};

			cv.wait(lock, [&]() {
				return exiting || empty || abort || havework();
			});
			if (empty || abort)
				return;
			if (gettask(task)) {
				results.emplace(std::get<0>(task).get_future(),
						std::get<1>(task));
				lock.unlock();
				std::get<0>(task)();
			} else if (exiting)
				return;
		} catch (...) {
			std::promise<T> p;
			try {
				results.emplace(p.get_future(),
						callback_t(task_helper<T>::fallback));
				// store anything thrown in the promise
				p.set_exception(std::current_exception());
				abort = true; // start abort ASAP
				return; // assume this thread should die
			} catch(...) { // set_exception() may throw too
				throw; // just throw from the thread
			}
		}
	}
	// This runs in the main thread.
	void result(bool wait) {
		std::unique_lock<std::mutex> lock{mtx};
		result_t result;

		if (getresult(result,wait) && !abort) {
			lock.unlock();
			try {
				// The internal get() may re-throw an exception.
				task_helper<T>::callback(result);
			} catch (...) {
				// try to abort all threads ASAP.
				abort = true;
				throw;
			}
		}
	}
	// Check if there is work (call with mtx locked)
	bool havework() const {
		return !tasks.empty();
	}
	// Check if there is work (call with mtx locked)
	bool haveresults() const {
		return !results.empty();
	}
	// try to get a task for our worker (call with mtx locked)
	bool gettask(task_t & task) {
		if (tasks.empty())
			return false;
		task = std::move(tasks.front());
		tasks.pop();
		return true;
	}
	// get a result and try to do the call back (call with mtx locked)
	bool getresult(result_t & result, bool wait) {
		if (results.empty())
			return false;
		const result_t & r = results.front();
		if (wait) {
			std::get<0>(r).wait();
		} else {
			std::future_status status =
					std::get<0>(r).wait_for(std::chrono::seconds(0));
			if (status != std::future_status::ready)
				return false;
		}
		result = std::move(results.front());
		results.pop();
		return true;
	}

	std::queue<task_t>       tasks;
	std::queue<result_t>     results;
	std::vector<std::thread> workers;
	std::condition_variable  cv;
	std::mutex               mtx;
	std::atomic<bool>        exiting{false};
	std::atomic<bool>        empty{false};
	std::atomic<bool>        abort{false};
	std::atomic<bool>        single_thread{false};
};

} // namespace OP

#endif // TASK_H
