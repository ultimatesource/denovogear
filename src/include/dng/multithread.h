/*
 * Copyright (c) 2016 Reed A. Cartwright
 * Authors:  Reed A. Cartwright <reed@cartwrig.ht>
 *
 * This file is part of DeNovoGear.
 *
 * DeNovoGear is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once
#ifndef DNG_MULTITHREAD_H
#define DNG_MULTITHREAD_H

#include <vector>
#include <queue>
#include <tuple>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future>
#include <memory>
#include <functional>
#include <type_traits>

#include <boost/fusion/functional/invocation/invoke.hpp>
#include <boost/fusion/adapted/std_tuple.hpp>

#include <dng/detail/unit_test.h>

namespace dng {
namespace multithread {

// basic thread pool supporting only a single type of job
template<typename... Args>
class BasicPool {
public:
    typedef std::tuple<typename std::remove_reference<Args>::type...> tuple_t;

    template<typename F>
    BasicPool(F f, size_t num_threads = std::thread::hardware_concurrency());

    virtual ~BasicPool();

    // add new work item to the pool
    template<typename... Args2>
    void Enqueue(Args2&&... args);

    // Erase all pending jobs
    void Clear();

    // do not copy, move, or assign
    BasicPool(const BasicPool&) = delete;
    BasicPool& operator=(const BasicPool&) = delete;
    BasicPool(BasicPool&&) = delete;
    BasicPool& operator=(BasicPool&&) = delete;

private:
    // need to keep track of threads so we can join them
    std::vector<std::thread> workers_;
    // the task queue
    std::queue<tuple_t> tasks_;

    // synchronization
    std::mutex mutex_;
    std::condition_variable condition_;
    bool stop_;    
};

template<typename... Args>
template<typename F>
BasicPool<Args...>::
BasicPool(F f, size_t num_threads) : stop_{false} {
    // if threads_n is 0, set it to 1
    if(num_threads == 0) {
        num_threads = 1;
    }

    auto lambda = [this,f]() {
        // wrap our callback in std::function so invoke will work correctly
        std::function<void(Args...)> g(f);
    	for(;;) {
            tuple_t args;
            {
                // wait until a new task is available 
                std::unique_lock<std::mutex> lock(mutex_);
                condition_.wait(lock,
                    [this]{ return stop_ || !tasks_.empty(); });
                // quit the worker if we have been told to stop and we are done.
                if(stop_ && tasks_.empty()) {
                    return;
                }
                args = std::move(tasks_.front());
                tasks_.pop();
            }
            boost::fusion::invoke(g,args);
        }
    };

    workers_.reserve(num_threads);
    for(; num_threads; --num_threads) {
        // construct a vector of workers to process data from the queue
        workers_.emplace_back(lambda);
    }
}

template<typename... Args>
BasicPool<Args...>::
~BasicPool() {
    {
        std::lock_guard<std::mutex> lock(mutex_);
        stop_ = true;
    }
    condition_.notify_all();
    for(auto && worker : workers_) {
        worker.join();
    }
}

// for perfect forwarding to work this must be a template function
template<typename... Args>
template<typename... Args2>
void BasicPool<Args...>::
Enqueue(Args2&&... args) {
    {
        std::lock_guard<std::mutex> lock(mutex_);
        tasks_.emplace(std::forward<Args2>(args)...);
    }
    condition_.notify_one();    
}

template<typename... Args>
void BasicPool<Args...>::
Clear() {
    std::lock_guard<std::mutex> lock(mutex_);
    tasks_.clear();
}

} // namespace dng::multithread
} // namespace dng

#endif // DNG_MULTITHREAD_H
