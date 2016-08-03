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

/* This file was originally downloaded from the following URL
 * https://raw.githubusercontent.com/Youka/ThreadPool/master/ThreadPool.hpp
 * Below is the copyright found in the original file.
 * 
 * Copyright (c) 2012 Jakob Progsch, Václav Zeman
 *
 * This software is provided 'as-is', without any express or implied
 * warranty. In no event will the authors be held liable for any damages
 * arising from the use of this software.
 *
 * Permission is granted to anyone to use this software for any purpose,
 * including commercial applications, and to alter it and redistribute it
 * freely, subject to the following restrictions:
 *
 *   1. The origin of this software must not be misrepresented; you must not
 *   claim that you wrote the original software. If you use this software
 *   in a product, an acknowledgment in the product documentation would be
 *   appreciated but is not required.
 *
 *   2. Altered source versions must be plainly marked as such, and must not be
 *   misrepresented as being the original software.
 *
 *   3. This notice may not be removed or altered from any source
 *   distribution.
 */

#ifndef DNG_THREAD_POOL_H
#define DNG_THREAD_POOL_H

// containers
#include <vector>
#include <queue>
// threading
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <future>
// utility wrappers
#include <memory>
#include <functional>
// exceptions
#include <stdexcept>

namespace dng {
// a pool of std::thread for resources recycling
class ThreadPool {
public:
    // the constructor just launches some amount of workers
    explicit ThreadPool(size_t threads_n);

    // deleted copy&move ctors&assignments
    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;
    ThreadPool(ThreadPool&&) = delete;
    ThreadPool& operator=(ThreadPool&&) = delete;

    // the destructor joins all threads
    virtual ~ThreadPool();

    // add new work item to the pool
    template<class F, class... Args>
    std::future<typename std::result_of<F(Args...)>::type> Enqueue(F&& f, Args&&... args);

 private:
    // need to keep track of threads so we can join them
    std::vector< std::thread > workers_;
    // the task queue
    std::queue< std::function<void()> > tasks_;

    // synchronization
    std::mutex queue_mutex_;
    std::condition_variable condition_;
    bool stop_;
};

inline
ThreadPool::ThreadPool(size_t threads_n = std::thread::hardware_concurrency()) : stop_{false}
{
    // if threads_n is 0, set it to 1
    if(threads_n == 0) {
        threads_n = 1;
    }

    this->workers_.reserve(threads_n);

    for(; threads_n; --threads_n) {
        this->workers_.emplace_back( [this] { for(;;) {
            std::function<void()> task;

            {
                // wait until more jobs is available 
                std::unique_lock<std::mutex> lock(this->queue_mutex_);
                this->condition_.wait(lock,
                    [this]{ return this->stop || !this->tasks_.empty(); });
                // quit the worker if we have been told to stop
                if(this->stop_) {
                    return;
                }
                task = std::move(this->tasks_.front());
                this->tasks_.pop();
            }

            task();
        }});
    }
}

// add new work item to the pool
template<class F, class... Args>
inline
std::future<typename std::result_of<F(Args...)>::type>
ThreadPool::Enqueue(F&& f, Args&&... args)
{
    using packaged_task_t = std::packaged_task<typename std::result_of<F(Args...)>::type ()>;

    std::shared_ptr<packaged_task_t> task(new packaged_task_t(
        std::bind(std::forward<F>(f), std::forward<Args>(args)...)
    ));
    auto res = task->get_future();
    {
        std::lock_guard<std::mutex> lock(this->queue_mutex_);
        this->tasks_.emplace([task](){ (*task)(); });
    }
    this->condition_.notify_one();
    return res;
}

// the destructor joins all threads
inline
virtual ThreadPool::~ThreadPool()
{
    {
        std::lock_guard<std::mutex> lock(this->queue_mutex);
        this->stop = true;
    }
    this->condition.notify_all();
    for(std::thread& worker : this->workers) {
        worker_.join();
    }
}

} // namespace dng

#endif // DNG_THREAD_POOL_H
