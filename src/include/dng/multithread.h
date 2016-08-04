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

namespace dng {
namespace multithread {

// basic thread pool supporting heterogeneous jobs
class TaskPool {
public:
    // the constructor just launches some amount of workers
    explicit TaskPool(size_t threads_n);

    // deleted copy&move ctors&assignments
    TaskPool(const TaskPool&) = delete;
    TaskPool& operator=(const TaskPool&) = delete;
    TaskPool(TaskPool&&) = delete;
    TaskPool& operator=(TaskPool&&) = delete;

    // the destructor joins all threads
    virtual ~TaskPool();

    // add new work item to the pool
    template<class F, class... Args>
    std::future<typename std::result_of<F(Args...)>::type> Enqueue(F&& f, Args&&... args);

 private:
    // need to keep track of threads so we can join them
    std::vector< std::thread > workers_;
    // the task queue
    std::queue< std::function<void()> > tasks_;

    // synchronization
    std::mutex mutex_;
    std::condition_variable condition_;
    bool stop_;
};

inline
TaskPool::TaskPool(size_t threads_n = std::thread::hardware_concurrency()) : stop_{false}
{
    // if threads_n is 0, set it to 1
    if(threads_n == 0) {
        threads_n = 1;
    }

    workers_.reserve(threads_n);

    for(; threads_n; --threads_n) {
        workers_.emplace_back( [this]() { for(;;) {
            std::function<void()> task;

            {
                // wait until a new task is available 
                std::unique_lock<std::mutex> lock(mutex_);
                condition_.wait(lock,
                    [this]{ return stop || !tasks_.empty(); });
                // quit the worker if we have been told to stop
                if(stop_) {
                    return;
                }
                task = std::move(tasks_.front());
                tasks_.pop();
            }

            task();
        }});
    }
}

// add new work item to the pool
template<class F, class... Args>
inline
std::future<typename std::result_of<F(Args...)>::type>
TaskPool::Enqueue(F&& f, Args&&... args)
{
    using packaged_task_t = std::packaged_task<typename std::result_of<F(Args...)>::type ()>;

    std::shared_ptr<packaged_task_t> task(new packaged_task_t(
        std::bind(std::forward<F>(f), std::forward<Args>(args)...)
    ));
    auto res = task->get_future();
    {
        std::lock_guard<std::mutex> lock(mutex_);
        tasks_.emplace([task](){ (*task)(); });
    }
    condition_.notify_one();
    return res;
}

// the destructor joins all threads
inline
virtual TaskPool::~TaskPool()
{
    {
        std::lock_guard<std::mutex> lock(mutex_);
        stop = true;
    }
    condition.notify_all();
    for(auto && worker : workers_) {
        worker.join();
    }
}

} // namespace dng::multithread
} // namespace dng

/* TaskPool derives from the ThreadPool class found here
 * https://raw.githubusercontent.com/Youka/ThreadPool/master/ThreadPool.hpp
 * Below is the copyright found in the original class.
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

#endif // DNG_MULTITHREAD_H
