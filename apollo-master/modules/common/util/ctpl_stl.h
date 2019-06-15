/*********************************************************
*
*  Copyright (C) 2014 by Vitaliy Vitsentiy
*
*  Modifications Copyright 2017 The Apollo Authors. All Rights Reserved.
*
*  Licensed under the Apache License, Version 2.0 (the "License");
*  you may not use this file except in compliance with the License.
*  You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
*  Unless required by applicable law or agreed to in writing, software
*  distributed under the License is distributed on an "AS IS" BASIS,
*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*  See the License for the specific language governing permissions and
*  limitations under the License.
*
*********************************************************/

// This file is a modification of
// https://github.com/vit-vit/CTPL/blob/master/ctpl_stl.h

#ifndef MODULES_COMMON_UTIL_CTPL_STL_H_
#define MODULES_COMMON_UTIL_CTPL_STL_H_

#include <atomic>
#include <chrono>
#include <exception>
#include <functional>
#include <future>
#include <memory>
#include <mutex>
#include <queue>
#include <thread>
#include <utility>
#include <vector>

/*
 * thread pool to run user's functions with signature
 *      ret func(int id, other_params)
 * where id is the index of the thread that runs the functions
 * ret is some return type
 */

namespace apollo {
namespace common {
namespace util {

namespace detail {
template <typename T>
class Queue {
 public:
  bool push(T const &value) {
    std::unique_lock<std::mutex> lock(mutex_);
    q_.push(value);
    return true;
  }
  // deletes the retrieved element, do not use for non integral types
  bool pop(T &v) {  // NOLINT
    std::unique_lock<std::mutex> lock(mutex_);
    if (q_.empty()) {
      return false;
    }
    v = q_.front();
    q_.pop();
    return true;
  }
  bool empty() {
    std::unique_lock<std::mutex> lock(mutex_);
    return q_.empty();
  }

 private:
  std::queue<T> q_;
  std::mutex mutex_;
};
}  // namespace detail

class ThreadPool {
 public:
  ThreadPool() {
    Init();
  }
  explicit ThreadPool(int n_threads) {
    Init();
    Resize(n_threads);
  }

  // the destructor waits for all the functions in the queue to be finished
  ~ThreadPool() {
    Stop(true);
  }

  // get the number of running threads in the pool
  int size() {
    return static_cast<int>(threads_.size());
  }

  // number of idle threads
  int NumIdle() {
    return n_waiting_;
  }
  std::thread &GetThread(const int i) {
    return *(threads_[i]);
  }

  // change the number of threads in the pool
  // should be called from one thread, otherwise be careful to not interleave,
  // also with stop()
  // n_threads must be >= 0
  void Resize(const int n_threads) {
    if (!is_stop_ && !is_done_) {
      int old_n_threads = static_cast<int>(threads_.size());

      if (old_n_threads <=
          n_threads) {  // if the number of threads is increased
        threads_.resize(n_threads);
        flags_.resize(n_threads);
        for (int i = old_n_threads; i < n_threads; ++i) {
          flags_[i] = std::make_shared<std::atomic<bool>>(false);
          SetThread(i);
        }
      } else {  // the number of threads is decreased
        for (int i = old_n_threads - 1; i >= n_threads; --i) {
          *(flags_[i]) = true;  // this thread will finish
          threads_[i]->detach();
        }
        {
          // stop the detached threads that were waiting
          std::unique_lock<std::mutex> lock(mutex_);
          cv_.notify_all();
        }
        // safe to delete because the threads are detached
        threads_.resize(n_threads);

        // safe to delete because the threads
        // have copies of shared_ptr of the
        // flags, not originals
        flags_.resize(n_threads);
      }
    }
  }

  // empty the queue
  void ClearQueue() {
    std::function<void(int id)> *_f;
    // empty the queue
    while (q_.pop(_f)) {
      delete _f;
    }
  }

  // pops a functional wrapper to the original function
  std::function<void(int)> Pop() {
    std::function<void(int id)> *_f = nullptr;
    q_.pop(_f);
    std::unique_ptr<std::function<void(int id)>> func(
        _f);  // at return, delete the function even if an exception occurred
    std::function<void(int)> f;
    if (_f) f = *_f;
    return f;
  }

  // wait for all computing threads to finish and join all threads
  void JoinAll() {
    std::unique_lock<std::mutex> lock(job_mutex_);
    while (n_jobs_to_do_ != 0) {
      job_cv_.wait(lock);
    }
  }

  // wait for all computing threads to finish and stop all threads
  // may be called asynchronously to not pause the calling thread while waiting
  // if is_wait == true, all the functions in the queue are run, otherwise the
  // queue is cleared without running the functions
  void Stop(bool is_wait = false) {
    if (!is_wait) {
      if (is_stop_) {
        return;
      }
      is_stop_ = true;
      for (int i = 0, n = size(); i < n; ++i) {
        *(flags_[i]) = true;  // command the threads to stop
      }
      ClearQueue();  // empty the queue
    } else {
      if (is_done_ || is_stop_) return;
      is_done_ = true;  // give the waiting threads a command to finish
    }
    {
      std::unique_lock<std::mutex> lock(mutex_);
      cv_.notify_all();  // stop all waiting threads
    }
    for (int i = 0; i < static_cast<int>(threads_.size());
         ++i) {  // wait for the computing threads to finish
      if (threads_[i]->joinable()) {
        threads_[i]->join();
      }
    }
    // if there were no threads in the pool but some functions in the queue, the
    // functions are not deleted by the threads
    // therefore delete them here
    ClearQueue();
    threads_.clear();
    flags_.clear();
  }

  template <typename F, typename... Rest>
  auto Push(F &&f, Rest &&... rest) -> std::future<decltype(f(0, rest...))> {
    auto pck =
        std::make_shared<std::packaged_task<decltype(f(0, rest...))(int)>>(
            std::bind(std::forward<F>(f), std::placeholders::_1,
                      std::forward<Rest>(rest)...));
    auto _f = new std::function<void(int id)>([pck](int id) { (*pck)(id); });
    q_.push(_f);
    ++n_jobs_to_do_;
    std::unique_lock<std::mutex> lock(mutex_);
    cv_.notify_one();
    return pck->get_future();
  }

  // run the user's function that excepts argument int - id of the running
  // thread. returned value is templatized
  // operator returns std::future, where the user can get the result and rethrow
  // the catched exceptins
  template <typename F>
  auto Push(F &&f) -> std::future<decltype(f(0))> {
    auto pck = std::make_shared<std::packaged_task<decltype(f(0))(int)>>(
        std::forward<F>(f));
    auto _f = new std::function<void(int id)>([pck](int id) { (*pck)(id); });
    q_.push(_f);
    ++n_jobs_to_do_;
    std::unique_lock<std::mutex> lock(mutex_);
    cv_.notify_one();
    return pck->get_future();
  }

 private:
  // deleted
  ThreadPool(const ThreadPool &);             // = delete;
  ThreadPool(ThreadPool &&);                  // = delete;
  ThreadPool &operator=(const ThreadPool &);  // = delete;
  ThreadPool &operator=(ThreadPool &&);       // = delete;

  void SetThread(int i) {
    std::shared_ptr<std::atomic<bool>> flag(
        flags_[i]);  // a copy of the shared ptr to the flag
    auto f = [this, i, flag /* a copy of the shared ptr to the flag */]() {
      std::atomic<bool> &_flag = *flag;
      std::function<void(int id)> *_f;
      bool is_pop_ = q_.pop(_f);
      while (true) {
        while (is_pop_) {  // if there is anything in the queue
          std::unique_ptr<std::function<void(int id)>> func(
              _f);  // at return, delete the function even if an exception
                    // occurred
          (*_f)(i);
          --n_jobs_to_do_;
          job_cv_.notify_one();
          if (_flag) {
            // the thread is wanted to stop, return even if the queue is not
            // empty yet
            return;
          } else {
            is_pop_ = q_.pop(_f);
          }
        }
        // the queue is empty here, wait for the next command
        std::unique_lock<std::mutex> lock(mutex_);
        ++n_waiting_;
        cv_.wait(lock, [this, &_f, &is_pop_, &_flag]() {
          is_pop_ = q_.pop(_f);
          return is_pop_ || is_done_ || _flag;
        });
        --n_waiting_;
        if (!is_pop_) {
          // if the queue is empty and is_done_ == true or *flag
          // then return
          return;
        }
      }
    };
    threads_[i].reset(
        new std::thread(f));  // compiler may not support std::make_unique()
  }

  void Init() {
    is_stop_ = false;
    is_done_ = false;
    n_waiting_ = 0;
    n_jobs_to_do_ = 0;
  }

  std::vector<std::unique_ptr<std::thread>> threads_;
  std::vector<std::shared_ptr<std::atomic<bool>>> flags_;
  detail::Queue<std::function<void(int id)> *> q_;
  std::atomic<bool> is_done_;
  std::atomic<bool> is_stop_;
  std::atomic<int> n_waiting_;  // how many threads are waiting

  std::mutex job_mutex_;
  std::condition_variable job_cv_;
  std::atomic<int> n_jobs_to_do_;  // how many jobs are NOT finished.

  std::mutex mutex_;
  std::condition_variable cv_;
};

}  // namespace util
}  // namespace common
}  // namespace apollo

#endif  // MODULES_COMMON_UTIL_CTPL_STL_H_
