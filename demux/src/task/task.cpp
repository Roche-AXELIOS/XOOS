#include "task/task.h"

namespace xoos::demux {
std::exception_ptr Task::task_exception;
std::mutex Task::task_mutex;
std::atomic<bool> Task::has_exception(false);

std::exception_ptr Task::GetException() { return Task::task_exception; }

bool Task::HasException() {
  // using memory_order_acquire that is sufficient to synchronize with the store in SetTaskException(), but
  // is faster than the default std::memory_order_seq_cst, which is guaranteed
  return has_exception.load(std::memory_order_acquire);
}

// As the demux module schedules tasks using silent_dependent_async() it is not
// possible to propagate exception from a task to the main thread using the
// std::future.get() mechanism. Instead tasks can call this method on exception
// which gets stored in a static variable (the first exception). Main thread
// can access it later and act on it.
void Task::SetTaskException(const std::exception_ptr& eptr) {
  // store the first exception
  // it is possible an exception AFTER the first exception gets here first, but that is okay, so long as we set a
  // TaskException
  if (!HasException()) {
    std::lock_guard lock(task_mutex);
    task_exception = eptr;
    // Guarantees that all preceding reads and writes in the current thread are visible to other threads before the
    // atomic store. Faster than the default std::memory_order_seq_cst, which is guaranteed, but sufficient to
    // synchronize with the load in HasException().
    has_exception.store(true, std::memory_order_release);
  }
}

void Task::ResetException() {
  const std::scoped_lock lock(Task::task_mutex);
  task_exception = nullptr;
  has_exception.store(false);
}
}  // namespace xoos::demux
