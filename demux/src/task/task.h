#pragma once
#include <atomic>
#include <exception>
#include <mutex>
#include <string>

namespace xoos::demux {
class FlowContext;  // forward declaration

// Used to specify the mode used for pairwise alignment.
enum class AlignmentMode { kPrefix, kInfix };

/// @brief Base class for all tasks - stores information that all tasks need
class Task {
 public:
  Task(FlowContext& exec, size_t batch_nr, const std::string& name)  // NOLINT
      : task_name(name), context(exec), batch_nr(batch_nr) {}

  virtual ~Task() = default;

  const std::string& Name() const { return task_name; }

  static std::exception_ptr GetException();
  static bool HasException();
  // For testing we need to be able to reset the exception state of the Task class, so we can test that exceptions are
  // properly propagated and handled.
  static void ResetException();

  // This has to be public since this need to be called outside of the task class, for example in the reader task when
  // we catch an exception in order to terminate the program gracefully.
  static void SetTaskException(const std::exception_ptr& eptr);

 protected:
  const std::string task_name;
  FlowContext& context;
  size_t batch_nr{0};

 private:
  static std::exception_ptr task_exception;
  static std::mutex task_mutex;
  static std::atomic<bool> has_exception;
};
}  // namespace xoos::demux
