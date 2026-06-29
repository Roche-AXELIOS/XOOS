#pragma once

#include <functional>
#include <memory>
#include <mutex>
#include <vector>

namespace xoos::concurrent {

/**
 * Intended to be used as a thread_local variable, this class enables
 * enumeration of all copies of the thread_local variable from a single thread
 * by tracking the thread_local in a collection.
 *
 * A thread_local copy is not automatically removed when the thread ends, and will
 * still be enumerated.
 */
template <class T>
class EnumerableThreadLocal {
 public:
  using Ptr = std::shared_ptr<T>;

  explicit EnumerableThreadLocal(Ptr local) : _local{local} {
    std::scoped_lock elements_lock(elements_mutex);
    elements.emplace_back(local);
  }

  Ptr Local() {
    return _local;
  }

  /**
   * To safely access the list of thread_local we must hold a lock,
   * this function ensures the lock is held while the list is accessed.
   */
  void ForEach(std::function<void(const T&)> func) const {
    std::scoped_lock elements_lock(elements_mutex);
    for (const auto& element : elements) {
      func(*element);
    }
  }

  /**
   * To safely access the list of thread_local we must hold a lock,
   * this function ensures the lock is held while the list is accessed.
   */
  void ForEachNonConst(std::function<void(T&)> func) {
    std::scoped_lock elements_lock(elements_mutex);
    for (const auto& element : elements) {
      func(*element);
    }
  }

  /**
   * Removes all thread_local objects from the collection. This should be called after a collection of objects has
   * been enumerated and a new set of objects is needed of the same type. Otherwise the old objects will still be
   * in the enumeration set and will interfere with the new one when enumerated.
   **/
  static void RemoveAll() {
    std::scoped_lock elements_lock(elements_mutex);
    elements.clear();
  }

 private:
  Ptr _local;

 private:
  static std::mutex elements_mutex;
  static std::vector<Ptr> elements;
};

template <class T>
std::mutex EnumerableThreadLocal<T>::elements_mutex;

template <class T>
std::vector<typename EnumerableThreadLocal<T>::Ptr> EnumerableThreadLocal<T>::elements;

}  // namespace xoos::concurrent
