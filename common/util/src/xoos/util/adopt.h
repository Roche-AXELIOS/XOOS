#pragma once

#include <memory>
#include <utility>

namespace xoos {

/// Adopt a unique_ptr from a C API into a std::shared_ptr, preserving its deleter.
///
/// Prefer this overload: the unique_ptr is the transfer token — no raw pointer escapes,
/// and if shared_ptr control-block allocation throws, the unique_ptr destructor fires
/// the deleter automatically.
template <typename T, typename Deleter>
std::shared_ptr<T> AdoptShared(std::unique_ptr<T, Deleter> ptr) {  // NOSONAR(M23_259,M23_286)
  return std::shared_ptr<T>(std::move(ptr));
}

/// Adopt a raw pointer from a C API into a std::shared_ptr with a custom deleter.
///
/// Use only when the deleter differs from the source unique_ptr's deleter (e.g., a
/// lambda). In all other cases, prefer the unique_ptr overload above.
///
/// Exception safety: the raw pointer is immediately wrapped in a unique_ptr guard so
/// the deleter fires on control-block allocation failure.
template <typename T, typename Deleter>
std::shared_ptr<T> AdoptShared(T* ptr, Deleter deleter) {  // NOSONAR(M23_259,M23_286)
  return std::shared_ptr<T>(std::unique_ptr<T, Deleter>(ptr, std::move(deleter)));
}

/// Adopt a raw pointer from a C API into a std::unique_ptr with a custom deleter.
///
/// make_unique cannot be used here: it allocates its own object. These pointers come
/// from C API functions that own their allocation.
template <typename T, typename Deleter>
std::unique_ptr<T, Deleter> AdoptUnique(T* ptr, Deleter deleter = Deleter{}) {  // NOSONAR(M23_259,M23_286)
  return std::unique_ptr<T, Deleter>(ptr, std::move(deleter));
}

}  // namespace xoos
