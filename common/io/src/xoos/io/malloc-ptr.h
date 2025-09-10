#pragma once

#include <gsl/gsl>
#include <memory>

namespace xoos::io {

struct MallocDeleter {
  void operator()(gsl::owner<void*> ptr);
};

/**
 * @brief A unique pointer that uses free() to deallocate memory instead of delete. This smart pointer should be used
 * when the memory was allocated using malloc(), calloc(), or realloc().
 */
template <typename T>
using MallocPtr = std::unique_ptr<T, MallocDeleter>;

}  // namespace xoos::io
