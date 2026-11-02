#include "xoos/error/error.h"

#include <stdexcept>
#include <string>

namespace xoos::error {

std::runtime_error Error(const std::string& fmt) {
  return std::runtime_error(fmt);
}

}  // namespace xoos::error
