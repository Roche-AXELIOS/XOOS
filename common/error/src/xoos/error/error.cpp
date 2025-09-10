#include "xoos/error/error.h"

namespace xoos::error {

std::runtime_error Error(const std::string& fmt) {
  return std::runtime_error(fmt);
}

}  // namespace xoos::error
