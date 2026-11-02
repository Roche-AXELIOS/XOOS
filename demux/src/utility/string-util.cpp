#include "utility/string-util.h"

namespace xoos::demux {
std::string_view Substr(const std::string_view& str,
                        const std::string::size_type start_pos,
                        const std::string::size_type end_pos) {
  return start_pos > end_pos ? "" : str.substr(start_pos, end_pos - start_pos);
}

std::string Reverse(const std::string_view& str) {
  return std::string{str.crbegin(), str.crend()};
}
}  // namespace xoos::demux
