#include "utility/string-util.h"

#include <algorithm>
#include <cstring>

namespace xoos::demux {
std::string_view Substr(const std::string_view& str, std::string::size_type start_pos, std::string::size_type end_pos) {
  return start_pos > end_pos ? "" : str.substr(start_pos, end_pos - start_pos);
}

std::string Reverse(const std::string_view& str) { return std::string{str.crbegin(), str.crend()}; }

std::optional<LociRange> FindExactMatch5p(const char* search_seq, uint search_seq_length, int max_wiggle_right,
                                          uint start_pos, const char* read_seq, uint seq_length) {
  auto end_search_pos = std::min(seq_length - search_seq_length, start_pos + max_wiggle_right);
  for (auto spos = start_pos; spos <= end_search_pos; spos += 1) {
    if (memcmp(search_seq, read_seq + spos, search_seq_length) == 0) {
      return std::make_optional<>(LociRange{spos, spos + search_seq_length});
    }
  }
  return std::nullopt;
}

std::optional<LociRange> FindExactMatch3p(const char* search_seq, uint search_seq_len, int max_wiggle_left,
                                          uint start_pos, const char* read_seq, uint seq_length) {
  int start_search_pos = static_cast<int>(start_pos) - static_cast<int>(search_seq_len);
  int end_search_pos = std::max(0, start_search_pos - max_wiggle_left);
  for (int spos = start_search_pos; spos >= end_search_pos; spos -= 1) {
    if (memcmp(search_seq, read_seq + spos, search_seq_len) == 0) {
      return std::make_optional<>(LociRange{static_cast<uint>(spos), static_cast<uint>(spos + search_seq_len)});
    }
  }
  return std::nullopt;
}
}  // namespace xoos::demux
